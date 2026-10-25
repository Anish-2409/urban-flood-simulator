"""Diagnostic hydrological and flood risk metrics.

Quantifies volume conservation, boundary losses, peak inundation depths,
and waterlogging duration hotspots.

Units:
- Depth: meters (m)
- Cell dimension dx: meters (m)
- Volume: cubic meters (m³)
- Discharge: cubic meters per second (m³/s)
- Time & duration: seconds (s)
"""

from typing import Dict, Tuple, Optional
import numpy as np
from floodsim.config import SimResult


# ---------------------------------------------------------------------------
# Primary metric functions
# ---------------------------------------------------------------------------

def stored_volume(h: np.ndarray, dx: float) -> float:
    """Calculate total stored water volume across the grid domain.

    Mathematical Formulation:
        Approximates the surface double integral of water depth h(x, y)
        over the rectangular spatial domain Omega using the midpoint / rectangle
        quadrature rule on a uniform Cartesian grid of cell area dx²:

            V = iint_Omega h(x, y) dx dy
              ≈ sum_{i=0}^{rows-1} sum_{j=0}^{cols-1} h_{i,j} * dx * dx   [m³]

        Each term h_{i,j} * dx² represents the water volume stored within
        cell (i, j) of footprint area dx × dx (m²) and average depth h_{i,j} (m).
        Summation is O(rows * cols) and exact to floating-point precision for
        uniform grid spacing.

    Discrete Approximation:
        Midpoint (rectangle) rule — constant value h_{i,j} assumed uniform
        within cell of area dx².  Error is O(dx²) relative to the true integral.

    Args:
        h (np.ndarray): 2D float64 array of shape (rows, cols) of water depth
            in meters (m).  All entries must be >= 0.
        dx (float): Grid cell edge length in meters (m).  dx > 0.

    Returns:
        float: Total stored water volume in cubic meters (m³).  Always >= 0.
    """
    return float(np.sum(h, dtype=np.float64) * dx * dx)


def boundary_discharge(
    h: np.ndarray,
    z: np.ndarray,
    dx: float,
    k_flow: float,
) -> float:
    """Calculate instantaneous volumetric water loss rate across open outer boundaries.

    Mathematical Formulation:
        Integrates the outward normal discharge flux along the closed boundary
        contour Gamma of the rectangular domain, approximated as a perimeter
        line integral (sum over boundary edges):

            Q_boundary = oint_Gamma max(0, q_n) dl
                       = sum_{b in boundary_edges} max(0, q_{out,b}) * dx   [m³/s]

        where dl = dx is the edge length of each boundary cell face, and the
        outward normal flux q_{out,b} [m²/s] is approximated using the diffusive
        wave (kinematic) approximation from solver.step:

            eta_b = h_b + z_b                         (water surface elevation)
            S_{out,b} = -grad(eta) · n_out            (outward surface slope)
            q_{out,b} = k_flow * h_b * S_{out,b}     (upwind, diffusive)

        where n_out is the unit outward normal at the respective boundary edge.
        The max(0, ·) condition ensures only genuine outflow is counted; negative
        values (inflow from outside domain) are clipped to zero.

    Discrete Approximation:
        First-order finite difference for the outward slope at each boundary:
            - Left edge   (j=0):      S_out = (eta_{i,0} - eta_{i,1})      / dx
            - Right edge  (j=cols-1): S_out = (eta_{i,cols-1} - eta_{i,cols-2}) / dx
            - Top edge    (i=0):      S_out = (eta_{0,j} - eta_{1,j})      / dx
            - Bottom edge (i=rows-1): S_out = (eta_{rows-1,j} - eta_{rows-2,j}) / dx

        Total discharge sums each face contribution * dx (m) to yield m³/s.

    Args:
        h (np.ndarray): 2D float64 array of water depth in meters (m),
            shape (rows, cols).
        z (np.ndarray): 2D float64 array of ground elevation in meters (m),
            shape (rows, cols).
        dx (float): Grid cell edge length in meters (m).
        k_flow (float): Flow routing conveyance coefficient (dimensionless).

    Returns:
        float: Total boundary discharge rate in cubic meters per second (m³/s).
            Always >= 0 (outflow only).
    """
    eta = h + z  # water surface elevation [m], shape (rows, cols)
    rows, cols = h.shape
    total_q: float = 0.0

    # ----- Left boundary (j = 0), outward normal points in -x direction -----
    if cols >= 2:
        S_left = (eta[:, 0] - eta[:, 1]) / dx          # outward slope [m/m]
        q_left = k_flow * h[:, 0] * S_left             # [m²/s]
        total_q += float(np.sum(np.maximum(0.0, q_left)) * dx)  # [m³/s]

    # ----- Right boundary (j = cols-1), outward normal points in +x direction -----
    if cols >= 2:
        S_right = (eta[:, -1] - eta[:, -2]) / dx
        q_right = k_flow * h[:, -1] * S_right
        total_q += float(np.sum(np.maximum(0.0, q_right)) * dx)

    # ----- Top boundary (i = 0), outward normal points in -y direction -----
    if rows >= 2:
        S_top = (eta[0, :] - eta[1, :]) / dx
        q_top = k_flow * h[0, :] * S_top
        total_q += float(np.sum(np.maximum(0.0, q_top)) * dx)

    # ----- Bottom boundary (i = rows-1), outward normal points in +y direction -----
    if rows >= 2:
        S_bot = (eta[-1, :] - eta[-2, :]) / dx
        q_bot = k_flow * h[-1, :] * S_bot
        total_q += float(np.sum(np.maximum(0.0, q_bot)) * dx)

    return total_q


def hotspot_duration(
    h_history: np.ndarray,
    times: np.ndarray,
    threshold: float,
) -> np.ndarray:
    """Compute total duration of inundation exceeding depth threshold for each cell.

    Mathematical Formulation:
        Calculates the time-integral of the indicator function for each cell (i, j):

            T_{i,j}(h >= threshold) = int_0^{t_end} 1_{h_{i,j}(t) >= threshold} dt   [s]

        where 1_{·} is 1 when the condition holds, 0 otherwise.

        Discretized using a piecewise-constant (rectangle) rule over the saved
        snapshot intervals.  For each consecutive pair of saved time steps
        (t_k, t_{k+1}), the interval width is Dt_k = t_{k+1} - t_k.  At the
        last saved step the interval is assumed zero (closed-right endpoint):

            T_{i,j} ≈ sum_{k=0}^{n_saved-2} 1_{h_{k,i,j} >= threshold} * (t_{k+1} - t_k)

        This is equivalent to the trapezoidal rule for a binary indicator
        (values 0 or 1) evaluated at left endpoints.

        Common threshold benchmarks:
            - 0.05 m: Pedestrian / sidewalk nuisance waterlogging
            - 0.15 m: Curb height / vehicle wheel rim danger limit
            - 0.30 m: Property damage / vehicular stalling threshold

    Discrete Approximation:
        Piecewise-constant (rectangle) rule on unequal intervals Delta_t_k.
        If the simulator saves every step, Delta_t_k = dt for all k.

    Args:
        h_history (np.ndarray): 3D float64 array of shape (n_saved, rows, cols)
            of water depth in m.
        times (np.ndarray): 1D float64 array of shape (n_saved,) of saved
            timestamps in seconds (s).  Must be monotonically non-decreasing.
        threshold (float): Critical depth threshold in meters (m).  threshold >= 0.

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) of inundation
            duration in seconds (s).  Each entry is in [0, t_end].
    """
    n_saved = len(times)
    if n_saved == 0:
        rows, cols = h_history.shape[1], h_history.shape[2]
        return np.zeros((rows, cols), dtype=np.float64)

    # Indicator: 1 where depth >= threshold, else 0.  Shape (n_saved, rows, cols)
    indicator = (h_history >= threshold).astype(np.float64)

    if n_saved == 1:
        # Only one snapshot; duration is 0 (no finite interval).
        return np.zeros(h_history.shape[1:], dtype=np.float64)

    # Time intervals between consecutive saved snapshots. Shape (n_saved-1,)
    dt_intervals = np.diff(times).astype(np.float64)  # [s]

    # Rectangle rule: weight each snapshot (except the last) by the forward interval.
    # indicator[:-1] has shape (n_saved-1, rows, cols); dt[:, None, None] broadcasts.
    duration = np.sum(
        indicator[:-1] * dt_intervals[:, np.newaxis, np.newaxis],
        axis=0,
        dtype=np.float64,
    )
    return duration


def peak_depth_map(h_history: np.ndarray) -> np.ndarray:
    """Extract spatial envelope of maximum water depth experienced throughout the simulation.

    Mathematical Formulation:
        Computes the pointwise supremum (element-wise maximum) across the
        temporal dimension (axis 0) of the depth history tensor:

            h_peak(i, j) = sup_{0 <= t <= t_end} h(t, i, j)
                         = max_{k=0,...,n_saved-1} h_{k,i,j}   [m]

        This identifies the worst-case inundation depth at every spatial
        location and is the primary indicator for flood-damage assessment.

    Discrete Approximation:
        Exact element-wise maximum over all n_saved snapshots.  No
        approximation error: the result equals the true maximum at saved
        instants (any peak between saves is not captured if save_every > 1).

    Args:
        h_history (np.ndarray): 3D float64 array of shape (n_saved, rows, cols)
            of water depth in m.

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) of peak water
            depth in meters (m).  All entries are >= 0.
    """
    return np.max(h_history, axis=0).astype(np.float64)


# ---------------------------------------------------------------------------
# Boundary volume helper (uses cum_boundary when fluxes are unavailable)
# ---------------------------------------------------------------------------

def boundary_volume_from_result(result: SimResult) -> np.ndarray:
    """Return the cumulative boundary outflow volume series from a SimResult.

    Mathematical Formulation:
        The boundary outflow volume accumulated up to each saved time t_k is:

            V_bound(t_k) = int_0^{t_k} Q_boundary(t) dt   [m³]

        In the simulator this is tracked directly as the field
        ``result.cum_boundary`` and is returned without recomputation:

            V_bound(t_k) = cum_boundary[k]   [m³]

        This function exists as a convenience wrapper for callers that only
        have a SimResult (not the raw flux arrays) and need cumulative boundary
        losses.  The underlying discrete summation in the solver is:

            V_bound(t_{k+1}) = V_bound(t_k) + Q_boundary(k) * dt   [m³]

    Args:
        result (SimResult): Completed simulation result object whose
            ``cum_boundary`` field contains the cumulative boundary discharge
            volume in cubic meters (m³) at each saved time step.

    Returns:
        np.ndarray: 1D float64 array of shape (n_saved,) of cumulative
            boundary outflow volume in cubic meters (m³).
    """
    return result.cum_boundary.astype(np.float64)


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def summary_statistics(result: SimResult, dx: float, threshold: float = 0.05) -> Dict[str, float]:
    """Compute a concise set of flood severity summary statistics from a SimResult.

    Mathematical Formulation:
        1. **Peak depth** (scalar supremum over space and time):
               h_max = max_{i,j,k} h_{k,i,j}   [m]

        2. **Maximum stored volume** (temporal maximum of V(t)):
               V_max = max_k V(t_k)   [m³]
               where V(t_k) = sum_{i,j} h_{k,i,j} * dx²   [m³]

        3. **Final stored volume** (value at last saved time):
               V_final = V(t_{n_saved-1})   [m³]

        4. **Peak volume time** (time when stored volume reaches V_max):
               t_peak = t_{k^*},  k^* = argmax_k V(t_k)   [s]

        5. **Hotspot area** (spatial extent flooded beyond threshold):
               A_hotspot = dx² * sum_{i,j} 1_{h_peak(i,j) >= threshold}   [m²]
               where h_peak(i,j) = max_k h_{k,i,j}.

    Args:
        result (SimResult): Simulation output with times, h_history, and volume.
        dx (float): Grid cell edge length in meters (m).  dx > 0.
        threshold (float): Depth threshold for hotspot classification in meters (m).
            Default 0.05 m (5 cm, pedestrian nuisance level).

    Returns:
        Dict[str, float]: Dictionary with keys:
            - ``"peak_depth_m"``       : max depth anywhere in domain [m]
            - ``"max_volume_m3"``      : maximum ponded volume [m³]
            - ``"final_volume_m3"``    : ponded volume at end of simulation [m³]
            - ``"peak_volume_time_s"`` : time of peak volume [s]
            - ``"hotspot_area_m2"``    : grid area with peak depth >= threshold [m²]
    """
    h_peak = peak_depth_map(result.h_history)

    peak_depth = float(np.max(h_peak))
    max_volume = float(np.max(result.volume))
    final_volume = float(result.volume[-1])
    peak_volume_time = float(result.times[int(np.argmax(result.volume))])
    hotspot_area = float(np.sum(h_peak >= threshold) * dx * dx)

    return {
        "peak_depth_m": peak_depth,
        "max_volume_m3": max_volume,
        "final_volume_m3": final_volume,
        "peak_volume_time_s": peak_volume_time,
        "hotspot_area_m2": hotspot_area,
    }
