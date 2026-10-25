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

import numpy as np


def stored_volume(h: np.ndarray, dx: float) -> float:
    """Calculate total stored water volume across the grid domain.

    Mathematical Formulation:
        Approximates surface integral of water depth over spatial domain Omega:
            V = iint_Omega h(x, y) dx dy approx sum_{i,j} (h_{i, j} * dx^2)   [m³]

    Args:
        h (np.ndarray): 2D float64 array of shape (rows, cols) of water depth in meters (m).
        dx (float): Grid cell edge length in meters (m).

    Returns:
        float: Total stored water volume in cubic meters (m³).

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("stored_volume is not yet implemented.")


def boundary_discharge(
    h: np.ndarray,
    z: np.ndarray,
    dx: float,
    k_flow: float,
) -> float:
    """Calculate instantaneous volumetric water loss rate across open outer boundaries.

    Mathematical Formulation:
        Integrates outward discharge flux along boundary contour Gamma:
            Q_boundary = oint_Gamma max(0, q_out) dl
                       = sum_{b in boundary edges} max(0, q_{out, b}) * dx   [m³/s]
        where outward normal flux q_out is computed from boundary water surface slope:
            eta_b = h_b + z_b
            q_{out, b} = k_flow * h_b * S_{out, b}

    Args:
        h (np.ndarray): 2D float64 array of water depth in meters (m).
        z (np.ndarray): 2D float64 array of ground elevation in meters (m).
        dx (float): Grid cell edge length in meters (m).
        k_flow (float): Flow routing conveyance coefficient.

    Returns:
        float: Total boundary discharge rate in cubic meters per second (m³/s).

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("boundary_discharge is not yet implemented.")


def hotspot_duration(
    h_history: np.ndarray,
    times: np.ndarray,
    threshold: float,
) -> np.ndarray:
    """Compute total duration of inundation exceeding depth threshold for each cell.

    Mathematical Formulation:
        Calculates time-integrated indicator function for waterlogging:
            T_{i, j}(h >= threshold) = int_0^{t_end} 1_{h_{i, j}(t) >= threshold} dt   [s]

        Using trapezoidal or piecewise-constant integration across saved time steps:
            T_{i, j} approx sum_{k} 1_{h_{k, i, j} >= threshold} * Delta t_k

        Common threshold benchmarks:
            - 0.05 m: Pedestrian / sidewalk nuisance waterlogging
            - 0.15 m: Curb height / vehicle wheel rim danger limit
            - 0.30 m: Property damage / vehicular stalling threshold

    Args:
        h_history (np.ndarray): 3D float64 array of shape (n_saved, rows, cols) of water depth in m.
        times (np.ndarray): 1D float64 array of shape (n_saved,) of saved timestamps in seconds (s).
        threshold (float): Critical depth threshold in meters (m).

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) of inundation duration in seconds (s).

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("hotspot_duration is not yet implemented.")


def peak_depth_map(h_history: np.ndarray) -> np.ndarray:
    """Extract spatial envelope of maximum water depth experienced throughout the simulation.

    Mathematical Formulation:
        Computes pointwise supremum across temporal dimension:
            h_peak(i, j) = max_{0 <= t <= t_end} h(t, i, j)   [m]

    Args:
        h_history (np.ndarray): 3D float64 array of shape (n_saved, rows, cols) of water depth in m.

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) of peak water depth in meters (m).

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("peak_depth_map is not yet implemented.")
