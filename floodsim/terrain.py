"""Terrain generation and surface boundary characteristics.

Provides synthetic digital elevation models (DEM), infiltration maps based on land cover,
and storm drain placement grids.

Units:
- Elevations z: meters (m), relative to reference datum.
- Infiltration capacity: meters per second (m/s).
- Drain capacity: meters per second (m/s) per cell.
"""

from typing import Literal, Sequence, Tuple
import numpy as np


def _smooth_field_2d(arr: np.ndarray, passes: int = 2) -> np.ndarray:
    """Apply separable 2D binomial/Gaussian smoothing filter to a 2D array.

    Mathematical Formulation:
        Uses discrete 5-point binomial weights K = [1, 4, 6, 4, 1] / 16
        applied separably along rows and columns with edge replication boundary padding:
            A_{i, j}^{(k+1)} = sum_{m=-2}^{2} sum_{n=-2}^{2} K_m K_n A_{i+m, j+n}^{(k)}

    Args:
        arr (np.ndarray): 2D float64 array of shape (rows, cols).
        passes (int): Number of smoothing passes (default 2).

    Returns:
        np.ndarray: Smoothed 2D float64 array of same shape (rows, cols).
    """
    out = arr.astype(np.float64, copy=True)
    rows, cols = out.shape

    if rows < 2 and cols < 2:
        return out

    for _ in range(passes):
        # Smooth along rows if rows >= 2
        if rows >= 2:
            p_rows = np.pad(out, ((2, 2), (0, 0)), mode="edge")
            out = (
                1.0 * p_rows[:-4, :]
                + 4.0 * p_rows[1:-3, :]
                + 6.0 * p_rows[2:-2, :]
                + 4.0 * p_rows[3:-1, :]
                + 1.0 * p_rows[4:, :]
            ) / 16.0

        # Smooth along cols if cols >= 2
        if cols >= 2:
            p_cols = np.pad(out, ((0, 0), (2, 2)), mode="edge")
            out = (
                1.0 * p_cols[:, :-4]
                + 4.0 * p_cols[:, 1:-3]
                + 6.0 * p_cols[:, 2:-2]
                + 4.0 * p_cols[:, 3:-1]
                + 1.0 * p_cols[:, 4:]
            ) / 16.0

    return out


def make_terrain(
    rows: int = 30,
    cols: int = 30,
    kind: Literal["bowl", "slope", "random"] = "bowl",
    seed: int = 0,
) -> np.ndarray:
    """Generate synthetic ground elevation topography z(x, y).

    Mathematical Formulation:
        - "bowl": Paraboloid depression centered at (r_0, c_0) with low center:
            dr = (r - r_0) / (rows / 2)
            dc = (c - c_0) / (cols / 2)
            z(r, c) = z_min + A * (dr² + dc²)
          where z_min = 0.0 m, A = 1.0 m, r_0 = (rows - 1)/2, c_0 = (cols - 1)/2.
          Water naturally ponds in the central depression where grad(z) -> 0.
        - "slope": Constant planar gradient from upstream (row 0) to downstream (row rows-1):
            z(r, c) = z_top * (rows - 1 - r) / (rows - 1)
          Water drains uniformly along the direction of steepest descent.
        - "random": Smooth continuous random field generated via Gaussian convolution
          on a uniform random grid seeded with `seed`, normalized to [0, 1.5] meters,
          creating realistic micro-topography.

    Args:
        rows (int): Number of grid cells in the y-direction (rows >= 1, default 30).
        cols (int): Number of grid cells in the x-direction (cols >= 1, default 30).
        kind (Literal["bowl", "slope", "random"]): Topographic profile archetype (default "bowl").
        seed (int): Pseudorandom seed for reproducible noise generation (default 0).

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) with ground elevation z in meters (m).

    Raises:
        ValueError: If rows or cols are not positive, or if kind is unrecognized.
    """
    if rows <= 0 or cols <= 0:
        raise ValueError(f"Grid dimensions must be positive integers, got rows={rows}, cols={cols}")

    if kind == "bowl":
        r_0 = (rows - 1) / 2.0
        c_0 = (cols - 1) / 2.0
        r_scale = max(rows / 2.0, 1.0)
        c_scale = max(cols / 2.0, 1.0)

        r_coords = np.arange(rows, dtype=np.float64)[:, None]
        c_coords = np.arange(cols, dtype=np.float64)[None, :]

        dr = (r_coords - r_0) / r_scale
        dc = (c_coords - c_0) / c_scale
        z = dr**2 + dc**2
        return z.astype(np.float64)

    elif kind == "slope":
        if rows > 1:
            z_row = np.linspace(1.0, 0.0, rows, dtype=np.float64)[:, None]
        else:
            z_row = np.zeros((1, 1), dtype=np.float64)
        z = np.tile(z_row, (1, cols))
        return z.astype(np.float64)

    elif kind == "random":
        rng = np.random.default_rng(seed)
        raw = rng.uniform(0.0, 1.0, size=(rows, cols))
        smoothed = _smooth_field_2d(raw, passes=2)
        s_min = float(smoothed.min())
        s_max = float(smoothed.max())
        if s_max > s_min:
            z = (smoothed - s_min) / (s_max - s_min) * 1.5
        else:
            z = np.zeros((rows, cols), dtype=np.float64)
        return z.astype(np.float64)

    else:
        raise ValueError(f"Unknown terrain kind: '{kind}'. Expected 'bowl', 'slope', or 'random'.")


def make_land_cover(
    shape: Tuple[int, int],
    paved_fraction: float,
    seed: int = 0,
) -> np.ndarray:
    """Generate spatial infiltration capacity field I_max(x, y) based on urban land cover.

    Mathematical Formulation:
        Urban areas are modeled as a spatial composite of permeable (green space)
        and impermeable (paved / asphalt / concrete) surfaces:
            I_max(i, j) = p(i, j) * I_paved(i, j) + (1 - p(i, j)) * I_green(i, j)
        where:
            - p(i, j) in {0, 1} is paved indicator with expectation E[p] = paved_fraction.
            - Paved cells have very low infiltration ~ 0 to 1e-7 m/s (asphalt/concrete).
            - Soil/green cells have higher infiltration ~ 1e-6 to 5e-6 m/s (turf, loam).

    Args:
        shape (Tuple[int, int]): Grid dimensions (rows, cols).
        paved_fraction (float): Fraction of cells that are impermeable (0.0 <= paved_fraction <= 1.0).
        seed (int): Pseudorandom seed for repeatable spatial assignment (default 0).

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) with infiltration capacity in m/s.

    Raises:
        ValueError: If shape dimensions are non-positive or paved_fraction is not in [0.0, 1.0].
    """
    if len(shape) != 2 or shape[0] <= 0 or shape[1] <= 0:
        raise ValueError(f"shape must be a 2-tuple of positive integers, got {shape}")
    if not (0.0 <= paved_fraction <= 1.0):
        raise ValueError(f"paved_fraction must be in range [0.0, 1.0], got {paved_fraction}")

    rng = np.random.default_rng(seed)
    is_paved = rng.random(size=shape) < paved_fraction
    paved_infil = rng.uniform(0.0, 1e-7, size=shape)
    green_infil = rng.uniform(1e-6, 5e-6, size=shape)

    infil = np.where(is_paved, paved_infil, green_infil)
    return infil.astype(np.float64)


def place_drains(
    shape: Tuple[int, int],
    positions: Sequence[Tuple[int, int]],
    capacity: float,
) -> np.ndarray:
    """Create a drainage inlet capacity grid D_max(x, y).

    Mathematical Formulation:
        Drain inlets act as localized point sinks removing ponded water at rate D <= D_max:
            D_max(i, j) = capacity (m/s)  if (i, j) in positions
            D_max(i, j) = 0.0 (m/s)       otherwise
        where capacity expresses the inlet removal rate normalized across cell area:
            capacity = Q_drain_capacity (m³/s) / (dx * dx) (m²) = [m/s].

    Args:
        shape (Tuple[int, int]): Grid dimensions (rows, cols).
        positions (Sequence[Tuple[int, int]]): Sequence of (row, col) coordinates with drain inlets.
        capacity (float): Maximum intake rate for each drain cell in meters per second (m/s).

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) with maximum drain intake in m/s.

    Raises:
        ValueError: If shape dimensions are non-positive or capacity is negative.
        IndexError: If any position is outside domain bounds.
    """
    if len(shape) != 2 or shape[0] <= 0 or shape[1] <= 0:
        raise ValueError(f"shape must be a 2-tuple of positive integers, got {shape}")
    if capacity < 0:
        raise ValueError(f"capacity must be non-negative, got {capacity}")

    drain_cap = np.zeros(shape, dtype=np.float64)
    rows, cols = shape
    for pos in positions:
        r, c = pos
        if not (0 <= r < rows and 0 <= c < cols):
            raise IndexError(f"Drain position {pos} is outside grid bounds {shape}")
        drain_cap[r, c] = float(capacity)

    return drain_cap
