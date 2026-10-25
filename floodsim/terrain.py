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


def make_terrain(
    rows: int,
    cols: int,
    kind: Literal["bowl", "slope", "random"] = "bowl",
    seed: int = 0,
) -> np.ndarray:
    """Generate synthetic ground elevation topography z(x, y).

    Mathematical Formulation:
        - "bowl": Paraboloid depression centered at (r_0, c_0):
            z(r, c) = z_min + A * ((r - r_0)² / (rows/2)² + (c - c_0)² / (cols/2)²)
          Water naturally ponds in the central depression where grad(z) -> 0.
        - "slope": Constant planar gradient from upstream to downstream:
            z(r, c) = z_top - S_0 * r * dx
          Water drains uniformly along the direction of steepest descent.
        - "random": Smooth random field generated via Gaussian convolution / Perlin-like noise
          on a uniform random grid seeded with `seed`, creating realistic micro-topography.

    Args:
        rows (int): Number of grid cells in the y-direction (rows >= 1).
        cols (int): Number of grid cells in the x-direction (cols >= 1).
        kind (Literal["bowl", "slope", "random"]): Topographic profile archetype.
        seed (int): Pseudorandom seed for reproducible noise generation.

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) with ground elevation z in meters (m).

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("make_terrain is not yet implemented.")


def make_land_cover(
    shape: Tuple[int, int],
    paved_fraction: float,
    seed: int = 0,
) -> np.ndarray:
    """Generate spatial infiltration capacity field I_max(x, y) based on urban land cover.

    Mathematical Formulation:
        Urban areas are modeled as a binary or mixed composite of permeable (green space)
        and impermeable (paved / asphalt / concrete) surfaces:
            I_max(i, j) = (1 - p(i, j)) * I_permeable + p(i, j) * I_impermeable
        where:
            - p(i, j) in {0, 1} is paved indicator with mean E[p] = paved_fraction.
            - Permeable infiltration capacity I_permeable ~ 10^-5 to 10^-4 m/s (sandy loam/turf).
            - Impermeable infiltration capacity I_impermeable ~ 0.0 m/s.

    Args:
        shape (Tuple[int, int]): Grid dimensions (rows, cols).
        paved_fraction (float): Fraction of cells that are impermeable (0.0 <= paved_fraction <= 1.0).
        seed (int): Pseudorandom seed for repeatable spatial assignment.

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) with infiltration capacity in m/s.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("make_land_cover is not yet implemented.")


def place_drains(
    shape: Tuple[int, int],
    positions: Sequence[Tuple[int, int]],
    capacity: float,
) -> np.ndarray:
    """Create a drainage inlet capacity grid D_max(x, y).

    Mathematical Formulation:
        Drain inlets act as point sinks removing water at rate D <= D_max:
            D_max(i, j) = capacity (m/s)  if (i, j) in positions
            D_max(i, j) = 0.0 (m/s)       otherwise
        where capacity expresses the inlet removal rate normalized across cell area:
            capacity = Q_drain_capacity (m³/s) / (dx * dx) (m²) = [m/s].

    Args:
        shape (Tuple[int, int]): Grid dimensions (rows, cols).
        positions (Sequence[Tuple[int, int]]): List of (row, col) coordinates with drain inlets.
        capacity (float): Maximum intake rate for each drain cell in meters per second (m/s).

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) with maximum drain intake in m/s.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("place_drains is not yet implemented.")
