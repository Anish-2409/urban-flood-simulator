"""Discrete spatial differential operators on 2D uniform Cartesian grids.

Implements finite-difference and finite-volume spatial operators:
- Spatial gradient of scalar fields (elevation z or water surface elevation eta = h + z)
- Flux divergence operator ensuring strict discrete mass conservation via shared cell-edge fluxes.

Units:
- Grid cell spacing dx: meters (m)
- Elevation / depth: meters (m)
- Fluxes qx, qy: meters squared per second (m²/s) [discharge per unit width]
- Divergence div(q): meters per second (m/s) [rate of depth change]
"""

from typing import Tuple
import numpy as np


def gradient(z: np.ndarray, dx: float) -> Tuple[np.ndarray, np.ndarray]:
    """Compute discrete spatial gradient of a 2D scalar field: grad(z) = (dz/dx, dz/dy).

    Mathematical Formulation:
        For scalar field z of shape (rows, cols) and cell resolution dx:
        Interior cells use second-order central differences:
            dz/dx(i, j) = (z(i, j+1) - z(i, j-1)) / (2 * dx)
            dz/dy(i, j) = (z(i+1, j) - z(i-1, j)) / (2 * dx)
        Boundaries use one-sided first-order differences:
            dz/dx(i, 0) = (z(i, 1) - z(i, 0)) / dx
            dz/dx(i, cols-1) = (z(i, cols-1) - z(i, cols-2)) / dx
        (similarly for dz/dy along rows).

    Args:
        z (np.ndarray): 2D float64 array of shape (rows, cols) representing scalar field (m).
        dx (float): Grid cell edge length in meters (m).

    Returns:
        Tuple[np.ndarray, np.ndarray]:
            dzdx: Spatial derivative along columns (x-direction, m/m, dimensionless).
            dzdy: Spatial derivative along rows (y-direction, m/m, dimensionless).

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("gradient is not yet implemented.")


def divergence(qx: np.ndarray, qy: np.ndarray, dx: float) -> np.ndarray:
    """Compute the discrete spatial divergence of 2D face fluxes: div(q) = dqx/dx + dqy/dy.

    Mathematical Formulation:
        Using finite volume balance over cell (i, j) of area dx²:
            div(q)_{i, j} = (q_{x, i, j+1/2} - q_{x, i, j-1/2}) / dx
                          + (q_{y, i+1/2, j} - q_{y, i-1/2, j}) / dx

        Conservation Property:
            Because cell-edge fluxes q_{x, i, j+1/2} and q_{y, i+1/2, j} are shared identically
            between adjacent cells with opposing normal vectors, the domain integral:
                sum_{i,j} div(q)_{i, j} * dx² = sum(boundary fluxes) * dx
            When closed boundary conditions are enforced, sum_{i,j} div(q)_{i, j} = 0
            identically to machine precision.

    Args:
        qx (np.ndarray): Net flux in x-direction across cell faces or cell centers (m²/s).
        qy (np.ndarray): Net flux in y-direction across cell faces or cell centers (m²/s).
        dx (float): Grid cell edge length in meters (m).

    Returns:
        np.ndarray: 2D float64 array of shape (rows, cols) representing div(q) in m/s.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("divergence is not yet implemented.")
