"""Configuration and result dataclasses for the urban flood simulator.

All quantities use SI units exclusively:
- Lengths, elevations, and water depths: meters (m)
- Spatial resolution: dx in meters (m)
- Time and duration: dt, t_end in seconds (s)
- Hydraulic rates (rain, infiltration, drainage): meters per second (m/s)
- Volumes: cubic meters (m³)
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class SimConfig:
    """Simulation configuration parameters for 2D urban waterlogging.

    Mathematical Context:
        Governs the discrete 2D water balance equation on a grid:
            h_new = h + dt * (R - I - D - div(q))
        under cell dimension dx (m) and discrete time step dt (s).

    Attributes:
        z (np.ndarray): 2D float64 array of shape (rows, cols) representing
            ground bed elevation above a datum in meters (m).
        dx (float): Grid cell edge length in meters (m). Cell area is dx * dx (m²).
        dt (float): Numerical time step size in seconds (s), constrained by CFL stability.
        t_end (float): Total simulation duration in seconds (s).
        rain (np.ndarray): 1D float64 array of rainfall rate R for each time step
            in meters per second (m/s). Length equals int(ceil(t_end / dt)).
        infil (np.ndarray): 2D float64 array of shape (rows, cols) denoting maximum
            soil infiltration capacity I_max in meters per second (m/s).
        drain_cap (np.ndarray): 2D float64 array of shape (rows, cols) denoting maximum
            drain removal capacity D_max in meters per second (m/s), 0 where no drain.
        k_flow (float): Flow conveyance / routing coefficient (dimensionless / scaling factor),
            governing flux speed across cell faces: q = k_flow * h * sqrt(g * h) * S or diffusive slope.
            Default is 0.25.
        closed_boundary (bool): If True, boundary cell faces have zero flux (closed basin / walls).
            If False, water may freely exit at domain edges. Default is True.
        h0 (Optional[np.ndarray]): Initial water depth in meters (m) as 2D float64 array
            of shape (rows, cols). If None, defaults to zeros everywhere.
        save_every (int): Frequency (in steps) at which state history is recorded.
            Default is 1.
    """

    z: np.ndarray
    dx: float
    dt: float
    t_end: float
    rain: np.ndarray
    infil: np.ndarray
    drain_cap: np.ndarray
    k_flow: float = 0.25
    closed_boundary: bool = True
    h0: Optional[np.ndarray] = None
    save_every: int = 1


@dataclass
class SimResult:
    """Simulation output results containing spatial state history and bulk water balance.

    Mathematical Context:
        Tracks conservation of mass across all saved time snapshots:
            V(t) = V(0) + V_rain(t) - V_infil(t) - V_drain(t) - V_boundary(t)
        where total stored volume V(t) = sum_{i,j}(h_{i,j}(t) * dx²).

    Attributes:
        times (np.ndarray): 1D float64 array of saved timestamps in seconds (s),
            length n_saved.
        h_history (np.ndarray): 3D float64 array of shape (n_saved, rows, cols)
            containing ponded water depth in meters (m) across the domain at each saved step.
        volume (np.ndarray): 1D float64 array of shape (n_saved,) containing total ponded
            water volume on the surface in cubic meters (m³).
        cum_rain (np.ndarray): 1D float64 array of shape (n_saved,) containing cumulative
            volume of rainfall added to domain in cubic meters (m³).
        cum_infil (np.ndarray): 1D float64 array of shape (n_saved,) containing cumulative
            volume of water infiltrated into the soil in cubic meters (m³).
        cum_drain (np.ndarray): 1D float64 array of shape (n_saved,) containing cumulative
            volume of water removed by drainage networks in cubic meters (m³).
        cum_boundary (np.ndarray): 1D float64 array of shape (n_saved,) containing cumulative
            volume of water discharged across outer boundaries in cubic meters (m³).
    """

    times: np.ndarray
    h_history: np.ndarray
    volume: np.ndarray
    cum_rain: np.ndarray
    cum_infil: np.ndarray
    cum_drain: np.ndarray
    cum_boundary: np.ndarray
