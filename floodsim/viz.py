"""Visualization utilities for 2D inundation maps and time-series hydrographs.

All plotting routines return clean, standalone matplotlib Figure instances
suitable for embedding in reports, notebooks, or Streamlit dashboards.
"""

from typing import Dict, Optional
import matplotlib.figure
import numpy as np
from floodsim.config import SimResult


def plot_depth_map(
    h: np.ndarray,
    z: Optional[np.ndarray] = None,
    title: str = "Water Depth Map",
) -> matplotlib.figure.Figure:
    """Render a 2D spatial colormap of water depth with optional topographic contour overlay.

    Mathematical Context:
        Displays scalar field h(x, y) [m] mapped to a perceptually uniform colormap
        (e.g., 'Blues' or 'YlGnBu') masked where h < 1e-4 m (dry ground).
        If bed elevation z(x, y) is provided, elevation contours are overlaid
        to visualize hydraulic flow paths down steepest topographic gradients.

    Args:
        h (np.ndarray): 2D float64 array of shape (rows, cols) of water depth in meters (m).
        z (Optional[np.ndarray]): Optional 2D float64 array of ground elevation in meters (m).
        title (str): Title for the figure.

    Returns:
        matplotlib.figure.Figure: Generated matplotlib Figure object.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("plot_depth_map is not yet implemented.")


def plot_timeseries(
    times: np.ndarray,
    volume: np.ndarray,
    cum_rain: np.ndarray,
    cum_infil: np.ndarray,
    cum_drain: np.ndarray,
) -> matplotlib.figure.Figure:
    """Plot global mass balance hydrograph and cumulative volume conservation components.

    Mathematical Context:
        Displays dynamic conservation of mass over time:
            V_stored(t) vs [V_rain(t) - V_infil(t) - V_drain(t)]
        Illustrates how surface ponding peaks and recedes as rainfall ceases
        and infiltration / drainage sinks absorb standing water.

    Args:
        times (np.ndarray): 1D float64 array of timestamps in seconds (s).
        volume (np.ndarray): 1D float64 array of total stored surface volume in m³.
        cum_rain (np.ndarray): 1D float64 array of cumulative rain volume in m³.
        cum_infil (np.ndarray): 1D float64 array of cumulative infiltration in m³.
        cum_drain (np.ndarray): 1D float64 array of cumulative drainage in m³.

    Returns:
        matplotlib.figure.Figure: Generated matplotlib Figure object.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("plot_timeseries is not yet implemented.")


def plot_scenarios(
    results: Dict[str, SimResult],
) -> matplotlib.figure.Figure:
    """Generate comparative hydrograph curves across multiple policy / intervention scenarios.

    Mathematical Context:
        Compares surface ponding volume V_stored(t) and maximum depth progression
        across different interventions (e.g., Baseline vs Improved Drains vs Green Infrastructure),
        enabling engineering trade-off evaluations of peak attenuation and drainage lag.

    Args:
        results (Dict[str, SimResult]): Mapping from scenario names to their respective SimResult.

    Returns:
        matplotlib.figure.Figure: Generated matplotlib Figure object.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("plot_scenarios is not yet implemented.")
