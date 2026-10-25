"""Sensitivity analysis and parametric sweeps for engineering design exploration.

Quantifies how system response (peak depth, total waterlogged duration, flood volume)
varies with physical and engineering control parameters:
- Soil infiltration capacity I_max (permeable pavement ratio)
- Storm drain conveyance capacity D_max
- Flow routing coefficient k_flow (hydraulic roughness / Manning's n analogue)
"""

from typing import Dict, Sequence
import pandas as pd
from floodsim.config import SimConfig, SimResult


def run_parameter_sweep(
    base_cfg: SimConfig,
    param_name: str,
    values: Sequence[float],
) -> Dict[float, SimResult]:
    """Execute multiple simulations sweeping a single scalar physical parameter.

    Mathematical Context:
        Evaluates system response gradient:
            d(Metrics) / d(theta)
        where theta in {k_flow, paved_fraction, drain_capacity_multiplier}.

    Args:
        base_cfg (SimConfig): Baseline template simulation configuration.
        param_name (str): Name of parameter to vary (e.g. 'k_flow').
        values (Sequence[float]): Numerical values to test.

    Returns:
        Dict[float, SimResult]: Mapping from tested parameter values to resulting SimResult.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("run_parameter_sweep is not yet implemented.")


def summarize_sensitivity(
    results: Dict[float, SimResult],
    dx: float,
) -> pd.DataFrame:
    """Tabulate key flood metrics across swept parameter levels into a pandas DataFrame.

    Mathematical Formulation:
        Extracts for each parameter level theta:
            - Peak ponded volume: V_peak = max_t V(t)   [m³]
            - Peak inundation depth: h_max = max_{t, i, j} h(t, i, j)   [m]
            - Inundation clearing time: t_clear = min {t > t_peak | V(t) <= 0.01 * V_peak}   [s]

    Args:
        results (Dict[float, SimResult]): Parametric simulation results.
        dx (float): Grid cell edge length in meters (m).

    Returns:
        pd.DataFrame: Table with columns ['param_value', 'peak_volume_m3', 'max_depth_m', 'clearing_time_s'].

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("summarize_sensitivity is not yet implemented.")
