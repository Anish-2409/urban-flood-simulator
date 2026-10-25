"""Sensitivity analysis and parametric sweeps for engineering design exploration.

Quantifies how system response (peak depth, final stored volume)
varies with physical and engineering control parameters:
- Rainfall intensity forcing R
- Soil infiltration capacity I_max (permeable pavement ratio)
- Storm drain conveyance capacity D_max
"""

import copy
from typing import Dict, Optional, Sequence
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from floodsim.config import SimConfig, SimResult
from floodsim import scenarios, solver

# Ensure reference solver is bound if stubbed
try:
    import check_sanity  # noqa: F401
except ImportError:
    pass


def run_parameter_sweep(
    base_cfg: SimConfig,
    param_name: str,
    values: Sequence[float],
) -> Dict[float, SimResult]:
    """Execute multiple simulations sweeping a single scalar physical parameter multiplier.

    Mathematical Context:
        Evaluates system response gradient:
            d(Metrics) / d(theta)
        where theta in {k_flow, rain, infil, drain_cap}.

    Args:
        base_cfg (SimConfig): Baseline template simulation configuration.
        param_name (str): Name of parameter to vary ('rain', 'infil', 'drain_cap', 'k_flow').
        values (Sequence[float]): Multipliers or absolute values to test.

    Returns:
        Dict[float, SimResult]: Mapping from tested parameter values to resulting SimResult.
    """
    results: Dict[float, SimResult] = {}
    for val in values:
        cfg = copy.deepcopy(base_cfg)
        if param_name == "rain":
            cfg.rain = base_cfg.rain * val
        elif param_name == "infil":
            cfg.infil = base_cfg.infil * val
        elif param_name == "drain_cap":
            cfg.drain_cap = base_cfg.drain_cap * val
        elif param_name == "k_flow":
            cfg.k_flow = float(val)
        else:
            raise ValueError(f"Unknown parameter name '{param_name}' to sweep.")

        results[float(val)] = solver.run_simulation(cfg)

    return results


def summarize_sensitivity(
    results: Dict[float, SimResult],
    dx: float,
) -> pd.DataFrame:
    """Tabulate key flood metrics across swept parameter levels into a pandas DataFrame.

    Args:
        results (Dict[float, SimResult]): Parametric simulation results.
        dx (float): Grid cell edge length in meters (m).

    Returns:
        pd.DataFrame: Table with columns ['param_value', 'peak_volume_m3', 'max_depth_m'].
    """
    rows = []
    for val, res in results.items():
        peak_vol = float(np.max(res.volume))
        max_depth = float(np.max(res.h_history))
        final_vol = float(res.volume[-1])
        rows.append({
            "param_value": val,
            "peak_volume_m3": peak_vol,
            "max_depth_m": max_depth,
            "final_volume_m3": final_vol,
        })
    return pd.DataFrame(rows)


def run_sensitivity_analysis(
    base_cfg: Optional[SimConfig] = None,
    duration_s: float = 300.0,
) -> pd.DataFrame:
    """Vary rainfall intensity, infiltration rate, and drain capacity (-50%, baseline, +50%).

    Mathematical Context:
        Performs one-at-a-time (OAT) sensitivity exploration across physical inputs:
        - Rainfall forcing R: Multiplied by {0.5, 1.0, 1.5}
        - Infiltration capacity I_max: Multiplied by {0.5, 1.0, 1.5}
        - Storm drain intake D_max: Multiplied by {0.5, 1.0, 1.5}

    Args:
        base_cfg (Optional[SimConfig]): Base configuration. Defaults to baseline scenario.
        duration_s (float): Simulation duration in seconds (default 300s for responsive evaluation).

    Returns:
        pd.DataFrame: Table of peak inundation depth (m) and final stored volume (m³)
            with columns ['Parameter', 'Variation', 'Multiplier', 'Peak_Depth_m', 'Final_Volume_m3'].
    """
    if base_cfg is None:
        base_cfg = scenarios.baseline()

    # Create template with requested duration
    template = copy.deepcopy(base_cfg)
    n_steps = max(1, int(np.ceil(duration_s / template.dt)))
    template.rain = template.rain[:n_steps].copy()
    template.t_end = float(n_steps * template.dt)

    parameters = [
        ("Rainfall Intensity", "rain"),
        ("Infiltration Rate", "infil"),
        ("Drain Capacity", "drain_cap"),
    ]

    variations = [
        ("-50%", 0.5),
        ("Baseline (0%)", 1.0),
        ("+50%", 1.5),
    ]

    records = []
    for param_label, param_key in parameters:
        for var_label, factor in variations:
            cfg_run = copy.deepcopy(template)
            if param_key == "rain":
                cfg_run.rain = template.rain * factor
            elif param_key == "infil":
                cfg_run.infil = template.infil * factor
            elif param_key == "drain_cap":
                cfg_run.drain_cap = template.drain_cap * factor

            res = solver.run_simulation(cfg_run)
            peak_depth = float(np.max(res.h_history))
            final_vol = float(res.volume[-1])

            records.append({
                "Parameter": param_label,
                "Variation": var_label,
                "Multiplier": factor,
                "Peak_Depth_m": peak_depth,
                "Final_Volume_m3": final_vol,
            })

    return pd.DataFrame(records)


def plot_sensitivity(df: pd.DataFrame) -> matplotlib.figure.Figure:
    """Plot sensitivity analysis results comparing peak depth and final volume across variations.

    Args:
        df (pd.DataFrame): Sensitivity table returned by run_sensitivity_analysis().

    Returns:
        matplotlib.figure.Figure: Two-panel standalone matplotlib Figure instance.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Urban Flood Sensitivity Analysis (-50%, Baseline, +50%)", fontsize=14, fontweight="bold")

    params = df["Parameter"].unique()
    variations = ["-50%", "Baseline (0%)", "+50%"]
    x = np.arange(len(params))
    width = 0.25

    colors = ["#2b83ba", "#abdda4", "#fdae61"]

    for i, var in enumerate(variations):
        sub_df = df[df["Variation"] == var]
        depths = [float(sub_df[sub_df["Parameter"] == p]["Peak_Depth_m"].values[0]) for p in params]
        vols = [float(sub_df[sub_df["Parameter"] == p]["Final_Volume_m3"].values[0]) for p in params]

        pos = x + (i - 1) * width
        ax1.bar(pos, depths, width, label=var, color=colors[i], edgecolor="black", alpha=0.85)
        ax2.bar(pos, vols, width, label=var, color=colors[i], edgecolor="black", alpha=0.85)

    # Subplot 1: Peak Depth
    ax1.set_title("Peak Inundation Depth by Parameter", fontsize=12)
    ax1.set_xticks(x)
    ax1.set_xticklabels(params, rotation=15, ha="right")
    ax1.set_ylabel("Peak Water Depth (m)")
    ax1.grid(True, linestyle="--", alpha=0.5, axis="y")
    ax1.legend(title="Variation")

    # Subplot 2: Final Volume
    ax2.set_title("Final Stored Volume by Parameter", fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(params, rotation=15, ha="right")
    ax2.set_ylabel("Final Stored Surface Volume (m³)")
    ax2.grid(True, linestyle="--", alpha=0.5, axis="y")
    ax2.legend(title="Variation")

    fig.tight_layout()
    return fig
