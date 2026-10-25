"""Visualization utilities for 2D inundation maps and time-series hydrographs.

All plotting routines return clean, standalone matplotlib Figure instances
suitable for embedding in reports, notebooks, or Streamlit dashboards.

The Agg backend is selected at import time when no interactive display is
available (i.e., when DISPLAY / WAYLAND_DISPLAY environment variables are
absent), so every function is safe to call in headless CI/server environments.

Units throughout:
- Depth h: meters (m)
- Elevation z: meters (m)
- Volume: cubic meters (m³)
- Time: seconds (s)
- Area: square meters (m²)
"""

import os
from typing import Dict, Optional

import numpy as np

# ------------------------------------------------------------------
# Select a safe backend before importing pyplot so that the module
# can be imported in headless / CI environments without an X display.
# ------------------------------------------------------------------
import matplotlib
_HAS_DISPLAY = bool(
    os.environ.get("DISPLAY", "") or os.environ.get("WAYLAND_DISPLAY", "")
)
if not _HAS_DISPLAY:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.figure
import matplotlib.colors as mcolors

from floodsim.config import SimResult


# ---------------------------------------------------------------------------
# Depth Map
# ---------------------------------------------------------------------------

def plot_depth_map(
    h: np.ndarray,
    z: Optional[np.ndarray] = None,
    title: str = "Water Depth Map",
) -> matplotlib.figure.Figure:
    """Render a 2D spatial colormap of water depth with optional topographic contour overlay.

    Mathematical Context:
        Displays the scalar field h(x, y) [m] as a discrete colormap on a
        uniform Cartesian grid.  Dry cells (h < 1e-4 m) are rendered as
        transparent / white using a masked array so the underlying topography
        remains visible.

        If bed elevation z(x, y) [m] is provided, contour lines of z are
        overlaid.  Flow paths follow the steepest descent of the total water
        surface eta = h + z, but iso-contours of z alone indicate the
        underlying topographic structure of the domain.

        Colorbar ticks are in meters (m) and the colormap 'Blues' encodes
        shallow depths as light blue and deep inundation as dark blue,
        consistent with standard hydrological convention.

    Args:
        h (np.ndarray): 2D float64 array of shape (rows, cols) of water depth
            in meters (m).
        z (Optional[np.ndarray]): Optional 2D float64 array of ground elevation
            in meters (m), same shape as h.  If provided, contour lines of z
            are drawn in grey.
        title (str): Title string for the figure panel.

    Returns:
        matplotlib.figure.Figure: Generated matplotlib Figure object (no
            plt.show() call is made; caller controls display / saving).
    """
    fig, ax = plt.subplots(figsize=(7, 5.5))

    # Mask dry cells so they render as white/transparent
    dry = h < 1e-4
    h_masked = np.ma.array(h, mask=dry)

    vmax = float(np.nanmax(h)) if float(np.nanmax(h)) > 0 else 1.0
    cmap = plt.cm.Blues.with_extremes(bad="whitesmoke")

    im = ax.imshow(
        h_masked,
        origin="upper",
        cmap=cmap,
        vmin=0.0,
        vmax=vmax,
        aspect="equal",
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Water Depth (m)", fontsize=10)

    # Optional elevation contours
    if z is not None:
        rows, cols = z.shape
        contour = ax.contour(
            np.arange(cols),
            np.arange(rows),
            z,
            levels=8,
            colors="gray",
            linewidths=0.6,
            alpha=0.55,
        )
        ax.clabel(contour, inline=True, fontsize=7, fmt="%.2f m")

    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Column index (j)", fontsize=10)
    ax.set_ylabel("Row index (i)", fontsize=10)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Hotspot Map
# ---------------------------------------------------------------------------

def plot_hotspot_map(
    hotspot: np.ndarray,
    title: str = "Flood Duration Hotspot Map",
) -> matplotlib.figure.Figure:
    """Render a 2D spatial map of per-cell flood inundation duration.

    Mathematical Context:
        Displays the inundation duration field T_{i,j} [s] computed by
        metrics.hotspot_duration.  For cell (i, j):

            T_{i,j} = int_0^{t_end} 1_{h_{i,j}(t) >= threshold} dt   [s]

        Cells with T_{i,j} = 0 (never flooded beyond threshold) are rendered
        white.  Warm colours (yellow → red) indicate progressively longer
        sustained flooding, which is proportional to waterlogging risk.

        Duration values are converted to minutes for readability by dividing
        by 60.

    Args:
        hotspot (np.ndarray): 2D float64 array of shape (rows, cols) of
            inundation duration in seconds (s).  Returned by hotspot_duration.
        title (str): Title string for the figure panel.

    Returns:
        matplotlib.figure.Figure: Generated matplotlib Figure object.
    """
    fig, ax = plt.subplots(figsize=(7, 5.5))

    hotspot_min = hotspot / 60.0  # convert s → min for readability
    vmax = float(np.max(hotspot_min)) if float(np.max(hotspot_min)) > 0 else 1.0

    cmap = plt.cm.YlOrRd.with_extremes(bad="whitesmoke")
    h_masked = np.ma.array(hotspot_min, mask=(hotspot_min <= 0))

    im = ax.imshow(
        h_masked,
        origin="upper",
        cmap=cmap,
        vmin=0.0,
        vmax=vmax,
        aspect="equal",
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Flood Duration (min)", fontsize=10)

    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Column index (j)", fontsize=10)
    ax.set_ylabel("Row index (i)", fontsize=10)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Time-Series Hydrograph
# ---------------------------------------------------------------------------

def plot_timeseries(
    times: np.ndarray,
    volume: np.ndarray,
    cum_rain: np.ndarray,
    cum_infil: np.ndarray,
    cum_drain: np.ndarray,
) -> matplotlib.figure.Figure:
    """Plot global mass balance hydrograph and cumulative volume conservation components.

    Mathematical Context:
        Displays the dynamic conservation of mass over time.  The stored
        surface volume V_stored(t) must satisfy:

            V_stored(t) = V_stored(0)
                        + V_rain(t) - V_infil(t) - V_drain(t) - V_bound(t)

        where:
            - V_stored(t)  = sum_{i,j} h_{i,j}(t) * dx²          [m³]
            - V_rain(t)    = cumulative rain volume added          [m³]
            - V_infil(t)   = cumulative infiltration removed       [m³]
            - V_drain(t)   = cumulative drain extraction removed   [m³]
            - V_bound(t)   = cumulative boundary outflow           [m³]

        The gap between the rainfall input curve and the sum of sink curves
        corresponds to ponded storage, allowing visual inspection of hydrograph
        timing: rising limb during storm, recession limb during drainage.

        Time axis is converted to minutes for readability (÷ 60).

    Args:
        times (np.ndarray): 1D float64 array of timestamps in seconds (s).
        volume (np.ndarray): 1D float64 array of total stored surface volume
            in cubic meters (m³).
        cum_rain (np.ndarray): 1D float64 array of cumulative rain volume
            in cubic meters (m³).
        cum_infil (np.ndarray): 1D float64 array of cumulative infiltration
            in cubic meters (m³).
        cum_drain (np.ndarray): 1D float64 array of cumulative drainage
            in cubic meters (m³).

    Returns:
        matplotlib.figure.Figure: Generated matplotlib Figure object.
    """
    times_min = times / 60.0  # s → min

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    # ---- Top panel: stored surface volume ----
    ax_top = axes[0]
    ax_top.plot(times_min, volume, color="steelblue", linewidth=2.0,
                label="Stored Surface Volume")
    ax_top.set_ylabel("Volume (m³)", fontsize=10)
    ax_top.set_title("Surface Ponded Volume over Time", fontsize=11, fontweight="bold")
    ax_top.legend(fontsize=9, loc="upper right")
    ax_top.grid(True, linestyle="--", alpha=0.4)
    ax_top.set_ylim(bottom=0)

    # ---- Bottom panel: cumulative mass components ----
    ax_bot = axes[1]
    ax_bot.plot(times_min, cum_rain, color="royalblue", linewidth=1.8,
                label="Cumul. Rainfall (m³)")
    ax_bot.plot(times_min, cum_infil, color="forestgreen", linewidth=1.8,
                linestyle="--", label="Cumul. Infiltration (m³)")
    ax_bot.plot(times_min, cum_drain, color="darkorange", linewidth=1.8,
                linestyle="-.", label="Cumul. Drainage (m³)")
    ax_bot.set_ylabel("Cumulative Volume (m³)", fontsize=10)
    ax_bot.set_xlabel("Time (min)", fontsize=10)
    ax_bot.set_title("Cumulative Mass Balance Components", fontsize=11, fontweight="bold")
    ax_bot.legend(fontsize=9, loc="upper left")
    ax_bot.grid(True, linestyle="--", alpha=0.4)
    ax_bot.set_ylim(bottom=0)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Scenario Comparison
# ---------------------------------------------------------------------------

def plot_scenarios(
    results: Dict[str, SimResult],
) -> matplotlib.figure.Figure:
    """Generate comparative hydrograph curves across multiple policy / intervention scenarios.

    Mathematical Context:
        Compares surface ponding volume V_stored(t) [m³] and peak depth
        progression h_max(t) = max_{i,j} h_{i,j}(t) [m] across different
        engineering interventions (e.g., Baseline vs Improved Drains vs
        Green Infrastructure).

        Metrics plotted for each scenario s:
            V_s(t)     = stored surface volume at time t   [m³]
            h_max_s(t) = max_{i,j} h_{i,j}(t)            [m]

        Comparing these curves quantifies:
            - Peak attenuation:  ΔV_max = V_max^{base} - V_max^{int}  [m³]
            - Timing shift:      Δt_peak = t_peak^{int} - t_peak^{base} [s]
            - Recession rate:    slope of V_s(t) after storm cessation

        Time axis is converted to minutes (÷ 60).

    Args:
        results (Dict[str, SimResult]): Mapping from scenario names (str) to
            their respective SimResult objects.  Expected keys include at least
            one of: "Baseline", "Improved Drains", "Green Infrastructure".

    Returns:
        matplotlib.figure.Figure: Generated matplotlib Figure object with
            two stacked panels (stored volume and max depth).
    """
    # Colour palette — cycle through a consistent set for up to 6 scenarios
    _PALETTE = [
        "#e41a1c", "#377eb8", "#4daf4a",
        "#984ea3", "#ff7f00", "#a65628",
    ]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    ax_vol = axes[0]
    ax_dep = axes[1]

    for idx, (name, result) in enumerate(results.items()):
        color = _PALETTE[idx % len(_PALETTE)]
        times_min = result.times / 60.0  # s → min

        # Stored volume curve
        ax_vol.plot(
            times_min, result.volume,
            label=name, color=color, linewidth=2.0,
        )

        # Max depth across domain at each saved step
        h_max_series = np.max(
            result.h_history.reshape(result.h_history.shape[0], -1),
            axis=1,
        )
        ax_dep.plot(
            times_min, h_max_series,
            label=name, color=color, linewidth=2.0,
        )

    ax_vol.set_ylabel("Stored Volume (m³)", fontsize=10)
    ax_vol.set_title("Scenario Comparison — Stored Surface Volume", fontsize=11, fontweight="bold")
    ax_vol.legend(fontsize=9, loc="upper right")
    ax_vol.grid(True, linestyle="--", alpha=0.4)
    ax_vol.set_ylim(bottom=0)

    ax_dep.set_ylabel("Max Depth (m)", fontsize=10)
    ax_dep.set_xlabel("Time (min)", fontsize=10)
    ax_dep.set_title("Scenario Comparison — Peak Inundation Depth", fontsize=11, fontweight="bold")
    ax_dep.legend(fontsize=9, loc="upper right")
    ax_dep.grid(True, linestyle="--", alpha=0.4)
    ax_dep.set_ylim(bottom=0)

    fig.tight_layout()
    return fig
