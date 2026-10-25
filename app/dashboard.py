"""Streamlit interactive web dashboard for the 2D urban flood simulator.

Educational tool for first-year engineering mathematics / numerical methods.
Run via:
    streamlit run app/dashboard.py
"""

import copy
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Fallback reference solver if solver is stubbed
try:
    import check_sanity  # noqa: F401
except ImportError:
    pass

from floodsim.config import SimConfig, SimResult
from floodsim import rainfall, scenarios, solver, terrain

MM_HR_TO_M_S = 1e-3 / 3600.0  # 1 mm/hr in m/s (~2.778e-7 m/s)


def run_dashboard_simulation(
    scenario_name: str,
    intensity_mm_hr: float,
    duration_min: float,
    drain_cap_mm_hr: float,
    paved_fraction: float,
    grid_size: int = 20,
    sim_duration_min: float = 30.0,
) -> tuple[SimConfig, SimResult]:
    """Build configuration from interactive controls and execute 2D simulation."""
    dx = 2.0  # meters
    dt = 1.0  # seconds
    t_end = sim_duration_min * 60.0  # seconds
    grid_shape = (grid_size, grid_size)

    # 1. Topography
    z = terrain.make_terrain(rows=grid_size, cols=grid_size, kind="bowl", seed=42)

    # 2. Rainfall forcing
    rain = rainfall.constant_storm(
        intensity_mm_hr=intensity_mm_hr,
        duration_min=duration_min,
        dt=dt,
        t_end=t_end,
    )

    # 3. Land cover infiltration
    infil = terrain.make_land_cover(shape=grid_shape, paved_fraction=paved_fraction, seed=42)

    # 4. Drainage placement
    r_mid, c_mid = grid_size // 2, grid_size // 2
    if scenario_name == "improved_drains":
        # 8 drain inlets surrounding central bowl
        drain_positions = [
            (r_mid - 1, c_mid - 1), (r_mid - 1, c_mid),
            (r_mid, c_mid - 1), (r_mid, c_mid),
            (r_mid - 2, c_mid - 1), (r_mid - 2, c_mid),
            (r_mid + 1, c_mid - 1), (r_mid + 1, c_mid),
        ]
        drain_cap_rate = (drain_cap_mm_hr * 1.5) * MM_HR_TO_M_S
    else:
        # Standard central 4 inlets
        drain_positions = [
            (r_mid - 1, c_mid - 1), (r_mid - 1, c_mid),
            (r_mid, c_mid - 1), (r_mid, c_mid),
        ]
        drain_cap_rate = drain_cap_mm_hr * MM_HR_TO_M_S

    drain_cap = terrain.place_drains(
        shape=grid_shape,
        positions=drain_positions,
        capacity=drain_cap_rate,
    )

    # 5. Sponge City green infrastructure modification
    if scenario_name == "green_infrastructure":
        rng = np.random.default_rng(101)
        paved_mask = infil <= 1e-7
        convert_mask = paved_mask & (rng.random(grid_shape) < 0.6)
        infil[convert_mask] = rng.uniform(2e-6, 6e-6, size=int(convert_mask.sum()))

    cfg = SimConfig(
        z=z,
        dx=dx,
        dt=dt,
        t_end=t_end,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
        k_flow=0.25,
        closed_boundary=True,
        save_every=max(1, int(t_end / 100)),  # ~100 time snapshots
    )

    res = solver.run_simulation(cfg)
    return cfg, res


def main() -> None:
    """Render the Streamlit dashboard."""
    st.set_page_config(
        page_title="Urban Flood Simulator 2D",
        page_icon="🌧️",
        layout="wide",
    )

    st.title("🌧️ Urban Flood Simulator — 2D Waterlogging Dashboard")
    st.markdown(
        "> **Educational Use Only**: 2D discrete cell-water balance model based on "
        r"$\frac{\partial h}{\partial t} = R - I - D - \nabla \cdot \mathbf{q}$. "
        "Outputs are scenario estimates for learning, not regulatory flood forecasts."
    )

    # Sidebar: Interactive Controls
    st.sidebar.header("🕹️ Simulation Controls")

    scenario_choice = st.sidebar.selectbox(
        "Scenario Selector",
        ["baseline", "improved_drains", "green_infrastructure"],
        format_func=lambda x: {
            "baseline": "Baseline (Impervious Urban)",
            "improved_drains": "Improved Drains (Expanded Grey Infrastructure)",
            "green_infrastructure": "Green Infrastructure (Sponge City)",
        }[x],
    )

    st.sidebar.subheader("🌧️ Rainfall Parameters")
    intensity_mm_hr = st.sidebar.slider(
        "Rainfall Intensity (mm/hr)",
        min_value=10.0,
        max_value=120.0,
        value=50.0,
        step=5.0,
    )
    duration_min = st.sidebar.slider(
        "Storm Duration (minutes)",
        min_value=5.0,
        max_value=60.0,
        value=20.0,
        step=5.0,
    )

    st.sidebar.subheader("🏗️ Drainage & Surface Controls")
    drain_cap_mm_hr = st.sidebar.slider(
        "Drain Capacity (mm/hr)",
        min_value=0.0,
        max_value=100.0,
        value=25.0,
        step=5.0,
    )
    paved_fraction = st.sidebar.slider(
        "Paved Fraction (Impervious)",
        min_value=0.0,
        max_value=1.0,
        value=0.80 if scenario_choice != "green_infrastructure" else 0.35,
        step=0.05,
    )

    st.sidebar.subheader("⚙️ Computational Domain")
    grid_res = st.sidebar.selectbox("Grid Resolution", [20, 30], index=0)
    sim_duration_min = st.sidebar.slider("Simulation Length (min)", 15.0, 60.0, 25.0, 5.0)

    # Execute simulation
    with st.spinner("Computing 2D overland water balance..."):
        cfg, res = run_dashboard_simulation(
            scenario_name=scenario_choice,
            intensity_mm_hr=intensity_mm_hr,
            duration_min=duration_min,
            drain_cap_mm_hr=drain_cap_mm_hr,
            paved_fraction=paved_fraction,
            grid_size=grid_res,
            sim_duration_min=sim_duration_min,
        )

    # Summary Metrics Row
    peak_depth = float(np.max(res.h_history))
    peak_volume = float(np.max(res.volume))
    total_rain_vol = float(res.cum_rain[-1])
    total_infil_vol = float(res.cum_infil[-1])
    total_drain_vol = float(res.cum_drain[-1])
    final_stored = float(res.volume[-1])
    mass_balance_error = abs(total_rain_vol - total_infil_vol - total_drain_vol - final_stored)

    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("Peak Depth", f"{peak_depth * 100:.1f} cm", help="Maximum ponded water depth across all cells")
    col2.metric("Peak Surface Volume", f"{peak_volume:.2f} m³", help="Maximum total ponded water stored")
    col3.metric("Rainfall Volume", f"{total_rain_vol:.2f} m³", help="Total precipitation deposited")
    col4.metric("Infiltrated Volume", f"{total_infil_vol:.2f} m³", help="Water absorbed into soil")
    col5.metric("Drained Volume", f"{total_drain_vol:.2f} m³", help="Water captured by storm inlets")

    st.divider()

    # Visualizations: 2 Columns
    viz_col1, viz_col2 = st.columns(2)

    with viz_col1:
        st.subheader("🗺️ Peak Water Depth Heatmap")
        fig_heat, ax_heat = plt.subplots(figsize=(6, 5))
        peak_map = np.max(res.h_history, axis=0)

        im = ax_heat.imshow(
            peak_map,
            cmap="Blues",
            origin="lower",
            extent=[0, cfg.dx * cfg.z.shape[1], 0, cfg.dx * cfg.z.shape[0]],
            vmin=0.0,
        )
        cbar = fig_heat.colorbar(im, ax=ax_heat)
        cbar.set_label("Depth h (m)", fontsize=10)

        # Topographic elevation contours
        cs = ax_heat.contour(
            cfg.z,
            levels=6,
            colors="black",
            origin="lower",
            extent=[0, cfg.dx * cfg.z.shape[1], 0, cfg.dx * cfg.z.shape[0]],
            linewidths=0.7,
            alpha=0.6,
        )
        ax_heat.clabel(cs, inline=True, fontsize=8, fmt="%.1fm")

        # Mark drain inlets
        drain_y, drain_x = np.where(cfg.drain_cap > 0)
        ax_heat.scatter(
            (drain_x + 0.5) * cfg.dx,
            (drain_y + 0.5) * cfg.dx,
            color="red",
            marker="x",
            s=40,
            label="Drain Inlets",
        )
        ax_heat.set_xlabel("X (meters)")
        ax_heat.set_ylabel("Y (meters)")
        ax_heat.set_title(f"Peak Inundation Map ({scenario_choice})", fontsize=11)
        ax_heat.legend(loc="upper right", fontsize=8)
        fig_heat.tight_layout()
        st.pyplot(fig_heat)

    with viz_col2:
        st.subheader("📈 Water Balance Hydrograph (Volume vs Time)")
        fig_ts, ax_ts = plt.subplots(figsize=(6, 5))
        t_minutes = res.times / 60.0

        ax_ts.plot(t_minutes, res.volume, label="Stored Surface Volume V(t)", color="#1f78b4", linewidth=2.0)
        ax_ts.plot(t_minutes, res.cum_rain, label="Cumulative Rain", color="#33a02c", linestyle="--")
        ax_ts.plot(t_minutes, res.cum_infil, label="Cumulative Infiltration", color="#ff7f00", linestyle=":")
        ax_ts.plot(t_minutes, res.cum_drain, label="Cumulative Drainage", color="#e31a1c", linestyle="-.")

        ax_ts.set_xlabel("Time (minutes)")
        ax_ts.set_ylabel("Water Volume (m³)")
        ax_ts.set_title("Global Catchment Mass Conservation", fontsize=11)
        ax_ts.grid(True, linestyle="--", alpha=0.5)
        ax_ts.legend(fontsize=8, loc="center right")
        fig_ts.tight_layout()
        st.pyplot(fig_ts)

    # Metrics Table
    st.subheader("📊 Catchment Water Balance Summary")
    summary_data = {
        "Metric": [
            "Total Rain Inflow (m³)",
            "Soil Infiltration (m³)",
            "Storm Drain Intake (m³)",
            "Boundary Discharge (m³)",
            "Final Surface Ponding (m³)",
            "Peak Surface Ponding (m³)",
            "Peak Depth Across Grid (cm)",
            "Continuity Residual Error (m³)",
        ],
        "Value": [
            f"{total_rain_vol:.4f}",
            f"{total_infil_vol:.4f}",
            f"{total_drain_vol:.4f}",
            f"{float(res.cum_boundary[-1]):.4f}",
            f"{final_stored:.4f}",
            f"{peak_volume:.4f}",
            f"{peak_depth * 100:.2f}",
            f"{mass_balance_error:.2e}",
        ],
        "Fraction of Rain (%)": [
            "100.0 %",
            f"{(total_infil_vol / max(total_rain_vol, 1e-6)) * 100:.1f} %",
            f"{(total_drain_vol / max(total_rain_vol, 1e-6)) * 100:.1f} %",
            f"{(float(res.cum_boundary[-1]) / max(total_rain_vol, 1e-6)) * 100:.1f} %",
            f"{(final_stored / max(total_rain_vol, 1e-6)) * 100:.1f} %",
            "-",
            "-",
            f"{(mass_balance_error / max(total_rain_vol, 1e-6)) * 100:.2e} %",
        ],
    }
    st.table(pd.DataFrame(summary_data))


if __name__ == "__main__":
    main()
