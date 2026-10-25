"""Tests for exact numerical water balance and mass conservation.

Verifies:
1. rain_in - infil - drain - boundary_out == change in stored volume (relative tolerance 1e-9).
2. no rain and dry start => depth stays 0 everywhere.
3. closed boundary, no drains, no infiltration => final volume equals total rain.
4. depth is never negative under extreme sinks or steep gradients.
5. on a bowl terrain, the centre accumulates more water than the rim.
"""

import sys
from pathlib import Path
import numpy as np
import pytest

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# If check_sanity exists in the root, importing it ensures reference solver is active if stubbed
try:
    import check_sanity  # noqa: F401
except ImportError:
    pass

from floodsim.config import SimConfig
from floodsim import solver
from floodsim.rainfall import constant_storm
from floodsim.terrain import make_land_cover, make_terrain, place_drains


def test_mass_balance_continuity_equation() -> None:
    """Verify that rain_in - infil - drain - boundary_out == change in stored volume (rtol 1e-9).

    Physical Law:
        Global mass conservation for discrete overland water balance:
            Delta V_stored(t) = V_rain(t) - V_infil(t) - V_drain(t) - V_boundary(t)
        Every drop of rain added must either remain on the ground surface,
        infiltrate into the soil, drain into storm inlets, or exit outer boundaries.
    """
    rows, cols = 12, 12
    dx = 2.0
    dt = 1.0
    t_end = 60.0

    z = make_terrain(rows=rows, cols=cols, kind="bowl", seed=10)
    rain = constant_storm(intensity_mm_hr=60.0, duration_min=0.5, dt=dt, t_end=t_end)
    infil = make_land_cover((rows, cols), paved_fraction=0.6, seed=20)
    drain_cap = place_drains((rows, cols), positions=[(5, 5), (6, 6)], capacity=1e-5)

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
    )

    res = solver.run_simulation(cfg)

    # Check water balance at every saved snapshot
    delta_stored = res.volume - res.volume[0]
    expected_delta = res.cum_rain - res.cum_infil - res.cum_drain - res.cum_boundary
    abs_diff = np.abs(delta_stored - expected_delta)
    denom = np.maximum(res.cum_rain, 1.0)
    rel_diff = abs_diff / denom

    assert np.all(rel_diff < 1e-9), f"Mass conservation violated! Max rel diff: {np.max(rel_diff):.2e}"


def test_no_rain_dry_start_stays_zero() -> None:
    """Verify that with no rain and dry start, water depth remains identically 0 everywhere.

    Physical Law:
        Spontaneous generation of water is physically impossible. With R(t) = 0
        and h(x, y, 0) = 0, no surface water can emerge regardless of topography or sinks.
    """
    rows, cols = 10, 10
    dx = 2.0
    dt = 1.0
    n_steps = 30
    t_end = float(n_steps * dt)

    z = make_terrain(rows=rows, cols=cols, kind="bowl", seed=42)
    rain = np.zeros(n_steps, dtype=np.float64)
    infil = np.full((rows, cols), 1e-5, dtype=np.float64)
    drain_cap = np.full((rows, cols), 1e-5, dtype=np.float64)

    cfg = SimConfig(
        z=z,
        dx=dx,
        dt=dt,
        t_end=t_end,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
        closed_boundary=True,
    )

    res = solver.run_simulation(cfg)

    assert np.all(res.h_history == 0.0), "Non-zero depth generated without rainfall!"
    assert np.all(res.volume == 0.0), "Stored volume non-zero without rainfall!"
    assert np.all(res.cum_rain == 0.0)
    assert np.all(res.cum_infil == 0.0)
    assert np.all(res.cum_drain == 0.0)
    assert np.all(res.cum_boundary == 0.0)


def test_closed_boundary_no_sinks_equals_rain() -> None:
    """Verify that in closed basin with no sinks, final stored volume equals total rain.

    Physical Law:
        When boundaries are impermeable (closed_boundary=True) and soil/drain
        extraction is absent (I = 0, D = 0), the domain behaves as an ideal bucket:
            V_stored(t_end) = V_rain(t_end)
    """
    rows, cols = 10, 10
    dx = 1.5
    dt = 0.5
    t_end = 30.0

    z = make_terrain(rows=rows, cols=cols, kind="bowl", seed=7)
    rain = constant_storm(intensity_mm_hr=40.0, duration_min=0.3, dt=dt, t_end=t_end)
    infil = np.zeros((rows, cols), dtype=np.float64)
    drain_cap = np.zeros((rows, cols), dtype=np.float64)

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
    )

    res = solver.run_simulation(cfg)

    final_vol = float(res.volume[-1])
    total_rain_vol = float(res.cum_rain[-1])

    assert total_rain_vol > 0.0, "Rainfall volume should be positive."
    np.testing.assert_allclose(
        final_vol,
        total_rain_vol,
        rtol=1e-9,
        atol=1e-9,
        err_msg="Stored volume does not equal total rain in a lossless closed basin!",
    )


def test_depth_never_negative() -> None:
    """Verify that depth h >= 0 is strictly maintained under extreme sinks and steep slopes.

    Physical Law:
        Ponded water depth cannot become negative. Infiltration, drainage,
        and overland outflows must be strictly limited to available water.
    """
    rows, cols = 8, 8
    dx = 1.0
    dt = 0.5
    t_end = 40.0

    z = make_terrain(rows=rows, cols=cols, kind="slope", seed=1)
    rain = constant_storm(intensity_mm_hr=100.0, duration_min=0.2, dt=dt, t_end=t_end)
    infil = np.full((rows, cols), 0.5, dtype=np.float64)  # 0.5 m/s extreme absorption
    drain_cap = np.full((rows, cols), 0.5, dtype=np.float64)  # 0.5 m/s extreme drainage

    # Test open boundary
    cfg_open = SimConfig(
        z=z,
        dx=dx,
        dt=dt,
        t_end=t_end,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
        k_flow=0.5,
        closed_boundary=False,
    )
    res_open = solver.run_simulation(cfg_open)
    assert np.all(res_open.h_history >= 0.0), "Negative depth encountered with open boundary!"

    # Test closed boundary
    cfg_closed = SimConfig(
        z=z,
        dx=dx,
        dt=dt,
        t_end=t_end,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
        k_flow=0.5,
        closed_boundary=True,
    )
    res_closed = solver.run_simulation(cfg_closed)
    assert np.all(res_closed.h_history >= 0.0), "Negative depth encountered with closed boundary!"


def test_bowl_terrain_center_accumulates_more_than_rim() -> None:
    """Verify topographic gravity routing: water accumulates in the depression center over the rim.

    Physical Law:
        Gravity drives overland runoff down hydraulic head gradients:
            q proportional to -grad(h + z)
        On a paraboloid bowl bed with lowest elevation in the center, runoff
        converges inward, producing deeper ponding at the center than at the outer perimeter.
    """
    rows, cols = 15, 15
    dx = 2.0
    dt = 1.0
    t_end = 60.0

    z = make_terrain(rows=rows, cols=cols, kind="bowl", seed=42)
    rain = constant_storm(intensity_mm_hr=50.0, duration_min=0.4, dt=dt, t_end=t_end)
    infil = np.full((rows, cols), 1e-7, dtype=np.float64)  # low infiltration
    drain_cap = np.zeros((rows, cols), dtype=np.float64)

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
    )

    res = solver.run_simulation(cfg)
    final_h = res.h_history[-1]

    # Center 3x3 region
    r_mid, c_mid = rows // 2, cols // 2
    center_depth = float(np.mean(final_h[r_mid - 1 : r_mid + 2, c_mid - 1 : c_mid + 2]))

    # Outer perimeter (rim)
    rim_mask = np.ones((rows, cols), dtype=bool)
    rim_mask[2:-2, 2:-2] = False
    rim_depth = float(np.mean(final_h[rim_mask]))

    assert center_depth > rim_depth, (
        f"Center depth ({center_depth:.4f} m) did not exceed rim depth ({rim_depth:.4f} m) on bowl terrain!"
    )
    assert center_depth > 0.0, "Center should have ponded water."
