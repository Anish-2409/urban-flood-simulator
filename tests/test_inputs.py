"""Tests for terrain generation, land cover, drains, rainfall, and scenarios.

Validates that:
- Array shapes match specifications (30x30 default grid, hyetograph time steps).
- Hydraulic rates (rain, infiltration, drain capacity) are strictly non-negative.
- Rainfall conversion factor from mm/hr to m/s is correct (1 mm/hr = 2.78e-7 m/s).
- Scenarios return valid SimConfig instances with identical shared terrain and rainfall.
"""

import sys
from pathlib import Path
import numpy as np
import pytest

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from floodsim.config import SimConfig
from floodsim.rainfall import MM_HR_TO_M_S, constant_storm, triangular_storm
from floodsim.scenarios import baseline, green_infrastructure, improved_drains
from floodsim.terrain import make_land_cover, make_terrain, place_drains


# ---------------------------------------------------------------------------
# Terrain Tests
# ---------------------------------------------------------------------------

def test_make_terrain_default_shape_and_dtype() -> None:
    """Verify default terrain generates a 30x30 float64 grid with finite values."""
    z = make_terrain()
    assert z.shape == (30, 30)
    assert z.dtype == np.float64
    assert np.all(np.isfinite(z))
    assert np.all(z >= 0.0)


def test_make_terrain_bowl_low_center() -> None:
    """Verify bowl terrain has its minimum in the center and slopes upwards to borders."""
    rows, cols = 30, 30
    z = make_terrain(rows=rows, cols=cols, kind="bowl")

    assert z.shape == (rows, cols)
    assert z.dtype == np.float64

    # Center cells should have the lowest elevation in the domain
    center_val = z[rows // 2, cols // 2]
    assert np.isclose(center_val, z.min(), atol=0.01)

    # Perimeter cells must be significantly higher than the center
    corner_vals = [z[0, 0], z[0, cols - 1], z[rows - 1, 0], z[rows - 1, cols - 1]]
    for corner in corner_vals:
        assert corner > center_val + 0.5


def test_make_terrain_slope_tilted_plane() -> None:
    """Verify slope terrain forms a planar tilt decreasing from upstream to downstream."""
    rows, cols = 30, 30
    z = make_terrain(rows=rows, cols=cols, kind="slope")

    assert z.shape == (rows, cols)
    # Upstream (row 0) should be higher than downstream (row rows - 1)
    assert z[0, 0] > z[rows - 1, 0]
    # Rows should decrease monotonically
    diffs = np.diff(z[:, 0])
    assert np.all(diffs <= 0.0)
    # Along columns within each row, elevation should be constant
    for r in range(rows):
        assert np.all(z[r, :] == z[r, 0])


def test_make_terrain_random_smoothed_and_seeded() -> None:
    """Verify random terrain is smoothed, non-negative, and reproducible with seed."""
    z1 = make_terrain(30, 30, kind="random", seed=42)
    z2 = make_terrain(30, 30, kind="random", seed=42)
    z3 = make_terrain(30, 30, kind="random", seed=99)

    assert z1.shape == (30, 30)
    assert np.all(z1 >= 0.0)
    # Identical seed must produce identical terrain
    assert np.array_equal(z1, z2)
    # Different seed should produce different terrain
    assert not np.array_equal(z1, z3)


def test_make_terrain_invalid_inputs() -> None:
    """Verify make_terrain raises ValueError for invalid shapes or kinds."""
    with pytest.raises(ValueError):
        make_terrain(0, 30)
    with pytest.raises(ValueError):
        make_terrain(30, -5)
    with pytest.raises(ValueError):
        make_terrain(30, 30, kind="unsupported_kind")  # type: ignore


# ---------------------------------------------------------------------------
# Land Cover & Infiltration Tests
# ---------------------------------------------------------------------------

def test_make_land_cover_shape_and_rates() -> None:
    """Verify land cover generates matching shape, non-negative rates, and expected ranges."""
    shape = (30, 30)
    infil = make_land_cover(shape, paved_fraction=0.8, seed=42)

    assert infil.shape == shape
    assert infil.dtype == np.float64
    # All infiltration rates must be non-negative
    assert np.all(infil >= 0.0)

    # 100% paved: every cell should be in paved range ~0 to 1e-7 m/s
    infil_paved = make_land_cover(shape, paved_fraction=1.0, seed=42)
    assert np.all(infil_paved >= 0.0)
    assert np.all(infil_paved <= 1e-7)

    # 0% paved (100% green/soil): every cell should be in green range ~1e-6 to 5e-6 m/s
    infil_green = make_land_cover(shape, paved_fraction=0.0, seed=42)
    assert np.all(infil_green >= 1e-6)
    assert np.all(infil_green <= 5e-6)


def test_make_land_cover_paved_fraction_scaling() -> None:
    """Verify higher paved fraction yields lower overall average infiltration capacity."""
    shape = (40, 40)
    low_paved = make_land_cover(shape, paved_fraction=0.2, seed=123)
    high_paved = make_land_cover(shape, paved_fraction=0.8, seed=123)

    assert low_paved.mean() > high_paved.mean()


def test_make_land_cover_invalid_inputs() -> None:
    """Verify make_land_cover validates shape and paved_fraction boundaries."""
    with pytest.raises(ValueError):
        make_land_cover((30, 30), paved_fraction=-0.1)
    with pytest.raises(ValueError):
        make_land_cover((30, 30), paved_fraction=1.5)
    with pytest.raises(ValueError):
        make_land_cover((-1, 30), paved_fraction=0.5)


# ---------------------------------------------------------------------------
# Drain Placement Tests
# ---------------------------------------------------------------------------

def test_place_drains_shape_capacity_and_non_negative() -> None:
    """Verify place_drains sets exact capacity at inlets and zero elsewhere."""
    shape = (30, 30)
    positions = [(5, 5), (10, 10), (14, 15)]
    capacity = 1.5e-4  # m/s

    drain_cap = place_drains(shape, positions, capacity=capacity)

    assert drain_cap.shape == shape
    assert drain_cap.dtype == np.float64
    # All rates non-negative
    assert np.all(drain_cap >= 0.0)

    # Inlets have exact specified capacity
    for r, c in positions:
        assert drain_cap[r, c] == capacity

    # Non-drain cells have zero capacity
    assert np.sum(drain_cap > 0.0) == len(positions)


def test_place_drains_invalid_inputs() -> None:
    """Verify place_drains checks bounds and negative capacity."""
    with pytest.raises(ValueError):
        place_drains((30, 30), positions=[(5, 5)], capacity=-1.0)
    with pytest.raises(IndexError):
        place_drains((30, 30), positions=[(35, 10)], capacity=1e-4)


# ---------------------------------------------------------------------------
# Rainfall Hyetograph & Conversion Tests
# ---------------------------------------------------------------------------

def test_rainfall_mm_hr_to_m_s_conversion_constant_storm() -> None:
    """Verify constant storm converts mm/hr correctly: 1 mm/hr = 2.78e-7 m/s."""
    dt = 1.0  # s
    t_end = 600.0  # s (10 minutes)
    duration_min = 5.0  # min (300 s)

    # 1 mm/hr test
    rain_1 = constant_storm(intensity_mm_hr=1.0, duration_min=duration_min, dt=dt, t_end=t_end)
    n_steps = int(np.ceil(t_end / dt))

    # Shape match: one value per time step
    assert rain_1.shape == (n_steps,)
    assert rain_1.dtype == np.float64
    # Rates non-negative
    assert np.all(rain_1 >= 0.0)

    # 1 mm/hr conversion check: during storm, rate equals 2.78e-7 m/s
    assert np.isclose(rain_1[0], 2.78e-7, atol=1e-12)
    assert np.isclose(rain_1[299], 2.78e-7, atol=1e-12)
    # After storm ends (t >= 300 s), rate is 0.0
    assert rain_1[300] == 0.0
    assert rain_1[599] == 0.0

    # Multi-value intensity test (50 mm/hr)
    rain_50 = constant_storm(intensity_mm_hr=50.0, duration_min=duration_min, dt=dt, t_end=t_end)
    expected_50 = 50.0 * 2.78e-7
    assert np.isclose(rain_50[0], expected_50, atol=1e-12)


def test_rainfall_triangular_storm_peak_and_shape() -> None:
    """Verify triangular storm peak conversion, shape matching, and temporal bounds."""
    dt = 1.0  # s
    t_end = 1200.0  # s (20 min)
    duration_min = 10.0  # min (600 s)
    peak_mm_hr = 60.0  # mm/hr

    rain = triangular_storm(peak_mm_hr=peak_mm_hr, duration_min=duration_min, dt=dt, t_end=t_end)
    n_steps = int(np.ceil(t_end / dt))

    # Shape match
    assert rain.shape == (n_steps,)
    assert rain.dtype == np.float64
    # Non-negative rates
    assert np.all(rain >= 0.0)

    # Conversion check at peak (t_half = 300 s -> index 300)
    expected_peak = 60.0 * 2.78e-7
    assert np.isclose(rain[300], expected_peak, atol=1e-12)

    # Starts at 0.0
    assert rain[0] == 0.0

    # Zero after storm duration (t >= 600 s)
    assert rain[600] == 0.0
    assert rain[1000] == 0.0


def test_rainfall_invalid_inputs() -> None:
    """Verify rainfall functions reject invalid negative or non-positive parameters."""
    with pytest.raises(ValueError):
        constant_storm(50.0, 30.0, dt=-1.0, t_end=100.0)
    with pytest.raises(ValueError):
        constant_storm(50.0, 30.0, dt=1.0, t_end=0.0)
    with pytest.raises(ValueError):
        constant_storm(-5.0, 30.0, dt=1.0, t_end=100.0)
    with pytest.raises(ValueError):
        triangular_storm(-10.0, 30.0, dt=1.0, t_end=100.0)


# ---------------------------------------------------------------------------
# Scenarios Tests
# ---------------------------------------------------------------------------

def test_scenarios_instantiation_and_shapes() -> None:
    """Verify baseline, improved_drains, and green_infrastructure return valid SimConfig with matching shapes."""
    cfg_base = baseline()
    cfg_imp = improved_drains()
    cfg_green = green_infrastructure()

    for name, cfg in [("baseline", cfg_base), ("improved_drains", cfg_imp), ("green_infrastructure", cfg_green)]:
        assert isinstance(cfg, SimConfig), f"{name} must return SimConfig"
        # 30x30 default grid shapes
        assert cfg.z.shape == (30, 30)
        assert cfg.infil.shape == (30, 30)
        assert cfg.drain_cap.shape == (30, 30)

        # Rainfall array length matches int(ceil(t_end / dt))
        expected_steps = int(np.ceil(cfg.t_end / cfg.dt))
        assert cfg.rain.shape == (expected_steps,)

        # All rates strictly non-negative
        assert np.all(cfg.rain >= 0.0), f"{name} rain contains negative rates"
        assert np.all(cfg.infil >= 0.0), f"{name} infil contains negative rates"
        assert np.all(cfg.drain_cap >= 0.0), f"{name} drain_cap contains negative rates"


def test_scenarios_fair_comparison_identical_terrain_and_rainfall() -> None:
    """Verify all scenarios use strictly identical topography and rainfall forcing for fair comparison."""
    cfg_base = baseline()
    cfg_imp = improved_drains()
    cfg_green = green_infrastructure()

    # Identical terrain bed elevation
    assert np.array_equal(cfg_base.z, cfg_imp.z)
    assert np.array_equal(cfg_base.z, cfg_green.z)

    # Identical rainfall forcing hyetograph
    assert np.array_equal(cfg_base.rain, cfg_imp.rain)
    assert np.array_equal(cfg_base.rain, cfg_green.rain)

    # Identical spatial and temporal discretization
    assert cfg_base.dx == cfg_imp.dx == cfg_green.dx
    assert cfg_base.dt == cfg_imp.dt == cfg_green.dt
    assert cfg_base.t_end == cfg_imp.t_end == cfg_green.t_end


def test_scenarios_improved_drains_intervention() -> None:
    """Verify improved_drains has 2x capacity and more drain inlets than baseline."""
    cfg_base = baseline()
    cfg_imp = improved_drains()

    # More active drain cells
    base_drain_count = int(np.sum(cfg_base.drain_cap > 0.0))
    imp_drain_count = int(np.sum(cfg_imp.drain_cap > 0.0))
    assert imp_drain_count > base_drain_count

    # 2x capacity at drain cells
    assert np.isclose(cfg_imp.drain_cap.max(), 2.0 * cfg_base.drain_cap.max())


def test_scenarios_green_infrastructure_intervention() -> None:
    """Verify green_infrastructure increases infiltration in previously paved cells."""
    cfg_base = baseline()
    cfg_green = green_infrastructure()

    # Overall infiltration capacity must be higher
    assert cfg_green.infil.mean() > cfg_base.infil.mean()

    # Increased infiltration occurred in previously paved cells (infil <= 1e-7 in baseline)
    upgraded_cells = cfg_green.infil > cfg_base.infil
    assert np.any(upgraded_cells)
    assert np.all(cfg_base.infil[upgraded_cells] <= 1e-7)
    # The upgraded cells now have green infiltration rates ~1e-6 to 5e-6 m/s
    assert np.all(cfg_green.infil[upgraded_cells] >= 1e-6)
    assert np.all(cfg_green.infil[upgraded_cells] <= 5e-6)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
