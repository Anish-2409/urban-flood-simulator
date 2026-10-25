"""Sanity test suite verifying physical behavior across policy scenarios.

Based on check_sanity.py in the project root.
Validates:
1. Water balance holds for baseline, improved_drains, and green_infrastructure (tol 1e-6).
2. Minimum water depth is strictly non-negative (h >= 0) for every scenario.
3. Under zero rainfall, all depths remain identically zero without altering the original configuration.
4. Peak surface ponding volume for improved_drains is less than or equal to baseline.

All tests include docstrings explaining the physical rule checked, and use shortened
simulation durations to execute within a few seconds.
"""

import copy
import sys
from pathlib import Path
from typing import Callable
import numpy as np
import pytest

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Base tests on check_sanity if available in the project root
try:
    import check_sanity  # noqa: F401
except ImportError:
    pass

from floodsim.config import SimConfig, SimResult
from floodsim import scenarios, solver

# Short duration for fast regression testing (seconds)
_FAST_DURATION_S: float = 120.0


def _load_fast_scenario(loader: Callable[[], SimConfig], duration: float = _FAST_DURATION_S) -> SimConfig:
    """Helper to create a deep copy of a scenario with shortened duration for fast testing."""
    cfg = copy.deepcopy(loader())
    n_steps = max(1, int(np.ceil(duration / cfg.dt)))
    cfg.rain = cfg.rain[:n_steps].copy()
    cfg.t_end = float(n_steps * cfg.dt)
    return cfg


@pytest.mark.parametrize("scenario_name", ["baseline", "improved_drains", "green_infrastructure"])
def test_water_balance_holds(scenario_name: str) -> None:
    """Verify global water volume balance conservation (tolerance 1e-6 m³).

    Physical Rule:
        Continuity / conservation of mass:
            V_rain(t) - V_infil(t) - V_drain(t) - V_boundary(t) - Delta V_stored(t) = 0
        In a closed urban catchment, all precipitation entering the system must be
        strictly accounted for by surface storage, soil infiltration, storm sewer
        intake, or boundary discharge.
    """
    loader = getattr(scenarios, scenario_name)
    cfg = _load_fast_scenario(loader)
    res: SimResult = solver.run_simulation(cfg)

    rain_vol = float(res.cum_rain[-1])
    infil_vol = float(res.cum_infil[-1])
    drain_vol = float(res.cum_drain[-1])
    boundary_vol = float(res.cum_boundary[-1])
    water_left = float(res.volume[-1]) - float(res.volume[0])

    diff = rain_vol - infil_vol - drain_vol - boundary_vol - water_left

    assert abs(diff) < 1e-6, (
        f"Water balance violated in '{scenario_name}'! "
        f"Residual difference: {diff:.6e} m³ (expected < 1e-6 m³)."
    )


@pytest.mark.parametrize("scenario_name", ["baseline", "improved_drains", "green_infrastructure"])
def test_minimum_depth_non_negative(scenario_name: str) -> None:
    """Verify that water depth h is strictly non-negative across all grid cells and times.

    Physical Rule:
        Physical state admissibility: Negative water depth is mathematically and physically
        impossible. Sinks (soil absorption, storm inlets) and convective outflows cannot
        extract more water than currently occupies a cell.
    """
    loader = getattr(scenarios, scenario_name)
    cfg = _load_fast_scenario(loader)
    res: SimResult = solver.run_simulation(cfg)

    min_depth = float(np.min(res.h_history))
    assert min_depth >= 0.0, (
        f"Negative depth detected in '{scenario_name}'! Min depth: {min_depth:.6e} m."
    )


@pytest.mark.parametrize("scenario_name", ["baseline", "improved_drains", "green_infrastructure"])
def test_zero_rainfall_depths_stay_zero(scenario_name: str) -> None:
    """Verify that with zero rainfall and dry initial state, water depth remains 0 everywhere.

    Physical Rule:
        Quiescent equilibrium: Without external meteorological forcing (rainfall rate R = 0)
        and starting from a dry bed (h0 = 0), the overland continuity equation admits only
        the trivial stationary solution h(x, y, t) = 0 for all t >= 0.
    """
    loader = getattr(scenarios, scenario_name)
    cfg = _load_fast_scenario(loader)

    # Copy config array without modifying original
    cfg.rain = np.zeros_like(cfg.rain)
    cfg.h0 = None

    res: SimResult = solver.run_simulation(cfg)

    assert np.all(res.h_history == 0.0), (
        f"Spontaneous water generated in zero-rain run of '{scenario_name}'! "
        f"Max depth: {np.max(res.h_history):.6e} m."
    )
    assert np.all(res.volume == 0.0)


def test_peak_stored_volume_improved_drains_vs_baseline() -> None:
    """Verify that peak surface ponding volume for improved_drains is <= baseline's peak volume.

    Physical Rule:
        Grey storm drainage intervention effectiveness: Doubling storm drain intake
        capacity and expanding inlet placement at local topographic minima accelerates
        water removal, thereby reducing or capping peak inundation volume compared to
        conventional drainage under identical rainfall forcing.
    """
    cfg_base = _load_fast_scenario(scenarios.baseline, duration=180.0)
    cfg_improved = _load_fast_scenario(scenarios.improved_drains, duration=180.0)

    res_base = solver.run_simulation(cfg_base)
    res_improved = solver.run_simulation(cfg_improved)

    peak_base = float(np.max(res_base.volume))
    peak_improved = float(np.max(res_improved.volume))

    assert peak_improved <= peak_base, (
        f"Improved drains peak volume ({peak_improved:.4f} m³) exceeded baseline ({peak_base:.4f} m³)!"
    )
