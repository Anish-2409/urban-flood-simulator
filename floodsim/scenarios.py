"""Pre-configured urban flood scenarios for policy and intervention comparisons.

Scenarios define benchmark setups to evaluate:
1. Baseline: High paved fraction (impervious urban area) with conventional drainage.
2. Improved Drains: Baseline land cover with expanded storm drain capacity (2x) and additional inlets.
3. Green Infrastructure: Reduced imperviousness (higher infiltration in previously paved cells).

All scenarios share identical bed topography and rainfall forcing to ensure fair comparisons.
Units: SI units exclusively (m, s, m/s, m³).
"""

from typing import Tuple
import numpy as np

from floodsim.config import SimConfig
from floodsim.rainfall import constant_storm
from floodsim.terrain import make_land_cover, make_terrain, place_drains

# Shared simulation domain geometry and meteorological parameters
_GRID_SHAPE: Tuple[int, int] = (30, 30)
_DX: float = 2.0  # meters
_DT: float = 1.0  # seconds
_T_END: float = 3600.0  # seconds (1 hour)
_STORM_INTENSITY_MM_HR: float = 50.0  # mm/hr
_STORM_DURATION_MIN: float = 30.0  # minutes

# Base drainage parameters
_BASE_DRAIN_POSITIONS = [(14, 14), (14, 15), (15, 14), (15, 15)]
_BASE_DRAIN_CAPACITY: float = 25.0 * 2.78e-7  # 25 mm/hr in m/s (~6.95e-6 m/s)


def _get_shared_terrain() -> np.ndarray:
    """Return shared 30x30 bowl topography for fair scenario comparison."""
    return make_terrain(rows=_GRID_SHAPE[0], cols=_GRID_SHAPE[1], kind="bowl", seed=42)


def _get_shared_rain() -> np.ndarray:
    """Return shared rainfall hyetograph for fair scenario comparison."""
    return constant_storm(
        intensity_mm_hr=_STORM_INTENSITY_MM_HR,
        duration_min=_STORM_DURATION_MIN,
        dt=_DT,
        t_end=_T_END,
    )


def baseline() -> SimConfig:
    """Construct the baseline urban waterlogging scenario.

    Mathematical Context:
        Represents a typical dense urban area before intervention:
        - Bowl depression topography prone to central water accumulation.
        - High paved fraction (~80% impervious with low infiltration I ~ 0 to 1e-7 m/s).
        - Conventional storm drainage capacity D_max ~ 25 mm/hr (6.95e-6 m/s) at central inlets.
        - Design storm event: 50 mm/hr for 30 minutes over a 1-hour simulation.

    Returns:
        SimConfig: Complete simulation configuration for baseline run.
    """
    z = _get_shared_terrain()
    rain = _get_shared_rain()
    infil = make_land_cover(_GRID_SHAPE, paved_fraction=0.8, seed=42)
    drain_cap = place_drains(_GRID_SHAPE, positions=_BASE_DRAIN_POSITIONS, capacity=_BASE_DRAIN_CAPACITY)

    return SimConfig(
        z=z,
        dx=_DX,
        dt=_DT,
        t_end=_T_END,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
    )


def improved_drains() -> SimConfig:
    """Construct scenario with upgraded grey storm drain infrastructure.

    Mathematical Context:
        Maintains identical bed topography and baseline land cover infiltration,
        but upgrades grey drainage infrastructure:
        - 2x drainage intake capacity (D_max_improved = 2 * D_max_base = 50 mm/hr in m/s).
        - Expanded drain network with additional inlets installed at adjacent low-lying cells.

    Returns:
        SimConfig: Complete simulation configuration with upgraded drainage.
    """
    z = _get_shared_terrain()
    rain = _get_shared_rain()
    infil = make_land_cover(_GRID_SHAPE, paved_fraction=0.8, seed=42)

    # 2x capacity + more drains surrounding central depression
    improved_positions = [
        (14, 14), (14, 15), (15, 14), (15, 15),
        (13, 14), (13, 15), (16, 14), (16, 15),
    ]
    improved_capacity = 2.0 * _BASE_DRAIN_CAPACITY
    drain_cap = place_drains(_GRID_SHAPE, positions=improved_positions, capacity=improved_capacity)

    return SimConfig(
        z=z,
        dx=_DX,
        dt=_DT,
        t_end=_T_END,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
    )


def green_infrastructure() -> SimConfig:
    """Construct scenario with nature-based green infrastructure (Sponge City measures).

    Mathematical Context:
        Maintains identical bed topography, rainfall forcing, and baseline drainage,
        but retrofits impermeable urban surfaces with permeable pavements, bio-retention cells,
        and rain gardens. This increases effective infiltration capacity in ~60% of previously
        paved cells to green soil rates (~1e-6 to 5e-6 m/s).

    Returns:
        SimConfig: Complete simulation configuration with green infrastructure.
    """
    z = _get_shared_terrain()
    rain = _get_shared_rain()
    base_infil = make_land_cover(_GRID_SHAPE, paved_fraction=0.8, seed=42)

    # Increase infiltration in a subset of previously paved cells
    rng = np.random.default_rng(101)
    paved_mask = base_infil <= 1e-7
    convert_mask = paved_mask & (rng.random(_GRID_SHAPE) < 0.6)

    infil = base_infil.copy()
    infil[convert_mask] = rng.uniform(1e-6, 5e-6, size=int(convert_mask.sum()))

    drain_cap = place_drains(_GRID_SHAPE, positions=_BASE_DRAIN_POSITIONS, capacity=_BASE_DRAIN_CAPACITY)

    return SimConfig(
        z=z,
        dx=_DX,
        dt=_DT,
        t_end=_T_END,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
    )
