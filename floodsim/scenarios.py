"""Pre-configured urban flood scenarios for policy and intervention comparisons.

Scenarios define benchmark setups to evaluate:
1. Baseline: High paved fraction (impervious urban area) with conventional drainage.
2. Improved Drains: Baseline land cover with expanded storm drain capacity / inlets.
3. Green Infrastructure: Reduced paved fraction (permeable pavements, rain gardens).
"""

from floodsim.config import SimConfig


def baseline() -> SimConfig:
    """Construct the baseline urban waterlogging scenario.

    Mathematical Context:
        Represents a heavily urbanized area:
        - High paved fraction (~80% impervious, low infiltration I ~ 0 m/s).
        - Conventional drainage capacity D_max (~10-25 mm/hr capacity).
        - Standard design storm event (e.g., 50 mm/hr for 30 minutes).
        - Bowl-shaped or low-lying depression topography prone to ponding.

    Returns:
        SimConfig: Complete simulation configuration for baseline run.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("baseline scenario is not yet implemented.")


def improved_drains() -> SimConfig:
    """Construct scenario with upgraded grey storm drain infrastructure.

    Mathematical Context:
        Maintains baseline topography and land cover imperviousness,
        but increases local drain intake capacity D_max (e.g., 2x to 4x baseline capacity)
        and/or places additional drain inlets at local topographic minima.

    Returns:
        SimConfig: Complete simulation configuration with upgraded drainage.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("improved_drains scenario is not yet implemented.")


def green_infrastructure() -> SimConfig:
    """Construct scenario with nature-based green infrastructure (Sponge City measures).

    Mathematical Context:
        Replaces impervious surfaces with permeable pavements, bio-retention cells,
        and urban parks, increasing effective infiltration capacity I_max
        (e.g., paved fraction reduced from 80% to 30%, increasing aggregate soil absorption).

    Returns:
        SimConfig: Complete simulation configuration with green infrastructure.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("green_infrastructure scenario is not yet implemented.")
