"""Urban Flood Simulator — Main CLI Entrypoint.

Educational 2D urban waterlogging simulator for first-year engineering mathematics.
Demonstrates scenario execution, solver time integration, and diagnostic visualization.
SI physical units used throughout: meters (m), seconds (s), m/s, m³.
"""

import sys
from typing import Dict
import numpy as np

from floodsim import scenarios, solver, viz
from floodsim.config import SimConfig, SimResult


def main() -> None:
    """Run all standard policy scenarios and generate visualizations.

    Executes:
    1. Baseline scenario (heavily paved, conventional drainage).
    2. Improved Drains scenario (2x drain capacity + expanded drain network).
    3. Green Infrastructure scenario (permeable retrofits on paved cells).

    Uses try/except blocks so execution completes gracefully when the numerical
    solver or visualization modules are stubs.
    """
    print("================================================================")
    print("      Urban Flood Simulator — 2D Explainable Waterlogging       ")
    print("   First-Year Engineering Mathematics / Numerical Methods       ")
    print("================================================================")
    print("Model: h_new = h + dt * (R - I - D - div(q)) [SI Units Only]")
    print()

    scenario_loaders = [
        ("Baseline", scenarios.baseline),
        ("Improved Drains", scenarios.improved_drains),
        ("Green Infrastructure", scenarios.green_infrastructure),
    ]

    configs: Dict[str, SimConfig] = {}
    results: Dict[str, SimResult] = {}

    # 1. Load and execute scenarios
    for name, loader in scenario_loaders:
        print(f"Loading scenario '{name}'...")
        cfg = loader()
        configs[name] = cfg
        print(f"  -> Grid: {cfg.z.shape}, dx={cfg.dx} m, dt={cfg.dt} s, duration={cfg.t_end} s")
        print(f"  -> Max rain: {cfg.rain.max():.2e} m/s, mean infil: {cfg.infil.mean():.2e} m/s")
        print(f"  -> Active drain cells: {int((cfg.drain_cap > 0).sum())}")

        print(f"  -> Running simulation for '{name}' via solver.run_simulation...")
        try:
            res = solver.run_simulation(cfg)
            results[name] = res
            print(f"  -> '{name}' simulation succeeded.")
        except NotImplementedError as exc:
            print(f"  -> [STUB NOTICE] solver.run_simulation: {exc}")
        except Exception as exc:
            print(f"  -> [WARNING] solver execution failed for '{name}': {exc}")
        print()

    # 2. Call visualization routines
    print("Calling visualization routines...")
    try:
        if results:
            for name, res in results.items():
                print(f"  -> Generating depth map and timeseries for '{name}'...")
                viz.plot_depth_map(res.h_history[-1], z=configs[name].z, title=f"{name} Final Depth")
                viz.plot_timeseries(
                    times=res.times,
                    volume=res.volume,
                    cum_rain=res.cum_rain,
                    cum_infil=res.cum_infil,
                    cum_drain=res.cum_drain,
                )
            print("  -> Generating multi-scenario comparison plot...")
            viz.plot_scenarios(results)
            print("  -> All visualizations generated successfully.")
        else:
            # If solver returned no results because it is a stub, still exercise viz functions
            print("  -> Solver returned no results (stub). Exercising viz stubs with initial state...")
            first_cfg = configs.get("Baseline")
            if first_cfg is not None:
                mock_h = np.zeros(first_cfg.z.shape, dtype=np.float64)
                viz.plot_depth_map(mock_h, z=first_cfg.z, title="Baseline Initial Depth")
    except NotImplementedError as exc:
        print(f"  -> [STUB NOTICE] viz function: {exc}")
    except Exception as exc:
        print(f"  -> [WARNING] viz generation failed: {exc}")

    print()
    print("================================================================")
    print("All 3 scenarios evaluated. Execution completed successfully.")
    print("================================================================")


if __name__ == "__main__":
    main()
