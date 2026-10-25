"""Urban Flood Simulator — Main CLI Entrypoint.

Educational 2D urban waterlogging simulator for first-year engineering mathematics.
Demonstrates scenario execution and contract status.
"""

import sys
import numpy as np
from floodsim.config import SimConfig


def main() -> None:
    """Print project status and verify basic environment setup."""
    print("================================================================")
    print("      Urban Flood Simulator — 2D Explainable Waterlogging       ")
    print("   First-Year Engineering Mathematics / Numerical Methods       ")
    print("================================================================")
    print("Status: Architectural contract initialized (stubs only).")
    print("Model:  h_new = h + dt * (R - I - D - div(q)) [SI Units Only]")
    print()
    print("Verifying SimConfig dataclass creation...")
    demo_cfg = SimConfig(
        z=np.zeros((10, 10), dtype=np.float64),
        dx=2.0,
        dt=1.0,
        t_end=60.0,
        rain=np.zeros(60, dtype=np.float64),
        infil=np.zeros((10, 10), dtype=np.float64),
        drain_cap=np.zeros((10, 10), dtype=np.float64),
    )
    print(f"-> Successfully created SimConfig for {demo_cfg.z.shape} grid.")
    print("-> Cell resolution: dx =", demo_cfg.dx, "m, dt =", demo_cfg.dt, "s")
    print()
    print("To run contract tests:")
    print("    pytest tests/")
    print()
    print("To launch Streamlit dashboard:")
    print("    streamlit run app/dashboard.py")
    print("================================================================")


if __name__ == "__main__":
    main()
