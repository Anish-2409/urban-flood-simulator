"""Test suite validating that SimConfig and SimResult dataclasses can be created.

Verifies instantiation, default values, and array attributes according to the architectural contract.
Compatible with pytest and standard python execution.
"""

import sys
from pathlib import Path
import numpy as np

# Ensure repository root is on sys.path for direct script execution and pytest
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from floodsim.config import SimConfig, SimResult


def test_sim_config_can_be_created() -> None:
    """Verify that SimConfig can be instantiated with valid mock arrays and parameters."""
    rows, cols = 10, 15
    z = np.zeros((rows, cols), dtype=np.float64)
    dx = 2.0  # meters
    dt = 1.0  # seconds
    t_end = 60.0  # seconds
    n_steps = int(np.ceil(t_end / dt))
    rain = np.full(n_steps, 1e-5, dtype=np.float64)  # m/s
    infil = np.full((rows, cols), 5e-6, dtype=np.float64)  # m/s
    drain_cap = np.zeros((rows, cols), dtype=np.float64)  # m/s

    cfg = SimConfig(
        z=z,
        dx=dx,
        dt=dt,
        t_end=t_end,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
    )

    assert cfg.z.shape == (rows, cols)
    assert cfg.dx == 2.0
    assert cfg.dt == 1.0
    assert cfg.t_end == 60.0
    assert cfg.rain.shape == (n_steps,)
    assert cfg.infil.shape == (rows, cols)
    assert cfg.drain_cap.shape == (rows, cols)

    # Check defaults
    assert cfg.k_flow == 0.25
    assert cfg.closed_boundary is True
    assert cfg.h0 is None
    assert cfg.save_every == 1


def test_sim_config_with_custom_options() -> None:
    """Verify SimConfig accepts optional initial depth and custom flags."""
    rows, cols = 5, 5
    z = np.ones((rows, cols), dtype=np.float64)
    h0 = np.full((rows, cols), 0.05, dtype=np.float64)
    rain = np.zeros(10, dtype=np.float64)
    infil = np.zeros((rows, cols), dtype=np.float64)
    drain_cap = np.zeros((rows, cols), dtype=np.float64)

    cfg = SimConfig(
        z=z,
        dx=1.0,
        dt=0.5,
        t_end=5.0,
        rain=rain,
        infil=infil,
        drain_cap=drain_cap,
        k_flow=0.3,
        closed_boundary=False,
        h0=h0,
        save_every=2,
    )

    assert cfg.k_flow == 0.3
    assert cfg.closed_boundary is False
    assert cfg.h0 is not None
    assert cfg.h0.shape == (rows, cols)
    assert cfg.save_every == 2


def test_sim_result_can_be_created() -> None:
    """Verify that SimResult can be instantiated with matching output series."""
    n_saved = 11
    rows, cols = 4, 6

    times = np.linspace(0.0, 10.0, n_saved, dtype=np.float64)
    h_history = np.zeros((n_saved, rows, cols), dtype=np.float64)
    volume = np.zeros(n_saved, dtype=np.float64)
    cum_rain = np.zeros(n_saved, dtype=np.float64)
    cum_infil = np.zeros(n_saved, dtype=np.float64)
    cum_drain = np.zeros(n_saved, dtype=np.float64)
    cum_boundary = np.zeros(n_saved, dtype=np.float64)

    result = SimResult(
        times=times,
        h_history=h_history,
        volume=volume,
        cum_rain=cum_rain,
        cum_infil=cum_infil,
        cum_drain=cum_drain,
        cum_boundary=cum_boundary,
    )

    assert result.times.shape == (n_saved,)
    assert result.h_history.shape == (n_saved, rows, cols)
    assert result.volume.shape == (n_saved,)
    assert result.cum_rain.shape == (n_saved,)
    assert result.cum_infil.shape == (n_saved,)
    assert result.cum_drain.shape == (n_saved,)
    assert result.cum_boundary.shape == (n_saved,)


if __name__ == "__main__":
    test_sim_config_can_be_created()
    test_sim_config_with_custom_options()
    test_sim_result_can_be_created()
    print("All contract dataclass tests passed successfully!")
