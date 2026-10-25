"""Rainfall hyetograph generators and meteorological forcing functions.

Provides synthetic storm rainfall profiles converted from standard hydrological units
(mm/hr, minutes) to SI physical simulation units (m/s, seconds).

Conversion:
    1 mm/hr = 10⁻³ m / 3600 s ≈ 2.777778 × 10⁻⁷ m/s
"""

import numpy as np


def constant_storm(
    intensity_mm_hr: float,
    duration_min: float,
    dt: float,
    t_end: float,
) -> np.ndarray:
    """Generate a constant intensity rectangular rainfall hyetograph R(t).

    Mathematical Formulation:
        Convert input intensity to SI units:
            R_0 = intensity_mm_hr / (1000 * 3600)   [m/s]
            t_dur = duration_min * 60                [s]

        The rainfall rate at time t is a step function:
            R(t) = R_0   if 0 <= t < t_dur
            R(t) = 0.0   if t >= t_dur

        Discretized into n_steps = int(ceil(t_end / dt)) time increments:
            R[k] = R(k * dt),   k = 0, ..., n_steps - 1.

    Args:
        intensity_mm_hr (float): Rainfall intensity in millimeters per hour (mm/hr).
        duration_min (float): Storm duration in minutes (min).
        dt (float): Simulation discrete time step in seconds (s).
        t_end (float): Total simulation duration in seconds (s).

    Returns:
        np.ndarray: 1D float64 array of shape (n_steps,) representing rainfall rate R in m/s.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("constant_storm is not yet implemented.")


def triangular_storm(
    peak_mm_hr: float,
    duration_min: float,
    dt: float,
    t_end: float,
) -> np.ndarray:
    """Generate a symmetric triangular design storm hyetograph R(t).

    Mathematical Formulation:
        Convert peak intensity and storm duration to SI units:
            R_peak = peak_mm_hr / (1000 * 3600)   [m/s]
            t_dur = duration_min * 60              [s]
            t_half = t_dur / 2                     [s]

        The rainfall rate forms an isosceles triangle:
            R(t) = R_peak * (t / t_half)                  if 0 <= t < t_half
            R(t) = R_peak * (t_dur - t) / t_half          if t_half <= t < t_dur
            R(t) = 0.0                                    if t >= t_dur

        Total rainfall depth deposited:
            P_total = 0.5 * R_peak * t_dur   [m]

        Discretized into n_steps = int(ceil(t_end / dt)) time increments:
            R[k] = R(k * dt),   k = 0, ..., n_steps - 1.

    Args:
        peak_mm_hr (float): Peak storm intensity at t_dur/2 in millimeters per hour (mm/hr).
        duration_min (float): Total rainfall duration in minutes (min).
        dt (float): Simulation discrete time step in seconds (s).
        t_end (float): Total simulation duration in seconds (s).

    Returns:
        np.ndarray: 1D float64 array of shape (n_steps,) representing rainfall rate R in m/s.

    Raises:
        NotImplementedError: Stub implementation for architectural contract.
    """
    raise NotImplementedError("triangular_storm is not yet implemented.")
