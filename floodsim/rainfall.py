"""Rainfall hyetograph generators and meteorological forcing functions.

Provides synthetic storm rainfall profiles converted from standard hydrological units
(mm/hr, minutes) to SI physical simulation units (m/s, seconds).

Conversion:
    1 mm/hr = 2.78 × 10⁻⁷ m/s
"""

import numpy as np

# Hydrological conversion factor: 1 mm/hr in meters per second (m/s)
MM_HR_TO_M_S: float = 2.78e-7


def constant_storm(
    intensity_mm_hr: float,
    duration_min: float,
    dt: float,
    t_end: float,
) -> np.ndarray:
    """Generate a constant intensity rectangular rainfall hyetograph R(t).

    Mathematical Formulation:
        Convert input intensity to SI units:
            R_0 = intensity_mm_hr * 2.78e-7          [m/s]
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
        ValueError: If dt <= 0, t_end <= 0, intensity_mm_hr < 0, or duration_min < 0.
    """
    if dt <= 0:
        raise ValueError(f"dt must be positive, got {dt}")
    if t_end <= 0:
        raise ValueError(f"t_end must be positive, got {t_end}")
    if intensity_mm_hr < 0:
        raise ValueError(f"intensity_mm_hr must be non-negative, got {intensity_mm_hr}")
    if duration_min < 0:
        raise ValueError(f"duration_min must be non-negative, got {duration_min}")

    n_steps = int(np.ceil(t_end / dt))
    rate_m_s = float(intensity_mm_hr * MM_HR_TO_M_S)
    t_dur = float(duration_min * 60.0)

    times = np.arange(n_steps, dtype=np.float64) * dt
    rain = np.where(times < t_dur, rate_m_s, 0.0)
    return rain.astype(np.float64)


def triangular_storm(
    peak_mm_hr: float,
    duration_min: float,
    dt: float,
    t_end: float,
) -> np.ndarray:
    """Generate a symmetric triangular design storm hyetograph R(t).

    Mathematical Formulation:
        Convert peak intensity and storm duration to SI units:
            R_peak = peak_mm_hr * 2.78e-7          [m/s]
            t_dur = duration_min * 60              [s]
            t_half = t_dur / 2                     [s]

        The rainfall rate forms an isosceles triangle peaking at t_half:
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
        ValueError: If dt <= 0, t_end <= 0, peak_mm_hr < 0, or duration_min < 0.
    """
    if dt <= 0:
        raise ValueError(f"dt must be positive, got {dt}")
    if t_end <= 0:
        raise ValueError(f"t_end must be positive, got {t_end}")
    if peak_mm_hr < 0:
        raise ValueError(f"peak_mm_hr must be non-negative, got {peak_mm_hr}")
    if duration_min < 0:
        raise ValueError(f"duration_min must be non-negative, got {duration_min}")

    n_steps = int(np.ceil(t_end / dt))
    peak_m_s = float(peak_mm_hr * MM_HR_TO_M_S)
    t_dur = float(duration_min * 60.0)

    rain = np.zeros(n_steps, dtype=np.float64)
    if t_dur > 0 and peak_m_s > 0:
        times = np.arange(n_steps, dtype=np.float64) * dt
        t_half = t_dur / 2.0
        rising = (times >= 0.0) & (times < t_half)
        falling = (times >= t_half) & (times < t_dur)
        rain[rising] = peak_m_s * (times[rising] / t_half)
        rain[falling] = peak_m_s * ((t_dur - times[falling]) / t_half)

    return rain.astype(np.float64)
