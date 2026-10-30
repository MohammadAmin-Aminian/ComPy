"""Water-column pressure/vertical relations and multitaper transfer estimates.

Pure functions for the local-P check (docs/local_p_method.md). SI units; z and
vertical motion positive upward; pressure is compression-positive excess
pressure. Gains follow the repository multiplier convention: g * p_nominal =
p_model, so g = |H_model| / |H_observed|.
"""

import numpy as np
from scipy.signal.windows import dpss

EARTH_RADIUS_KM = 6371.0
KM_PER_DEGREE = np.pi * EARTH_RADIUS_KM / 180.0
# Repository reference settings (Rayleigh/microseism checks; Deng et al. 2022).
DENSITY = 1028.0
SOUND_SPEED = 1500.0
# Sensitivity ranges propagated as model intervals, never fitted.
DENSITY_RANGE = (1020.0, 1050.0)
SOUND_SPEED_RANGE = (1450.0, 1540.0)
P_VELOCITY_MIN = 5000.0  # slowest P velocity assumed along a P ray, m/s
RESONANCE_MARGIN = 0.9  # theta <= 0.9*pi/2, as rayleigh_calibration.water_ratio
STEADY_STATE_ROUND_TRIPS = 10


def slowness_from_ray_parameter(ray_param_s_per_rad):
    """Horizontal slowness at the surface (s/m) from a TauP ray parameter."""
    value = float(ray_param_s_per_rad)
    if not np.isfinite(value) or value < 0:
        raise ValueError("Ray parameter must be finite and nonnegative")
    return value / (EARTH_RADIUS_KM * 1000.0)


def slowness_from_s_per_degree(ray_param_s_per_deg):
    """Horizontal slowness (s/m) from s/deg, e.g. Deng's q < 10 s/deg."""
    return float(ray_param_s_per_deg) / (KM_PER_DEGREE * 1000.0)


def vertical_cosine(q, sound_speed=SOUND_SPEED):
    """eta = sqrt(1 - q^2 c0^2) = cos(i2); reject evanescent waves in water."""
    q = np.asarray(q, dtype=float)
    if np.any(~np.isfinite(q)) or np.any(q < 0):
        raise ValueError("Slowness must be finite and nonnegative")
    if not np.isfinite(sound_speed) or sound_speed <= 0:
        raise ValueError("Sound speed must be positive and finite")
    x = (q * sound_speed) ** 2
    if np.any(x >= 1):
        raise ValueError("Wave is evanescent in water (q * c0 >= 1)")
    return np.sqrt(1 - x)


def eta_min_for_p(sound_speed=SOUND_SPEED, p_velocity_min=P_VELOCITY_MIN):
    """Smallest vertical cosine for any P-type wave with q <= 1/v_P,min."""
    return float(vertical_cosine(1 / p_velocity_min, sound_speed))


def _positive(**values):
    for name, value in values.items():
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")


def round_trip_time(depth, q=0.0, sound_speed=SOUND_SPEED):
    """Delay of the first sea-surface multiple at the seafloor: 2 h eta / c0."""
    _positive(depth=depth)
    return 2 * depth * float(vertical_cosine(q, sound_speed)) / sound_speed


def resonance_frequencies(depth, q=0.0, sound_speed=SOUND_SPEED, count=1):
    """Poles of p/a: f_n = (2n+1) c0 / (4 h eta)."""
    _positive(depth=depth)
    eta = float(vertical_cosine(q, sound_speed))
    n = np.arange(count)
    return (2 * n + 1) * sound_speed / (4 * depth * eta)


def water_phase(f, depth, q=0.0, sound_speed=SOUND_SPEED):
    """theta = 2 pi f h eta / c0."""
    _positive(depth=depth)
    f = np.asarray(f, dtype=float)
    if np.any(~np.isfinite(f)) or np.any(f < 0):
        raise ValueError("Frequencies must be finite and nonnegative")
    return 2 * np.pi * f * depth * vertical_cosine(q, sound_speed) / sound_speed


def tan_ratio(theta):
    """F = tan(theta)/theta with its series near zero; no pole guard."""
    theta = np.asarray(theta, dtype=float)
    out = np.ones_like(theta)
    small = np.abs(theta) < 1e-4
    out[small] = 1 + theta[small] ** 2 / 3 + 2 * theta[small] ** 4 / 15
    big = ~small
    out[big] = np.tan(theta[big]) / theta[big]
    return out


def steady_state_p_over_a(
    f,
    depth,
    q=0.0,
    density=DENSITY,
    sound_speed=SOUND_SPEED,
    margin=RESONANCE_MARGIN,
):
    """Deng Eq. 5 / Zha Eq. A5: p / a_up = rho h tan(theta)/theta (Pa per m/s^2).

    Steady state of the reverberating water column for one plane wave. Raises
    when any frequency reaches margin * (pi/2), the first-pole guard.
    """
    _positive(density=density)
    theta = water_phase(f, depth, q, sound_speed)
    if np.any(theta > margin * np.pi / 2):
        raise ValueError("Frequency reaches the first water-column resonance margin")
    return density * depth * tan_ratio(theta)


def first_transmission_p_over_v(q=0.0, density=DENSITY, sound_speed=SOUND_SPEED):
    """Upgoing plane wave before the first surface multiple: p / v_up = rho c0/eta.

    Valid only for the segment t - t_onset < round_trip_time (method note Eq. B).
    """
    _positive(density=density)
    return density * sound_speed / float(vertical_cosine(q, sound_speed))


def classify_regime(
    window_seconds,
    depth_bounds,
    eta_bounds,
    sound_speed_bounds,
    round_trips=STEADY_STATE_ROUND_TRIPS,
):
    """Return 'steady_state', 'pre_multiple' or 'mixed' for a window after onset.

    Pre-multiple requires the window to end before the earliest possible
    multiple (shallow depth, smallest eta, fastest c0). Steady state requires
    round_trips times the latest possible round trip.
    """
    h_lo, h_hi = depth_bounds
    e_lo, e_hi = eta_bounds
    c_lo, c_hi = sound_speed_bounds
    _positive(window=window_seconds, h_lo=h_lo, e_lo=e_lo, c_lo=c_lo)
    tau_min = 2 * h_lo * e_lo / c_hi
    tau_max = 2 * h_hi * e_hi / c_lo
    if window_seconds <= tau_min:
        return "pre_multiple", tau_min, tau_max
    if window_seconds >= round_trips * tau_max:
        return "steady_state", tau_min, tau_max
    return "mixed", tau_min, tau_max


def model_variation(f, half_bandwidth, depth, q=0.0, sound_speed=SOUND_SPEED):
    """Largest |F(f +- W)/F(f) - 1|; inf where f + W reaches the first pole."""
    f = np.asarray(f, dtype=float)
    centre = tan_ratio(water_phase(f, depth, q, sound_speed))
    worst = np.zeros_like(f)
    for shift in (-half_bandwidth, half_bandwidth):
        g = np.clip(f + shift, 0, None)
        theta = water_phase(g, depth, q, sound_speed)
        below = theta < np.pi / 2
        side = np.where(below, tan_ratio(np.where(below, theta, 0.0)), np.inf)
        worst = np.maximum(worst, np.abs(side / centre - 1))
    return worst


def coherence_zero_level(tapers, alpha=0.05):
    """Squared coherence exceeded by chance with probability alpha (K tapers)."""
    if tapers < 2:
        raise ValueError("At least two independent tapers are required")
    return 1 - alpha ** (1 / (tapers - 1))


def multitaper_spectra(x, y, sample_rate, nw=2.5, tapers=None, nfft=None):
    """Eigencoefficients of two equal-length series with DPSS tapers.

    Returns f and complex arrays (K, nf) for x and y, scaled so that
    mean(|X_k|^2) is a one-sided PSD. Averaging over K tapers gives 2K degrees
    of freedom; there is no single-periodogram coherence.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be equal-length 1-D arrays")
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("Nonfinite samples")
    n = len(x)
    k = int(tapers or (2 * nw - 1))
    if k < 2 or k > 2 * nw - 1 or n < 4 * nw:
        raise ValueError("Need 2 <= K <= 2NW-1 tapers and a longer window")
    nfft = int(nfft or n)
    if nfft < n:
        raise ValueError("nfft must not truncate the window")
    windows = dpss(n, nw, Kmax=k)
    scale = np.sqrt(2.0 / sample_rate)

    def coefficients(series):
        series = series - np.mean(series)
        trend = np.polyfit(np.arange(n), series, 1)
        series = series - np.polyval(trend, np.arange(n))
        return np.fft.rfft(windows * series, n=nfft, axis=1) * scale

    f = np.fft.rfftfreq(nfft, 1 / sample_rate)
    return f, coefficients(x), coefficients(y)


def cross_spectral_products(pk, vk):
    """Taper-averaged auto/cross spectra, H1, H2 and squared coherence."""
    spp = np.mean(np.abs(pk) ** 2, axis=0)
    svv = np.mean(np.abs(vk) ** 2, axis=0)
    spv = np.mean(pk * np.conj(vk), axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        h1 = spv / svv
        h2 = spp / np.conj(spv)
        coherence = np.abs(spv) ** 2 / (spp * svv)
    return {
        "spp": spp,
        "svv": svv,
        "spv": spv,
        "h1": h1,
        "h2": h2,
        "coherence": np.clip(np.nan_to_num(coherence), 0, 1),
    }


def noise_corrected_ratio(spp, svv, npp, nvv):
    """sqrt((S_pp - N_pp)/(S_vv - N_vv)); nan where noise dominates."""
    num, den = spp - npp, svv - nvv
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((num > 0) & (den > 0), np.sqrt(num / den), np.nan)


def polarity_class(phase_deg, tolerance_deg=45.0):
    """+1 near 0 deg, -1 near +-180 deg, 0 otherwise (indeterminate)."""
    wrapped = (np.asarray(phase_deg, dtype=float) + 180.0) % 360.0 - 180.0
    out = np.zeros_like(wrapped, dtype=int)
    out[np.abs(wrapped) <= tolerance_deg] = 1
    out[np.abs(wrapped) >= 180.0 - tolerance_deg] = -1
    return out


def circular_mean_deg(phase_deg):
    phase = np.deg2rad(np.asarray(phase_deg, dtype=float))
    if not phase.size:
        return float("nan")
    return float(np.rad2deg(np.angle(np.mean(np.exp(1j * phase)))))


def aic_onset(series, minimum=10):
    """Maeda (1985) AIC onset index of a 1-D series; None if undefined."""
    x = np.asarray(series, dtype=float)
    n = len(x)
    if n < 2 * minimum + 2 or not np.isfinite(x).all():
        return None
    c1 = np.cumsum(x)
    c2 = np.cumsum(x**2)
    k = np.arange(minimum, n - minimum)
    left_var = c2[k - 1] / k - (c1[k - 1] / k) ** 2
    right_n = n - k
    right_var = (c2[-1] - c2[k - 1]) / right_n - ((c1[-1] - c1[k - 1]) / right_n) ** 2
    with np.errstate(divide="ignore"):
        aic = k * np.log(np.maximum(left_var, 1e-300)) + (right_n - 1) * np.log(
            np.maximum(right_var, 1e-300)
        )
    return int(k[np.argmin(aic)])
