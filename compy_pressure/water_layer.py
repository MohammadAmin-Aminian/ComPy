"""Water-layer pressure/acceleration ratio, its sensitivities and gain-fit bias.

ObsPy-free companion of ``compy_pressure.rayleigh_calibration`` for audits:
the same Zha & Webb (2016, eq. A5) ratio as ``rayleigh_calibration.water_ratio``
and ``microseism_check.admittance_height`` (no gravity), with

* :func:`log_sensitivity`: d ln(P/a)/d ln(sound speed), d ln(P/a)/d ln(phase
  velocity) and d ln(P/a)/d ln(depth), analytic;
* :func:`fit_gain`: the least-squares pressure multiplier of
  ``rayleigh_calibration.fit_multiplier`` (``weights="ls"``), or a fit of
  log ratios with equal weight per bin (``weights="log"``);
* :func:`gain_bias`: the multiplier a calibration returns when the true water
  column or crust differ from the model;
* :func:`column_mean_density`: depth-mean in-situ seawater density.

The ratio is P/a = rho H tan(x)/x with x = 2 pi f H sqrt(1/alpha² - 1/c²)
(tanh for c < alpha). Near the first water-column resonance (x -> pi/2) it
grows without bound, and so does its sensitivity to alpha, c and H.
"""

import numpy as np

GUARD = 0.9 * np.pi / 2  # rayleigh_calibration rejects bins at or beyond this x


def _x(f, phase_velocity, depth, sound_speed):
    f, c = np.broadcast_arrays(np.asarray(f, float), np.asarray(phase_velocity, float))
    q2 = sound_speed**-2.0 - c**-2.0
    return 2 * np.pi * f * depth * np.sqrt(np.abs(q2)), q2


def ratio(f, phase_velocity, depth, density=1028.0, sound_speed=1500.0):
    """P / upward acceleration [Pa s²/m] of a Rayleigh wave under a water layer."""
    x, q2 = _x(f, phase_velocity, depth, sound_speed)
    with np.errstate(divide="ignore", invalid="ignore"):
        shape = np.where(q2 >= 0, np.tan(x) / x, np.tanh(x) / x)
    shape = np.where(x < 1e-5, 1 + np.sign(q2) * x**2 / 3, shape)
    return density * depth * shape


def usable(f, phase_velocity, depth, sound_speed=1500.0):
    """Bins below the resonance guard of ``rayleigh_calibration.water_ratio``."""
    x, q2 = _x(f, phase_velocity, depth, sound_speed)
    return ~((q2 >= 0) & (x >= GUARD))


def log_sensitivity(f, phase_velocity, depth, sound_speed=1500.0):
    """Analytic d ln(P/a) / d ln(parameter) for sound speed, phase velocity, depth.

    With T = tan(x)/x, d ln T/d ln x = x/(sin x cos x) - 1 (x/(sinh x cosh x)
    - 1 for c < alpha), d ln x/d ln H = 1, d ln x/d ln alpha = -alpha^-2/q2
    and d ln x/d ln c = c^-2/q2, with q2 = alpha^-2 - c^-2.

    Returns:
        dict: ``sound_speed``, ``phase_velocity``, ``depth`` (arrays) and ``x``.
    """
    x, q2 = _x(f, phase_velocity, depth, sound_speed)
    c = np.broadcast_to(np.asarray(phase_velocity, float), x.shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        dlt = np.where(
            q2 >= 0, x / (np.sin(x) * np.cos(x)) - 1, x / (np.sinh(x) * np.cosh(x)) - 1
        )
    dlt = np.where(x < 1e-5, 2 * np.sign(q2) * x**2 / 3, dlt)
    return dict(
        x=x,
        depth=1 + dlt,  # P/a = rho H T(x): one power of H outside T
        sound_speed=dlt * (-(sound_speed**-2.0) / q2),
        phase_velocity=dlt * (c**-2.0 / q2),
    )


def fit_gain(observed, expected, weights="ls"):
    """Pressure multiplier g with g * observed ≈ expected.

    ``weights="ls"`` is ``rayleigh_calibration.fit_multiplier``: least squares
    on the ratios, so bins with the largest expected ratio (near resonance)
    dominate. ``weights="log"`` fits ln g = mean(ln expected - ln observed),
    one weight per bin.
    """
    r, t = np.asarray(observed, float), np.asarray(expected, float)
    if weights == "ls":
        return float(np.dot(r, t) / np.dot(r, r))
    if weights == "log":
        return float(np.exp(np.mean(np.log(t) - np.log(r))))
    raise ValueError("weights must be 'ls' or 'log'")


def gain_bias(f, phase_velocity, depth, *, model_density=1028.0, model_sound_speed=1500.0,
              true_density=None, true_sound_speed=None, true_phase_velocity=None,
              weights="ls"):  # fmt: skip
    """Multiplier returned for a perfect gauge when the model is wrong.

    The observed ratio is the true P/a (g = 1); the expected ratio uses the
    model water column and phase velocity. Bins past the guard are dropped
    as in ``rayleigh_calibration``. A value of 1.02 means the calibration
    reports a pressure multiplier 2 % too high.
    """
    c = np.broadcast_to(np.asarray(phase_velocity, float), np.shape(f))
    true_c = (
        c
        if true_phase_velocity is None
        else np.broadcast_to(np.asarray(true_phase_velocity, float), c.shape)
    )
    rho_t = model_density if true_density is None else true_density
    alpha_t = model_sound_speed if true_sound_speed is None else true_sound_speed
    keep = usable(f, c, depth, model_sound_speed) & usable(f, true_c, depth, alpha_t)
    if keep.sum() < 3:
        raise ValueError("fewer than three bins below the resonance guard")
    observed = ratio(np.asarray(f)[keep], true_c[keep], depth, rho_t, alpha_t)
    expected = ratio(
        np.asarray(f)[keep], c[keep], depth, model_density, model_sound_speed
    )
    return fit_gain(observed, expected, weights)


def column_mean_density(depth, *, deep_potential_density=1027.8, compressibility=4.4e-10,
                        gravity=9.81):  # fmt: skip
    """Depth-mean in-situ density [kg/m³] of a column of uniform potential density.

    In-situ density grows linearly with pressure, rho(z) = rho0 (1 + beta rho0
    g z), so the column mean is rho0 (1 + beta rho0 g H / 2). The warm upper
    ocean lowers it by a few tenths of kg/m³, which is neglected.
    """
    return deep_potential_density * (
        1
        + 0.5
        * compressibility
        * deep_potential_density
        * gravity
        * np.asarray(depth, float)
    )
