"""Validated numerical primitives shared by ComPy's public modules.

All physical quantities use SI units. Coherence arguments are explicitly
amplitude coherence, rather than SciPy's magnitude-squared coherence.
"""

import operator

import numpy as np
from scipy.optimize import brentq


def positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def integer(value, name, minimum=1):
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if isinstance(value, (bool, np.bool_)) or result < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return result


def residual(d, m, s=1):
    d, m = np.asarray(d, dtype=float), np.asarray(m, dtype=float)
    s = np.asarray(s, dtype=float)
    if d.shape != m.shape or d.size == 0:
        raise ValueError("data and model must have equal, nonempty shapes")
    if not np.all(np.isfinite(d)) or not np.all(np.isfinite(m)):
        raise ValueError("data and model must be finite")
    if not np.all(np.isfinite(s)) or np.any(s <= 0):
        raise ValueError("uncertainty must be finite and positive")
    try:
        s = np.broadcast_to(s, d.shape)
    except ValueError as exc:
        raise ValueError("uncertainty must broadcast to the data shape") from exc
    return (d - m) / s


def log_likelihood(d, m, s=1):
    """Unnormalized independent Gaussian log likelihood for fixed uncertainty."""
    r = residual(d, m, s)
    return -0.5 * np.sum(r * r)


def roughness(vs, order=2):
    """Squared derivative energy on a uniform, unit-spaced depth grid."""
    order = integer(order, "order")
    values = np.asarray(vs, dtype=float).reshape(-1)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("profile must be nonempty and finite")
    if values.size == 1:
        return 0.0
    for _ in range(order):
        values = np.gradient(values)
    return float(np.sum(values * values))


def wavenumber(omega, depth, gravity=9.79329):
    """Solve omega**2 = g*k*tanh(k*depth), including DC and scalar inputs.

    A dimensionless bracketed solve avoids the previous first-bin assumption
    and unbounded secant iteration. Returns an array with the input shape.
    """
    depth, gravity = positive(depth, "depth"), positive(gravity, "gravity")
    omega = np.asarray(omega, dtype=float)
    if np.any(~np.isfinite(omega)) or np.any(omega < 0):
        raise ValueError("angular frequencies must be finite and nonnegative")
    q = omega.ravel() ** 2 * depth / gravity
    roots = np.zeros_like(q)
    for i, value in enumerate(q):
        if value == 0:
            continue
        if value >= 20:
            roots[i] = value  # tanh(value) differs from 1 below double precision
        else:
            scale = max(np.sqrt(value), value)
            # Scale the unknown so tiny positive frequencies retain accuracy.
            roots[i] = scale * brentq(
                lambda y: y * np.tanh(scale * y) / (value / scale) - 1,
                0,
                2,
                xtol=1e-14,
                rtol=1e-14,
            )
    return (roots / depth).reshape(omega.shape)


def sliding_window(a, ws, ss=None, hann=True):
    """Copy complete windows; taper integer input in floating point."""
    ws = integer(ws, "window size")
    ss = ws if ss is None else integer(ss, "step size")
    a = np.asarray(a)
    if a.ndim != 1 or len(a) < ws:
        raise ValueError("input must be a 1D array at least as long as the window")
    out = np.lib.stride_tricks.sliding_window_view(a, ws)[::ss].copy()
    if hann:
        out = out * np.hanning(ws)
    return out, len(out)


def compliance_uncertainty(compliance, coherence, n_windows):
    """Amplitude-coherence error estimate; zero coherence gives infinity.

    n_windows is the independent spectral-average count. Overlapping windows
    must not be interpreted as independent observations without correction.
    """
    n_windows = positive(n_windows, "n_windows")
    com, coh = np.broadcast_arrays(
        np.asarray(compliance, float), np.asarray(coherence, float)
    )
    if (
        np.any(~np.isfinite(com))
        or np.any(~np.isfinite(coh))
        or np.any((coh < 0) | (coh > 1))
    ):
        raise ValueError(
            "compliance must be finite; amplitude coherence must be in [0, 1]"
        )
    result = np.full(com.shape, np.inf)
    np.divide(
        np.abs(com) * np.sqrt(1 - coh**2),
        coh * np.sqrt(2 * n_windows),
        out=result,
        where=coh > 0,
    )
    return result


def plot_smooth(data, window_length, polyorder, **kwargs):
    """Savitzky-Golay display smoothing bounded by the available axis length.

    Return unchanged values when there are too few samples for the requested
    polynomial. This helper is for plots, not inversion or scientific estimators.
    """
    from scipy.signal import savgol_filter

    values = np.asarray(data)
    axis = kwargs.get("axis", -1)
    width = min(int(window_length), values.shape[axis])
    if width % 2 == 0:
        width -= 1
    if width <= polyorder:
        return values.copy()
    return savgol_filter(values, width, polyorder, **kwargs)
