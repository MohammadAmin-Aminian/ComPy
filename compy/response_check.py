"""Check a pressure gauge's low-frequency response with two independent observables.

A differential pressure gauge (DPG) is a high-pass sensor: its leak (and, in
some metadata, a coupling stage) removes the very low frequencies. If the
metadata corner differs from the true one, deconvolving the metadata leaves

* a **phase** error in the compliance band: the Z/P phase of an elastic
  seafloor is real (0 or 180 degrees), so any frequency-dependent residual
  that is absent with absolute gauges points to the pressure response; and
* an **amplitude** error that grows toward low frequency, largest at tidal
  periods, where the tide itself is a known signal (Doran et al. 2019).

This module fits a factor ``k`` on the metadata poles below ``leak_hz`` to
the residual phase (:func:`fit_leak_factor`), predicts from it the tidal
amplitude ratio (:func:`predicted_ratio`), measures the tide in the raw counts
(:func:`tidal_harmonics`) and gives the compliance bias the fitted response
implies once a gain has been fitted in the Rayleigh-wave band
(:func:`compliance_bias`).

Responses are lists of poles-and-zeros stages ``(zeros, poles, radians)``;
:func:`paz_stages` builds them from an ObsPy channel, everything else is
NumPy only. In the whole-archive audit (docs/audit_2026-11-data/REPORT.md),
the tide measured in 145 DPG records agrees with the phase prediction
(Spearman 0.48), with a median offset of +3.4 dB that a one-factor model
does not explain: treat the response correction as a model, not a fact.
"""

from __future__ import annotations

import numpy as np

#: Angular speeds of the main constituents [deg/h].
TIDAL_SPEED_DEG_H = {"M2": 28.9841042, "S2": 30.0, "N2": 28.4397295,
                     "K1": 15.0410686, "O1": 13.9430356}  # fmt: skip
#: M2 frequency [Hz].
F_M2 = TIDAL_SPEED_DEG_H["M2"] / 360.0 / 3600.0


def paz_stages(channel):
    """Poles-and-zeros stages of an ObsPy channel as ``(zeros, poles, radians)``."""
    return [(np.array(s.zeros, complex), np.array(s.poles, complex),
             "RADIANS" in s.pz_transfer_function_type.upper())
            for s in channel.response.response_stages
            if s.__class__.__name__ == "PolesZerosResponseStage"]  # fmt: skip


def _pole_hz(p, rad):
    return np.abs(p) / (2 * np.pi) if rad else np.abs(p)


def analog_response(stages, f, k=1.0, *, mode="lowest", leak_hz=0.05):
    """Analog response of ``stages`` at ``f`` with low poles scaled by ``k``.

    Args:
        stages (list): ``(zeros, poles, radians)`` per stage.
        f (array): Frequencies [Hz].
        k (float): Factor on the scaled poles (k < 1: lower corner, longer
            time constant).
        mode (str): ``"lowest"`` scales only the lowest pole (the leak);
            ``"all"`` scales every pole below ``leak_hz``.
        leak_hz (float): Upper limit of the poles scaled in mode ``"all"``.

    Returns:
        ndarray: Complex response (no gain normalisation).
    """
    if mode not in ("lowest", "all"):
        raise ValueError("mode must be 'lowest' or 'all'")
    f = np.atleast_1d(np.asarray(f, float))
    poles_hz = [_pole_hz(p, rad) for _, ps, rad in stages for p in ps if abs(p) > 0]
    lowest = min(poles_hz, default=np.inf)
    h = np.ones(len(f), complex)
    for z, p, rad in stages:
        s = 1j * (2 * np.pi if rad else 1.0) * f
        hz = _pole_hz(p, rad)
        scale = (
            np.isclose(hz, lowest) if mode == "lowest" else (hz > 0) & (hz < leak_hz)
        )
        p = np.where(scale, p * k, p)
        h *= np.prod(s[:, None] - z[None, :], axis=1) / np.prod(
            s[:, None] - p[None, :], axis=1
        )
    return h


def predicted_residual(stages, f, k, **kw):
    """Residual Z/P phase [deg] if the true response has its poles scaled by ``k``.

    Deconvolving the metadata H_m from a record made with H_t leaves
    P_est = P H_t/H_m, so the Z/P phase is offset by phase(H_m) - phase(H_t).
    """
    hm, ht = analog_response(stages, f, **kw), analog_response(stages, f, k, **kw)
    return np.degrees(np.angle(hm / ht))


def fit_leak_factor(stages, f, residual_deg, *, ks=None, tolerance=1.5, **kw):
    """Least-squares ``k`` for the residual Z/P phase of a compliance stack.

    Args:
        stages (list): Metadata response of the pressure channel.
        f (array): Frequencies of the compliance band [Hz].
        residual_deg (array): Phase relative to 0/180 degrees
            (``compy.anomaly.residual_phase``).
        ks (array or None): Trial factors (default 0.02-20, log spaced).
        tolerance (float): Factors with RMS within ``tolerance`` times the
            best (or best + 0.5 degree) are acceptable.

    Returns:
        dict: ``k``, ``rms_deg``, ``k_low``, ``k_high`` (acceptable range),
        ``at_grid_edge`` (True if the best k is the first or last trial).
    """
    ks = np.geomspace(0.02, 20, 901) if ks is None else np.asarray(ks, float)
    res = np.asarray(residual_deg, float)
    rms = np.array(
        [
            np.sqrt(np.mean((res - predicted_residual(stages, f, k, **kw)) ** 2))
            for k in ks
        ]
    )
    i = int(np.argmin(rms))
    ok = rms <= max(tolerance * rms[i], rms[i] + 0.5)
    return {"k": float(ks[i]), "rms_deg": float(rms[i]), "k_low": float(ks[ok].min()),
            "k_high": float(ks[ok].max()), "at_grid_edge": i in (0, len(ks) - 1)}  # fmt: skip


def predicted_ratio(stages, k, f=F_M2, **kw):
    """|H_true / H_meta| at ``f``: the tidal amplitude ratio the fitted k predicts."""
    return float(
        np.abs(analog_response(stages, f, k, **kw) / analog_response(stages, f, **kw))[
            0
        ]
    )


def compliance_bias(stages, f, k, *, calibration_band=(0.02, 0.06), **kw):
    """Compliance bias [%] implied by ``k`` after a gain fitted in ``calibration_band``.

    Compliance scales as 1/|P_est| = |H_meta/H_true|; a Rayleigh-wave gain
    removes its mean over the calibration band, the rest stays.
    """
    fcal = np.geomspace(*calibration_band, 9)

    def r(x):
        return np.abs(
            analog_response(stages, x, **kw) / analog_response(stages, x, k, **kw)
        )

    return 100.0 * (r(f) / np.mean(r(fcal)) - 1.0)


def _spline_basis(h, spacing):
    """Cubic B-spline design matrix [len(h), n] with knots every ``spacing`` hours."""
    from scipy.interpolate import BSpline

    lo, hi = float(h.min()), float(h.max())
    inner = np.arange(lo, hi + spacing, spacing)
    if len(inner) < 2:
        inner = np.array([lo, hi + 1.0])
    knots = np.r_[[inner[0]] * 3, inner, [inner[-1]] * 3]
    m = BSpline.design_matrix(np.clip(h, inner[0], inner[-1]), knots, 3)
    m = m.toarray()
    return m[:, m.any(axis=0)]  # drop knots that fall in gaps


def tidal_harmonics(t, x, *, constituents=("M2", "S2", "N2", "K1", "O1"), detrend_hours=73,
                    n_controls=60, clip_mad=8.0):  # fmt: skip
    """Amplitudes of tidal constituents in an hourly series, with a noise level.

    Outliers are found against a running median of ``detrend_hours``: hours
    beyond ``clip_mad`` robust standard deviations, and hours without a
    running median (next to gaps), are excluded explicitly. The drift is
    fitted jointly with the tides (cubic B-spline, knots every
    ``detrend_hours``), by least squares (no nodal correction). The noise level is the median amplitude fitted alone at
    ``n_controls`` non-tidal periods between 8 and 40 h.

    Args:
        t (array): POSIX times [s] of hourly values (gaps allowed).
        x (array): Values (e.g. hourly mean counts).

    Returns:
        dict: ``amplitude`` {name: value}, ``noise``, ``n_used``, ``days``.
    """
    from compy.anomaly import control_periods

    t, x = np.asarray(t, float), np.asarray(x, float)
    order = np.argsort(t)
    t, x = t[order], x[order]
    half = detrend_hours * 3600 / 2
    base = np.full(len(t), np.nan)
    lo = np.searchsorted(t, t - half)
    hi = np.searchsorted(t, t + half, side="right")
    for i in range(len(t)):
        seg = x[lo[i] : hi[i]]
        seg = seg[np.isfinite(seg)]
        if len(seg) >= detrend_hours // 3:
            base[i] = np.median(seg)
    r = x - base
    fin = np.isfinite(r)
    if fin.sum() < 48:
        raise ValueError("too few hours with a running median")
    med = np.median(r[fin])
    mad = 1.4826 * np.median(np.abs(r[fin] - med))
    good = fin & (np.abs(np.where(fin, r, med) - med) <= clip_mad * mad + 1e-12)
    # The running median only finds outliers. A median of summed sinusoids is
    # not linear (it biased M2 by ~4 % in tests), so the drift is fitted with
    # the tides: a cubic B-spline with knots every ``detrend_hours``.
    h = (t - t[0]) / 3600.0
    drift = _spline_basis(h, detrend_hours)
    n_drift = drift.shape[1]
    cols = [drift]
    for c in constituents:
        w = np.radians(TIDAL_SPEED_DEG_H[c] * h)
        cols.append(np.column_stack([np.cos(w), np.sin(w)]))
    coef, *_ = np.linalg.lstsq(np.hstack(cols)[good], x[good], rcond=None)
    amp = {c: float(np.hypot(coef[n_drift + 2 * i], coef[n_drift + 2 * i + 1]))
           for i, c in enumerate(constituents)}  # fmt: skip
    noise = []
    for p in control_periods(
        [360.0 / TIDAL_SPEED_DEG_H[c] for c in constituents], n=n_controls
    ):
        B = np.column_stack(
            [drift, np.cos(2 * np.pi * h / p), np.sin(2 * np.pi * h / p)]
        )
        c2, *_ = np.linalg.lstsq(B[good], x[good], rcond=None)
        noise.append(np.hypot(c2[n_drift], c2[n_drift + 1]))
    return {"amplitude": amp, "noise": float(np.median(noise)), "n_used": int(good.sum()),
            "days": float((t[-1] - t[0]) / 86400.0)}  # fmt: skip


def pressure_sign_from_phase(phase_deg, *, max_residual=45.0):
    """Pressure sign for the Rayleigh calibration from the stacked Z/P phase.

    The compliance Z/P of an elastic seafloor is real; with metadata that give
    upward Z and positive compression it is near 180 degrees (sign +1 in
    ``compy_pressure``'s calibration), and near 0 degrees when the pressure
    polarity is inverted (sign -1). In the whole-archive run this rule matched
    the sign of every resolved stage D calibration but one (107 of 108; the
    exception has a 53 degree residual), and the local-P polarity of 61 of 64
    channels.

    Args:
        phase_deg (array): Phase of S_pz in the compliance band [deg].
        max_residual (float): Largest |phase - (0 or 180)| accepted [deg].

    Returns:
        int or None: -1, +1, or None when the phase is not near 0 or 180.
    """
    from compy.anomaly import residual_phase

    res, polarity = residual_phase(phase_deg)
    if not np.all(np.isfinite(res)) or np.max(np.abs(res)) > max_residual:
        return None
    return 1 if polarity == 180.0 else -1


def correct_compliance(
    compliance, stages, f, k, *, calibration_band=(0.02, 0.06), **kw
):
    """Compliance corrected for a pressure response whose low poles are scaled by ``k``.

    Equivalent to deconvolving the corrected response and refitting the
    Rayleigh gain in ``calibration_band``: compliance / (1 + bias), with the
    bias of :func:`compliance_bias`. Apply the same factor to uncertainties.
    """
    bias = compliance_bias(
        stages, np.asarray(f, float), k, calibration_band=calibration_band, **kw
    )
    return np.asarray(compliance, float) / (1.0 + bias / 100.0)


def scale_leak_poles(channel, k, *, mode="lowest", leak_hz=0.05):
    """Scale the low poles of an ObsPy channel's response in place (for stage A).

    The poles scaled are those of :func:`analog_response` (``mode``). The
    normalisation factor of each modified stage is recomputed so that the
    stage keeps unit gain at its normalisation frequency.

    Returns:
        int: Number of poles changed.
    """
    stages = channel.response.response_stages
    pz = [s for s in stages if s.__class__.__name__ == "PolesZerosResponseStage"]
    hz = [(_pole_hz(p, "RADIANS" in s.pz_transfer_function_type.upper()), s, i)
          for s in pz for i, p in enumerate(s.poles) if abs(p) > 0]  # fmt: skip
    lowest = min((h for h, _, _ in hz), default=np.inf)
    changed = 0
    for h, s, i in hz:
        if (mode == "lowest" and np.isclose(h, lowest)) or (
            mode == "all" and h < leak_hz
        ):
            s.poles[i] = s.poles[i] * k
            changed += 1
            rad = "RADIANS" in s.pz_transfer_function_type.upper()
            fn = s.normalization_frequency
            sn = 1j * (2 * np.pi if rad else 1.0) * fn
            h_n = np.prod([sn - z for z in s.zeros]) / np.prod(
                [sn - p for p in s.poles]
            )
            s.normalization_factor = float(1.0 / abs(h_n))
    return changed
