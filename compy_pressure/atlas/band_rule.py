"""Two confident limits for every 3-hour window of a pressure channel ("rule W").

The atlas keeps, for every channel, the part of the spectrum that can be trusted.
Rule W states that part as a low and a high frequency, from the instrument and
from the window itself, so that each limit has a physical reason and can be
checked against the metadata, the data or the phase of the pressure against the
vertical seismometer.

High limit ``f_hi`` (the lowest of)
    sampling      ``nyquist_fraction * fs / 2**bin_half_octave``: the upper edge
                  of the 1/4-octave bin stays under 0.8 Nyquist (anti-alias
                  filter of the digitiser; StationXML ``fs``)
    response      the frequency where the StationXML response of the epoch of
                  the window is ``response_drop_db`` (3 dB) under the nominal
                  sensitivity, on the short-period side
    self-noise    a cap by sensor type: the highest frequency up to which the
                  median spectrum of that type stays within ``tol_db`` (6 dB) of
                  the hydrophones (``empirical_caps``); on the archive of
                  November 2026 APG 0.745 Hz, DPG 13 Hz, hydrophone none

Low limit ``f_lo`` (the highest of)
    estimator     the 1/4-octave band must hold one Welch bin of the 2700-s
                  segment: ``f * (2**h - 2**-h) >= 1 / segment`` (h = 1/8), a
                  period of 468 s
    context       a window without the 6 h of context on each side has no level
                  beyond 300 s (``native.LONG_CONTEXT_PERIOD_S``)
    instrument    hydrophone: the cut-off validated on the pressure-vertical
                  phase, else the response (-3 dB); DPG: the response (-3 dB)
                  and a period of at most 300 s (the time constant implied by
                  the data is 241 s against 80 s in the metadata, so the levels
                  beyond 300 s are too high); APG: none

The digitiser was also tested: the response-corrected level of every channel is
at least 30 dB above the quantisation floor ``-10 log10(6 fs) - 10 log10 |R|^2``
(one-sided PSD of a uniform error of 1 count), so it never limits the band and
is not part of the rule.

Rule "V" (validated limits, 24 Nov 2026, after a four-reviewer review) is rule W with the changes that the evidence supports
    estimator     468 s for every sensor, as in rule W. A third-Welch-bin limit (900 s, ``apg_deep_depth_m``) is available but off: the
                  production estimator has a bias of +1 to +2 dB there on spectra as steep as the real ones, and DPG levels beyond
                  about 500 s are amplified noise (coherence with APGs falls from 0.97 to 0.46 between 486 and 817 s)
    DPG           unverified responses keep rule W (response -3 dB, at most 300 s). The SIO DPGs of the verified networks (column
                  ``dpg_validated``: 7D, response corner 29-31 s; the corrected level agrees with 12 depth-matched APGs, offsets
                  within the 2.5 dB tolerance to 490 s) are divided out to the estimator limit as long as every pole of the response
                  (``dpg_poles_s``: 29.5 s and a leak near 600 s) stays within ``dpg_tol_db`` for a time constant wrong by a factor
                  ``k`` (1.2 for the short pole, 3 for the leak, one-sided because the data give longer time constants than the
                  metadata). The high end is where the level is 6 dB above the gauge's own noise floor (``f_hi_data_hz``, from the
                  median spectrum), not a type-wide cap
    hydrophone    the cut-off validated on the pressure-vertical phase, but never beyond the period where the response correction
                  reaches ``hydrophone_max_correction_db`` (the phase is blind to a corner error once the correction is large)
    each limit records its evidence class in ``band_rule_lo_evidence`` and ``band_rule_hi_evidence``

A bin of a window is kept when its centre frequency lies in ``[f_lo, f_hi]``.
Frequencies are in Hz, periods in s, levels in dB re 1 unit^2/Hz.

Effect on the archive of November 2026 (381 gauges): the median spectrum moves
by less than 0.5 dB between 0.3 and 20 s; gauges covering 100-170 s fall from
270 to 175, those covering 3-10 s stay at 325, and nothing is left beyond 468 s.
See docs/PRESSURE_ATLAS.md.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from . import native

DEFAULTS = {
    "response_drop_db": 3.0,
    "nyquist_fraction": 0.4,
    "bin_half_octave": 0.125,
    "dpg_period_limit_s": 300.0,
    "context_period_s": native.LONG_CONTEXT_PERIOD_S,
    "self_noise_cap_hz": {"DPG": 13.0, "APG": 0.745, "Hydrophone": None},
    "self_noise_cap_hz_v": {"DPG": None, "APG": 0.745, "Hydrophone": None},
    "caps_tol_db": 6.0,
    "caps_min_gauges": 8,
    "caps_min_depth_m": 1000.0,
    "caps_start_period_s": 2.0,
    # rule "V" (validated limits); "W" keeps the behaviour of 24 Nov 2026
    "rule": "W",
    "estimator_period_s": None,  # V: None = one Welch bin in the 1/4-octave band (468 s)
    "apg_deep_depth_m": None,  # V: set (with apg_deep_period_s) to let deep absolute gauges reach 900 s; provisional (bias +1 to +2 dB)
    "apg_deep_period_s": 900.0,
    "dpg_tol_db": 2.5,
    "dpg_k_unvalidated": 3.0,
    "dpg_k_validated": 1.2,
    "hydrophone_max_correction_db": 10.0,
    "dpg_fit_from_period_s": 3.0,
}


def _cfg(cfg):
    out = dict(DEFAULTS)
    out.update(cfg or {})
    return out


def welch_limit_hz(segment_s=native.SEGMENT_S, half_octave=0.125):
    """Lowest centre frequency whose 1/4-octave band holds one Welch bin [Hz].

    The band is ``[f / 2**h, f * 2**h]`` with ``h = half_octave``; one Welch bin
    of a segment of ``segment_s`` is ``1 / segment_s`` wide. Below this
    frequency the smoothing is the Welch resolution, not 1/4 octave.
    """
    return 1.0 / segment_s / (2.0**half_octave - 2.0**-half_octave)


def response_band_hz(response_db, periods=None, drop_db=3.0):
    """Frequencies ``(f_lo, f_hi)`` [Hz] where a response leaves its pass band.

    Args:
        response_db: ``10 log10 |R|^2 / S^2`` at ``periods`` (response over the
            nominal sensitivity ``S``), as stored by the native estimator.
        periods: ascending periods [s] of the grid (default ``native.GRID``).
        drop_db: how far under the nominal sensitivity the band ends [dB].

    The pass band is the longest run of nodes whose level is at least
    ``-drop_db``. A crossing is placed by linear interpolation of the level in
    dB against log period between the last node inside and the first node
    outside. NaN means the response stays inside the band to the end of the
    grid on that side (no limit); both NaN when no node is inside.
    """
    P = np.asarray(native.GRID if periods is None else periods, float)
    db = np.asarray(response_db, float)
    ok = np.isfinite(db) & (db >= -drop_db)
    if not ok.any():
        return np.nan, np.nan
    runs, i = [], 0
    while i < len(ok):
        if ok[i]:
            j = i
            while j + 1 < len(ok) and ok[j + 1]:
                j += 1
            runs.append((i, j))
            i = j + 1
        else:
            i += 1
    a, b = max(runs, key=lambda r: r[1] - r[0])
    f_lo = f_hi = np.nan
    if a > 0:
        t = (-drop_db - db[a - 1]) / (db[a] - db[a - 1])
        f_hi = 1.0 / np.exp(np.log(P[a - 1]) + t * (np.log(P[a]) - np.log(P[a - 1])))
    if b < len(ok) - 1:
        t = (-drop_db - db[b]) / (db[b + 1] - db[b])
        f_lo = 1.0 / np.exp(np.log(P[b]) + t * (np.log(P[b + 1]) - np.log(P[b])))
    return f_lo, f_hi


def empirical_caps(
    median_db, periods, sensor_type, depth_m, reference="Hydrophone", cfg=None
):
    """Self-noise cap of each sensor type from the median spectra of the gauges [Hz].

    Args:
        median_db: ``(n_gauges, n_periods)`` median spectrum of each gauge, NaN
            outside its usable band.
        periods: ascending periods [s] of the columns.
        sensor_type: type of each gauge (``DPG``, ``APG``, ``Hydrophone``).
        depth_m: water depth of each gauge [m].
        reference: the type taken as free of self-noise at short periods.

    For each other type, the median over its gauges deeper than
    ``caps_min_depth_m`` is compared with the median of the reference type at
    each period, from ``caps_start_period_s`` toward short periods, where both
    have at least ``caps_min_gauges`` gauges. The cap is the highest frequency
    reached before the difference exceeds ``caps_tol_db`` or the data end.
    The reference gets ``inf``; a type that already differs at the start gets
    NaN.
    """
    c = _cfg(cfg)
    P = np.asarray(periods, float)
    M = np.asarray(median_db, float)
    typ = np.asarray(sensor_type)
    deep = np.asarray(depth_m, float) >= c["caps_min_depth_m"]

    def med(t):
        X = M[deep & (typ == t)]
        with np.errstate(all="ignore"):
            return np.nanmedian(X, axis=0), np.isfinite(X).sum(axis=0)

    ref, n_ref = med(reference)
    out = {reference: np.inf}
    i0 = int(np.argmin(np.abs(P - c["caps_start_period_s"])))
    for t in sorted(set(typ.tolist()) - {reference}):
        m, n = med(t)
        ok = (n >= c["caps_min_gauges"]) & (n_ref >= c["caps_min_gauges"])
        cap = None
        for i in range(i0, -1, -1):
            if not ok[i] or abs(m[i] - ref[i]) > c["caps_tol_db"]:
                break
            cap = P[i]
        out[t] = np.nan if cap is None else 1.0 / cap
    return out


def correction_db(period_s, corner_s):
    """Correction ``-10 log10 |H|^2`` of a first-order high-pass with its -3 dB corner at ``corner_s`` [dB]."""
    return 10.0 * np.log10(1.0 + (np.asarray(period_s, float) / corner_s) ** 2)


def dpg_poles_s(response_db, periods=None, fit_from_s=3.0):
    """Corner periods ``(T_short, T_long)`` [s] of the first-order high-pass poles that make up a DPG response.

    The response (dB under the nominal sensitivity) at periods of ``fit_from_s`` and longer is fitted by one pole, or by two when that
    lowers the rms misfit by more than half and by more than 0.3 dB (the Cascadia DPGs: 29.5 s and a leak at 595 s). ``T_long`` is
    NaN for one pole; both are NaN when there is no long-period roll-off.
    """
    P = np.asarray(native.GRID if periods is None else periods, float)
    y = -np.asarray(response_db, float)
    sel = (P >= fit_from_s) & np.isfinite(y)
    if sel.sum() < 10 or np.nanmax(y[sel]) < 1.0:
        return np.nan, np.nan
    Pc, yc = P[sel], y[sel]
    grid = np.logspace(np.log10(3.0), 5.0, 121)
    c = {t: correction_db(Pc, t) for t in grid}
    rms1 = {t: float(np.sqrt(np.mean((yc - c[t]) ** 2))) for t in grid}
    t1 = min(rms1, key=rms1.get)
    best = (rms1[t1], t1, np.nan)
    for i, ta in enumerate(grid[::2]):
        for tb in grid[::2]:
            if tb < ta:
                continue
            r = float(np.sqrt(np.mean((yc - c[ta] - c[tb]) ** 2)))
            if r < best[0] - 1e-12 and r < 0.5 * rms1[t1] and rms1[t1] - r > 0.3:
                best = (r, ta, tb)
    return (best[1], best[2]) if np.isfinite(best[2]) else (t1, np.nan)


def tolerance_period_s(corner_s, tol_db=2.5, k=3.0, two_sided=False):
    """Longest period [s] at which a pole with its corner at ``corner_s`` is divided out within ``tol_db`` of uncertainty.

    The time constant may be wrong by a factor ``k`` (longer: ``corner * k``; ``two_sided`` also shorter). The error of the correction
    at period T is ``correction_db(T, corner) - correction_db(T, k * corner)``. Infinity when ``20 log10 k`` stays under ``tol_db``.
    """
    T = corner_s * np.logspace(-1.0, 3.0, 4001)
    err = correction_db(T, corner_s) - correction_db(T, k * corner_s)
    if two_sided:
        err = np.maximum(err, np.abs(correction_db(T, corner_s) - correction_db(T, corner_s / k)))
    bad = np.where(err > tol_db)[0]
    if not len(bad):
        return np.inf
    i = int(bad[0])
    if i == 0:
        return float(T[0])
    t = (tol_db - err[i - 1]) / (err[i] - err[i - 1])
    return float(np.exp(np.log(T[i - 1]) + t * (np.log(T[i]) - np.log(T[i - 1]))))


def dpg_tolerance_limit_hz(response_db, periods=None, tol_db=2.5, k=3.0, k_short=None, fit_from_s=3.0):
    """Lowest frequency [Hz] at which every pole of a DPG response is divided out within ``tol_db``; 0 when there is no limit.

    ``k`` applies to every pole (one-sided); ``k_short`` (two-sided) replaces it for the short-period pole of a two-pole response
    when the response is verified (see the module docstring). Returns ``(f_lo_hz, T_short, T_long)``.
    """
    t_short, t_long = dpg_poles_s(response_db, periods, fit_from_s)
    if not np.isfinite(t_short):
        return 0.0, np.nan, np.nan
    lims = []
    if np.isfinite(t_long):
        lims.append(tolerance_period_s(t_short, tol_db, k_short or k, two_sided=k_short is not None))
        lims.append(tolerance_period_s(t_long, tol_db, k))
    else:
        lims.append(tolerance_period_s(t_short, tol_db, k))
    T = min(lims)
    return (0.0 if not np.isfinite(T) else 1.0 / T), t_short, t_long


class Limits:
    """Limits of one channel for each response epoch [Hz], and what sets them."""

    def __init__(self, f_lo_hz, f_hi_hz, lo_by, hi_by, lo_evidence=None, hi_evidence=None):
        self.f_lo_hz = np.asarray(f_lo_hz, float)
        self.f_hi_hz = np.asarray(f_hi_hz, float)
        self.lo_by = list(lo_by)
        self.hi_by = list(hi_by)
        self.lo_evidence = list(lo_evidence) if lo_evidence is not None else [""] * len(self.lo_by)
        self.hi_evidence = list(hi_evidence) if hi_evidence is not None else [""] * len(self.hi_by)


def channel_limits(
    sensor_type,
    fs,
    response_db,
    periods=None,
    validated_pz_hz=np.nan,
    cfg=None,
    depth_m=np.nan,
    dpg_validated=False,
    f_hi_data_hz=np.nan,
):
    """Low and high limit of a channel for each of its response epochs.

    Args:
        sensor_type: ``DPG``, ``APG`` or ``Hydrophone``.
        fs: sampling rate [Hz].
        response_db: ``(n_epochs, n_periods)`` or ``(n_periods,)``, see
            ``response_band_hz``.
        periods: ascending periods [s] (default ``native.GRID``).
        validated_pz_hz: low cut-off validated on the pressure-vertical phase
            [Hz] (hydrophones), NaN when there is none.
        cfg: overrides of ``DEFAULTS``.
        depth_m: water depth [m] (rule V: deep absolute gauges reach 900 s).
        dpg_validated: the DPG response is verified against APGs (rule V).
        f_hi_data_hz: high limit from the gauge's own noise floor [Hz] (rule V, DPG).
    """
    c = _cfg(cfg)
    rule_v = c["rule"] == "V"
    caps = {
        k: (np.inf if v is None else float(v))
        for k, v in c["self_noise_cap_hz_v" if rule_v else "self_noise_cap_hz"].items()
    }
    db = np.atleast_2d(np.asarray(response_db, float))
    cap = caps.get(sensor_type, np.inf)
    h1 = c["nyquist_fraction"] * fs / 2.0 ** c["bin_half_octave"]
    welch = welch_limit_hz(half_octave=c["bin_half_octave"])
    if rule_v and c["estimator_period_s"]:
        welch = 1.0 / float(c["estimator_period_s"])
    if rule_v and sensor_type == "APG" and c["apg_deep_depth_m"] and np.isfinite(depth_m) and depth_m >= c["apg_deep_depth_m"]:
        welch = 1.0 / float(c["apg_deep_period_s"])
    f_lo, f_hi, lo_by, hi_by, lo_ev, hi_ev = [], [], [], [], [], []
    for x in db:
        r_lo, r_hi = response_band_hz(x, periods, c["response_drop_db"])
        cands = {
            "sampling": h1,
            "response": r_hi if np.isfinite(r_hi) else np.inf,
            "self-noise": cap,
        }
        if rule_v and sensor_type == "DPG" and np.isfinite(f_hi_data_hz):
            cands["noise floor"] = float(f_hi_data_hz)
        k = min(cands, key=cands.get)
        f_hi.append(cands[k])
        hi_by.append(k)
        hi_ev.append(
            {
                "sampling": "physical bound",
                "response": "metadata",
                "self-noise": "intercomparison with hydrophones",
                "noise floor": "level 6 dB above the gauge noise floor",
            }[k]
        )
        r_lo = r_lo if np.isfinite(r_lo) else 0.0
        ev = "metadata"
        if sensor_type == "Hydrophone":
            if np.isfinite(validated_pz_hz):
                inst, why, ev = float(validated_pz_hz), "validated P-Z phase", "P-Z phase"
                if rule_v:
                    f_c, _ = response_band_hz(x, periods, c["hydrophone_max_correction_db"])
                    if np.isfinite(f_c) and f_c > inst:
                        inst, why, ev = f_c, "response correction cap", "metadata (10 dB correction cap)"
            else:
                inst, why, ev = r_lo, "response", "metadata (-3 dB)"
        elif sensor_type == "DPG" and rule_v and dpg_validated:
            inst, _t_s, _t_l = dpg_tolerance_limit_hz(
                x, periods, c["dpg_tol_db"], c["dpg_k_unvalidated"], c["dpg_k_validated"], c["dpg_fit_from_period_s"]
            )
            why, ev = "DPG response, verified", "APG intercomparison + tolerance"
        elif sensor_type == "DPG":
            tc = 1.0 / c["dpg_period_limit_s"]
            inst, why = (r_lo, "response") if r_lo >= tc else (tc, "DPG time constant")
        else:
            inst, why = 0.0, "none"
        if welch >= inst:
            f_lo.append(welch)
            lo_by.append("estimator")
            lo_ev.append("physical bound")
        else:
            f_lo.append(inst)
            lo_by.append(why)
            lo_ev.append(ev)
    return Limits(f_lo, f_hi, lo_by, hi_by, lo_ev, hi_ev)


def window_limits(limits, epoch_index, long_context_ok, cfg=None):
    """Limits ``(f_lo, f_hi)`` [Hz] of each window.

    The limits of the response epoch of the window, and, for a window without
    the context that periods beyond ``context_period_s`` need, a low limit of at
    least ``1 / context_period_s``.
    """
    c = _cfg(cfg)
    ep = np.asarray(epoch_index, int)
    ctx = np.asarray(long_context_ok, bool)
    f_lo = np.maximum(
        limits.f_lo_hz[ep], np.where(ctx, 0.0, 1.0 / c["context_period_s"])
    )
    return f_lo, limits.f_hi_hz[ep]


def band_mask(f_lo_hz, f_hi_hz, periods=None):
    """``(n_windows, n_periods)`` mask: the bins whose centre lies in the band."""
    f = 1.0 / np.asarray(native.GRID if periods is None else periods, float)
    return (f[None, :] >= np.asarray(f_lo_hz)[:, None]) & (
        f[None, :] <= np.asarray(f_hi_hz)[:, None]
    )


def add_limits(table, cfg=None, native_column="product"):
    # optional columns of the table: depth_m, dpg_validated (bool), f_hi_data_hz (rule V)
    """Channel table with the limits of rule W in ``f_lo_hz`` and ``f_hi_hz``.

    Reads the response of each channel from its native-rate file (column
    ``native_column``, path of the ``.npz`` of ``native.process_channel``). For a
    channel with several response epochs the most restrictive limits are kept.
    A hydrophone whose ``f_lo_source`` starts with ``validated`` keeps that
    ``f_lo_hz`` as the cut-off validated on the pressure-vertical phase.
    The earlier limits are kept in ``f_lo_hz_before`` and ``f_hi_hz_before``;
    ``band_rule_lo_by`` and ``band_rule_hi_by`` say what set each limit.
    """
    t = table.copy()
    if native_column not in t:
        raise ValueError(
            f"rule W needs the column {native_column!r} (native-rate file of each channel)"
        )
    t["f_lo_hz_before"] = t.get("f_lo_hz", np.nan)
    t["f_hi_hz_before"] = t.get("f_hi_hz", np.nan)
    lo, hi, lo_by, hi_by, lo_ev, hi_ev = [], [], [], [], [], []
    for r in t.itertuples():
        path = Path(str(getattr(r, native_column)))
        if not path.exists():
            lo.append(np.nan)
            hi.append(np.nan)
            lo_by.append("no native file")
            hi_by.append("no native file")
            lo_ev.append("")
            hi_ev.append("")
            continue
        with np.load(path, allow_pickle=True) as z:
            db = np.atleast_2d(np.asarray(z["response_db"], float))
            periods = np.asarray(z["periods"], float) if "periods" in z.files else None
        src = str(getattr(r, "f_lo_source", ""))
        pz = float(r.f_lo_hz) if src.startswith("validated") else np.nan
        lim = channel_limits(
            r.type,
            float(r.fs),
            db,
            periods,
            pz,
            cfg,
            depth_m=float(getattr(r, "depth_m", np.nan)),
            dpg_validated=bool(getattr(r, "dpg_validated", False)),
            f_hi_data_hz=float(getattr(r, "f_hi_data_hz", np.nan)),
        )
        k_lo, k_hi = int(np.argmax(lim.f_lo_hz)), int(np.argmin(lim.f_hi_hz))
        lo.append(float(lim.f_lo_hz[k_lo]))
        hi.append(float(lim.f_hi_hz[k_hi]))
        lo_by.append(lim.lo_by[k_lo])
        hi_by.append(lim.hi_by[k_hi])
        lo_ev.append(lim.lo_evidence[k_lo])
        hi_ev.append(lim.hi_evidence[k_hi])
    t["f_lo_hz"], t["f_hi_hz"] = lo, hi
    t["band_rule_lo_by"], t["band_rule_hi_by"] = lo_by, hi_by
    t["band_rule_lo_evidence"], t["band_rule_hi_evidence"] = lo_ev, hi_ev
    return t


def octaves(f_lo_hz, f_hi_hz):
    """Width of a band in octaves."""
    return math.log2(f_hi_hz / f_lo_hz)
