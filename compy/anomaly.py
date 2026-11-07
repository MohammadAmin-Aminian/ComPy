"""Automatic detection of anomalies and unexpected correlations in OBS products.

The detectors work on what stage C archives (``compy.results``): per-window
transfer-function and PSD levels, accepted-window times, monthly stacks and
the stacked spectra. They use robust statistics so that a few storms or
glitches do not set the scale:

* :func:`robust_z`: median/MAD z-scores.
* :func:`changepoints`: binary segmentation for shifts in the median of a
  series, with a penalty inflated for serial correlation.
* :func:`periodicity`: amplitude of a sinusoid at tidal and solar periods,
  with an empirical p-value against many nearby non-tidal periods, which
  keeps red (storm-driven) variability from looking periodic.
* :func:`rayleigh_periodicity`: the same for event times (e.g. accepted
  windows), with the Rayleigh statistic of circular uniformity.
* :func:`level_dependence`: does the transfer function depend on the
  pressure level (pressure-gauge noise or nonlinearity)?
* :func:`kmeans`: clusters of windows with distinct behaviour.
* :func:`phase_models`: fits instrument-response hypotheses to the residual
  Z/P phase and converts each into the compliance bias it implies.
* :func:`common_mode`: cross-station correlation of monthly deviations.
* :func:`simulate_gate_bias`: Monte Carlo of the bias that a per-window
  coherence gate gives the stacked H1 estimate.

:func:`scan_station` and :func:`scan_network` run them and return
:class:`Anomaly` records; :func:`rank` orders them by a transparent score
(statistical strength times the effect on compliance) and
:func:`to_markdown` / :func:`to_csv` write the table.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass, field

import numpy as np
from scipy import stats

#: Tidal and solar periods [h] tested by default (Doodson/Darwin names).
TIDAL_PERIODS_H = {
    "M2": 12.4206012,
    "S2": 12.0,
    "N2": 12.65834751,
    "K1": 23.93447213,
    "O1": 25.81933871,
    "M4": 6.210300601,
    "S1": 24.0,
}

#: Period at which the stated compliance effect of a phase model is referenced.
CALIBRATION_BAND_HZ = (0.02, 0.06)


# --------------------------------------------------------------------------
# Records and ranking
# --------------------------------------------------------------------------
@dataclass
class Anomaly:
    """One flagged behaviour of a station or of the network.

    Attributes:
        station (str): Station code, or "network".
        category (str): "instrument", "noise", "temporal", "spectral",
            "selection" or "network".
        kind (str): Detector name, e.g. "changepoint", "tidal_acceptance".
        summary (str): One-line description with the numbers.
        significance (float): Gaussian-equivalent strength (|z| or the
            normal quantile of an empirical p-value), capped at 10.
        effect_pct (float): Size of the effect on compliance [%] where it
            applies, otherwise on the tested quantity; NaN if unknown.
        time (str): ISO time of the event, if any.
        evidence (dict): Numbers behind the record.
    """

    station: str
    category: str
    kind: str
    summary: str
    significance: float
    effect_pct: float = float("nan")
    time: str = ""
    evidence: dict = field(default_factory=dict)

    @property
    def score(self):
        """Significance (capped at 8) times 1 + |effect|/5 % (effect capped at 20 %)."""
        effect = 0.0 if not np.isfinite(self.effect_pct) else abs(self.effect_pct)
        return min(self.significance, 8.0) * (1.0 + min(effect, 20.0) / 5.0)


def p_to_sigma(p):
    """Two-sided Gaussian-equivalent significance of a p-value (capped at 10)."""
    p = float(np.clip(p, 1e-23, 1.0))
    return float(min(stats.norm.isf(p / 2), 10.0))


def rank(anomalies, min_significance=2.0):
    """Anomalies at or above ``min_significance``, highest score first."""
    keep = [a for a in anomalies if a.significance >= min_significance]
    return sorted(keep, key=lambda a: (-a.score, a.station, a.kind))


def to_rows(anomalies):
    """Flat dictionaries (rank, score and fields) for tables."""
    rows = []
    for i, a in enumerate(anomalies, 1):
        row = {"rank": i, "score": round(a.score, 2), **asdict(a)}
        row["significance"] = round(a.significance, 2)
        row["effect_pct"] = (
            None if not np.isfinite(a.effect_pct) else round(a.effect_pct, 2)
        )
        row["evidence"] = {k: _short(v) for k, v in a.evidence.items()}
        rows.append(row)
    return rows


def _short(value):
    if isinstance(value, float):
        return float(f"{value:.4g}")
    if isinstance(value, (list, tuple)) and len(value) > 12:
        return f"[{len(value)} values]"
    return value


def to_markdown(anomalies, max_rows=None):
    """Markdown table of ranked anomalies."""
    rows = to_rows(anomalies)[:max_rows]
    lines = [
        "| # | Station | Category | Detector | Finding | σ | Effect [%] | Score |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        effect = "–" if r["effect_pct"] is None else f"{r['effect_pct']:+.2f}"
        lines.append(
            f"| {r['rank']} | {r['station']} | {r['category']} | {r['kind']} | "
            f"{r['summary']} | {r['significance']:.1f} | {effect} | {r['score']:.1f} |"
        )
    return "\n".join(lines)


def to_csv(anomalies, path):
    """Write ranked anomalies to ``path`` as CSV."""
    rows = to_rows(anomalies)
    names = ["rank", "score", "station", "category", "kind", "summary",
             "significance", "effect_pct", "time", "evidence"]  # fmt: skip
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=names)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r[k] for k in names})


# --------------------------------------------------------------------------
# Robust statistics
# --------------------------------------------------------------------------
def mad_sigma(x):
    """1.4826 x median absolute deviation (Gaussian-consistent scale)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not x.size:
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def robust_z(x):
    """(x - median) / (1.4826 MAD); falls back to the standard deviation."""
    x = np.asarray(x, float)
    scale = mad_sigma(x)
    if not np.isfinite(scale) or scale == 0:
        scale = float(np.nanstd(x)) or 1.0
    return (x - np.nanmedian(x)) / scale


def lag1_autocorrelation(x):
    """Lag-1 autocorrelation of a series (finite values, mean removed)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)] - np.nanmean(x)
    if len(x) < 3 or not np.any(x):
        return 0.0
    return float(np.clip(np.dot(x[:-1], x[1:]) / np.dot(x, x), -0.99, 0.99))


# --------------------------------------------------------------------------
# Changepoints
# --------------------------------------------------------------------------
@dataclass
class Changepoint:
    """A shift in the median of a series."""

    index: int
    time: float
    shift: float
    before: float
    after: float
    significance: float
    n_before: int
    n_after: int


def _best_split(x, min_size):
    n = len(x)
    if n < 2 * min_size:
        return None, 0.0
    c1 = np.cumsum(x)
    c2 = np.cumsum(x * x)
    k = np.arange(min_size, n - min_size + 1)
    left = c2[k - 1] - c1[k - 1] ** 2 / k
    right = (c2[-1] - c2[k - 1]) - (c1[-1] - c1[k - 1]) ** 2 / (n - k)
    total = c2[-1] - c1[-1] ** 2 / n
    gain = total - (left + right)
    i = int(np.argmax(gain))
    return int(k[i]), float(gain[i])


def changepoints(t, x, *, min_size=50, max_changes=6, threshold=5.0):
    """Shifts in the level of ``x`` by binary segmentation.

    Each candidate split maximizes the reduction of the squared residuals
    (least-squares change in mean) on robustly clipped data. It is kept
    when the robust two-sample statistic of the medians on either side,
    deflated for lag-1 autocorrelation rho by sqrt((1 - rho)/(1 + rho)),
    exceeds ``threshold`` (a Gaussian-equivalent significance).

    Args:
        t (array): Times (sorted).
        x (array): Values (e.g. log10 transfer function per window).
        min_size (int): Minimum segment length.
        max_changes (int): Largest number of changepoints.
        threshold (float): Required significance.

    Returns:
        list of Changepoint: Sorted by time.
    """
    t, x = np.asarray(t, float), np.asarray(x, float)
    ok = np.isfinite(t) & np.isfinite(x)
    t, x = t[ok], x[ok]
    order = np.argsort(t)
    t, x = t[order], x[order]
    if len(x) < 2 * min_size:
        return []
    med, scale = np.median(x), mad_sigma(x) or float(np.std(x)) or 1.0
    clipped = np.clip(x, med - 5 * scale, med + 5 * scale)
    segments, found = [(0, len(x))], []
    while len(found) < max_changes:
        best = None
        for a, b in segments:
            k, gain = _best_split(clipped[a:b], min_size)
            if k is None:
                continue
            k += a
            left, right = x[a:k], x[k:b]
            rho = lag1_autocorrelation(np.concatenate([left - np.median(left),
                                                       right - np.median(right)]))  # fmt: skip
            s = mad_sigma(
                np.concatenate([left - np.median(left), right - np.median(right)])
            )
            s = s or float(np.std(x)) or 1.0
            se = (
                s
                * math.sqrt(1 / len(left) + 1 / len(right))
                * math.sqrt((1 + rho) / (1 - rho))
            )
            z = abs(np.median(right) - np.median(left)) / (1.2533 * se)
            if best is None or z > best[0]:
                best = (z, k, a, b)
        if best is None or best[0] < threshold:
            break
        z, k, a, b = best
        segments.remove((a, b))
        segments += [(a, k), (k, b)]
        found.append((k, z))
    out = []
    bounds = sorted({0, len(x), *[k for k, _ in found]})
    for k, z in sorted(found):
        i = bounds.index(k)
        left, right = x[bounds[i - 1] : k], x[k : bounds[i + 1]]
        out.append(Changepoint(index=int(order[k]), time=float(t[k]),
                               shift=float(np.median(right) - np.median(left)),
                               before=float(np.median(left)), after=float(np.median(right)),
                               significance=float(min(z, 10.0)),
                               n_before=len(left), n_after=len(right)))  # fmt: skip
    return out


def segment_residuals(t, x, cps):
    """``x`` minus the median of its segment between changepoints."""
    t, x = np.asarray(t, float), np.asarray(x, float)
    out = np.full_like(x, np.nan)
    edges = [-np.inf, *sorted(c.time for c in cps), np.inf]
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (t >= a) & (t < b) & np.isfinite(x)
        if sel.any():
            out[sel] = x[sel] - np.median(x[sel])
    return out


# --------------------------------------------------------------------------
# Periodicity
# --------------------------------------------------------------------------
def _sinusoid_amplitude(t_h, x, period_h):
    w = 2 * np.pi / period_h
    design = np.column_stack([np.ones_like(t_h), np.cos(w * t_h), np.sin(w * t_h)])
    coef, *_ = np.linalg.lstsq(design, x, rcond=None)
    return float(np.hypot(coef[1], coef[2])), float(
        np.degrees(np.arctan2(coef[2], coef[1]))
    )


def control_periods(targets_h, span=(8.0, 40.0), n=300, exclude_frac=0.03):
    """Non-tidal periods [h] in ``span``, away from targets and their halves."""
    grid = np.geomspace(*span, n)
    lines = list(targets_h) + [p / 2 for p in targets_h] + [2 * p for p in targets_h]
    keep = np.ones(len(grid), bool)
    for p in lines:
        keep &= np.abs(grid / p - 1) > exclude_frac
    return grid[keep]


def periodicity(t, x, periods_h=None, *, controls_h=None):
    """Sinusoidal modulation of a series at given periods, with empirical p-values.

    The least-squares amplitude at each target period is compared with the
    amplitudes at many control periods (:func:`control_periods`); the
    p-value is the fraction of controls at least as large. Storm-driven
    (red) variability raises the controls too, so it does not register as
    periodic, unlike with an F-test that assumes white residuals.

    Args:
        t (array): Times [s].
        x (array): Values.
        periods_h (dict): ``{name: period [h]}``; :data:`TIDAL_PERIODS_H`.
        controls_h (array): Control periods [h].

    Returns:
        dict: ``{name: {"period_h", "amplitude", "phase_deg", "p_value",
        "p_tail", "null_median"}}``. ``p_value`` is the empirical fraction of
        controls (it cannot go below 1/(n_controls + 1)); ``p_tail`` extends
        it with an exponential law for amplitude² scaled to the controls.
    """
    periods_h = TIDAL_PERIODS_H if periods_h is None else periods_h
    t, x = np.asarray(t, float), np.asarray(x, float)
    ok = np.isfinite(t) & np.isfinite(x)
    t_h, x = (t[ok] - t[ok].min()) / 3600.0, x[ok] - np.mean(x[ok])
    controls_h = (
        control_periods(periods_h.values()) if controls_h is None else controls_h
    )
    null = np.array([_sinusoid_amplitude(t_h, x, p)[0] for p in controls_h])
    scale = np.median(null**2) / np.log(2)  # exponential null for amplitude**2
    out = {}
    for name, p in periods_h.items():
        amp, ph = _sinusoid_amplitude(t_h, x, p)
        p_emp = (1 + np.sum(null >= amp)) / (1 + len(null))
        out[name] = dict(period_h=p, amplitude=amp, phase_deg=ph, p_value=float(p_emp),
                         p_tail=float(np.exp(-amp**2 / scale)),
                         null_median=float(np.median(null)))  # fmt: skip
    return out


def rayleigh_periodicity(times, periods_h=None, *, controls_h=None):
    """Clustering of event times at a phase of given periods.

    For each period the Rayleigh resultant length R = |mean(exp(i 2 pi t/T))|
    and Z = n R² are computed, with an empirical p-value against control
    periods. A significant Z means events (e.g. accepted windows) favour
    one part of the cycle.

    Returns:
        dict: ``{name: {"period_h", "R", "Z", "phase_deg", "p_value",
        "p_tail", "p_rayleigh"}}``; ``p_tail`` scales the exponential law of
        R² to the controls (overdispersion from clustered gaps), and
        ``p_rayleigh`` assumes independent events.
    """
    periods_h = TIDAL_PERIODS_H if periods_h is None else periods_h
    t_h = np.asarray(times, float) / 3600.0
    t_h = t_h[np.isfinite(t_h)]
    n = len(t_h)
    controls_h = (
        control_periods(periods_h.values()) if controls_h is None else controls_h
    )

    def resultant(p):
        z = np.mean(np.exp(2j * np.pi * t_h / p))
        return abs(z), float(np.degrees(np.angle(z)))

    null = np.array([resultant(p)[0] for p in controls_h])
    scale = np.median(null**2) / np.log(2)  # overdispersed Rayleigh null
    out = {}
    for name, p in periods_h.items():
        r, ph = resultant(p)
        out[name] = dict(period_h=p, R=r, Z=n * r * r, phase_deg=ph,
                         p_value=float((1 + np.sum(null >= r)) / (1 + len(null))),
                         p_tail=float(np.exp(-r * r / scale)),
                         p_rayleigh=float(np.exp(-n * r * r)),
                         null_median_R=float(np.median(null)))  # fmt: skip
    return out


# --------------------------------------------------------------------------
# Level dependence and clustering
# --------------------------------------------------------------------------
def level_dependence(level, value, n_bins=8, min_count=20):
    """Binned medians of ``value`` against ``level`` and the upper-half slope.

    The slope uses only the upper half of the bins, where the per-window
    estimator bias from incoherent noise (larger at low signal) is small.

    Returns:
        dict: ``centres``, ``medians``, ``counts``, ``slope`` (value units
        per level unit), ``spearman`` and its ``p_value`` over all windows.
    """
    level, value = np.asarray(level, float), np.asarray(value, float)
    ok = np.isfinite(level) & np.isfinite(value)
    level, value = level[ok], value[ok]
    edges = np.quantile(level, np.linspace(0, 1, n_bins + 1))
    idx = np.clip(np.searchsorted(edges, level, side="right") - 1, 0, n_bins - 1)
    centres, medians, counts = [], [], []
    for b in range(n_bins):
        sel = idx == b
        if sel.sum() >= min_count:
            centres.append(float(np.median(level[sel])))
            medians.append(float(np.median(value[sel])))
            counts.append(int(sel.sum()))
    slope = float("nan")
    if len(centres) >= 4:
        h = len(centres) // 2
        slope = float(np.polyfit(centres[h:], medians[h:], 1)[0])
    rho, p = stats.spearmanr(level, value) if len(level) > 3 else (np.nan, np.nan)
    return dict(centres=centres, medians=medians, counts=counts, slope=slope,
                spearman=float(rho), p_value=float(p))  # fmt: skip


def kmeans(features, k=3, *, seed=0, n_init=10, n_iter=100):
    """K-means on robust-scaled features (k-means++ starts, best of n_init).

    Returns:
        tuple: (labels, centres in scaled units, inertia).
    """
    x = np.asarray(features, float)
    x = np.column_stack([robust_z(col) for col in x.T])
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_init):
        centres = [x[rng.integers(len(x))]]
        for _ in range(1, k):
            d2 = np.min([np.sum((x - c) ** 2, axis=1) for c in centres], axis=0)
            centres.append(x[rng.choice(len(x), p=d2 / d2.sum())])
        centres = np.array(centres)
        for _ in range(n_iter):
            labels = np.argmin(((x[:, None, :] - centres[None]) ** 2).sum(-1), axis=1)
            new = np.array([x[labels == j].mean(0) if np.any(labels == j) else centres[j]
                            for j in range(k)])  # fmt: skip
            if np.allclose(new, centres):
                break
            centres = new
        inertia = float(((x - centres[labels]) ** 2).sum())
        if best is None or inertia < best[2]:
            best = (labels, centres, inertia)
    return best


# --------------------------------------------------------------------------
# Instrument-response hypotheses from the Z/P phase
# --------------------------------------------------------------------------
def _hp1(f, fc):
    """First-order high-pass i(f/fc)/(1 + i f/fc) (DPG leak, Doran et al. 2019)."""
    x = 1j * np.asarray(f, float) / fc
    return x / (1 + x)


def _hp2(f, f0, damping):
    """Second-order high-pass of a feedback seismometer (velocity flat above f0)."""
    r = np.asarray(f, float) / f0
    return -(r**2) / (1 - r**2 + 2j * damping * r)


def residual_phase(phase_deg):
    """Phase relative to the nearest of 0 and 180 degrees, and that polarity."""
    phase = np.asarray(phase_deg, float)
    polarity = 180.0 if np.median(np.abs(phase)) > 90 else 0.0
    res = (phase - polarity + 180.0) % 360.0 - 180.0
    return res, polarity


def phase_models(f, phase_deg, *, seismometer_period=240.0, seismometer_damping=0.707,
                 dpg_corners_hz=np.geomspace(2e-4, 8e-3, 25), tolerance=1.5,
                 calibration_band=CALIBRATION_BAND_HZ, check_band=(0.03, 0.07)):  # fmt: skip
    """Fit instrument hypotheses to the residual Z/P phase of a stack.

    The Z/P transfer function of an elastic seafloor is real: its phase is
    0 or 180 degrees, depending on the pressure polarity convention. A
    residual that varies with frequency points to a response that differs
    from the one removed. Four hypotheses are fitted by least squares:

    * ``constant``: a frequency-independent offset (null model);
    * ``timing``: a clock offset dt between channels, phase = 360 f dt;
    * ``dpg_corner``: the DPG leak corner differs from the metadata one,
      phase = atan(f/f_true) - atan(f/f_meta) (first-order high-pass with
      a zero at DC; Doran et al. 2019); f_meta on a grid, f_true fitted;
    * ``seismometer_corner``: the seismometer corner (nominal
      ``seismometer_period`` and ``seismometer_damping``) is off by a
      factor, with damping 0.6, nominal or 0.8.

    The phase constrains a corner *difference*, not the amplitude error,
    which depends on where the corners sit. So every parameter set whose
    RMS misfit is within ``tolerance`` times the best of its family
    (or 0.5 degree, whichever is larger) is kept, and the range of the
    compliance bias they imply is returned. The bias is what remains after
    a Rayleigh-wave gain fitted in ``calibration_band`` absorbed the error
    there. ``rayleigh_slope`` is each model's ln(r(f2)/r(f1)) over
    ``check_band`` for r = measured/expected P/a, comparable with the
    calibration's saved ``median_log_ratio_slope``.

    Args:
        f (array): Frequencies [Hz] of the compliance band.
        phase_deg (array): Phase of S_pz at ``f`` [deg].

    Returns:
        dict: ``polarity`` (0 or 180), ``residual_deg``, ``freqs`` and per
        model ``params`` (best fit), ``rms_deg``, ``prediction_deg`` and,
        for response models, ``compliance_bias_pct`` (best fit) with
        ``bias_low_pct``/``bias_high_pct`` (envelope of acceptable fits per
        frequency) and ``rayleigh_slope`` (range over acceptable fits).
    """
    f = np.asarray(f, float)
    res, polarity = residual_phase(phase_deg)
    fcal = np.geomspace(*calibration_band, 9)
    out = {"polarity": polarity, "residual_deg": res, "freqs": f}

    def rms(pred):
        return np.sqrt(np.mean((res - pred) ** 2, axis=-1))

    c = float(np.mean(res))
    out["constant"] = dict(params={"offset_deg": c}, rms_deg=float(rms(np.full_like(f, c))),
                           prediction_deg=np.full_like(f, c))  # fmt: skip
    dt = float(np.dot(360 * f, res) / np.dot(360 * f, 360 * f))
    out["timing"] = dict(params={"dt_s": dt}, rms_deg=float(rms(360 * f * dt)),
                         prediction_deg=360 * f * dt)  # fmt: skip

    def summarize(name, fits, amplitude, params):
        # fits: (rms, key, prediction); amplitude(key, freqs): compliance factor.
        # r = measured/expected P/a scales as 1/amplitude in both families.
        best = min(fits, key=lambda x: x[0])
        limit = max(tolerance * best[0], best[0] + 0.5)
        ok = [x for x in fits if x[0] <= limit]
        biases = np.array(
            [amplitude(k, f) / np.mean(amplitude(k, fcal)) for _, k, _ in ok]
        )
        slopes = [float(np.log(amplitude(k, check_band[0]) / amplitude(k, check_band[1])))
                  for _, k, _ in ok]  # fmt: skip
        best_bias = amplitude(best[1], f) / np.mean(amplitude(best[1], fcal))
        out[name] = dict(params=params(best[1]), rms_deg=float(best[0]),
                         prediction_deg=best[2], compliance_bias_pct=100 * (best_bias - 1),
                         bias_low_pct=100 * (biases.min(axis=0) - 1),
                         bias_high_pct=100 * (biases.max(axis=0) - 1),
                         rayleigh_slope=(min(slopes), max(slopes)), n_acceptable=len(ok))  # fmt: skip

    # DPG: compliance ~ 1/|P_corrected|, P_corrected = P |H_t/H_m|.
    trial_t = np.geomspace(2e-5, 2e-2, 400)
    fits = []
    for fm in dpg_corners_hz:
        pred = np.degrees(
            np.arctan(f[None, :] / trial_t[:, None]) - np.arctan(f / fm)[None, :]
        )
        err = rms(pred)
        for i in range(len(trial_t)):
            fits.append((float(err[i]), (float(fm), float(trial_t[i])), pred[i]))

    def dpg_amp(key, freqs):
        fm, ft = key
        return np.abs(_hp1(freqs, fm)) / np.abs(_hp1(freqs, ft))

    summarize("dpg_corner", fits, dpg_amp,
              lambda k: {"f_meta_hz": k[0], "f_true_hz": k[1],
                         "tau_meta_s": 1 / (2 * np.pi * k[0]),
                         "tau_true_s": 1 / (2 * np.pi * k[1])})  # fmt: skip

    # Seismometer: compliance ~ |Z_corrected| = |Z| |H_t/H_m|; r = P/a ~ 1/|H_t/H_m|.
    f0 = 1.0 / seismometer_period
    fits = []
    for scale in np.geomspace(0.5, 2.0, 301):
        for damp in (seismometer_damping, 0.6, 0.8):
            ratio = _hp2(f, f0 * scale, damp) / _hp2(f, f0, seismometer_damping)
            pred = np.degrees(np.angle(ratio))
            fits.append((float(rms(pred)), (float(scale), damp), pred))

    def seis_amp(key, freqs):
        scale, damp = key
        return np.abs(
            _hp2(freqs, f0 * scale, damp) / _hp2(freqs, f0, seismometer_damping)
        )

    summarize("seismometer_corner", fits, seis_amp,
              lambda k: {"period_true_s": seismometer_period / k[0], "damping_true": k[1],
                         "period_meta_s": seismometer_period})  # fmt: skip
    return out


# --------------------------------------------------------------------------
# Network
# --------------------------------------------------------------------------
def common_mode(table, min_overlap=6):
    """Correlation of series shared by several stations.

    Args:
        table (dict): ``{station: {key: value}}``, e.g. monthly relative
            deviations of the band-mean compliance.
        min_overlap (int): Common keys needed for a pair.

    Returns:
        dict: ``stations``, ``keys``, ``matrix`` (stations x keys, NaN where
        missing), pairwise Spearman ``corr``, ``mean_corr`` per station
        (average with the others) and ``network_mean_corr``.
    """
    stations = sorted(table)
    keys = sorted({k for v in table.values() for k in v})
    m = np.full((len(stations), len(keys)), np.nan)
    for i, s in enumerate(stations):
        for j, k in enumerate(keys):
            if k in table[s]:
                m[i, j] = table[s][k]
    n = len(stations)
    corr = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            ok = np.isfinite(m[i]) & np.isfinite(m[j])
            if ok.sum() >= min_overlap:
                corr[i, j] = stats.spearmanr(m[i, ok], m[j, ok])[0]
    off = corr.copy()
    np.fill_diagonal(off, np.nan)
    with np.errstate(invalid="ignore"):
        mean_corr = np.nanmean(off, axis=1) if n > 1 else np.array([np.nan])
    return dict(stations=stations, keys=keys, matrix=m, corr=corr,
                mean_corr=mean_corr, network_mean_corr=float(np.nanmean(off)) if n > 1 else np.nan)  # fmt: skip


# --------------------------------------------------------------------------
# Selection bias of coherence-gated stacks
# --------------------------------------------------------------------------
def simulate_gate_bias(coherence, *, seastate_sigma=0.5, gate=0.9, n_windows=2000,
                       window=2625, nseg=1024, n_bins=9, exclusive=False, seed=0):  # fmt: skip
    """Monte Carlo bias of a stacked H1 when windows are gated on coherence.

    Pressure is white Gaussian noise whose amplitude varies between windows
    by a log-normal sea-state factor (``seastate_sigma`` in natural log);
    the vertical is the pressure plus constant incoherent noise, set so the
    coherence is ``coherence`` at the median sea state. Windows are kept
    when the median amplitude coherence over ``n_bins`` bins reaches
    ``gate`` (``exclusive``: per bin, ignoring that bin and its neighbours),
    then S_pz and S_pp are stacked as in ``stacked_compliance``. The true
    transfer function is 1.

    Returns:
        dict: ``acceptance`` and ``bias_pct`` (stacked |S_pz|/S_pp - 1, in %).
    """
    from scipy import signal

    rng = np.random.default_rng(seed)
    a = np.exp(seastate_sigma * rng.standard_normal(n_windows))[:, None]
    noise = math.sqrt(1 / coherence**2 - 1)
    p = a * rng.standard_normal((n_windows, window))
    z = p + noise * rng.standard_normal((n_windows, window))
    kw = dict(
        nperseg=nseg, noverlap=nseg // 2, window="hann", detrend="constant", axis=-1
    )
    _, zz = signal.welch(z, **kw)
    _, pp = signal.welch(p, **kw)
    _, pz = signal.csd(p, z, **kw)
    first = max(2, (zz.shape[1] - n_bins) // 2)
    idx = np.arange(first, first + n_bins)
    coh = np.abs(pz[:, idx]) / np.sqrt(zz[:, idx] * pp[:, idx])
    if exclusive:
        mask = np.zeros_like(coh)
        for j in range(n_bins):
            others = [k for k in range(n_bins) if abs(k - j) >= 2]
            mask[:, j] = np.median(coh[:, others], axis=1) >= gate
    else:
        mask = np.repeat((np.median(coh, axis=1) >= gate)[:, None], n_bins, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        h1 = np.abs((pz[:, idx] * mask).sum(0)) / (pp[:, idx] * mask).sum(0)
    return dict(
        acceptance=float(mask.mean()), bias_pct=float(100 * (np.nanmean(h1) - 1))
    )


def expected_gate_bias(
    acceptance, *, gate=0.9, seastate_sigma=0.5, seed=0, n_windows=1500
):
    """Selection bias of a coherence-gated stack, matched to an acceptance rate.

    Finds by bisection the median coherence for which
    :func:`simulate_gate_bias` accepts the observed fraction of windows,
    and returns the stacked-H1 bias there. ``seastate_sigma`` (log-normal
    spread of the window pressure amplitude) controls how much the gate
    selects on signal strength rather than on noise; the bias falls
    quickly as it grows.

    Returns:
        dict: ``coherence`` (matched), ``acceptance`` and ``bias_pct``.
    """
    lo, hi = 0.5, 0.995
    sim = None
    for _ in range(9):
        mid = 0.5 * (lo + hi)
        sim = simulate_gate_bias(mid, seastate_sigma=seastate_sigma, gate=gate,
                                 n_windows=n_windows, seed=seed)  # fmt: skip
        if sim["acceptance"] < acceptance:
            lo = mid
        else:
            hi = mid
    return dict(coherence=mid, **sim)


# --------------------------------------------------------------------------
# Scans
# --------------------------------------------------------------------------
def _iso(t):
    from datetime import datetime, timezone

    return datetime.fromtimestamp(float(t), tz=timezone.utc).strftime("%Y-%m-%dT%H:%M")


def _pct_from_log10(x):
    return 100.0 * (10.0 ** float(x) - 1.0)


def scan_phase(s):
    """Residual Z/P phase of a station's stack and the instrument models it fits."""
    sel = s.in_band
    if sel.sum() < 3:
        return [], None
    m = phase_models(s.freqs[sel], s.phase_deg[sel])
    res = m["residual_deg"]
    out = []
    if np.max(np.abs(res)) > 45:
        out.append(Anomaly(
            s.station, "instrument", "phase_inconsistent",
            f"Z/P phase {res[0]:+.0f}° at {1e3 * s.freqs[sel][0]:.1f} mHz, far from 0/180°: "
            "the Z or P response (or polarity) is not the one removed; |compliance| is unreliable",
            10.0, float("nan"), evidence={"residual_deg": res.round(1).tolist(), "label": s.label}))  # fmt: skip
        return out, m
    n_avg = float(s.summary.get("results", {}).get(s.label, {}).get("n_averages", 1e4))
    coh = np.clip(s.coherence[sel], 1e-3, 0.999)
    sigma_deg = np.degrees(np.sqrt(1 - coh**2) / (coh * np.sqrt(2 * n_avg)))
    sigma_deg = np.maximum(sigma_deg, 0.5)  # floor for systematics
    z = float(np.max(np.abs(res) / sigma_deg))
    best = (
        "dpg_corner"
        if m["dpg_corner"]["rms_deg"] <= m["seismometer_corner"]["rms_deg"]
        else "seismometer_corner"
    )
    lo = min(
        m["dpg_corner"]["bias_low_pct"][0], m["seismometer_corner"]["bias_low_pct"][0]
    )
    hi = max(
        m["dpg_corner"]["bias_high_pct"][0], m["seismometer_corner"]["bias_high_pct"][0]
    )
    if np.max(np.abs(res)) >= 3.0:
        out.append(Anomaly(
            s.station, "instrument", "phase_residual",
            f"Z/P phase residual {res[0]:+.1f}° at {1e3 * s.freqs[sel][0]:.1f} mHz → {res[-1]:+.1f}° at "
            f"{1e3 * s.freqs[sel][-1]:.1f} mHz; timing model rms {m['timing']['rms_deg']:.1f}° vs response "
            f"models {m['dpg_corner']['rms_deg']:.1f}/{m['seismometer_corner']['rms_deg']:.1f}°; implied compliance "
            f"bias at band bottom {lo:+.0f} to {hi:+.0f} %",
            min(z, 10.0), 0.5 * (lo + hi),
            evidence={"residual_deg": res.round(2).tolist(), "best_model": best,
                      "dpg": m["dpg_corner"]["params"], "seismometer": m["seismometer_corner"]["params"],
                      "rms_constant": m["constant"]["rms_deg"], "rms_timing": m["timing"]["rms_deg"],
                      "rms_dpg": m["dpg_corner"]["rms_deg"], "rms_seismometer": m["seismometer_corner"]["rms_deg"],
                      "bias_range_band_bottom_pct": [float(lo), float(hi)]}))  # fmt: skip
    return out, m


def stack_shift(t, log_tf, p_db, cps, end=None):
    """Change of a stacked |Z/P| caused by windows after the first changepoint.

    The stack of ``stacked_compliance`` weights each window by its pressure
    power, so the stacked level is approximated by sum(w 10**log_tf)/sum(w)
    with w = 10**(p_db/10). Returns the stack of all windows before ``end``
    relative to the stack of the windows before the first changepoint [%].
    """
    t, r, p_db = (np.asarray(v, float) for v in (t, log_tf, p_db))
    if not cps:
        return 0.0
    first = min(c.time for c in cps)
    w, tf = 10 ** (p_db / 10), 10**r
    use = t < (np.inf if end is None else end)
    ref = t < first
    if not ref.any() or not use.any():
        return float("nan")
    return float(100 * (np.sum(w[use] * tf[use]) / np.sum(w[use])
                        / (np.sum(w[ref] * tf[ref]) / np.sum(w[ref])) - 1))  # fmt: skip


def scan_levels(s, *, cp_threshold=5.0, p_periodic=0.01):
    """Changepoints, periodicity, level dependence and clusters of per-window levels."""
    out, info = [], {}
    if s.window_levels is None or len(s.window_levels) < 200:
        return out, info
    t, r, p_db, a_db = np.asarray(s.window_levels, float).T
    names = {"transfer function": r, "pressure level": p_db, "vertical level": a_db}
    cps = {
        k: changepoints(t, v, min_size=100, threshold=cp_threshold)
        for k, v in names.items()
    }
    info["changepoints"] = cps
    for c in cps["transfer function"]:
        near = [f"{k} {c2.shift:+.1f} dB" for k in ("pressure level", "vertical level")
                for c2 in cps[k] if abs(c2.time - c.time) < 5 * 86400]  # fmt: skip
        effect = _pct_from_log10(c.shift)
        out.append(Anomaly(
            s.station, "temporal", "changepoint",
            f"|Z/P| level steps by {effect:+.1f} % at {_iso(c.time)} ({c.n_before} windows before, "
            f"{c.n_after} after)" + (f"; with {', '.join(near)}" if near else "; no concurrent P/Z level step"),
            c.significance, effect, _iso(c.time),
            evidence={"shift_log10": c.shift, "n_before": c.n_before, "n_after": c.n_after,
                      "concurrent": near}))  # fmt: skip
    shift = stack_shift(t, r, p_db, cps["transfer function"])
    info["stack_shift_pct"] = shift
    if cps["transfer function"] and abs(shift) >= 0.3:
        first = min(c.time for c in cps["transfer function"])
        out.append(Anomaly(
            s.station, "temporal", "epoch_mixing",
            f"the full-period stack is {shift:+.2f} % off the stack before the first |Z/P| step "
            f"({_iso(first)[:10]}; {100 * np.mean(t >= first):.0f} % of windows after it): one gain "
            "for the whole record mixes epochs",
            max(c.significance for c in cps["transfer function"]), shift, _iso(first),
            {"windows_after_pct": float(100 * np.mean(t >= first))}))  # fmt: skip
    resid = segment_residuals(t, r, cps["transfer function"])
    info["residual"] = resid
    zr = robust_z(resid)
    frac = float(np.mean(np.abs(zr) > 5))
    info["outlier_fraction"] = frac
    if frac > 0.01:
        out.append(Anomaly(s.station, "noise", "outliers",
                           f"{100 * frac:.1f} % of windows have |robust z| > 5 in |Z/P| level",
                           min(10.0, 2 + 100 * frac), float("nan"),
                           evidence={"fraction": frac}))  # fmt: skip
    per = periodicity(t, resid)
    info["periodicity"] = per
    for name, v in per.items():
        if v["p_value"] <= p_periodic and v["p_tail"] <= p_periodic:
            amp = 100 * np.log(10) * v["amplitude"]
            out.append(Anomaly(
                s.station, "temporal", "tidal_transfer_function",
                f"|Z/P| level modulated at {name} ({v['period_h']:.2f} h): ±{amp:.2f} % "
                f"(tail p = {v['p_tail']:.1e}, controls' median ±{100 * np.log(10) * v['null_median']:.2f} %)",
                p_to_sigma(v["p_tail"]), amp, evidence=v))  # fmt: skip
    dep = level_dependence(p_db, resid)
    info["level_dependence"] = dep
    if np.isfinite(dep["slope"]):
        effect = 100 * (10 ** (10 * dep["slope"]) - 1)  # per 10 dB of pressure
        if abs(effect) >= 1.0 and dep["p_value"] < 1e-3:
            out.append(Anomaly(
                s.station, "noise", "pressure_level_dependence",
                f"|Z/P| level changes by {effect:+.1f} % per 10 dB of pressure among strong windows "
                f"(Spearman {dep['spearman']:+.2f})",
                p_to_sigma(dep["p_value"]), effect, evidence=dep))  # fmt: skip
    labels, centres, _ = kmeans(np.column_stack([resid, p_db, a_db]), k=3)
    info["clusters"] = labels
    for j in range(3):
        sel = labels == j
        if sel.mean() < 0.03 or sel.all():
            continue
        shift = _pct_from_log10(np.median(resid[sel]))
        # A regime, not the tail of one distribution: its windows are
        # concentrated in time compared with the others.
        ks = stats.ks_2samp(t[sel], t[~sel])
        if abs(shift) >= 3.0 and ks.statistic >= 0.3:
            out.append(Anomaly(
                s.station, "noise", "window_regime",
                f"{sel.sum()} windows ({100 * sel.mean():.0f} %) form a regime in time "
                f"(KS {ks.statistic:.2f}, mostly {_iso(np.quantile(t[sel], 0.1))[:10]}–"
                f"{_iso(np.quantile(t[sel], 0.9))[:10]}) with |Z/P| {shift:+.1f} % "
                f"(P {np.median(p_db[sel]) - np.median(p_db):+.1f} dB, Z {np.median(a_db[sel]) - np.median(a_db):+.1f} dB)",
                p_to_sigma(ks.pvalue), shift, evidence={"fraction": float(sel.mean()),
                                                         "ks": float(ks.statistic)}))  # fmt: skip
    return out, info


def scan_acceptance(s, *, p_periodic=0.01):
    """Tidal clustering of accepted windows and coherence-gate risk."""
    out, info = [], {}
    t = np.asarray(s.window_starts, float)
    if len(t) < 200:
        return out, info
    half = 0.5 * (s.window_s or 0.0)
    ray = rayleigh_periodicity(t + half)
    info["rayleigh"] = ray
    for name, v in ray.items():
        if v["p_value"] <= p_periodic and v["p_tail"] <= p_periodic:
            out.append(Anomaly(
                s.station, "temporal", "tidal_acceptance",
                f"accepted windows cluster at one phase of {name} ({v['period_h']:.2f} h): modulation "
                f"{200 * v['R']:.1f} % (controls' median {200 * v['null_median_R']:.1f} %, tail p = {v['p_tail']:.1e})",
                p_to_sigma(v["p_tail"]), 200 * v["R"], evidence=v))  # fmt: skip
    acc = s.acceptance
    info["acceptance"] = acc
    return out, info


def scan_monthly(s, *, min_windows=30):
    """Excursions and drift of monthly compliance stacks."""
    out, info = [], {}
    keys, counts, values = s.monthly_matrix(min_windows)
    if len(keys) < 4:
        return out, info
    full = s.compliance[s.in_band]
    dev = 100 * (values / full[None, :] - 1)  # %, months x freqs
    mean_dev = dev.mean(axis=1)
    info.update(keys=keys, counts=counts, deviation_pct=dev, band_mean_pct=mean_dev)
    scale = max(mad_sigma(mean_dev), 0.5)
    for k, d, n in zip(keys, mean_dev, counts):
        z = abs(d - np.median(mean_dev)) / scale
        if z >= 3 and abs(d) >= 3:
            out.append(Anomaly(s.station, "temporal", "monthly_excursion",
                               f"{k}: monthly stack {d:+.1f} % from the full stack ({n} windows)",
                               min(z, 10.0), float(d), k, {"n_windows": int(n)}))  # fmt: skip
    rho, p = stats.spearmanr(np.arange(len(keys)), mean_dev)
    info["trend"] = (float(rho), float(p))
    if p < 0.01:
        slope = np.polyfit(np.arange(len(keys)), mean_dev, 1)[0]
        out.append(Anomaly(s.station, "temporal", "monthly_drift",
                           f"monthly compliance drifts {slope:+.2f} %/month (Spearman {rho:+.2f})",
                           p_to_sigma(p), float(slope * len(keys)),
                           evidence={"spearman": float(rho), "p_value": float(p)}))  # fmt: skip
    # frequency dependence of the scatter: low/high band ratio
    spread = dev.std(axis=0)
    info["spread_by_freq_pct"] = spread
    return out, info


def scan_spectra(s, *, notch_depth=0.03, line_db=8.0):
    """Coherence notches in the band and narrow lines in the stacked vertical PSD."""
    out, info = [], {}
    f, coh = s.all_freqs, s.mean_coherence_all
    lo, hi = s.band
    idx = np.where((f >= lo) & (f <= hi))[0]
    for i in idx[1:-1]:
        depth = 0.5 * (coh[i - 1] + coh[i + 1]) - coh[i]
        if depth >= notch_depth:
            out.append(Anomaly(s.station, "spectral", "coherence_notch",
                               f"coherence notch of {depth:.2f} at {1e3 * f[i]:.1f} mHz",
                               min(10.0, depth / 0.01), float("nan"),
                               evidence={"freq_hz": float(f[i]), "depth": float(depth)}))  # fmt: skip
    with np.errstate(divide="ignore"):
        acc = 10 * np.log10(s.S_zz * (2 * np.pi * f) ** 4)
    sel = (f > 0.02) & (f < 0.55) & np.isfinite(acc)
    if sel.sum() > 20:
        from scipy.ndimage import median_filter

        base = median_filter(acc[sel], size=15, mode="nearest")
        excess = acc[sel] - base
        peaks = np.where((excess > line_db) & (excess == median_filter(excess, 3)))[0]
        info["lines_hz"] = f[sel][peaks].tolist()
        for i in peaks:
            out.append(Anomaly(s.station, "spectral", "vertical_line",
                               f"narrow line {excess[i]:.0f} dB above background at {f[sel][i]:.3f} Hz",
                               min(10.0, excess[i] / 2), float("nan"),
                               evidence={"freq_hz": float(f[sel][i])}))  # fmt: skip
    return out, info


def scan_inversion(s, variant="final"):
    """Sampled-gain factor and fit quality of a stage E record."""
    out = []
    rec = s.inversion.get(variant)
    if not rec:
        return out
    par = rec.get("parameters", {})
    gains = {k: v for k, v in par.items() if k.startswith("gain")}
    for name, v in gains.items():
        med, lo, hi = v["median"], v["p05"], v["p95"]
        if not lo <= 1 <= hi and abs(med - 1) >= 0.02:
            out.append(Anomaly(s.station, "instrument", "gain_pull",
                               f"inversion pulls the {name} factor to {med:.3f} (90 % {lo:.3f}–{hi:.3f}), "
                               "outside 1: data and calibration disagree",
                               min(10.0, abs(med - 1) / max((hi - lo) / 3.29, 1e-3)),
                               100 * (med - 1), evidence=v))  # fmt: skip
    chi2, n = rec.get("chi2", {}).get("median"), rec.get("n_data")
    if chi2 and n and chi2 / n > 2:
        out.append(Anomaly(s.station, "instrument", "poor_fit",
                           f"χ²/N = {chi2:.1f}/{n}: the 1-D model does not fit the compliance",
                           min(10.0, (chi2 - n) / np.sqrt(2 * n)), float("nan"),
                           evidence={"chi2": chi2, "n": n}))  # fmt: skip
    return out


def scan_station(s):
    """All single-station detectors; returns (anomalies, info dict)."""
    found, info = [], {}
    for name, scan in (("phase", scan_phase), ("levels", scan_levels),
                       ("acceptance", scan_acceptance), ("monthly", scan_monthly),
                       ("spectra", scan_spectra)):  # fmt: skip
        a, i = scan(s)
        found += a
        info[name] = i
    info["phase_models"] = info.pop("phase")
    found += scan_inversion(s)
    return found, info


def scan_network(stations, infos, *, n_permutations=500, seed=0, name="network"):
    """Cross-station detectors on the results of :func:`scan_station`.

    Args:
        stations (dict): ``{code: StationProducts}``.
        infos (dict): ``{code: info}`` from :func:`scan_station`.

    Returns:
        tuple: (anomalies, info).
    """
    out, info = [], {}
    table = {}
    for code, inf in infos.items():
        m = inf.get("monthly", {})
        if "keys" in m:
            table[code] = dict(zip(m["keys"], m["band_mean_pct"]))
    if len(table) >= 3:
        cm = common_mode(table)
        rng = np.random.default_rng(seed)
        null = []
        for _ in range(n_permutations):
            shuffled = {k: dict(zip(v.keys(), rng.permutation(list(v.values()))))
                        for k, v in table.items()}  # fmt: skip
            null.append(common_mode(shuffled)["network_mean_corr"])
        p = (1 + np.sum(np.array(null) >= cm["network_mean_corr"])) / (
            1 + n_permutations
        )
        cm["p_value"] = float(p)
        info["common_mode"] = cm
        if p <= 0.01:
            out.append(Anomaly(name, "network", "common_mode_monthly",
                               f"monthly compliance deviations are correlated across {len(table)} stations "
                               f"(mean Spearman {cm['network_mean_corr']:+.2f}, permutation p = {p:.3f})",
                               p_to_sigma(p), float("nan"), evidence={"p_value": float(p)}))  # fmt: skip
    steps = sorted((c.time, code, c.shift) for code, inf in infos.items()
                   for c in inf.get("levels", {}).get("changepoints", {}).get("transfer function", []))  # fmt: skip
    groups, current = [], []
    for item in steps:
        if current and item[0] - current[0][0] > 15 * 86400:
            groups.append(current)
            current = []
        current.append(item)
    if current:
        groups.append(current)
    info["synchronous_steps"] = groups
    for g in groups:
        codes = sorted({c for _, c, _ in g})
        if len(codes) >= 3:
            effects = [_pct_from_log10(x) for _, _, x in g]
            out.append(Anomaly(
                name, "network", "synchronous_steps",
                f"|Z/P| steps at {len(codes)} stations within {(g[-1][0] - g[0][0]) / 86400:.0f} days "
                f"({_iso(g[0][0])[:10]} – {_iso(g[-1][0])[:10]}: "
                + ", ".join(f"{c} {e:+.1f} %" for (_, c, _), e in zip(g, effects)) + ")",
                p_to_sigma(stats.binom.sf(len(codes) - 1, len(infos), 30 / 400)), float(np.median(effects)),
                _iso(g[0][0]), {"stations": codes}))  # fmt: skip
    bands = {"compliance band": None, "microseism band": (0.1, 0.3)}
    for kind, spectrum in (("vertical", "S_zz"), ("pressure", "S_pp")):
        for label, band in bands.items():
            levels = {}
            for code, st in stations.items():
                f = st.all_freqs
                lo, hi = st.band if band is None else band
                sel = (f >= lo) & (f <= hi)
                psd = getattr(st, spectrum)[sel]
                if kind == "vertical":
                    psd = psd * (2 * np.pi * f[sel]) ** 4
                levels[code] = 10 * np.log10(np.mean(psd))
            values = np.array(list(levels.values()))
            scale = max(mad_sigma(values), 1.0)
            for code, v in levels.items():
                z = (v - np.median(values)) / scale
                if abs(z) >= 4 and abs(v - np.median(values)) >= 6:
                    out.append(Anomaly(code, "instrument", f"{kind}_level_outlier",
                                       f"stacked {kind} PSD in the {label} {v - np.median(values):+.0f} dB from the "
                                       f"{name} median: gain, units or a failed channel",
                                       min(abs(z), 10.0), float("nan"),
                                       evidence={"level_db": float(v), "network_median_db": float(np.median(values))}))  # fmt: skip
            info[f"{kind}_{label.split()[0]}_levels_db"] = levels
    signs = [np.sign(inf["phase_models"]["residual_deg"][0]) for inf in infos.values()
             if inf.get("phase_models") is not None
             and np.max(np.abs(inf["phase_models"]["residual_deg"])) < 45]  # fmt: skip
    if len(signs) >= 5:
        k = int(np.sum(np.array(signs) > 0))
        p = stats.binomtest(k, len(signs), 0.5).pvalue
        info["phase_sign"] = (k, len(signs), float(p))
        if p <= 0.01:
            out.append(Anomaly(name, "network", "systematic_phase",
                               f"{k} of {len(signs)} stations have a positive Z/P phase residual at the band "
                               f"bottom (binomial p = {p:.1e}): a shared instrument-response mismatch",
                               p_to_sigma(p), float("nan"), evidence={"k": k, "n": len(signs)}))  # fmt: skip
    pairs = [(s.sediment_thickness(), infos[c]["phase_models"]["residual_deg"][0])
             for c, s in stations.items()
             if s.sediment_thickness() is not None and infos[c].get("phase_models") is not None
             and np.max(np.abs(infos[c]["phase_models"]["residual_deg"])) < 45]  # fmt: skip
    if len(pairs) >= 6:
        x, y = np.array(pairs).T
        rho, p = stats.spearmanr(x, y)
        info["phase_vs_sediment"] = (float(rho), float(p), len(pairs))
    return out, info
