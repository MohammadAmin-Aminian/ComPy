"""Diagnostic and insight figures for compliance QC (ComPy 3 audit).

Every function takes plain arrays or ``compy.results.StationProducts`` and
returns a Matplotlib Figure; nothing here needs ObsPy. What each figure is
for and which anomaly it exposes:

=================================  ==============================================
Figure                             Reveals
=================================  ==============================================
:func:`plot_spectrogram`           storms, data gaps, glitch lines, IG band power
:func:`plot_coherogram`            when Z-P coherence (compliance) or Z-H (tilt)
                                   dominates; coherence drops after events
:func:`plot_compliance_confidence` the error budget: formal B&P vs month-to-month
                                   scatter vs model floor
:func:`plot_compliance_stability`  steps, drifts and seasonal changes of compliance
:func:`plot_tilt`                  tilt orientation/strength vs time and tide
:func:`plot_noise_models`          PSD levels against Peterson / Webb references
:func:`plot_multistation`          network panels of compliance and coherence
:func:`plot_phase_diagnostic`      residual Z/P phase and the response models
                                   that explain it, with implied compliance bias
:func:`plot_window_levels`         per-window |Z/P|, P and Z levels with steps
:func:`plot_periodicity`           tidal/diurnal modulation by station
:func:`plot_acceptance_phase`      accepted windows folded on tidal periods
:func:`plot_common_mode`           month x station deviations and correlation
:func:`plot_anomaly_ranking`       the ranked anomaly table as a chart
=================================  ==============================================

Colours follow one rule per job: a fixed categorical order for identity,
one blue ramp for magnitude, blue <-> red through grey for signed deviations.
"""

from __future__ import annotations

from datetime import datetime, timezone

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

#: Categorical hues in fixed order (identity); never cycled past eight.
CATEGORICAL = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#6250d6", "#e34948")  # fmt: skip
INK, MUTED, GRID = "#1f1f1d", "#6b6a63", "#e1e0d9"
#: One-hue magnitude ramp (light to dark blue).
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "compy_blues", ["#f4f8fd", "#b7d3f6", "#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"]
)
#: Signed deviation: blue (negative) - grey (zero) - red (positive).
DIVERGING = LinearSegmentedColormap.from_list(
    "compy_diverging", ["#104281", "#6da7ec", "#f0efec", "#ef8f8e", "#a32b2a"]
)
CATEGORY_COLOR = {"instrument": CATEGORICAL[0], "temporal": CATEGORICAL[1],
                  "noise": CATEGORICAL[2], "spectral": CATEGORICAL[3],
                  "network": CATEGORICAL[6], "selection": CATEGORICAL[4]}  # fmt: skip


def _style(ax):
    ax.grid(True, color=GRID, lw=0.6, ls="-")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=INK, labelsize=8)
    return ax


def _date_axis(ax):
    locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))


def _dates(t):
    return np.array(
        [datetime.fromtimestamp(float(x), tz=timezone.utc) for x in np.atleast_1d(t)]
    )


def _band(ax, band, vertical=True):
    if band is None:
        return
    for edge in band:
        (ax.axvline if vertical else ax.axhline)(edge, color=INK, lw=0.8, ls="--")


# --------------------------------------------------------------------------
# Time-frequency
# --------------------------------------------------------------------------
def plot_spectrogram(times, f, psd, *, band=None, fmin=1e-3, title=None, units=None):
    """Spectrograms of several channels over a deployment.

    Purpose: see storms (broadband pressure), infragravity power in the
    compliance band, glitches (horizontal lines), gaps and sensor changes.

    Args:
        times (array): Window centres [POSIX s].
        f (array): Frequencies [Hz].
        psd (dict): ``{label: (n_windows, n_freqs) PSD}`` in linear units.
        band (tuple): Compliance band, drawn as dashed lines.
        units (dict): ``{label: unit string}`` for colour bars.
    """
    n = len(psd)
    fig, axes = plt.subplots(n, 1, figsize=(10, 2.2 * n + 0.6), sharex=True,
                             constrained_layout=True, squeeze=False)  # fmt: skip
    sel = f >= fmin
    x = _dates(times)
    for ax, (label, values) in zip(axes[:, 0], psd.items()):
        with np.errstate(divide="ignore"):
            db = 10 * np.log10(np.asarray(values)[:, sel].T)
        lo, hi = np.nanpercentile(db, [2, 98])
        mesh = ax.pcolormesh(
            x, f[sel], db, cmap=SEQUENTIAL, vmin=lo, vmax=hi, shading="auto"
        )
        ax.set_yscale("log")
        ax.set_ylabel(f"{label}\nFrequency [Hz]", fontsize=8)
        _band(ax, band, vertical=False)
        fig.colorbar(mesh, ax=ax, label=(units or {}).get(label, "dB"), pad=0.01)
    _date_axis(axes[-1, 0])
    if title:
        axes[0, 0].set_title(title, fontsize=10)
    return fig


def plot_coherogram(
    times, f, coherences, *, band=None, threshold=None, fmin=1e-3, title=None
):
    """Coherence of channel pairs against time and frequency.

    Purpose: Z-P coherence marks when compliance dominates the vertical;
    Z-H1/Z-H2 coherence marks tilt noise; P-H coherence marks horizontals
    that carry pressure-related signal (which coherent-noise removal would
    then also remove from Z). A drop of Z-P coherence after an event or
    a re-levelling shows up as a vertical stripe.
    """
    n = len(coherences)
    fig, axes = plt.subplots(n, 1, figsize=(10, 2.0 * n + 0.6), sharex=True,
                             constrained_layout=True, squeeze=False)  # fmt: skip
    sel = f >= fmin
    x = _dates(times)
    for ax, (label, values) in zip(axes[:, 0], coherences.items()):
        mesh = ax.pcolormesh(x, f[sel], np.asarray(values)[:, sel].T, cmap=SEQUENTIAL,
                             vmin=0, vmax=1, shading="auto")  # fmt: skip
        ax.set_yscale("log")
        ax.set_ylabel(f"{label}\nFrequency [Hz]", fontsize=8)
        _band(ax, band, vertical=False)
        cb = fig.colorbar(mesh, ax=ax, pad=0.01, label="Coherence")
        if threshold is not None:
            cb.ax.axhline(threshold, color=INK, lw=1)
    _date_axis(axes[-1, 0])
    if title:
        axes[0, 0].set_title(title, fontsize=10)
    return fig


# --------------------------------------------------------------------------
# Compliance error budget and stability
# --------------------------------------------------------------------------
def error_budget(s, *, floor_rel=0.01, shrinkage=0.5, min_windows=30):
    """One-sigma components of a station's compliance error [1/Pa].

    ``random`` is the Bendat & Piersol error, ``months`` the standard error
    of the monthly stacks (their covariance divided by their number, shrunk
    to its diagonal by ``shrinkage`` as in ``ComplianceData.build``),
    ``floor`` the relative model floor and ``total`` their quadrature sum.
    """
    sel = s.in_band
    eta = s.compliance[sel]
    out = {"freqs": s.freqs[sel], "eta": eta, "random": s.uncertainty[sel],
           "floor": floor_rel * eta}  # fmt: skip
    _, _, months = s.monthly_matrix(min_windows)
    if len(months) >= 3:
        c = np.cov(months, rowvar=False) / len(months)
        c = (1 - shrinkage) * c + shrinkage * np.diag(np.diag(c))
        out["months"] = np.sqrt(np.diag(c))
        out["month_curves"] = months
    else:
        out["months"] = np.zeros_like(eta)
    out["total"] = np.sqrt(out["random"] ** 2 + out["months"] ** 2 + out["floor"] ** 2)
    return out


def plot_compliance_confidence(s, *, floor_rel=0.01, title=None):
    """Compliance with nested confidence bands and its error budget.

    Purpose: shows how small the formal (Bendat & Piersol) error is next to
    the month-to-month scatter and the 1 % model floor that the inversion
    adds, and whether monthly stacks stay inside the total band.
    """
    b = error_budget(s, floor_rel=floor_rel)
    f, eta = b["freqs"], b["eta"]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]},
                                  constrained_layout=True)  # fmt: skip
    for curve in b.get("month_curves", []):
        ax.plot(f, curve, color=MUTED, lw=0.6, alpha=0.5)
    ax.fill_between(f, eta - 2 * b["total"], eta + 2 * b["total"], color="#cde2fb",
                    label="±2σ total", lw=0)  # fmt: skip
    ax.fill_between(f, eta - b["total"], eta + b["total"], color="#86b6ef",
                    label="±1σ total (random + months + floor)", lw=0)  # fmt: skip
    ax.fill_between(f, eta - b["random"], eta + b["random"], color="#1c5cab",
                    label="±1σ Bendat & Piersol", lw=0)  # fmt: skip
    ax.plot(f, eta, color=INK, lw=1.6, marker="o", ms=3, label="stacked compliance")
    ax.set_ylabel("Normalized compliance [1/Pa]", fontsize=9)
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    _style(ax)
    for key, color, label in (("random", CATEGORICAL[0], "Bendat & Piersol"),
                              ("months", CATEGORICAL[1], "monthly SEM"),
                              ("floor", CATEGORICAL[2], f"floor {100 * floor_rel:.0f} %"),
                              ("total", INK, "total")):  # fmt: skip
        ax2.plot(f, 100 * b[key] / eta, color=color, lw=1.6 if key == "total" else 1.2,
                 marker="o", ms=3, label=label)  # fmt: skip
    ax2.set_xlabel("Frequency [Hz]", fontsize=9)
    ax2.set_ylabel("1σ [% of η]", fontsize=9)
    ax2.legend(fontsize=7, ncol=4, frameon=False, loc="upper left")
    _style(ax2)
    ax.set_title(
        title or f"{s.station}: compliance and error budget ({s.label} gain)",
        fontsize=10,
    )
    return fig


def plot_compliance_stability(s, *, min_windows=30, steps=None, title=None):
    """Monthly compliance relative to the full stack.

    Purpose: steps (sensor sensitivity changes), drifts and seasonal
    behaviour; frequency-dependent changes (a low-frequency-only change
    points to tilt or the DPG leak, a flat one to gain).

    Args:
        s (StationProducts): Station products.
        steps (list of float): POSIX times of detected steps, drawn as lines.
    """
    keys, counts, months = s.monthly_matrix(min_windows)
    if not keys:
        raise ValueError("no monthly stacks")
    f = s.freqs[s.in_band]
    dev = 100 * (months / s.compliance[s.in_band][None, :] - 1)
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 6), constrained_layout=True,
                                  gridspec_kw={"height_ratios": [1.3, 1]})  # fmt: skip
    lim = max(2.0, np.nanpercentile(np.abs(dev), 98))
    mesh = ax.pcolormesh(np.arange(len(keys) + 1), np.append(f, f[-1] + np.diff(f)[-1]) - np.diff(f)[0] / 2,
                         dev.T, cmap=DIVERGING, norm=TwoSlopeNorm(0, -lim, lim))  # fmt: skip
    ax.set_xticks(np.arange(len(keys)) + 0.5, keys, rotation=60, fontsize=7)
    ax.set_ylabel("Frequency [Hz]", fontsize=9)
    fig.colorbar(mesh, ax=ax, label="Month − full stack [%]", pad=0.01)
    x = np.arange(len(keys)) + 0.5
    sem = 100 * np.array([np.nanstd(d / 100) for d in dev]) / np.sqrt(dev.shape[1])
    ax2.axhline(0, color=MUTED, lw=0.8)
    ax2.errorbar(x, dev.mean(axis=1), yerr=sem, color=CATEGORICAL[0], marker="o", ms=4,
                 lw=1.2, capsize=2, label="band mean")  # fmt: skip
    ax2.plot(
        x,
        dev[:, 0],
        color=CATEGORICAL[1],
        lw=1,
        marker="s",
        ms=3,
        label=f"{1e3 * f[0]:.1f} mHz",
    )
    ax2.plot(
        x,
        dev[:, -1],
        color=CATEGORICAL[2],
        lw=1,
        marker="^",
        ms=3,
        label=f"{1e3 * f[-1]:.1f} mHz",
    )
    for k, n in zip(x, counts):
        ax2.annotate(str(n), (k, 0), xytext=(0, -12), textcoords="offset points",
                     ha="center", fontsize=6, color=MUTED)  # fmt: skip
    if steps:
        from compy.results import month_start

        starts = np.array([month_start(k) for k in keys])
        for t in steps:
            pos = np.interp(t, starts, np.arange(len(keys)))
            for a in (ax, ax2):
                a.axvline(pos, color=CATEGORICAL[7], lw=1.2, ls=":")
    ax2.set_xticks(x, keys, rotation=60, fontsize=7)
    ax2.set_xlim(0, len(keys))
    ax2.set_ylabel("Month − full stack [%]", fontsize=9)
    ax2.legend(fontsize=7, frameon=False, ncol=3)
    _style(ax2)
    ax.set_title(
        title or f"{s.station}: monthly stacks (window counts under the dots)",
        fontsize=10,
    )
    return fig


# --------------------------------------------------------------------------
# Tilt
# --------------------------------------------------------------------------
def plot_tilt(
    times,
    azimuth,
    angle,
    *,
    strength=None,
    period_h=12.4206012,
    angle_label="Tilt angle [°]",
    strength_label="strength",
    title=None,
):
    """Tilt orientation and strength against time and tidal phase.

    Purpose: a stable azimuth means a fixed instrument tilt; azimuth that
    flips with the tide means current-driven tilt (Corela et al. 2023);
    jumps in angle mean re-levelling or sinking. The right panel folds
    the azimuth on ``period_h`` (M2 by default).

    Args:
        times (array): POSIX times [s].
        azimuth (array): Tilt azimuth [deg].
        angle (array): Tilt angle [deg] (or any tilt-strength proxy).
        strength (array): Colour (e.g. variance reduction or coherence).
    """
    times, az, ang = (np.asarray(v, float) for v in (times, azimuth, angle))
    fig = plt.figure(figsize=(11, 5.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 3)
    ax1, ax2 = fig.add_subplot(grid[0, :2]), fig.add_subplot(grid[1, :2])
    ax3 = fig.add_subplot(grid[:, 2])
    x = _dates(times)
    c = strength if strength is not None else np.full_like(az, 0.5)
    sc = ax1.scatter(x, az % 360, c=c, s=3, cmap=SEQUENTIAL, rasterized=True)
    ax1.set_ylabel("Azimuth [°]", fontsize=9)
    ax1.set_ylim(0, 360)
    ax2.scatter(x, ang, c=c, s=3, cmap=SEQUENTIAL, rasterized=True)
    ax2.set_ylabel(angle_label, fontsize=9)
    phase = (times / 3600.0 / period_h) % 1.0
    h = ax3.hist2d(phase, az % 360, bins=(24, 36), cmap=SEQUENTIAL)
    ax3.set_xlabel(f"Phase of {period_h:.2f} h cycle", fontsize=9)
    ax3.set_ylabel("Azimuth [°]", fontsize=9)
    fig.colorbar(h[3], ax=ax3, label="Hours")
    if strength is not None:
        fig.colorbar(sc, ax=[ax1, ax2], label=strength_label, pad=0.01)
    for a in (ax1, ax2, ax3):
        _style(a)
    _date_axis(ax1)
    _date_axis(ax2)
    ax1.set_title(title or "Tilt orientation and strength", fontsize=10)
    return fig


# --------------------------------------------------------------------------
# Noise models and network panels
# --------------------------------------------------------------------------
def plot_noise_models(stations, *, kind="vertical", highlight=(), title=None):
    """Stacked PSDs of accepted windows against reference noise levels.

    Purpose: whether a station's vertical (as acceleration) sits between the
    Peterson NLNM/NHNM, or its pressure inside Webb's (1998) seafloor
    bounds; outliers point to gain or response errors. The stacks hold
    accepted windows only (high Z-P coherence), so they are not noise
    surveys of the whole record.
    """
    from compy import noise_models as nm

    fig, ax = plt.subplots(figsize=(7.5, 5), constrained_layout=True)
    labelled = False
    for i, (code, s) in enumerate(sorted(stations.items())):
        f = s.all_freqs[1:]
        if kind == "vertical":
            y = nm.displacement_to_acceleration_db(f, s.S_zz[1:])
        else:
            with np.errstate(divide="ignore"):
                y = 10 * np.log10(s.S_pp[1:])
        if code in highlight:
            k = list(highlight).index(code)
            ax.semilogx(f, y, color=CATEGORICAL[k % 8], lw=1.6, label=code)
        else:
            ax.semilogx(f, y, color=MUTED, lw=0.6, alpha=0.6,
                        label=None if labelled else "other stations")  # fmt: skip
            labelled = True
    fx = np.geomspace(1e-3, 0.6, 200)
    if kind == "vertical":
        ax.semilogx(
            fx, nm.nlnm(fx), color=INK, lw=1.2, ls="--", label="Peterson NLNM / NHNM"
        )
        ax.semilogx(fx, nm.nhnm(fx), color=INK, lw=1.2, ls="--")
        ax.set_ylabel("Vertical acceleration PSD [dB re 1 (m/s²)²/Hz]", fontsize=9)
    else:
        webb = nm.webb_pressure_bounds()
        if webb is not None:
            ax.vlines(
                webb[0], webb[1], webb[2], color=INK, lw=3, label="Webb (1998) range"
            )
        ax.set_ylabel("Pressure PSD [dB re 1 Pa²/Hz]", fontsize=9)
    ax.set_xlabel("Frequency [Hz]", fontsize=9)
    ax.set_xlim(1e-3, 0.6)
    ax.legend(fontsize=7, frameon=False)
    _style(ax)
    ax.set_title(title or f"Stacked {kind} PSD of accepted windows", fontsize=10)
    return fig


def plot_multistation(stations, *, ncols=6, floor_rel=0.01, title=None):
    """Small multiples of compliance (with ±1σ total error) and coherence.

    Purpose: compare shape and level across a network; a station whose
    curve bends differently or whose coherence collapses stands out.
    """
    codes = sorted(stations)
    nrows = int(np.ceil(len(codes) / ncols))
    fig, axes = plt.subplots(2 * nrows, ncols, figsize=(2.2 * ncols, 3.0 * nrows),
                             sharex=True, constrained_layout=True, squeeze=False,
                             gridspec_kw={"height_ratios": [2, 1] * nrows})  # fmt: skip
    for i, code in enumerate(codes):
        s = stations[code]
        r, c = divmod(i, ncols)
        a, b = axes[2 * r, c], axes[2 * r + 1, c]
        bud = error_budget(s, floor_rel=floor_rel)
        f = bud["freqs"]
        a.fill_between(f, 1e11 * (bud["eta"] - bud["total"]), 1e11 * (bud["eta"] + bud["total"]),
                       color="#86b6ef", lw=0)  # fmt: skip
        a.plot(f, 1e11 * bud["eta"], color=CATEGORICAL[0], lw=1.4)
        a.set_title(f"{code} ({s.label})", fontsize=8)
        b.plot(s.freqs, s.coherence, color=CATEGORICAL[1], lw=1.2)
        b.set_ylim(0.5, 1.0)
        _band(b, s.band)
        for ax in (a, b):
            _style(ax)
            ax.tick_params(labelsize=6)
        if c == 0:
            a.set_ylabel("η [1e-11/Pa]", fontsize=7)
            b.set_ylabel("Coherence", fontsize=7)
    for j in range(len(codes), nrows * ncols):
        r, c = divmod(j, ncols)
        axes[2 * r, c].axis("off")
        axes[2 * r + 1, c].axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("Frequency [Hz]", fontsize=7)
    fig.suptitle(title or "Compliance (±1σ total) and Z–P coherence", fontsize=10)
    return fig


def plot_phase_diagnostic(models, *, groups=None, title=None):
    """Residual Z/P phase per station and the compliance bias it may imply.

    Purpose: an elastic seafloor gives a real Z/P transfer function; a
    frequency-dependent residual phase exposes a response mismatch that a
    scalar Rayleigh-wave gain does not remove. Left: residual phase (dots)
    with the best response-model fit (line). Right: the range of
    compliance bias at each station's lowest band frequency allowed by the
    acceptable fits of the DPG-corner and seismometer-corner models.

    Args:
        models (dict): ``{station: compy.anomaly.phase_models(...)}``.
        groups (dict): ``{station: group label}`` (e.g. network) for colour.
    """
    codes = [
        c for c in sorted(models) if np.max(np.abs(models[c]["residual_deg"])) < 45
    ]
    labels = sorted({(groups or {}).get(c, "") for c in codes})
    color = {g: CATEGORICAL[i % 8] for i, g in enumerate(labels)}
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True,
                                  gridspec_kw={"width_ratios": [1.2, 1]})  # fmt: skip
    seen = set()
    for code in codes:
        m = models[code]
        g = (groups or {}).get(code, "")
        best = min(("dpg_corner", "seismometer_corner"), key=lambda k: m[k]["rms_deg"])
        ax.plot(1e3 * m["freqs"], m["residual_deg"], "o", ms=3, color=color[g],
                label=None if g in seen else (g or "stations"))  # fmt: skip
        ax.plot(
            1e3 * m["freqs"],
            m[best]["prediction_deg"],
            "-",
            lw=0.8,
            color=color[g],
            alpha=0.7,
        )
        seen.add(g)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("Frequency [mHz]", fontsize=9)
    ax.set_ylabel("Z/P phase − nearest of 0°/180° [°]", fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    y = np.arange(len(codes))
    for i, code in enumerate(codes):
        m = models[code]
        g = (groups or {}).get(code, "")
        for k, (key, dy) in enumerate(
            (("dpg_corner", -0.15), ("seismometer_corner", 0.15))
        ):
            lo, hi = m[key]["bias_low_pct"][0], m[key]["bias_high_pct"][0]
            ax2.plot([lo, hi], [i + dy, i + dy], color=CATEGORICAL[0] if k == 0 else CATEGORICAL[1],
                     lw=3, solid_capstyle="butt",
                     label=None if i else ("DPG-corner fits" if k == 0 else "seismometer-corner fits"))  # fmt: skip
    ax2.axvline(0, color=MUTED, lw=0.8)
    ax2.set_xlim(-25, 15)
    ax2.set_yticks(y, codes, fontsize=7)
    ax2.invert_yaxis()
    ax2.set_xlabel("Implied compliance bias at band bottom [%]", fontsize=9)
    ax2.legend(fontsize=7, frameon=True, loc="upper right")
    _style(ax2)
    fig.suptitle(
        title or "Residual Z/P phase and the response mismatch it implies", fontsize=10
    )
    return fig


# --------------------------------------------------------------------------
# Per-window levels, periodicity and anomalies
# --------------------------------------------------------------------------
def plot_window_levels(s, info=None, *, title=None):
    """Per-window |Z/P|, pressure and vertical levels with detected steps.

    Purpose: a |Z/P| step without a pressure-level step is a sensitivity
    change of one sensor; a step with both is a change of sea state; a
    slow trend is drift. ``info`` is ``scan_station(...)[1]["levels"]``.
    """
    t, r, p_db, a_db = np.asarray(s.window_levels, float).T
    x = _dates(t)
    fig, axes = plt.subplots(
        3, 1, figsize=(10, 6.5), sharex=True, constrained_layout=True
    )
    series = ((r - np.median(r)) * np.log(10) * 100, p_db, a_db)
    labels = ("|Z/P| − median [%]", "Pressure [dB re 1 Pa²/Hz,\nnominal]",
              "Vertical acc. [dB re 1 (m/s²)²/Hz]")  # fmt: skip
    keys = ("transfer function", "pressure level", "vertical level")
    cps = (info or {}).get("changepoints", {})
    for ax, y, label, key, color in zip(axes, series, labels, keys, CATEGORICAL[:3]):
        ax.plot(x, y, ".", ms=1.5, color=color, rasterized=True)
        edges = [t.min(), *sorted(c.time for c in cps.get(key, [])), t.max() + 1]
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (t >= a) & (t < b)
            if sel.any():
                ax.plot(_dates([a, b]), [np.median(y[sel])] * 2, color=INK, lw=1.5)
        for c in cps.get("transfer function", []):
            ax.axvline(_dates([c.time])[0], color=CATEGORICAL[7], lw=1, ls=":")
        ax.set_ylabel(label, fontsize=8)
        _style(ax)
    _date_axis(axes[-1])
    axes[0].set_title(title or f"{s.station}: accepted-window levels; black: segment medians; "
                      "dotted: |Z/P| steps", fontsize=10)  # fmt: skip
    return fig


def plot_periodicity(results, *, value="modulation", title=None):
    """Station x constituent map of tidal/diurnal modulation.

    Args:
        results (dict): ``{station: rayleigh_periodicity(...) or periodicity(...)}``.
        value (str): "modulation" (2R, accepted-window clustering, %) or
            "amplitude" (log10 amplitude, shown in %).
    """
    codes = sorted(results)
    names = list(next(iter(results.values())))
    grid = np.zeros((len(codes), len(names)))
    sig = np.zeros_like(grid, bool)
    for i, c in enumerate(codes):
        for j, n in enumerate(names):
            v = results[c][n]
            grid[i, j] = 200 * v["R"] if "R" in v else 100 * np.log(10) * v["amplitude"]
            sig[i, j] = v["p_value"] <= 0.01 and v.get("p_tail", 0) <= 0.01
    fig, ax = plt.subplots(figsize=(1.0 * len(names) + 2.5, 0.28 * len(codes) + 1.4),
                           constrained_layout=True)  # fmt: skip
    mesh = ax.pcolormesh(grid, cmap=SEQUENTIAL, vmin=0)
    for i, j in zip(*np.nonzero(sig)):
        ax.text(j + 0.5, i + 0.5, "●", ha="center", va="center", fontsize=7,
                color="white" if grid[i, j] > 0.6 * grid.max() else INK)  # fmt: skip
    ax.set_xticks(np.arange(len(names)) + 0.5, names, fontsize=8)
    ax.set_yticks(np.arange(len(codes)) + 0.5, codes, fontsize=7)
    ax.invert_yaxis()
    fig.colorbar(
        mesh,
        ax=ax,
        label="Modulation [%]" if value == "modulation" else "Amplitude [%]",
    )
    ax.set_title(
        title or "Periodic modulation (● empirical and tail p ≤ 0.01)", fontsize=9
    )
    return fig


def plot_acceptance_phase(
    stations, *, period_h=12.4206012, bins=24, ncols=6, title=None
):
    """Accepted windows folded on a period, one panel per station.

    Purpose: with no tidal effect the histogram is flat; a trough at one
    phase means windows there fail the coherence gate (e.g. tilt noise at
    peak tidal flow). Rate relative to the mean.
    """
    codes = sorted(stations)
    nrows = int(np.ceil(len(codes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.0 * ncols, 1.6 * nrows + 0.6),
                             sharex=True, sharey=True, constrained_layout=True, squeeze=False)  # fmt: skip
    edges = np.linspace(0, 1, bins + 1)
    for i, code in enumerate(codes):
        ax = axes.flat[i]
        t = np.asarray(stations[code].window_starts, float)
        h, _ = np.histogram((t / 3600.0 / period_h) % 1.0, edges)
        ax.bar(edges[:-1], h / h.mean(), width=1 / bins, align="edge", color=CATEGORICAL[0],
               edgecolor="white", lw=0.5)  # fmt: skip
        ax.axhline(1, color=INK, lw=0.8)
        ax.set_title(code, fontsize=7)
        _style(ax)
        ax.tick_params(labelsize=6)
    for j in range(len(codes), nrows * ncols):
        axes.flat[j].axis("off")
    fig.suptitle(title or f"Accepted windows folded on {period_h:.2f} h (rate / mean; phase 0 = epoch)",
                 fontsize=9)  # fmt: skip
    return fig


def plot_common_mode(cm, *, title=None):
    """Monthly deviations by station and their cross-station correlation.

    Args:
        cm (dict): Output of ``compy.anomaly.common_mode``.
    """
    m = cm["matrix"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 0.25 * len(cm["stations"]) + 2.2),
                                  constrained_layout=True, gridspec_kw={"width_ratios": [1.5, 1]})  # fmt: skip
    lim = max(2.0, np.nanpercentile(np.abs(m), 97))
    mesh = ax.pcolormesh(
        np.ma.masked_invalid(m), cmap=DIVERGING, norm=TwoSlopeNorm(0, -lim, lim)
    )
    ax.set_xticks(np.arange(len(cm["keys"])) + 0.5, cm["keys"], rotation=70, fontsize=6)
    ax.set_yticks(np.arange(len(cm["stations"])) + 0.5, cm["stations"], fontsize=7)
    ax.invert_yaxis()
    fig.colorbar(mesh, ax=ax, label="Band-mean monthly − full [%]", pad=0.01)
    mesh2 = ax2.pcolormesh(
        np.ma.masked_invalid(cm["corr"]), cmap=DIVERGING, norm=TwoSlopeNorm(0, -1, 1)
    )
    ax2.set_xticks(
        np.arange(len(cm["stations"])) + 0.5, cm["stations"], rotation=90, fontsize=6
    )
    ax2.set_yticks(np.arange(len(cm["stations"])) + 0.5, cm["stations"], fontsize=6)
    ax2.invert_yaxis()
    fig.colorbar(mesh2, ax=ax2, label="Spearman ρ", pad=0.01)
    ax.set_title(
        title or "Monthly compliance deviations (left) and their correlation (right)",
        fontsize=9,
    )
    return fig


def plot_anomaly_ranking(anomalies, *, top=25, title=None):
    """Horizontal bars of the highest-scoring anomalies, coloured by category."""
    items = list(anomalies)[:top]
    fig, ax = plt.subplots(
        figsize=(9, 0.32 * len(items) + 1.2), constrained_layout=True
    )
    y = np.arange(len(items))
    ax.barh(y, [a.score for a in items], color=[CATEGORY_COLOR.get(a.category, MUTED) for a in items],
            edgecolor="white", height=0.75)  # fmt: skip
    ax.set_yticks(y, [f"{a.station} · {a.kind}" for a in items], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Score = min(σ, 8) × (1 + min(|effect %|, 20)/5)", fontsize=8)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in CATEGORY_COLOR.values()]
    ax.legend(
        handles, list(CATEGORY_COLOR), fontsize=7, frameon=False, loc="lower right"
    )
    _style(ax)
    ax.set_title(title or "Ranked anomalies", fontsize=10)
    return fig
