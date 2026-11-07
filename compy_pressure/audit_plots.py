"""Diagnostic figures of the pressure-gauge audit (ObsPy-free).

=====================================  =========================================
Figure                                 Reveals
=====================================  =========================================
:func:`plot_gain_consistency`          which calibration methods agree, per
                                       channel, and where one method breaks
:func:`plot_gain_depth`                whether gains depend on depth (physics)
                                       or on the network (instruments/metadata)
:func:`plot_water_layer_sensitivity`   how much the Rayleigh/microseism gain
                                       depends on the assumed water column and
                                       crust, by depth and frequency
:func:`plot_infragravity_residuals`    Godin residuals by basin, depth, sensor
:func:`plot_network_scores`            the six QC/calibration scores by network
=====================================  =========================================
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from compy.plotting.diagnostics import (
    CATEGORICAL,
    DIVERGING,
    GRID,
    INK,
    MUTED,
    SEQUENTIAL,
    _style,
)
from compy_pressure import water_layer
from compy_pressure.table_audit import gains

TYPES = ("DPG", "Hydrophone", "APG (absolute)")
TYPE_COLOR = dict(zip(TYPES, CATEGORICAL[:3]))


def _db(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 20 * np.log10(x)


def plot_gain_consistency(table, *, limit_db=30.0):
    """Per-channel multipliers of pairs of methods, with ±1.5 and ±6 dB bands.

    Purpose: the provisional score passes |20 log10 g| <= 1.5 dB and fails
    >= 6 dB; points off the diagonal are channels where one method failed or
    the methods measure different things.
    """
    g = gains(table)
    pairs = (
        ("microseism", "rayleigh"),
        ("local_p", "rayleigh"),
        ("peer_offset", "microseism"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)
    for ax, (a, b) in zip(axes, pairs):
        both = np.isfinite(g[a]) & np.isfinite(g[b])
        x, y = (
            np.clip(_db(g[b][both]), -limit_db, limit_db),
            np.clip(_db(g[a][both]), -limit_db, limit_db),
        )
        line = np.array([-limit_db, limit_db])
        ax.fill_between(line, line - 6, line + 6, color="#eeeeea", lw=0, label="±6 dB")
        ax.fill_between(
            line, line - 1.5, line + 1.5, color="#cde2fb", lw=0, label="±1.5 dB"
        )
        ax.plot(line, line, color=INK, lw=0.8)
        for kind in TYPES:
            sel = table["type"][both] == kind
            if sel.any():
                ax.plot(
                    x[sel],
                    y[sel],
                    "o",
                    ms=4,
                    color=TYPE_COLOR[kind],
                    label=f"{kind} ({sel.sum()})",
                )
        d = _db(g[a][both] / g[b][both])
        ax.set_title(f"{a} vs {b}: n = {both.sum()}, median {np.median(d):+.1f} dB, "
                     f"{np.mean(np.abs(d) <= 6):.0%} within 6 dB", fontsize=8)  # fmt: skip
        ax.set_xlabel(f"{b} g [dB]", fontsize=9)
        ax.set_ylabel(f"{a} g [dB]", fontsize=9)
        ax.set_xlim(-limit_db, limit_db)
        ax.set_ylim(-limit_db, limit_db)
        ax.legend(fontsize=6, frameon=False, loc="upper left")
        _style(ax)
    fig.suptitle(f"Pressure multipliers of the same channel by two methods (clipped at ±{limit_db:.0f} dB)",
                 fontsize=10)  # fmt: skip
    return fig


def plot_gain_depth(table, info=None):
    """Gains against depth across the table, and within network x sensor groups.

    Purpose: a depth trend that survives within groups points to physics
    (water-column model, pressure-dependent sensitivity); one that vanishes
    is the networks' own instruments or metadata.
    """
    g = gains(table)
    depth = table["Profondeur"] / 1000
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    ok = np.isfinite(g["microseism"]) & np.isfinite(depth)
    for kind in TYPES:
        sel = ok & (table["type"] == kind)
        ax.plot(depth[sel], np.clip(_db(g["microseism"][sel]), -20, 20), "o", ms=3.5,
                color=TYPE_COLOR[kind], label=kind)  # fmt: skip
    sel = ok & np.isfinite(g["rayleigh"])
    ax.plot(depth[np.isfinite(g["rayleigh"])], np.clip(_db(g["rayleigh"][np.isfinite(g["rayleigh"])]), -20, 20),
            "x", ms=5, color=INK, label="Rayleigh (all types)")  # fmt: skip
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("Water depth [km]", fontsize=9)
    ax.set_ylabel("Pressure multiplier g [dB] (clipped ±20)", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    _style(ax)
    xs, ys, labels = [], [], []
    for net in np.unique(table["network"][ok]):
        for kind in TYPES:
            grp = ok & (table["network"] == net) & (table["type"] == kind)
            if grp.sum() >= 4:
                xs.append(depth[grp] - depth[grp].mean())
                ys.append(_db(g["microseism"][grp]) - _db(g["microseism"][grp]).mean())
                labels.append(f"{net} {kind.split()[0]}")
    for i, (x, y, lab) in enumerate(zip(xs, ys, labels)):
        ax2.plot(x, np.clip(y, -15, 15), "o", ms=3.5, color=CATEGORICAL[i % 8] if i < 8 else MUTED,
                 label=lab if i < 8 else ("other groups" if i == 8 else None))  # fmt: skip
    if xs:
        allx, ally = np.concatenate(xs), np.concatenate(ys)
        slope = np.polyfit(allx, ally, 1)[0]
        line = np.array([allx.min(), allx.max()])
        ax2.plot(
            line,
            slope * line,
            color=INK,
            lw=1.4,
            label=f"pooled slope {slope:+.1f} dB/km",
        )
    ax2.axhline(0, color=MUTED, lw=0.8)
    ax2.set_xlabel("Depth − group mean [km]", fontsize=9)
    ax2.set_ylabel("Microseism g − group mean [dB]", fontsize=9)
    ax2.legend(fontsize=6, frameon=False, ncol=2)
    _style(ax2)
    ax.set_title("Across the table", fontsize=9)
    ax2.set_title("Within network × sensor groups (≥ 4 channels)", fontsize=9)
    return fig


def plot_water_layer_sensitivity(dispersion_f, dispersion_c, *, depths=np.linspace(500, 6000, 56),
                                 freqs=np.linspace(0.02, 0.08, 61)):  # fmt: skip
    """Model sensitivity of the water-layer P/a ratio and the gain bias it causes.

    Left: d ln(P/a)/d ln(sound speed) over depth and frequency (white: past
    the resonance guard). Right: the multiplier a perfect gauge receives,
    30-70 mHz, when the true column or crust differ from the model.
    """
    c = np.interp(freqs, dispersion_f, dispersion_c)
    sens = np.full((len(depths), len(freqs)), np.nan)
    for i, h in enumerate(depths):
        s = water_layer.log_sensitivity(freqs, c, h)
        ok = water_layer.usable(freqs, c, h)
        sens[i, ok] = s["sound_speed"][ok]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    mesh = ax.pcolormesh(
        freqs * 1e3,
        depths / 1000,
        -sens,
        cmap=SEQUENTIAL,
        vmin=0,
        vmax=10,
        shading="auto",
    )
    ax.axvspan(30, 70, color="none", ec=INK, lw=1.2, ls="--")
    fig.colorbar(mesh, ax=ax, label="−d ln(P/a) / d ln(sound speed)")
    ax.set_xlabel("Frequency [mHz] (dashed: Rayleigh band)", fontsize=9)
    ax.set_ylabel("Water depth [km]", fontsize=9)
    ax.set_title(
        "Sensitivity to the assumed sound speed (white: rejected near resonance)",
        fontsize=9,
    )
    f = np.arange(0.03, 0.0701, 0.0025)
    cf = np.interp(f, dispersion_f, dispersion_c)
    cases = (("sound speed +1 % (1515 m/s)", dict(true_sound_speed=1515.0), CATEGORICAL[0], "-"),
             ("crust phase velocity −5 %", dict(true_phase_velocity=0.95 * cf), CATEGORICAL[1], "-"),
             ("crust phase velocity +5 %", dict(true_phase_velocity=1.05 * cf), CATEGORICAL[2], "-"),
             ("in-situ column-mean density", "density", CATEGORICAL[3], "-"),
             ("sound speed +1 %, log-ratio fit", dict(true_sound_speed=1515.0, weights="log"), CATEGORICAL[0], "--"))  # fmt: skip
    d = np.linspace(500, 5800, 54)
    for label, kw, color, style in cases:
        bias = []
        for h in d:
            try:
                if kw == "density":
                    b = water_layer.gain_bias(
                        f, cf, h, true_density=float(water_layer.column_mean_density(h))
                    )
                else:
                    b = water_layer.gain_bias(f, cf, h, **kw)
            except ValueError:
                b = np.nan
            bias.append(100 * (b - 1))
        ax2.plot(d / 1000, bias, style, color=color, lw=1.6, label=label)
    ax2.axhline(0, color=MUTED, lw=0.8)
    ax2.set_xlabel("Water depth [km]", fontsize=9)
    ax2.set_ylabel("Multiplier returned − 1 [%]", fontsize=9)
    ax2.legend(fontsize=7, frameon=False)
    ax2.set_title("Bias of a 30–70 mHz Rayleigh gain for a perfect gauge", fontsize=9)
    _style(ax2)
    return fig


def plot_infragravity_residuals(table):
    """Godin (2013) Eq. 6 residuals by depth, basin and sensor type.

    Purpose: a basin-wide excess is ocean climate (the model is a single
    average spectrum); a sensor-type deficit at one network is a response
    that does not reach the infragravity band.
    """
    res = table["infragravity_median_residual_db"]
    depth = table["Profondeur"] / 1000
    basins = [
        b
        for b in ("Pacific", "Atlantic", "Indian", "Mediterranean")
        if np.any((table["basin"] == b) & np.isfinite(res))
    ]
    fig, axes = plt.subplots(1, len(basins), figsize=(3.6 * len(basins), 4.2), sharey=True,
                             constrained_layout=True, squeeze=False)  # fmt: skip
    for ax, b in zip(axes[0], basins):
        for kind in TYPES:
            sel = np.isfinite(res) & (table["basin"] == b) & (table["type"] == kind)
            if sel.any():
                ax.plot(depth[sel], np.clip(res[sel], -50, 30), "o", ms=3, color=TYPE_COLOR[kind],
                        label=f"{kind} ({sel.sum()})")  # fmt: skip
        sel = np.isfinite(res) & (table["basin"] == b)
        ax.axhspan(-10, 10, color="#cde2fb", lw=0)
        ax.axhline(
            np.median(res[sel]),
            color=INK,
            lw=1.2,
            label=f"median {np.median(res[sel]):+.1f} dB",
        )
        ax.set_title(b, fontsize=9)
        ax.set_xlabel("Water depth [km]", fontsize=9)
        ax.legend(fontsize=6, frameon=False, loc="lower left")
        _style(ax)
    axes[0, 0].set_ylabel("Residual to Godin Eq. 6 [dB] (blue: score 1)", fontsize=9)
    return fig


def plot_network_scores(table):
    """Network x check map of mean scores (0 fail, 1 pass; blank: not scored)."""
    cols = ("score_expected_levels", "score_nearby", "score_godin", "score_rayleigh",
            "score_microseism", "score_local_p")  # fmt: skip
    names = ("Expected", "Nearby", "Godin", "Rayleigh", "Microseism", "Local P")
    nets = sorted(np.unique(table["network"]))
    grid = np.full((len(nets), len(cols)), np.nan)
    counts = np.zeros_like(grid, int)
    for i, n in enumerate(nets):
        sel = table["network"] == n
        for j, c in enumerate(cols):
            v = table[c][sel]
            v = v[np.isfinite(v)]
            if v.size:
                grid[i, j], counts[i, j] = v.mean(), v.size
    fig, ax = plt.subplots(figsize=(7, 0.26 * len(nets) + 1.5), constrained_layout=True)
    mesh = ax.pcolormesh(np.ma.masked_invalid(grid), cmap=DIVERGING.reversed(), norm=TwoSlopeNorm(0.5, 0, 1),
                         edgecolors=GRID, lw=0.5)  # fmt: skip
    for i, j in zip(*np.nonzero(counts)):
        ax.text(
            j + 0.5,
            i + 0.5,
            str(counts[i, j]),
            ha="center",
            va="center",
            fontsize=6,
            color=INK,
        )
    ax.set_xticks(np.arange(len(cols)) + 0.5, names, fontsize=8)
    ax.set_yticks(np.arange(len(nets)) + 0.5, nets, fontsize=7)
    ax.invert_yaxis()
    fig.colorbar(mesh, ax=ax, label="Mean score (numbers: channels scored)")
    ax.set_title("QC and calibration scores by network", fontsize=9)
    return fig
