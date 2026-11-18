"""Standard figures of the pressure atlas (Matplotlib, Agg-safe).

Every figure uses frequency and period axes on a log scale, levels in
dB re 1 Pa^2/Hz, and legends inside the axes. Colours are fixed per category so
that a sensor type or depth class has the same colour in every figure.
"""

from __future__ import annotations

import json
from importlib import resources

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from . import native

GRID = native.GRID
TYPE_COLOURS = {"DPG": "#1f5fa8", "APG": "#c0392b", "Hydrophone": "#2e8b57", "Other": "#7f7f7f"}
DEPTH_CMAP = "viridis_r"
BANDS_SHADE = {"secondary microseism": (3, 10), "notch": (15, 50), "infragravity": (100, 300)}


def brown_2014():
    """Digitized Brown et al. (2014) IDC hydroacoustic 5 % and 95 % curves: (frequency Hz, low dB, high dB)."""
    ref = resources.files("compy_pressure.reference_data")
    meta = json.loads((ref / "pressure_references.json").read_text())["brown"]
    rows = np.genfromtxt(str(ref / meta["data_file"]), delimiter=",", names=True)
    return rows["frequency_hz"], rows["low_db_pa2_hz"], rows["high_db_pa2_hz"]


def period_axis(ax, lo=0.02, hi=1300.0, label=True):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    if label:
        ax.set_xlabel("Period (s)")
    ax.grid(True, which="major", alpha=0.3)
    ax.grid(True, which="minor", alpha=0.08)


def shade_bands(ax, alpha=0.06):
    for name, (a, b) in BANDS_SHADE.items():
        ax.axvspan(a, b, color="k", alpha=alpha, lw=0)


def density_panel(ax, density, edges, title="", curves=None, gauges=None, vmin=1e-3, vmax=0.3, brown=True):
    """Probability density (period x dB) with optional percentile curves {label: (y, style)}."""
    P = np.ma.masked_invalid(np.where(density > 0, density, np.nan))
    centres = 0.5 * (np.asarray(edges)[:-1] + np.asarray(edges)[1:])
    pc = ax.pcolormesh(GRID, centres, P.T, shading="nearest", cmap="magma_r", norm=LogNorm(vmin, vmax), rasterized=True)
    for lab, (y, kw) in (curves or {}).items():
        ax.plot(GRID, y, label=lab, **kw)
    if brown:
        f, lo, hi = brown_2014()
        for y in (lo, hi):
            s = f >= 0.2
            ax.plot(1 / f[s], y[s], color="#0b7285", lw=1.3)
            ax.plot(1 / f[~s], y[~s], color="#0b7285", lw=1.1, ls=(0, (4, 2.5)))
    period_axis(ax)
    ax.set_ylabel("Pressure PSD (dB re 1 Pa$^2$/Hz)")
    if title:
        ax.set_title(title, fontsize=10, loc="left")
    if gauges is not None:
        ax2 = ax.twinx()
        ax2.plot(GRID, gauges, color="0.4", lw=0.8, ls=":")
        ax2.set_ylabel("gauges at full weight", color="0.4", fontsize=8)
        ax2.tick_params(labelsize=7, colors="0.4")
    return pc


def model_lines(ax, tables, colours, which=("p5_db", "p50_db", "p95_db"), min_gauges=8):
    """Median (solid) and 5/95 % (thin) of several model tables {label: DataFrame}; dashed where fewer than 8 gauges."""
    for lab, t in tables.items():
        c = colours.get(lab, None)
        thick = t.gauges.values >= min_gauges
        for col, lw in zip(which, (0.9, 2.0, 0.9)):
            y = t[col].values
            ax.plot(GRID, np.where(thick, y, np.nan), color=c, lw=lw, label=f"{lab} ({int(t.gauges.max())})" if col == "p50_db" else None)
            ax.plot(GRID, np.where(~thick, y, np.nan), color=c, lw=lw * 0.8, ls="--")
    period_axis(ax)
    ax.set_ylabel("Pressure PSD (dB re 1 Pa$^2$/Hz)")


def save(fig, path, dpi=160):
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
