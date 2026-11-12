"""Report figures 1-3: what compliance measures, and what its errors are.

Fig. 1  Thin soft sediment acts through its P modulus; relative sensitivity to Vs, Vp, density.
Fig. 2  Errors of a stacked compliance curve on 172 stations.
Fig. 3  Old and new covariance at one station (YV.RR52), and two independent gain calibrations.

Run: PYTHONPATH=<worktree> python g01_figures_method.py
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compy.inversion import dram
from compy.numerics import wavenumber
from compy.plotting.auto import GRID, INK, MUTED, POSTERIOR, PRIOR, REFERENCE, STYLE

OUT = Path(__file__).resolve().parents[1] / "figures"
RUN = Path.home() / "Desktop/ComPy/14_inversion_overhaul"
BLUE, ORANGE, AQUA = POSTERIOR, PRIOR, REFERENCE


def fig1():
    depth = 3000.0
    rock = np.array(
        [
            [400, 2500, 4600, 2400],
            [1500, 2750, 6000, 3300],
            [4500, 2950, 6900, 3850],
            [0, 3300, 8000, 4500],
        ],
        float,
    )
    sed = np.array([[100, 1700, 1560, 200], [200, 1900, 1750, 450]], float)
    f = np.linspace(0.004, 0.02, 60)
    fw = dram.Forward(depth, f)
    k = wavenumber(2 * np.pi * f, depth)
    e_rock = fw(rock)
    full = np.vstack([sed, rock])
    e_sed = fw(full)
    half, double = full.copy(), full.copy()
    half[:2, 3] *= 0.5
    double[:2, 3] *= 2
    uni = e_rock + k * np.sum(sed[:, 0] / (sed[:, 1] * sed[:, 2] ** 2))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), gridspec_kw=dict(wspace=0.28))
    ax = axes[0]
    s = 1e11
    ax.fill_between(
        1e3 * f,
        s * fw(double),
        s * fw(half),
        color=BLUE,
        alpha=0.2,
        lw=0,
        label="sediment Vs halved to doubled",
    )
    ax.plot(1e3 * f, s * e_sed, color=BLUE, lw=2, label="rock + 300 m of sediment")
    ax.plot(
        1e3 * f, s * uni, color=ORANGE, lw=1.6, ls="--", label="rock + k h / (ρ Vp²)"
    )
    ax.plot(1e3 * f, s * e_rock, color=INK, lw=1.6, label="bare rock")
    ax.set(
        xlabel="Frequency (mHz)",
        ylabel="Normalized compliance (10$^{-11}$ Pa$^{-1}$)",
        title="300 m of sediment under 3000 m of water",
        ylim=(0, 12),
        xlim=(4, 20),
    )
    ax.legend(loc="upper left")
    # sensitivities of bare rock
    ax = axes[1]
    fs = np.linspace(0.004, 0.02, 17)
    fw2 = dram.Forward(depth, fs)
    base = np.log(fw2(rock))
    for col, label, color in (
        (3, "shear velocity Vs", BLUE),
        (2, "compressional velocity Vp", ORANGE),
        (1, "density", AQUA),
    ):
        trial = rock.copy()
        trial[:, col] *= 1.01
        ax.plot(
            1e3 * fs, (np.log(fw2(trial)) - base) / 0.01, color=color, lw=2, label=label
        )
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set(
        xlabel="Frequency (mHz)",
        ylabel="Change of compliance per unit relative change",
        title="Sensitivity of compliance (bare rock, all layers together)",
        ylim=(-1.4, 0.1),
        xlim=(4, 20),
    )
    ax.legend(loc="lower right")
    fig.savefig(OUT / "fig01_physics.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def fig2():
    with open(RUN / "e02_errors/error_model.csv") as fh:
        rows = list(csv.DictReader(fh))
    with open(
        Path.home() / "Desktop/ComPy/12_inversion_all_stations/tables/station_plan.csv"
    ) as fh:
        tier = {r["sid"]: r["tier"] for r in csv.DictReader(fh)}
    rows = [r for r in rows if tier.get(r["sid"]) == "1"]

    def col(name):
        return np.array([float(r[name]) for r in rows if r[name] not in ("", "nan")])

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(wspace=0.32))
    ax = axes[0]
    series = (
        ("bp_med", "Bendat & Piersol formula", INK),
        ("iid_med", "windows resampled", BLUE),
        ("b30_med", "30-day blocks resampled", ORANGE),
    )
    for name, label, color in series:
        v = np.sort(100 * col(name))
        ax.plot(v, (np.arange(len(v)) + 0.5) / len(v), color=color, lw=2, label=label)
    ax.axvline(1.0, color=MUTED, ls=":", lw=1)
    ax.text(
        1.03,
        0.05,
        "1 % floor of ComPy 2.x",
        rotation=90,
        color=MUTED,
        fontsize=8,
        va="bottom",
    )
    ax.set(
        xscale="log",
        xlim=(0.05, 10),
        xlabel="Relative error of the stacked curve (%)",
        ylabel="Share of 172 stations",
        title="Scatter of a stacked compliance curve",
    )
    ax.legend(loc="lower right")
    ax = axes[1]
    lags = ("iid_corr_adjacent", "iid_corr_lag2", "b30_corr_adjacent", "b30_corr_ends")
    labels = (
        "windows\n1 bin apart",
        "windows\n2 bins apart",
        "blocks\n1 bin apart",
        "blocks\nband ends",
    )
    data = [col(n) for n in lags]
    parts = ax.boxplot(
        data,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color=INK, lw=1.5),
    )
    for patch, color in zip(parts["boxes"], (BLUE, BLUE, ORANGE, ORANGE)):
        patch.set(facecolor=color, alpha=0.35, edgecolor=color)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xticks(range(1, 5), labels, fontsize=8)
    ax.set(
        ylabel="Correlation between frequency bins",
        title="Errors are correlated between bins\n(resampling windows, or 30-day blocks)",
        ylim=(-0.6, 1.05),
    )
    ax = axes[2]
    for name, label, color in (
        ("half_level", "level", BLUE),
        ("half_tilt", "tilt across the band", ORANGE),
    ):
        v = np.sort(100 * np.abs(col(name)))
        ax.plot(v, (np.arange(len(v)) + 0.5) / len(v), color=color, lw=2, label=label)
    ax.set(
        xscale="log",
        xlim=(0.02, 30),
        xlabel="First half against second half of the record (%)",
        ylabel="Share of 172 stations",
        title="What changes over a deployment",
    )
    ax.legend(loc="upper left")
    fig.savefig(OUT / "fig02_errors.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def fig3():
    old = np.load(
        Path.home()
        / "Desktop/ComPy/12_inversion_all_stations/YV/results/RR52/inversion/RR52_final.npz"
    )
    new_path = RUN / "t03_pilot/YV/RR52/RR52_auto.npz"
    new = np.load(new_path)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(wspace=0.32))
    ax = axes[0]
    for d, label, color in (
        (old, "ComPy 2.8 (1 % independent floor)", ORANGE),
        (new, "ComPy 3 (measured + tilt + smooth)", BLUE),
    ):
        sig = np.sqrt(np.diag(d["covariance"])) / d["eta"]
        ax.plot(1e3 * d["freqs"], 100 * sig, "o-", color=color, lw=2, ms=4, label=label)
    ax.set(
        xlabel="Frequency (mHz)",
        ylabel="One-sigma error (%)",
        title="YV.RR52: size of the errors",
        ylim=(0, 2),
    )
    ax.legend(loc="lower left")
    ax = axes[1]
    for d, label, color in ((old, "ComPy 2.8", ORANGE), (new, "ComPy 3", BLUE)):
        c = d["covariance"]
        sd = np.sqrt(np.diag(c))
        corr = c / np.outer(sd, sd)
        n = len(sd)
        lag = np.arange(1, n)
        ax.plot(
            lag,
            [np.mean(np.diag(corr, k)) for k in lag],
            "o-",
            color=color,
            lw=2,
            ms=4,
            label=label,
        )
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set(
        xlabel="Separation (frequency bins of 1.22 mHz)",
        ylabel="Mean correlation of the errors",
        title="YV.RR52: correlation between bins",
        ylim=(-1.05, 1.05),
    )
    ax.legend(loc="lower left")
    ax = axes[2]
    with open(
        Path.home() / "Desktop/ComPy/05_report/tables/compliance_results.csv"
    ) as fh:
        rows = [r for r in csv.DictReader(fh) if r["gain_bruitfm"] and r["gain_compy"]]
    lr = np.array(
        [100 * np.log(float(r["gain_compy"]) / float(r["gain_bruitfm"])) for r in rows]
    )
    ax.hist(
        np.clip(lr, -30, 30),
        bins=np.arange(-30, 31, 2.5),
        color=BLUE,
        alpha=0.75,
        edgecolor="white",
    )
    x = np.linspace(-30, 30, 300)
    ax.plot(
        x,
        len(lr)
        * 2.5
        * np.exp(-0.5 * (x / (5 * np.sqrt(2))) ** 2)
        / (5 * np.sqrt(2) * np.sqrt(2 * np.pi)),
        color=ORANGE,
        lw=2,
        label="two gains with 5 % error each",
    )
    ax.set(
        xlabel="Rayleigh gain against microseism gain (%)",
        ylabel="Stations",
        title=f"Two calibrations of {len(lr)} gauges",
    )
    ax.legend(loc="upper left")
    fig.savefig(OUT / "fig03_covariance_gain.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    (OUT / "fig03_numbers.json").write_text(
        json.dumps(
            dict(
                n=len(lr),
                median=float(np.median(lr)),
                robust_sd=float(1.4826 * np.median(np.abs(lr - np.median(lr)))),
                p16=float(np.quantile(lr, 0.16)),
                p84=float(np.quantile(lr, 0.84)),
            ),
            indent=1,
        )
    )


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    with matplotlib.rc_context(STYLE):
        fig1()
        fig2()
        fig3()
    print("figures in", OUT)
