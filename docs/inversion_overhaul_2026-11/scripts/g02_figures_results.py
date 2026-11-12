"""Report figures 4-8: results of the archive run with ComPy 3.

Fig. 4  The eight YV stations: new Vs with its prior, and the ComPy 2.8 values.
Fig. 5  Archive: misfit, change of Vs against 2.8, share of stations resolved with depth.
Fig. 6  Sediment thickness from compliance against GlobSed.
Fig. 7  Sensitivity runs on the benchmark set.
Fig. 8  Calibration: coverage of the 90 % intervals in the synthetic tests.
Also draws the standard figure of every station (--stations) into its folder.

Run: PYTHONPATH=<worktree> python g02_figures_results.py [--stations]
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from compy.plotting.auto import GRID, INK, MUTED, POSTERIOR, PRIOR, REFERENCE, STYLE

OUT = Path(__file__).resolve().parents[1] / "figures"
RUN = Path.home() / "Desktop/ComPy/14_inversion_overhaul"
ARCHIVE = RUN / "r01_archive"
BLUE, ORANGE, AQUA = POSTERIOR, PRIOR, REFERENCE
DEPTHS = (500, 1000, 2000, 3000, 5000)


def table(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def num(rows, key):
    out = []
    for r in rows:
        try:
            out.append(float(r[key]))
        except (ValueError, KeyError, TypeError):
            out.append(np.nan)
    return np.array(out)


def fig4(rows):
    yv = [f"YV.RR{n}" for n in (28, 29, 34, 36, 38, 40, 50, 52)]
    by = {r["sid"]: r for r in rows}
    fig, axes = plt.subplots(
        2,
        4,
        figsize=(13, 8.6),
        sharex=True,
        sharey=True,
        gridspec_kw=dict(wspace=0.08, hspace=0.16),
    )
    for ax, sid in zip(axes.ravel(), yv):
        sta = sid.split(".")[1]
        p = np.genfromtxt(
            ARCHIVE / "YV" / sta / f"{sta}_auto_profile.csv", delimiter=",", names=True
        )
        z = p["depth_m"]
        keep = z >= 10
        ax.fill_betweenx(
            z[keep],
            p["vs_p05"][keep] / 1e3,
            p["vs_p95"][keep] / 1e3,
            color=BLUE,
            alpha=0.18,
            lw=0,
        )
        ax.fill_betweenx(
            z[keep],
            p["vs_p16"][keep] / 1e3,
            p["vs_p84"][keep] / 1e3,
            color=BLUE,
            alpha=0.38,
            lw=0,
        )
        ax.plot(
            p["vs_p50"][keep] / 1e3,
            z[keep],
            color=BLUE,
            lw=1.8,
            label="ComPy 3: median, 16-84 %, 5-95 %",
        )
        for key in ("prior_p16", "prior_p84"):
            ax.plot(
                p[key][keep] / 1e3,
                z[keep],
                color=ORANGE,
                lw=0.9,
                ls="--",
                label="prior, 16-84 %" if key == "prior_p16" else None,
            )
        r = by[sid]
        zs = np.array(DEPTHS, float)
        old = np.array([float(r[f"old_vs_{d}"]) for d in DEPTHS]) / 1e3
        ax.plot(old, zs, "s", color=INK, ms=5, label="ComPy 2.8 median")
        res = p["resolved"] > 0
        if res.any():
            ax.plot(
                [0.02, 0.02],
                [z[res].min(), z[res].max()],
                transform=ax.get_yaxis_transform(),
                color=INK,
                lw=4,
                solid_capstyle="butt",
                label="resolved by the data",
            )
        ax.set_yscale("log")
        ax.set_ylim(9000, 30)
        ax.set_xlim(0, 5)
        ax.set_title(
            f"{sid}: {r['quality']}, {float(r['depth_m']):.0f} m water", fontsize=9
        )
    for ax in axes[1]:
        ax.set_xlabel("Shear velocity (km/s)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Depth below seafloor (m)")
    handles, labels = axes[0, 3].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        fontsize=8,
        bbox_to_anchor=(0.5, 0.0),
    )
    fig.subplots_adjust(bottom=0.11)
    fig.savefig(OUT / "fig04_yv_profiles.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def fig5(rows):
    t1 = [r for r in rows if r["tier"] == "1"]
    ok = [r for r in t1 if r["quality"] != "poor"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(wspace=0.3))
    ax = axes[0]
    new = num(t1, "chi2_median") / num(t1, "n_data")
    old = num(t1, "old_chi2_over_n")
    bins = np.geomspace(0.05, 20, 28)
    ax.hist(
        old[np.isfinite(old)],
        bins,
        color=ORANGE,
        alpha=0.6,
        label="ComPy 2.8",
        edgecolor="white",
    )
    ax.hist(
        new[np.isfinite(new)],
        bins,
        color=BLUE,
        alpha=0.6,
        label="ComPy 3",
        edgecolor="white",
    )
    ax.axvline(1, color=MUTED, lw=0.8, ls=":")
    ax.set(
        xscale="log",
        xlabel="Median chi-square / number of data",
        ylabel="Stations (tier 1)",
        title="Misfit against the stated errors",
    )
    ax.legend(loc="upper right")
    ax = axes[1]
    data = []
    for d in DEPTHS:
        r = num(ok, f"vs_{d}") / num(ok, f"old_vs_{d}")
        data.append(100 * (r[np.isfinite(r)] - 1))
    parts = ax.boxplot(
        data,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color=INK, lw=1.5),
    )
    for patch in parts["boxes"]:
        patch.set(facecolor=BLUE, alpha=0.35, edgecolor=BLUE)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xticks(range(1, len(DEPTHS) + 1), [f"{d / 1000:g}" for d in DEPTHS])
    ax.set(
        xlabel="Depth below seafloor (km)",
        ylabel="Vs of ComPy 3 against ComPy 2.8 (%)",
        title=f"Change of the median Vs ({len(ok)} stations)",
    )
    ax = axes[2]
    zs = np.array([100, 300, 500, 1000, 2000, 3000, 5000])
    share = [np.nanmean(num(ok, f"resolved_{z}")) for z in zs]
    ax.plot(100 * np.array(share), zs, "o-", color=BLUE, lw=2)
    ax.set_yscale("log")
    ax.set_ylim(7000, 70)
    ax.set(
        xlabel="Stations with this depth resolved (%)",
        ylabel="Depth below seafloor (m)",
        title="Where the data constrain Vs",
        xlim=(0, 100),
    )
    fig.savefig(OUT / "fig05_archive.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    return dict(
        chi2_over_n_new=float(np.nanmedian(new)),
        chi2_over_n_old=float(np.nanmedian(old)),
        change_pct={
            str(d): [
                float(np.median(x)),
                float(np.quantile(x, 0.16)),
                float(np.quantile(x, 0.84)),
            ]
            for d, x in zip(DEPTHS, data)
        },
        resolved_share={str(int(z)): float(s) for z, s in zip(zs, share)},
        n=len(ok),
    )


def fig6(rows):
    glob = {
        r["sid"]: r["globsed_m"]
        for r in table(RUN / "reference/globsed_at_stations.csv")
    }
    ok = [
        r
        for r in rows
        if r["quality"] != "poor" and r["tier"] == "1" and glob.get(r["sid"])
    ]
    x = np.array([float(glob[r["sid"]]) for r in ok])
    y, lo, hi = (
        num(ok, "sediment_thickness_median"),
        num(ok, "sediment_thickness_p16"),
        num(ok, "sediment_thickness_p84"),
    )
    thick = np.array([float(r["crust1_sediment"]) >= 500 for r in ok])
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    for sel, color, label in (
        (~thick, BLUE, "soft unit over igneous reference"),
        (thick, ORANGE, "CRUST1.0 sediment >= 500 m: soft unit only"),
    ):
        ax.errorbar(
            np.maximum(x[sel], 3),
            np.maximum(y[sel], 3),
            yerr=[np.maximum(y[sel] - lo[sel], 0), hi[sel] - y[sel]],
            fmt="o",
            ms=4,
            color=color,
            elinewidth=0.6,
            alpha=0.75,
            label=label,
        )
    ax.plot([3, 8000], [3, 8000], color=MUTED, lw=0.8)
    ax.set(
        xscale="log",
        yscale="log",
        xlim=(3, 8000),
        ylim=(3, 8000),
        xlabel="Sediment thickness, GlobSed v3 (m)",
        ylabel="Soft-sediment thickness from compliance (m)",
        title=f"Sediment against an independent map ({len(ok)} stations)",
    )
    ax.legend(loc="upper left", fontsize=7)
    fig.savefig(OUT / "fig06_sediment_globsed.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    from scipy.stats import spearmanr

    thin = ~thick & (x > 0) & (y > 0)
    rho = spearmanr(x[thin], y[thin])
    return dict(
        n=len(ok),
        n_thin=int(thin.sum()),
        spearman_thin=float(rho.statistic),
        p_thin=float(rho.pvalue),
        median_ratio_thin=float(np.median(y[thin] / x[thin])),
        within_factor_3_thin=float(
            np.mean(np.abs(np.log(y[thin] / x[thin])) < np.log(3))
        ),
    )


def changed(tag, base_row):
    """Whether the variant differs from the standard run at this station.

    Stations with thick CRUST1.0 sediment already use a prior width of 0.5, and
    stations without a response correction already have the 4 % tilt allowance.
    """
    if tag == "sigma05":
        return float(base_row["crust1_sediment"]) < 500.0
    if tag == "unchecked":
        return base_row["response"] == "corrected"
    return True


def fig7(rows):
    base = {r["sid"]: r for r in rows}
    variants = (
        ("sigma05", "prior sigma 0.5\ninstead of 0.3"),
        ("corr09", "layer correlation 0.9\ninstead of 0.7"),
        ("unchecked", "no response correction,\ntilt allowance 4 %"),
        ("vp15", "Vp factor sigma 15 %\ninstead of 7 %"),
        ("fc06", "band cut at 0.6 f_c\ninstead of 0.8"),
    )
    found = [
        (tag, label)
        for tag, label in variants
        if (RUN / "r02_sensitivity/tables" / f"results_{tag}.csv").exists()
    ]
    if not found:
        return {}
    fig, axes = plt.subplots(
        1,
        len(found),
        figsize=(2.9 * len(found), 4.2),
        sharey=True,
        gridspec_kw=dict(wspace=0.08),
    )
    numbers = {}
    for ax, (tag, label) in zip(np.atleast_1d(axes), found):
        var = table(RUN / "r02_sensitivity/tables" / f"results_{tag}.csv")
        data, width = [], []
        for d in DEPTHS:
            pairs = [
                (
                    float(v[f"vs_{d}"]),
                    float(base[v["sid"]][f"vs_{d}"]),
                    float(v[f"vs_{d}_p84"]) - float(v[f"vs_{d}_p16"]),
                    float(base[v["sid"]][f"vs_{d}_p84"])
                    - float(base[v["sid"]][f"vs_{d}_p16"]),
                )
                for v in var
                if v["sid"] in base
                and v["quality"] != "poor"
                and base[v["sid"]]["quality"] != "poor"
                and changed(tag, base[v["sid"]])
            ]
            a = np.array(pairs)
            data.append(100 * (a[:, 0] / a[:, 1] - 1))
            width.append(float(np.median(a[:, 2] / a[:, 3])))
        parts = ax.boxplot(
            data,
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            vert=False,
            medianprops=dict(color=INK, lw=1.5),
        )
        for patch in parts["boxes"]:
            patch.set(facecolor=BLUE, alpha=0.35, edgecolor=BLUE)
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.set_title(f"{label}\n({len(data[0])} stations)", fontsize=8)
        ax.set_xlim(-25, 25)
        ax.set_xlabel("Change of median Vs (%)")
        numbers[tag] = dict(
            n=len(data[0]),
            median_change={str(d): float(np.median(x)) for d, x in zip(DEPTHS, data)},
            spread={
                str(d): [float(np.quantile(x, 0.16)), float(np.quantile(x, 0.84))]
                for d, x in zip(DEPTHS, data)
            },
            interval_width_ratio={str(d): w for d, w in zip(DEPTHS, width)},
            largest_change={
                str(d): float(np.max(np.abs(x))) for d, x in zip(DEPTHS, data)
            },
        )
    ax0 = np.atleast_1d(axes)[0]
    ax0.set_yticks(range(1, len(DEPTHS) + 1), [f"{d / 1000:g} km" for d in DEPTHS])
    ax0.set_ylabel("Depth below seafloor")
    ax0.invert_yaxis()
    fig.savefig(OUT / "fig07_sensitivity.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    return numbers


def fig8():
    summary = json.loads((RUN / "e03_sbc/summary.json").read_text())
    cases = (
        ("YV.RR52_prior", "truth from the prior\n(deep ocean)"),
        ("7D.J09B_prior", "truth from the prior\n(shelf)"),
        ("YV.RR52_gradient", "smooth-gradient Earth"),
        ("YV.RR52_prior_tilt4", "slope error 4 %"),
        ("YV.RR52_prior_gain10", "gain error 10 %"),
    )
    by = {s["name"]: s for s in summary}
    fig, axes = plt.subplots(
        1, len(cases), figsize=(13, 4.0), sharey=True, gridspec_kw=dict(wspace=0.08)
    )
    for ax, (name, label) in zip(axes, cases):
        s = by[name]
        z = np.array([float(k) for k in s["vs"]])
        c90 = np.array([v["cover90"] for v in s["vs"].values()])
        c50 = np.array([v["cover50"] for v in s["vs"].values()])
        ax.axvspan(
            100 * s["expected90"][0],
            100 * s["expected90"][1],
            color=BLUE,
            alpha=0.15,
            lw=0,
        )
        ax.axvspan(
            100 * s["expected50"][0],
            100 * s["expected50"][1],
            color=ORANGE,
            alpha=0.15,
            lw=0,
        )
        ax.plot(100 * c90, z, "o-", color=BLUE, lw=1.6, ms=4, label="90 % intervals")
        ax.plot(100 * c50, z, "s-", color=ORANGE, lw=1.6, ms=4, label="50 % intervals")
        ax.set_yscale("log")
        ax.set_ylim(9000, 35)
        ax.set_xlim(25, 102)
        ax.set_title(label, fontsize=9)
        ax.set_xlabel("Intervals holding the truth (%)")
    axes[0].set_ylabel("Depth below seafloor (m)")
    axes[1].legend(loc="center", bbox_to_anchor=(0.5, 0.5), fontsize=7)
    fig.savefig(OUT / "fig08_calibration.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def station_figures():
    from compy.inversion import auto
    from compy.plotting.auto import station_figure

    done = 0
    for path in sorted(ARCHIVE.glob("*/*/*_auto.json")):
        target = path.with_suffix(".png")
        if target.exists():
            continue
        station_figure(auto.load_result(path.with_suffix("")), target)
        done += 1
    print(done, "station figures drawn")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    rows = table(ARCHIVE / "tables/results_auto.csv")
    with matplotlib.rc_context(STYLE):
        fig4(rows)
        numbers = dict(archive=fig5(rows), sediment=fig6(rows), sensitivity=fig7(rows))
        fig8()
    (OUT / "results_numbers.json").write_text(json.dumps(numbers, indent=1))
    print(json.dumps(numbers, indent=1)[:2500])
    if "--stations" in sys.argv:
        station_figures()
