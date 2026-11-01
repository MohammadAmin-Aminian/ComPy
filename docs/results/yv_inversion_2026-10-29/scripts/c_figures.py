"""Stage C: figures and tables of the Vs inversion.

Writes figures/<STA>_inversion.png, the small multiples vs_profiles.png,
sensitivity.png, sediment.png and compliance_fit.png, and
results/tables/inversion_summary.{json,csv}.
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from a_invert import EPOCH, ROOT, build_model, load_data  # noqa: E402

import compy_inversion as ci  # noqa: E402

STATIONS = ("RR28", "RR29", "RR34", "RR36", "RR38", "RR40", "RR50", "RR52")
SLOTS = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#6250d6",
    "#e34948",
)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d3"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 150,
})  # fmt: skip
Z = np.linspace(0.0, 8000.0, 401)
DEPTHS_KM = (0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0)
# Vs-depth descriptions in Aminian et al. (2025, section 5.2), read from the
# text; approximate (depth m, Vs m/s) vertices, not a digitized Figure 9.
GJI = {
    "RR28": [(0, 400), (500, 400), (500, 2700), (1200, 3500), (4300, 3500), (7500, 4000), (8000, 4200)],
    "RR29": [(0, 2000), (4000, 4000), (8000, 4000)],
    "RR34": [(0, 900), (1000, 3500), (3000, 4000), (7500, 4000), (7500, 4200), (8000, 4200)],
    "RR36": [(0, 2000), (2000, 2400), (4200, 4100), (8000, 4100)],
    "RR38": [(0, 2200), (2000, 2200), (4200, 4100), (8000, 4100)],
    "RR40": [(0, 900), (1000, 3400), (3500, 3500), (4000, 4050), (8000, 4050)],
    "RR50": [(0, 2000), (2000, 2400), (4200, 4100), (8000, 4100)],
    "RR52": [(0, 2100), (1000, 2100), (1500, 3800), (3000, 3800), (3500, 4250), (8000, 4250)],
}  # fmt: skip
VARIANT_LABELS = {"main": "Main (Vp/Vs 1.85)", "brocher": "Brocher Vp", "sedscale": "Sediment Vs free",
                  "bruitfm": "gains.csv gain", "full": "Whole record", "moho": "Moho ±1 km",
                  "fixedvp": "Vp, ρ fixed at CRUST1.0", "vpvs": "Crustal Vp/Vs free"}  # fmt: skip


def load(station, variant):
    """Posterior of a saved run, rebuilt with its data and model."""
    d = np.load(ROOT / "results" / station / f"{station}_{variant}.npz")
    spec = dict(main={}, prior={}, brocher=dict(vp_vs=None), sedscale=dict(sediment_scale=(0.5, 2.0)),
                bruitfm={}, full={}, moho=dict(moho_shift=(-1000.0, 1000.0)),
                fixedvp=dict(reference="crust1"), vpvs=dict(vp_vs_range=(1.7, 2.4)))[variant]  # fmt: skip
    model, column = build_model(station, **spec)
    data = ci.ComplianceData(
        d["freqs"], d["eta"], float(d["depth"]), d["covariance"], station=station
    )
    summary = json.loads(
        (ROOT / "results" / station / f"{station}_{variant}.json").read_text()
    )
    post = ci.Posterior(data=data, model=model, prior_only=variant == "prior", theta=d["theta"],
                        alpha=d["alpha"], loglik=d["loglik"], roughness=d["roughness"], pred=d["pred"],
                        acceptance=np.array(summary["acceptance"]),
                        acceptance_dr=np.array(summary["acceptance_dr"]), settings=summary["settings"])  # fmt: skip
    return post, column, summary


def steps(thick, vs, bottom=8000.0):
    """(depth, vs) vertices of a layered profile; the half-space runs to ``bottom``."""
    tops = np.concatenate([[0.0], np.cumsum(thick[:-1])])
    keep = tops < bottom
    z = np.repeat(np.append(tops[keep], bottom), 2)[1:-1]
    return z, np.repeat(np.asarray(vs)[keep], 2)


def step(column):
    """CRUST1.0 Vs as (depth, vs) step vertices to 8 km."""
    m = ci.crust1_model(column)
    return steps(m[:, 0], m[:, 3])


def station_figure(station, post, prior, column, occam, summary):
    fig, axes = plt.subplots(
        1, 4, figsize=(13, 4.6), gridspec_kw=dict(width_ratios=[1.25, 1.2, 0.55, 0.8])
    )
    data = post.data
    ax = axes[0]
    pred = post.pred.reshape(-1, len(data.eta)) * 1e9
    lo, med, hi = np.quantile(pred, [0.05, 0.5, 0.95], axis=0)
    f = data.freqs * 1e3
    ax.fill_between(
        f, lo, hi, color=SLOTS[0], alpha=0.25, lw=0, label="Posterior 90 % prediction"
    )
    ax.plot(f, med, color=SLOTS[0], lw=2, label="Posterior median prediction")
    ax.errorbar(f, data.eta * 1e9, yerr=data.sigma * 1e9, fmt="o", color=INK, ms=4, capsize=2,
                label="Measured (±1σ, gain error included)")  # fmt: skip
    ax.set_xlabel("Frequency [mHz]")
    ax.set_ylabel("Normalized compliance [1/GPa]")
    ax.set_title(
        f"(a) Data fit, χ² = {summary['chi2']['median']:.1f} for N = {len(data.eta)}",
        loc="left",
    )
    ax.legend(fontsize=8, loc="upper left")
    ax = axes[1]
    vs = post.vs_draws(Z, max_draws=6000)
    bins = np.linspace(0, 5000, 101)
    hist = np.array([np.histogram(v, bins)[0] for v in vs.T]).astype(float)
    hist /= hist.max(axis=1, keepdims=True)
    ax.pcolormesh(
        0.5 * (bins[1:] + bins[:-1]),
        Z / 1e3,
        hist,
        cmap="Blues",
        vmin=0,
        vmax=1.15,
        shading="nearest",
        rasterized=True,
    )
    q = np.quantile(vs, [0.05, 0.5, 0.95], axis=0)
    pq = np.quantile(prior.vs_draws(Z, max_draws=6000), [0.05, 0.95], axis=0)
    ax.plot(q[1], Z / 1e3, color=INK, lw=2, label="Posterior median")
    ax.plot(q[0], Z / 1e3, color=INK, lw=0.8, ls="--", label="Posterior 5–95 %")
    ax.plot(q[2], Z / 1e3, color=INK, lw=0.8, ls="--")
    ax.plot(pq[0], Z / 1e3, color=MUTED, lw=0.8, ls=":", label="Prior 5–95 %")
    ax.plot(pq[1], Z / 1e3, color=MUTED, lw=0.8, ls=":")
    cz, cv = step(column)
    ax.plot(cv, cz / 1e3, color=SLOTS[1], lw=1.4, label="CRUST1.0")
    g = np.array(GJI[station])
    ax.plot(
        g[:, 1],
        g[:, 0] / 1e3,
        color=SLOTS[2],
        lw=1.4,
        ls="--",
        label="GJI 2025 (text, approx.)",
    )
    oz, ov = steps(np.array(occam["thickness"]), np.array(occam["vs"]))
    ax.plot(ov, oz / 1e3, color=SLOTS[6], lw=1.2, label="Occam (from CRUST1.0)")
    ax.set_ylim(8, 0)
    ax.set_xlim(0, 5000)
    ax.set_xlabel("Vs [m/s]")
    ax.set_ylabel("Depth below seafloor [km]")
    ax.set_title("(b) Vs posterior", loc="left")
    ax.legend(fontsize=7.5, loc="lower left")
    ax = axes[2]
    median = np.median(post.theta.reshape(-1, post.theta.shape[-1]), axis=0)
    zk, kern = ci.vs_kernels(data.depth, data.freqs[[0, -1]], post.model.layers(median),
                             dz=50.0, max_depth=30000.0)  # fmt: skip
    for k, color, label in zip(
        kern, (SLOTS[0], SLOTS[1]), (f"{f[0]:.1f} mHz", f"{f[-1]:.1f} mHz")
    ):
        w = np.abs(k) * np.gradient(zk)
        ax.plot(np.cumsum(w) / w.sum(), zk / 1e3, color=color, lw=1.4, label=label)
    sens = {q: ci.sensitivity_depth(zk, kern[0], q) for q in (0.5, 0.9)}
    for q, ls in ((0.5, "--"), (0.9, ":")):
        ax.axhline(sens[q] / 1e3, color=INK, lw=0.9, ls=ls)
        ax.text(0.02, sens[q] / 1e3 - 0.12, f"{int(q * 100)} % above {sens[q] / 1e3:.1f} km",
                color=INK, fontsize=7.5)  # fmt: skip
    ax.set_ylim(8, 0)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Cumulative fraction of\n|d ln η / d ln Vs|")
    ax.set_title("(c) Sensitivity", loc="left")
    ax.legend(fontsize=7.5, loc="lower left")
    ax = axes[3]
    i = post.names.index("sediment_thickness")
    h = post.theta[:, :, i].ravel()
    ax.hist(
        h,
        bins=np.linspace(0, 1000, 41),
        color=SLOTS[0],
        alpha=0.8,
        density=True,
        label="Posterior",
    )
    ax.axhline(1 / 1000, color=MUTED, ls=":", label="Prior (uniform)")
    ax.axvline(column["sediment_thickness"], color=SLOTS[1], lw=1.4, label="CRUST1.0")
    ax.set_xlabel("Sediment thickness [m]")
    ax.set_yticks([])
    ax.set_title("(d) Sediment", loc="left")
    ax.legend(fontsize=7.5)
    epoch = " (stable epoch)" if station in EPOCH else ""
    fig.suptitle(f"YV.{station}{epoch}: water {data.depth:.0f} m, band {f[0]:.1f}–{f[-1]:.1f} mHz, "
                 f"max R-hat {summary['max_rhat']:.3f}, min ESS {summary['min_ess']:.0f}", x=0.01, ha="left")  # fmt: skip
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / f"{station}_inversion.png")
    plt.close(fig)
    return sens


def small_multiples(store):
    fig, axes = plt.subplots(2, 4, figsize=(13, 8), sharey=True, sharex=True)
    for ax, station in zip(axes.ravel(), STATIONS):
        s = store[station]
        q = s["vs_q"]
        ax.fill_betweenx(
            Z / 1e3,
            q[0],
            q[4],
            color=SLOTS[0],
            alpha=0.18,
            lw=0,
            label="Posterior 90 %",
        )
        ax.fill_betweenx(
            Z / 1e3,
            q[1],
            q[3],
            color=SLOTS[0],
            alpha=0.35,
            lw=0,
            label="Posterior 68 %",
        )
        ax.plot(q[2], Z / 1e3, color=INK, lw=1.8, label="Posterior median")
        ax.plot(
            s["prior_q"][0], Z / 1e3, color=MUTED, lw=0.8, ls=":", label="Prior 90 %"
        )
        ax.plot(s["prior_q"][1], Z / 1e3, color=MUTED, lw=0.8, ls=":")
        cz, cv = step(s["column"])
        ax.plot(cv, cz / 1e3, color=SLOTS[1], lw=1.2, label="CRUST1.0")
        g = np.array(GJI[station])
        ax.plot(
            g[:, 1],
            g[:, 0] / 1e3,
            color=SLOTS[2],
            lw=1.2,
            ls="--",
            label="GJI 2025 (text, approx.)",
        )
        if s["resolved_m"] < 8000:
            ax.axhspan(s["resolved_m"] / 1e3, 8, color=GRID, alpha=0.45, lw=0)
            ax.text(150, 7.75, "< 10 % of sensitivity", color=MUTED, fontsize=7.5)
        flag = (
            "  ⚠ misfit"
            if s["summary"]["chi2"]["median"] > 2 * s["summary"]["n_data"]
            else ""
        )
        ax.set_title(f"{station}{flag}", loc="left")
        ax.set_xlim(0, 5000)
        ax.set_ylim(8, 0)
    for ax in axes[:, 0]:
        ax.set_ylabel("Depth below seafloor [km]")
    for ax in axes[1]:
        ax.set_xlabel("Vs [m/s]")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        ncol=6,
        loc="upper right",
        fontsize=8,
        bbox_to_anchor=(1.0, 0.965),
    )
    fig.suptitle("Vs posteriors (main run). Grey: below the depth holding 90 % of the lowest-frequency "
                 "Vs sensitivity (deeper than 8 km at RR29 and RR50).", x=0.01, ha="left")  # fmt: skip
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(ROOT / "figures" / "vs_profiles.png")
    plt.close(fig)


def sensitivity_figure():
    fig, axes = plt.subplots(2, 4, figsize=(13, 8), sharey=True, sharex=True)
    for ax, station in zip(axes.ravel(), STATIONS):
        for color, variant in zip(SLOTS, VARIANT_LABELS):
            path = ROOT / "results" / station / f"{station}_{variant}.npz"
            if not path.exists():
                continue
            post, _, _ = load(station, variant)
            med = np.median(post.vs_draws(Z, max_draws=3000), axis=0)
            ax.plot(med, Z / 1e3, color=color, lw=2.2 if variant == "main" else 1.2,
                    label=VARIANT_LABELS[variant])  # fmt: skip
        ax.set_title(station, loc="left")
        ax.set_xlim(0, 5000)
        ax.set_ylim(8, 0)
    for ax in axes[:, 0]:
        ax.set_ylabel("Depth below seafloor [km]")
    for ax in axes[1]:
        ax.set_xlabel("Median Vs [m/s]")
    handles = [
        plt.Line2D([], [], color=c, lw=2.2 if v == "main" else 1.2)
        for c, v in zip(SLOTS, VARIANT_LABELS)
    ]
    fig.legend(
        handles, list(VARIANT_LABELS.values()), ncol=4, loc="upper right", fontsize=8.5
    )
    fig.suptitle(
        "Sensitivity of the posterior median to modelling choices", x=0.01, ha="left"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(ROOT / "figures" / "sensitivity.png")
    plt.close(fig)


def sediment_figure(store):
    fig, ax = plt.subplots(figsize=(8, 3.6))
    x = np.arange(len(STATIONS))
    for k, (variant, color, dx) in enumerate(
        (("main", SLOTS[0], -0.12), ("sedscale", SLOTS[6], 0.12))
    ):
        for j, station in enumerate(STATIONS):
            p = store[station]["sediment"][variant]
            ax.plot([j + dx] * 2, [p["p05"], p["p95"]], color=color, lw=1.5)
            ax.plot(j + dx, p["median"], "o", color=color, ms=6,
                    label=(VARIANT_LABELS[variant] if j == 0 else None))  # fmt: skip
    ax.plot(x, [store[s]["column"]["sediment_thickness"] for s in STATIONS], "s", color=SLOTS[1],
            ms=6, mfc="none", label="CRUST1.0")  # fmt: skip
    gji = {"RR28": 145, "RR29": 110, "RR34": 400}
    ax.plot([STATIONS.index(s) for s in gji], list(gji.values()), "D", color=SLOTS[2], ms=6, mfc="none",
            label="GlobSed, as quoted by GJI 2025")  # fmt: skip
    ax.set_xticks(x, STATIONS)
    ax.set_ylabel("Sediment thickness [m]")
    ax.set_title(
        "Sediment (or very slow uppermost layer) thickness: median and 90 % interval",
        loc="left",
    )
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "sediment.png")
    plt.close(fig)


def compliance_figure(store):
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.5), sharex=True)
    for ax, station in zip(axes.ravel(), STATIONS):
        post = store[station]["post"]
        d = post.data
        pred = post.pred.reshape(-1, len(d.eta)) * 1e9
        lo, med, hi = np.quantile(pred, [0.05, 0.5, 0.95], axis=0)
        f = d.freqs * 1e3
        ax.fill_between(
            f,
            lo,
            hi,
            color=SLOTS[0],
            alpha=0.25,
            lw=0,
            label="Posterior 90 % prediction",
        )
        ax.plot(f, med, color=SLOTS[0], lw=1.8, label="Posterior median prediction")
        ax.errorbar(f, d.eta * 1e9, yerr=d.sigma * 1e9, fmt="o", color=INK, ms=3.5, capsize=2,
                    label="Measured ±1σ")  # fmt: skip
        ax.set_title(f"{station}: χ² {store[station]['summary']['chi2']['median']:.1f} / N {len(d.eta)}",
                     loc="left")  # fmt: skip
    for ax in axes[1]:
        ax.set_xlabel("Frequency [mHz]")
    for ax in axes[:, 0]:
        ax.set_ylabel("Compliance [1/GPa]")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "compliance_fit.png")
    plt.close(fig)


def main():
    (ROOT / "figures").mkdir(exist_ok=True)
    (ROOT / "results" / "tables").mkdir(exist_ok=True)
    store, rows = {}, []
    for station in STATIONS:
        post, column, summary = load(station, "main")
        prior, _, _ = load(station, "prior")
        occam = json.loads(
            (ROOT / "results" / station / f"{station}_occam.json").read_text()
        )
        sens = station_figure(station, post, prior, column, occam, summary)
        zr = sens[0.9]
        vs = post.vs_draws(Z, max_draws=6000)
        sed = {}
        for variant in ("main", "sedscale"):
            s = json.loads(
                (ROOT / "results" / station / f"{station}_{variant}.json").read_text()
            )
            sed[variant] = s["parameters"]["sediment_thickness"]
        store[station] = dict(post=post, column=column, summary=summary, resolved_m=zr, sediment=sed,
                              vs_q=np.quantile(vs, [0.05, 0.16, 0.5, 0.84, 0.95], axis=0),
                              prior_q=np.quantile(prior.vs_draws(Z, max_draws=6000), [0.05, 0.95], axis=0))  # fmt: skip
        at = np.quantile(
            post.vs_draws(np.array(DEPTHS_KM) * 1e3, max_draws=6000),
            [0.05, 0.5, 0.95],
            axis=0,
        )
        _, info = load_data(station)
        row = dict(station=station, epoch=station in EPOCH, water_depth_m=post.data.depth,
                   band_mhz=f"{post.data.freqs[0] * 1e3:.1f}-{post.data.freqs[-1] * 1e3:.1f}",
                   n_data=len(post.data.eta), months=info["months"],
                   chi2_median=round(summary["chi2"]["median"], 2), max_rhat=round(summary["max_rhat"], 3),
                   min_ess=round(summary["min_ess"]), sens50_km=round(sens[0.5] / 1e3, 2),
                   sens90_km=round(sens[0.9] / 1e3, 2),
                   sediment_m=f"{sed['main']['median']:.0f} ({sed['main']['p05']:.0f}-{sed['main']['p95']:.0f})",
                   crust1_sediment_m=column["sediment_thickness"], moho_m=column["moho_depth"],
                   occam_chi2=round(occam["chi2"], 2))  # fmt: skip
        for k, z in enumerate(DEPTHS_KM):
            row[f"vs_{z:g}km"] = f"{at[1, k]:.0f} ({at[0, k]:.0f}-{at[2, k]:.0f})"
        rows.append(row)
        print(station, row, flush=True)
    small_multiples(store)
    sensitivity_figure()
    sediment_figure(store)
    compliance_figure(store)
    tables = ROOT / "results" / "tables"
    (tables / "inversion_summary.json").write_text(
        json.dumps(rows, indent=1, default=float)
    )
    with open(tables / "inversion_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
