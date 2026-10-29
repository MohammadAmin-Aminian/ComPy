"""Stage I: figures in the style of Aminian et al. (2025) and the thesis.

From the 'final' runs (free crustal Vp/Vs, gain sampled per epoch, after-jump
data joined):

  paper_fig9_heatmaps.png    Fig. 9: Vs posterior heat maps, median (red) and
                             CRUST1.0 (white dashed), custom viridis map with
                             each depth row normalized to its maximum.
  paper_fig_misfit.png       Thesis Fig. 3.19b / supplement Fig. 5: misfit
                             against iteration for every chain and station.
  <STA>_paper.png            Thesis Fig. 3.19: (a) measured compliance and
                             simulated compliance of posterior models (green),
                             (b) misfit function, (c) Vs heat map.
  paper_fig8_compliance.png  Fig. 8: log compliance by tectonic setting.
  paper_fig10_profiles.png   Fig. 10: medians by setting and serpentinization.

Differences from the published plotting code, on purpose:
  * The depth axis is in true metres. The published heat maps have 100-m
    rows (a 1-m profile averaged over 100 samples) but tick labels placing
    row 12 at -2000 m, i.e. 167-m rows, while the CRUST1.0 lines use 167-m
    rows. Their posterior depths therefore read 5/3 too deep.
  * No misfit threshold: every post-burn-in draw is a posterior sample; the
    misfit panel shows chi2/N = 1, the expected value, and the burn-in end.
  * Serpentinization uses the relation of the Fig. 10 caption (Vs 4.2 km/s =
    0 %, 2.2 km/s = 100 %). The published code, (9000 - 2 Vs) / 30, gives 0 %
    at 4.5 km/s and 100 % at 3.0 km/s.
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from a_invert import EPOCH, ROOT, build_model  # noqa: E402

import compy_inversion as ci  # noqa: E402

STATIONS = ("RR28", "RR29", "RR34", "RR36", "RR38", "RR40", "RR50", "RR52")
TITLES = {  # Aminian et al. (2025), Fig. 9
    "RR28": "RR28 (Somali 22°S - 53°E)", "RR29": "RR29 (Somali 24°S - 51°E)",
    "RR34": "RR34 (Somali 32°S - 52°E)", "RR36": "RR36 (SWIR 33°S - 55°E)",
    "RR38": "RR38 (SWIR 30°S - 59°E)", "RR40": "RR40 (SWIR 28°S - 63°E)",
    "RR50": "RR50 (RTJ 25°S - 70°E)", "RR52": "RR52 (CIR 20°S - 68°E)",
}  # fmt: skip
GROUPS = {"Somali Plate": ("RR28", "RR29", "RR34"), "SWIR": ("RR36", "RR38", "RR40"),
          "RTJ and CIR": ("RR50", "RR52")}  # fmt: skip
STYLES = ("-", "--", "-.", ":")
VIRIDIS = plt.get_cmap("viridis")
CMAP = LinearSegmentedColormap.from_list(
    "aminian2025", list(zip([0, 0.4, 1], [VIRIDIS(0.0), VIRIDIS(0.5), VIRIDIS(1.0)])))  # fmt: skip
Z = np.arange(0.0, 8000.0 + 1, 50.0)
VS_BINS = np.linspace(0, 5000, 101)
plt.rcParams.update(
    {"font.family": "serif", "font.size": 12, "axes.grid": True, "grid.alpha": 0.35}
)


def load(station, variant="final"):
    d = np.load(ROOT / "results" / station / f"{station}_{variant}.npz")
    model, column = build_model(station, vp_vs_range=(1.7, 2.4))
    data = ci.ComplianceData(d["freqs"], d["eta"], float(d["depth"]), d["covariance"], station=station,
                             groups=d["groups"], group_names=tuple(d["group_names"]))  # fmt: skip
    names = list(d["names"])
    post = ci.Posterior(data=data, model=model, prior_only=False, theta=d["theta"], alpha=d["alpha"],
                        loglik=d["loglik"], roughness=d["roughness"], pred=d["pred"],
                        acceptance=np.zeros(1), acceptance_dr=np.zeros(1), settings={},
                        trace_chi2=d["trace_chi2"], names_extra=names[len(model.names):])  # fmt: skip
    return post, column


def heat(post):
    vs = post.vs_draws(Z, max_draws=8000)
    h = np.array([np.histogram(v, VS_BINS)[0] for v in vs.T], float)
    h /= np.maximum(h.max(axis=1, keepdims=True), 1)
    return h, np.median(vs, axis=0)


def crust1_steps(column):
    m = ci.crust1_model(column)
    tops = np.concatenate([[0.0], np.cumsum(m[:-1, 0])])
    z = np.repeat(np.append(tops, 8000.0), 2)[1:-1]
    return np.repeat(m[:, 3], 2), z


def crust1_vs_at(column, z):
    m = ci.crust1_model(column)
    tops = np.concatenate([[0.0], np.cumsum(m[:-1, 0])])
    return m[np.searchsorted(tops, z, side="right") - 1, 3]


def heat_axis(ax, post, column, cbar_label=True, legend=False):
    h, med = heat(post)
    im = ax.imshow(h, extent=(0, 5000, -Z[-1], 0), aspect="auto", cmap=CMAP, vmin=0, vmax=1,
                   origin="upper", interpolation="nearest")  # fmt: skip
    cv, cz = crust1_steps(column)
    ax.plot(cv, -cz, color="white", ls="--", lw=2, label="Crust1")
    ax.plot(med, -Z, color="red", lw=2.2, label="Median")
    ax.set_xlim(0, 5000)
    ax.set_ylim(-8000, 0)
    ax.set_xticks(range(0, 5001, 1000))
    ax.grid(color="white", alpha=0.25)
    cb = plt.colorbar(im, ax=ax, fraction=0.06, pad=0.02)
    if cbar_label:
        cb.set_label("Probability")
    if legend:
        ax.legend(
            loc="lower left", facecolor="0.15", labelcolor="white", framealpha=0.6
        )
    return med


def fig9(store):
    fig, axes = plt.subplots(2, 4, figsize=(22, 12), sharey=True)
    for k, (ax, s) in enumerate(zip(axes.ravel(), STATIONS)):
        post, column = store[s]
        heat_axis(ax, post, column, cbar_label=k in (3, 7), legend=k == 7)
        ax.set_title(TITLES[s] + (" *" if s in EPOCH else ""))
        ax.text(0.04, 0.97, f"({'abcdefgh'[k]})", transform=ax.transAxes, color="white",
                fontweight="bold", va="top", fontsize=14)  # fmt: skip
        if k % 4 == 0:
            ax.set_ylabel("Depth [m]")
        if k >= 4:
            ax.set_xlabel("Shear Velocity [m/s]")
    fig.suptitle("Vs posterior (ComPy 2.3, free crustal Vp/Vs, gain sampled; * stable and after-jump epochs "
                 "inverted together)", fontsize=13)  # fmt: skip
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "paper_fig9_heatmaps.png", dpi=150)
    plt.close(fig)


def misfit_axis(ax, post, title=None):
    n = len(post.data.eta)
    burn = 30000
    it = np.arange(1, post.trace_chi2.shape[1] + 1)
    keep = (
        np.unique(np.geomspace(1, len(it), 3000).astype(int)) - 1
    )  # even on a log axis
    for c, trace in enumerate(post.trace_chi2):
        ax.plot(
            it[keep],
            trace[keep] / n,
            lw=0.6,
            alpha=0.8,
            label="Misfit" if c == 0 else None,
        )
    ax.set_xscale("log")
    ax.axvline(burn, color="red", ls="--", lw=1.5, label="Burn-in end")
    ax.axhline(1.0, color="red", ls=":", lw=1.5, label="χ²/N = 1")
    ax.set_yscale("log")
    ax.set_xlabel("Iteration [i]")
    ax.set_ylabel("χ² / N")
    if title:
        ax.set_title(title)


def data_axis(ax, post):
    d = post.data
    pred = post.pred.reshape(-1, len(d.eta))
    pick = np.random.default_rng(0).choice(len(pred), 300, replace=False)
    markers = ("o", "^")
    for g, name in enumerate(d.group_names):
        sel = d.groups == g
        f = d.freqs[sel]
        order = np.argsort(f)
        for i in pick:
            ax.plot(
                f[order], pred[i, sel][order] * 1e9, color="green", lw=0.25, alpha=0.4
            )
        label = {
            "stable": "Measured (stable epoch)",
            "after": "Measured (after jump)",
        }.get(name, "Measured")
        ax.errorbar(f[order], d.eta[sel][order] * 1e9, yerr=d.sigma[sel][order] * 1e9, fmt=markers[g],
                    color="black" if g == 0 else "0.45", ms=5, capsize=4, lw=1.5, label=label)  # fmt: skip
    ax.plot([], [], color="green", lw=1, label="Simulated (posterior models)")
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Normalized compliance [1/GPa]")
    ax.legend(loc="upper left", fontsize=10)


def station_figure(s, post, column):
    fig, axes = plt.subplots(
        1, 3, figsize=(20, 6), gridspec_kw=dict(width_ratios=[1, 1.1, 0.9])
    )
    data_axis(axes[0], post)
    axes[0].set_title("(a) Measured and simulated compliance")
    misfit_axis(axes[1], post, "(b) Misfit function")
    axes[1].legend(loc="upper right", fontsize=10)
    heat_axis(axes[2], post, column, legend=True)
    axes[2].set_title("(c) " + TITLES[s])
    axes[2].set_xlabel("Shear Velocity [m/s]")
    axes[2].set_ylabel("Depth [m]")
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / f"{s}_paper.png", dpi=130)
    plt.close(fig)


def misfit_grid(store):
    fig, axes = plt.subplots(2, 4, figsize=(22, 9), sharex=True)
    for ax, s in zip(axes.ravel(), STATIONS):
        misfit_axis(ax, store[s][0], TITLES[s])
    axes[0, 0].legend(loc="upper right", fontsize=10)
    fig.suptitle(
        "Misfit function over iterations, all chains (burn-in: first 30 000 iterations)",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "paper_fig_misfit.png", dpi=130)
    plt.close(fig)


def fig8(store):
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5), sharey=True)
    for ax, (setting, stations) in zip(axes, GROUPS.items()):
        for ls, s in zip(STYLES, stations):
            d = store[s][0].data
            sel = d.groups == 0
            f, eta, sig = d.freqs[sel], d.eta[sel] * 1e9, d.sigma[sel] * 1e9
            line = ax.plot(f, eta, ls=ls, lw=2.5, label=s)[0]
            ax.fill_between(
                f, eta - sig, eta + sig, color=line.get_color(), alpha=0.25, lw=0
            )
        ax.set_yscale("log")
        ax.set_xlim(0.005, 0.02)
        ax.set_title(setting)
        ax.set_xlabel("Frequency [Hz]")
        ax.legend(loc="lower right")
    axes[0].set_ylabel("Log Compliance [1/GPa]")
    fig.suptitle(
        "Normalized compliance (stable epoch; shading ±1σ including the 5 % gain error)",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "paper_fig8_compliance.png", dpi=130)
    plt.close(fig)


def fig10(store):
    fig, axes = plt.subplots(1, 3, figsize=(20, 9), sharey=True)
    for ax, stations, title in ((axes[0], GROUPS["Somali Plate"], "Oceanic Plate"),
                                (axes[1], GROUPS["SWIR"] + GROUPS["RTJ and CIR"], "Mid-Ocean Ridge")):  # fmt: skip
        for s in stations:
            post, column = store[s]
            ax.plot(
                np.median(post.vs_draws(Z, 6000), axis=0),
                -Z,
                lw=2.2,
                label=f"{s} {TITLES[s][5:]}",
            )
        mean_col = np.mean([crust1_vs_at(store[s][1], Z) for s in stations], axis=0)
        ax.plot(mean_col, -Z, color="black", ls="--", lw=2, label="Avg Crust1")
        ax.set_title(title)
        ax.set_xlabel("Shear Velocity [m/s]")
        ax.set_xlim(0, 5000)
        ax.legend(loc="lower left", fontsize=10)
    axes[1].axvline(4500, color="grey", ls="--", lw=2)
    axes[1].text(4600, -2500, "Fresh peridotite", rotation=90, color="grey")
    axes[0].set_ylabel("Depth [m]")
    ax = axes[2]
    for s in GROUPS["SWIR"] + GROUPS["RTJ and CIR"]:
        med = np.median(store[s][0].vs_draws(Z, 6000), axis=0)
        ax.plot(np.clip((4200 - med) / 2000 * 100, 0, 100), -Z, lw=2.2, label=s)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Serpentinization (%)")
    ax.set_title("Serpentinization if all low Vs were serpentinite")
    ax.legend(loc="lower right", fontsize=10)
    axes[0].set_ylim(-8000, 0)
    fig.suptitle("Posterior median Vs; serpentinization from Vs 4.2 km/s = 0 % to 2.2 km/s = 100 % "
                 "(Miller & Christensen 1997), as in the Fig. 10 caption", fontsize=13)  # fmt: skip
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "paper_fig10_profiles.png", dpi=130)
    plt.close(fig)


def main():
    variant = sys.argv[1] if len(sys.argv) > 1 else "final"
    store = {s: load(s, variant) for s in STATIONS}
    fig9(store)
    misfit_grid(store)
    fig8(store)
    fig10(store)
    for s in STATIONS:
        station_figure(s, *store[s])
    print("done")


if __name__ == "__main__":
    main()
