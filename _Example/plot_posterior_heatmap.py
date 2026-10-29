"""Posterior Vs-depth heat map, data fit, marginals and chain traces.

Usage: python _Example/plot_posterior_heatmap.py <run directory of converged_inversion.py>
"""

import sys
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

D = sys.argv[1]
z = np.load(f"{D}/chains.npz")
S = json.load(open(f"{D}/summary.json"))
models, preds, misfit = z["models"], z["preds"], z["misfit"]
truth, start, f, obs, clean, unc = (
    z[k] for k in ("truth", "start", "f", "obs", "clean", "unc")
)
BURN, ADAPT = int(z["burnin"]), int(z["adapt"])
C, _, _, N = models.shape
P = S["parameters"]

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.edgecolor": "#52514e",
        "axes.linewidth": 0.8,
        "xtick.color": "#52514e",
        "ytick.color": "#52514e",
        "axes.labelcolor": "#0b0b0b",
        "xtick.direction": "out",
        "ytick.direction": "out",
    }
)
CHAIN = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # fixed categorical order
MED, CI, TRUE, START_C, PRIOR = "#e34948", "#ff9f40", "#ffffff", "#c3c2b7", "#7fd3ff"

# ---- retained, thinned posterior models --------------------------------------
thin = 5
R = models[:, :, :, BURN::thin]  # chains x layers x 4 x draws
R = R.transpose(0, 3, 1, 2).reshape(-1, models.shape[1], 4)  # draws x layers x 4
iface1 = R[:, 0, 0]  # depth of first interface
iface2 = R[:, 0, 0] + R[:, 1, 0]  # 1100 m, conserved


def profiles(zgrid):
    """Vs(z) for every retained draw on a depth grid (layer above z belongs to z)."""
    v = np.where(
        zgrid[None, :] < iface1[:, None],
        R[:, 0, 3][:, None],
        np.where(
            zgrid[None, :] < iface2[:, None], R[:, 1, 3][:, None], R[:, 2, 3][:, None]
        ),
    )
    return v


def density(zgrid, vbins):
    v = profiles(zgrid)
    nb = len(vbins) - 1
    idx = np.clip(np.digitize(v, vbins) - 1, 0, nb - 1)
    rows = np.repeat(np.arange(len(zgrid))[None, :], len(v), 0)
    H = np.bincount((rows * nb + idx).ravel(), minlength=len(zgrid) * nb).reshape(
        len(zgrid), nb
    )
    H = H / H.max(axis=1, keepdims=True)  # per-depth normalized PDF
    return H, np.median(v, 0), np.quantile(v, 0.05, 0), np.quantile(v, 0.95, 0)


def step_profile(m, zmax):
    zz = [0.0]
    vv = [m[0, 3]]
    top = 0.0
    for k in range(len(m) - 1):
        top += m[k, 0]
        zz += [top, top]
        vv += [m[k, 3], m[k + 1, 3]]
    zz.append(zmax)
    vv.append(m[-1, 3])
    return np.array(vv), np.array(zz)


def heat_panel(ax, zmax, dz, dv, title, prior=False):
    zgrid = np.arange(dz / 2, zmax, dz)
    vbins = np.arange(0, 5000 + dv, dv)
    H, med, lo, hi = density(zgrid, vbins)
    im = ax.pcolormesh(
        vbins,
        np.arange(0, zmax + dz / 2, dz)[: len(zgrid) + 1],
        H,
        cmap="viridis",
        vmin=0,
        vmax=1,
        shading="flat",
        rasterized=True,
    )
    ax.plot(lo, zgrid, color=CI, lw=1.1, ls=(0, (4, 2)))
    ax.plot(hi, zgrid, color=CI, lw=1.1, ls=(0, (4, 2)))
    ax.plot(med, zgrid, color=MED, lw=2.4, ls=(0, (5, 2.5)))
    tv, tz = step_profile(truth, zmax)
    ax.plot(tv, tz, color=TRUE, lw=1.6)
    sv, sz = step_profile(start, zmax)
    ax.plot(sv, sz, color=START_C, lw=1.2, ls="-.")
    if prior:
        top = truth[0, 0]
        for b in (0.15 * start[0, 3], 1.2 * start[0, 3]):
            ax.plot([b, b], [0, top], color=PRIOR, lw=1.2, ls=":")
    ax.set_xlim(0, 5000)
    ax.set_ylim(zmax, 0)
    ax.set_xlabel("Shear-wave velocity, Vs (m/s)")
    ax.set_ylabel("Depth below seafloor (m)")
    ax.set_title(title, loc="left")
    ax.grid(color="white", alpha=0.18, lw=0.6)
    return im


fig = plt.figure(figsize=(15, 9.2), dpi=170)
gs = GridSpec(
    4,
    5,
    figure=fig,
    width_ratios=[1.05, 1.05, 0.045, 0.22, 1.0],
    height_ratios=[1.5, 1, 1, 1],
    wspace=0.34,
    hspace=0.75,
)
ax_full = fig.add_subplot(gs[:, 0])
ax_zoom = fig.add_subplot(gs[:, 1])
im = heat_panel(ax_full, 1600, 4, 100, "a) Posterior Vs–depth density")
heat_panel(ax_zoom, 250, 1, 100, "b) Zoom: top 250 m (sediment layer)", prior=True)
cb = fig.colorbar(im, cax=fig.add_subplot(gs[:, 2]))
cb.set_label("Normalized probability density (per depth)")
handles = [
    Line2D([], [], color=TRUE, lw=1.6, label="True model"),
    Line2D([], [], color=MED, lw=2.4, ls=(0, (5, 2.5)), label="Posterior median"),
    Line2D([], [], color=CI, lw=1.1, ls=(0, (4, 2)), label="5–95 % credible bounds"),
    Line2D([], [], color=START_C, lw=1.2, ls="-.", label="Prior-defining model"),
    Line2D([], [], color=PRIOR, lw=1.2, ls=":", label="Vs₁ prior bounds"),
]
leg = ax_full.legend(
    handles=handles[:4],
    loc="lower left",
    fontsize=8.5,
    framealpha=0.85,
    facecolor="#1a1a19",
    edgecolor="none",
    labelcolor="white",
)
ax_zoom.legend(
    handles=handles,
    loc="lower left",
    fontsize=8.5,
    framealpha=0.85,
    facecolor="#1a1a19",
    edgecolor="none",
    labelcolor="white",
)

# ---- c) compliance fit --------------------------------------------------------
ax = fig.add_subplot(gs[0, 4])
ret = preds[:, BURN::thin].reshape(-1, len(f))
lo, md, hi = np.quantile(ret, [0.05, 0.5, 0.95], axis=0)
fm = f * 1e3
ax.fill_between(
    fm,
    lo * 1e11,
    hi * 1e11,
    color="#2a78d6",
    alpha=0.25,
    lw=0,
    label="Posterior predictive 5–95 %",
)
ax.plot(fm, md * 1e11, color="#2a78d6", lw=2, label="Posterior predictive median")
ax.plot(fm, clean * 1e11, color="#0b0b0b", lw=1.2, ls=":", label="True compliance")
ax.errorbar(
    fm,
    obs * 1e11,
    yerr=unc * 1e11,
    fmt="o",
    ms=4.5,
    color="#0b0b0b",
    mfc="white",
    mew=1.2,
    elinewidth=1,
    capsize=2.5,
    label="Observed ±1σ (3 %)",
)
ax.set_xlabel("Frequency (mHz)")
ax.set_ylabel("Normalized compliance\n(×10⁻¹¹ Pa⁻¹)")
ax.set_title("c) Data fit (water depth 4000 m)", loc="left")
ax.grid(color="#e6e5e0", lw=0.6)
ax.set_axisbelow(True)
ax.legend(fontsize=7.8, frameon=False, loc="upper left")


# ---- d–f) marginals ------------------------------------------------------------
def marginal(ax, key, values, true, prior_rng, label, title, bins):
    v = values[:, BURN:].ravel()
    if prior_rng is not None:
        ax.axvspan(*prior_rng, color="#e6e5e0", lw=0, label="Prior support")
    ax.hist(
        v,
        bins=bins,
        density=True,
        color="#2a78d6",
        alpha=0.85,
        edgecolor="white",
        linewidth=0.4,
        label="Posterior",
    )
    p = P[key]
    for q in (p["p05"], p["p95"]):
        ax.axvline(q, color="#52514e", lw=1, ls=(0, (4, 2)))
    ax.axvline(true, color="#e34948", lw=2, label="Truth")
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel(label)
    ax.set_title(
        f"{title}   R̂ = {p['rhat']:.3f}, ESS = {p['ess_bulk']:.0f}",
        loc="left",
        fontsize=9.5,
    )


ax_d = fig.add_subplot(gs[1, 4])
marginal(
    ax_d,
    "Vs1",
    models[:, 0, 3],
    1500,
    (0.15 * 1800, 1.2 * 1800),
    "Vs₁ (m/s)",
    "d) Layer 1 Vs",
    np.arange(250, 2200, 50),
)
marginal(
    fig.add_subplot(gs[2, 4]),
    "Vs2",
    models[:, 1, 3],
    3500,
    (0.9 * 3500, 1.1 * 3500),
    "Vs₂ (m/s)",
    "e) Layer 2 Vs",
    np.arange(3140, 3860, 20),
)
ax_f = fig.add_subplot(gs[3, 4])
marginal(
    ax_f,
    "h1",
    models[:, 0, 0],
    100,
    None,
    "h₁ (m)",
    "f) Layer 1 thickness",
    np.arange(60, 150, 2),
)
from matplotlib.patches import Patch

ax_f.legend(
    handles=[
        Patch(color="#2a78d6", alpha=0.85, label="Posterior"),
        Patch(color="#e6e5e0", label="Prior support"),
        Line2D([], [], color="#e34948", lw=2, label="Truth"),
        Line2D([], [], color="#52514e", lw=1, ls=(0, (4, 2)), label="5–95 %"),
    ],
    fontsize=7.5,
    frameon=False,
    loc="upper right",
    handlelength=1.6,
)

fig.suptitle(
    "ComPy compliance inversion — synthetic test, 4 converged chains "
    f"({C} × {N:,} iterations, burn-in {BURN:,})",
    fontsize=13,
    x=0.02,
    ha="left",
)
fig.savefig(f"{D}/posterior_heatmap.png", bbox_inches="tight")
plt.close(fig)

# ---- convergence figure ---------------------------------------------------------
fig, axs = plt.subplots(
    4, 1, figsize=(11, 9), dpi=150, sharex=True, gridspec_kw=dict(hspace=0.35)
)
series = [
    ("Vs1", models[:, 0, 3], "Vs₁ (m/s)", 1500),
    ("Vs2", models[:, 1, 3], "Vs₂ (m/s)", 3500),
    ("h1", models[:, 0, 0], "h₁ (m)", 100),
    ("misfit", misfit, "Weighted residual norm", None),
]
it = np.arange(N)
for ax, (key, v, lab, tr) in zip(axs, series):
    ax.axvspan(0, ADAPT, color="#e6e5e0", lw=0)
    ax.axvspan(ADAPT, BURN, color="#f2f1ed", lw=0)
    for c in range(C):
        ax.plot(
            it[::10],
            v[c, ::10],
            color=CHAIN[c],
            lw=0.6,
            alpha=0.85,
            label=f"Chain {c + 1}",
        )
    if tr is not None:
        ax.axhline(tr, color="#0b0b0b", lw=1.2, ls=":", label="Truth")
    p = P[key]
    ax.set_ylabel(lab)
    ax.set_title(
        f"{lab}: R̂ = {p['rhat']:.3f}, bulk ESS = {p['ess_bulk']:.0f}",
        loc="left",
        fontsize=10,
    )
    ax.grid(color="#e6e5e0", lw=0.5, axis="y")
    ax.set_axisbelow(True)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)
axs[0].legend(
    ncol=5, fontsize=8, frameon=False, loc="lower right", bbox_to_anchor=(1, 1.12)
)
axs[-1].text(
    ADAPT / 2,
    0.93,
    "step-size\nadaptation",
    transform=axs[-1].get_xaxis_transform(),
    ha="center",
    va="top",
    fontsize=8,
    color="#52514e",
)
axs[-1].text(
    (ADAPT + BURN) / 2,
    0.93,
    "burn-in",
    transform=axs[-1].get_xaxis_transform(),
    ha="center",
    va="top",
    fontsize=8,
    color="#52514e",
)
axs[-1].set_ylim(top=misfit.max() * 1.6)
axs[-1].set_xlabel("Iteration")
fig.suptitle(
    "Chain convergence — dispersed starts, frozen tuned steps after adaptation",
    fontsize=12,
    x=0.02,
    ha="left",
    y=0.995,
)
fig.savefig(f"{D}/convergence.png", bbox_inches="tight")
print("saved")
