"""Seismology-style posterior figures for _Example/deep_crust_inversion.py.

Usage: python _Example/plot_deep_posterior.py <output directory of the run>

Writes deep_posterior.png (Vs-depth posterior density in the convention of the
RR52 figure, sensitivity kernels, data fit, prior-versus-posterior widths) and
deep_convergence.png (chain traces, sampled smoothness, misfit).
"""
import sys
import json
from types import SimpleNamespace

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

D = sys.argv[1]
z = np.load(f"{D}/chains.npz")
S = json.load(open(f"{D}/summary.json"))
d = SimpleNamespace(TRUTH=z["truth"], EDGES=z["edges"], F=z["freq"], H=float(z["water_depth"]))
d.N_LAYER = len(d.EDGES) - 1
d.MID = 0.5 * (d.EDGES[:-1] + d.EDGES[1:])
d.sample_truth = lambda column, depth: d.TRUTH[
    np.searchsorted(np.cumsum(d.TRUTH[:-1, 0]), depth, side="right"), column]
kern = {"z": z["kernel_z"], "K": z["kernels"], "f": d.F}
prior = z["prior_vs"]
BURN = S["settings"]["burnin"]
vs = z["vs"][:, :d.N_LAYER, BURN::5].astype(float) / 1000          # km/s, thinned
post = vs.transpose(1, 0, 2).reshape(d.N_LAYER, -1)                # layers x draws
edges = d.EDGES / 1000
ZMAX = 14.0
HS_VS = d.TRUTH[-1, 3] / 1000
lam = 2 * np.pi / __import__("inv_compy").gravd(2 * np.pi * d.F, d.H)
depth_limit = lam.max() / 6 / 1000                                  # Crawford et al. (1991)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
    "axes.edgecolor": "#52514e", "axes.linewidth": 0.8, "xtick.color": "#52514e",
    "ytick.color": "#52514e", "axes.labelcolor": "#0b0b0b",
})
MED, CI, TRUE, OCC, LIM = "#e34948", "#ff9f40", "#ffffff", "#7fd3ff", "#c3c2b7"


def steps(values, bottom_val=HS_VS):
    """Step profile (x, z) through the inversion layers and the half-space."""
    zz, xx = [0.0], [values[0]]
    for k in range(d.N_LAYER):
        zz += [edges[k + 1], edges[k + 1]]
        xx += [values[k], values[k + 1] if k + 1 < d.N_LAYER else bottom_val]
    zz.append(ZMAX); xx.append(bottom_val)
    return np.array(xx), np.array(zz)


def truth_steps():
    t = d.TRUTH
    zz, xx, top = [0.0], [t[0, 3]], 0.0
    for k in range(len(t) - 1):
        top += t[k, 0]
        zz += [top, top]; xx += [t[k, 3], t[k + 1, 3]]
    zz.append(ZMAX * 1000); xx.append(t[-1, 3])
    return np.array(xx) / 1000, np.array(zz) / 1000


fig = plt.figure(figsize=(15.5, 9.6), dpi=170)
gs = GridSpec(3, 4, figure=fig, width_ratios=[1.35, 0.045, 0.16, 1.0],
              height_ratios=[1, 1, 1], wspace=0.25, hspace=0.62)

# ---- a) posterior Vs-depth density -----------------------------------------
ax = fig.add_subplot(gs[:, 0])
vbins = np.arange(0, 5.0001, 0.05)
H = np.zeros((d.N_LAYER + 1, len(vbins) - 1))
for k in range(d.N_LAYER):
    h, _ = np.histogram(post[k], bins=vbins)
    H[k] = h / h.max()
H[-1, np.searchsorted(vbins, HS_VS) - 1] = 1.0                    # fixed half-space
yedges = np.append(edges, ZMAX)
im = ax.pcolormesh(vbins, yedges, H, cmap="viridis", vmin=0, vmax=1,
                   shading="flat", rasterized=True)
q05, q50, q95 = np.quantile(post, [0.05, 0.5, 0.95], axis=1)
for qq in (q05, q95):
    x, y = steps(qq); ax.plot(x, y, color=CI, lw=1.1, ls=(0, (4, 2)))
x, y = steps(q50); ax.plot(x, y, color=MED, lw=2.6, ls=(0, (5, 2.5)))
x, y = steps(z["occam"][:d.N_LAYER, 3] / 1000); ax.plot(x, y, color=OCC, lw=1.4)
x, y = truth_steps(); ax.plot(x, y, color=TRUE, lw=1.8)
ax.axhline(depth_limit, color=LIM, lw=1, ls=":")
ax.text(0.08, depth_limit - 0.12, f"λmax/6 = {depth_limit:.1f} km: deepest sensed structure "
        "(Crawford et al. 1991)", color=LIM, fontsize=8, va="bottom")
ax.axhline(edges[-1], color=LIM, lw=0.8, ls="-", alpha=0.6)
ax.text(0.08, edges[-1] + 0.1, "fixed mantle half-space", color=LIM, fontsize=8, va="top")
ax.set_xlim(0, 5); ax.set_ylim(ZMAX, 0)
ax.set_xlabel("Shear-wave velocity, Vs (km/s)")
ax.set_ylabel("Depth below seafloor (km)")
ax.set_title("a) Posterior Vs–depth density", loc="left")
if "sigma_jump" in S:
    sj = S["sigma_jump"]
    ax.text(0.03, 0.235, f"Smoothness learned from the data (hierarchical Bayes):\n"
            f"jump scale σ = {sj['median']:.0f} m/s (90 %: {sj['p05']:.0f}–{sj['p95']:.0f})",
            transform=ax.transAxes, color="white", fontsize=8, va="top")
ax.grid(color="white", alpha=0.15, lw=0.6)
ax.legend(handles=[
    Line2D([], [], color=TRUE, lw=1.8, label="True model"),
    Line2D([], [], color=MED, lw=2.6, ls=(0, (5, 2.5)), label="Posterior median"),
    Line2D([], [], color=CI, lw=1.1, ls=(0, (4, 2)), label="5–95 % credible bounds"),
    Line2D([], [], color=OCC, lw=1.4, label="Occam minimum-structure model"),
], loc="center left", bbox_to_anchor=(0.0, 0.52), fontsize=8.5, framealpha=0.85, facecolor="#1a1a19",
    edgecolor="none", labelcolor="white")
cb = fig.colorbar(im, cax=fig.add_subplot(gs[:, 1]))
cb.set_label("Normalized probability density (per depth)")

# ---- b) sensitivity kernels -------------------------------------------------
axk = fig.add_subplot(gs[0, 3])
zk, K, fk = kern["z"] / 1000, kern["K"], kern["f"]
sel = [0, 3, 6, 9, 12, 15]
cmap = plt.get_cmap("viridis")
for n, i in enumerate(sel):
    axk.plot(-K[i] / np.abs(K[i]).max(), zk, color=cmap(n / (len(sel) - 1) * 0.9), lw=1.6,
             label=f"{fk[i] * 1e3:.1f} mHz")
axk.axvline(0, color="#52514e", lw=0.6)
axk.axhline(depth_limit, color="#52514e", lw=0.8, ls=":")
axk.set_ylim(ZMAX, 0)
axk.set_xlabel("−∂ln η/∂ln Vs, normalized per frequency")
axk.set_ylabel("Depth (km)")
axk.set_title("b) Vs sensitivity kernels (true model)", loc="left")
axk.legend(fontsize=7.5, frameon=False, ncol=2, loc="lower right")
axk.grid(color="#e6e5e0", lw=0.5); axk.set_axisbelow(True)

# ---- c) data fit -------------------------------------------------------------
axd = fig.add_subplot(gs[1, 3])
pr = z["preds"][:, BURN::5].reshape(-1, len(d.F)).astype(float)
lo, md, hi = np.quantile(pr, [0.05, 0.5, 0.95], axis=0)
fm, sc = d.F * 1e3, 1e11
import inv_compy as inv
occ_pred = inv.calc_norm_compliance(d.H, d.F, z["occam"])
axd.fill_between(fm, lo * sc, hi * sc, color="#2a78d6", alpha=0.25, lw=0, label="Posterior predictive 5–95 %")
axd.plot(fm, md * sc, color="#2a78d6", lw=2, label="Posterior predictive median")
axd.plot(fm, occ_pred * sc, color="#1baf7a", lw=1.3, ls="--", label="Occam model")
axd.plot(fm, z["clean"] * sc, color="#0b0b0b", lw=1.1, ls=":", label="True compliance")
axd.errorbar(fm, z["obs"] * sc, yerr=z["unc"] * sc, fmt="o", ms=4, color="#0b0b0b",
             mfc="white", mew=1.1, elinewidth=1, capsize=2, label="Observed ±1σ (3 %)")
axd.set_xlabel("Frequency (mHz)")
axd.set_ylabel("Normalized compliance\n(×10⁻¹¹ Pa⁻¹)")
axd.set_title(f"c) Data fit (water depth {d.H:.0f} m)", loc="left")
axd.legend(fontsize=7.2, frameon=False, loc="upper left")
axd.grid(color="#e6e5e0", lw=0.5); axd.set_axisbelow(True)

# ---- d) resolution: posterior vs prior width --------------------------------
axr = fig.add_subplot(gs[2, 3])
mid = d.MID / 1000
pw = (q95 - q05)
axr.plot(pw, mid, "o-", color="#2a78d6", ms=4, lw=1.6, label="Posterior 90 % width")
if prior is not None:
    pp = prior[:, :d.N_LAYER].transpose(1, 0, 2).reshape(d.N_LAYER, -1).astype(float) / 1000
    prw = np.quantile(pp, 0.95, axis=1) - np.quantile(pp, 0.05, axis=1)
    axr.plot(prw, mid, "s--", color="#52514e", ms=3.5, lw=1.2, label="Prior 90 % width")
axr.axhline(depth_limit, color="#52514e", lw=0.8, ls=":")
axr.set_ylim(ZMAX, 0); axr.set_xlim(left=0)
axr.set_xlabel("90 % interval width (km/s)")
axr.set_ylabel("Depth (km)")
rh = max(v["rhat"] for v in S["layers"].values())
es = min(v["ess_bulk"] for v in S["layers"].values())
axr.set_title(f"d) What the data resolve  (max R̂ {rh:.3f}, min ESS {es:.0f})", loc="left", fontsize=10)
axr.legend(fontsize=7.5, frameon=False, loc="center", bbox_to_anchor=(0.62, 0.5))
axr.grid(color="#e6e5e0", lw=0.5); axr.set_axisbelow(True)

fig.suptitle("ComPy compliance inversion — deep-crust synthetic test "
             f"(4 chains × {S['settings']['iterations']:,} iterations, 14 layers to 13 km)",
             fontsize=13, x=0.02, ha="left")
fig.savefig(f"{D}/deep_posterior.png", bbox_inches="tight")
plt.close(fig)

# ---- convergence -------------------------------------------------------------
CH = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
show = [0, 4, 7, 10]
hier = "alpha" in z.files
fig, axs = plt.subplots(len(show) + 1 + hier, 1, figsize=(11, 10 + 2 * hier), dpi=140,
                        sharex=True, gridspec_kw=dict(hspace=0.5))
full = z["vs"]
it = np.arange(full.shape[-1])
names = S["settings"]["starts"]
for ax, j in zip(axs, show):
    ax.axvspan(0, S["settings"]["adapt_until"], color="#e6e5e0", lw=0)
    ax.axvspan(S["settings"]["adapt_until"], BURN, color="#f2f1ed", lw=0)
    for c in range(full.shape[0]):
        ax.plot(it[::20], full[c, j, ::20] / 1000, color=CH[c], lw=0.5, alpha=0.85,
                label=f"Chain {c + 1} ({names[c]} start)")
    tr = d.sample_truth(3, d.MID[j]) / 1000
    ax.axhline(tr, color="#0b0b0b", lw=1.1, ls=":", label="Truth")
    L = S["layers"][f"layer_{j + 1:02d}"]
    ax.set_title(f"Layer {j + 1} ({L['top_m'] / 1000:.2f}–{L['bottom_m'] / 1000:.2f} km): "
                 f"R̂ = {L['rhat']:.3f}, bulk ESS = {L['ess_bulk']:.0f}", loc="left", fontsize=10)
    ax.set_ylabel("Vs (km/s)")
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
axs[0].legend(ncol=3, fontsize=7.5, frameon=False, loc="lower right", bbox_to_anchor=(1, 1.18))
if hier:
    ax = axs[len(show)]
    ax.axvspan(0, S["settings"]["adapt_until"], color="#e6e5e0", lw=0)
    ax.axvspan(S["settings"]["adapt_until"], BURN, color="#f2f1ed", lw=0)
    for c in range(full.shape[0]):
        ax.plot(it[::20], np.sqrt(1 / (0.25 * z["alpha"][c, ::20].astype(float))),
                color=CH[c], lw=0.5)
    ax.set_yscale("log"); ax.set_ylabel("σ (m/s)")
    sj = S["sigma_jump"]
    ax.set_title(f"Smoothness jump scale σ (sampled): R̂ = {sj['rhat']:.3f}, "
                 f"bulk ESS = {sj['ess_bulk']:.0f}", loc="left", fontsize=10)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
ax = axs[-1]
ax.axvspan(0, S["settings"]["adapt_until"], color="#e6e5e0", lw=0)
ax.axvspan(S["settings"]["adapt_until"], BURN, color="#f2f1ed", lw=0)
for c in range(full.shape[0]):
    ax.plot(it[::20], z["misfit"][c, ::20] ** 2, color=CH[c], lw=0.5)
ax.axhline(len(d.F), color="#0b0b0b", lw=1, ls="--")
ax.set_yscale("log")
ax.set_ylabel("χ²")
ax.set_title(f"Data misfit χ² (dashed: N = {len(d.F)})", loc="left", fontsize=10)
ax.set_xlabel("Iteration")
ax.text(S["settings"]["adapt_until"] / 2, 0.92, "adaptation", transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=8, color="#52514e")
ax.text((S["settings"]["adapt_until"] + BURN) / 2, 0.92, "burn-in",
        transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8, color="#52514e")
for sp in ("right", "top"):
    ax.spines[sp].set_visible(False)
fig.suptitle("Chain convergence from dispersed starts", fontsize=12, x=0.02, ha="left", y=0.995)
fig.savefig(f"{D}/deep_convergence.png", bbox_inches="tight")
print("saved")
