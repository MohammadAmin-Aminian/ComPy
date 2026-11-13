"""Standard figure of one automatic inversion (``compy.inversion.auto``).

``station_figure`` draws four panels: the data with the fit, the residuals in
units of the data error, shear velocity against depth with the prior and the
reference, and how much the data narrowed each parameter. Colours are one
fixed assignment (posterior blue, prior orange, reference green) from a
colour-blind-checked palette; no panel has two y-scales. The figure uses its
own local style and leaves Matplotlib's global settings alone.
"""

from __future__ import annotations

import numpy as np

POSTERIOR, PRIOR, REFERENCE = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#dcdbd6", "#fcfcfb"
STYLE = {
    "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.linewidth": 0.6, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.5, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "legend.frameon": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.spines.top": False, "axes.spines.right": False,
}  # fmt: skip


def _short(name):
    return (name.replace("sediment_", "sed. ").replace("vs_crust_", "Vs layer ")
            .replace("vs_mantle", "Vs mantle").replace("moho_depth", "Moho depth")
            .replace("vp_upper", "Vp factor, upper").replace("vp_lower", "Vp factor, lower")
            .replace("density", "density factor").replace("vs_scale", "Vs scale")
            .replace("modulus", "P modulus").replace("thickness", "thickness"))  # fmt: skip


def station_figure(result, path=None, title=None, max_depth=None):
    """Four-panel summary of an ``InversionResult``.

    Args:
        result (InversionResult): From ``compy.inversion.auto.invert``.
        path (str or Path): File to write (PNG or PDF); None returns the figure only.
        title (str): Figure title; default the station label and its quality.
        max_depth (float): Deepest depth drawn [m]; default the deepest Moho allowed.

    Returns:
        matplotlib.figure.Figure
    """
    import matplotlib

    with matplotlib.rc_context(STYLE):
        import matplotlib.pyplot as plt

        data, model, post = result.data, result.model, result.posterior
        quality = result.quality()
        fig = plt.figure(figsize=(11.0, 7.4))
        grid = fig.add_gridspec(2, 3, height_ratios=[2.2, 1.0], width_ratios=[1.25, 1.0, 0.9],
                                hspace=0.32, wspace=0.34, left=0.07, right=0.985, top=0.9,
                                bottom=0.08)  # fmt: skip
        f = 1e3 * data.freqs
        scale = 1e11
        pred = post.pred.reshape(-1, len(f))
        lo, med, hi = np.quantile(pred, [0.05, 0.5, 0.95], axis=0)

        ax = fig.add_subplot(grid[0, 0])
        ax.fill_between(f, scale * lo, scale * hi, color=POSTERIOR, alpha=0.25, lw=0,
                        label="model, 5-95 %")  # fmt: skip
        ax.plot(f, scale * med, color=POSTERIOR, lw=1.6, label="model, median")
        ax.errorbar(f, scale * data.eta, yerr=scale * data.sigma, fmt="o", ms=4, color=INK,
                    elinewidth=0.9, capsize=2, label="measured, 1 sigma")  # fmt: skip
        ax.set(ylabel="Normalized compliance (10$^{-11}$ Pa$^{-1}$)",
               title="Data and fit")  # fmt: skip
        ax.tick_params(labelbottom=False)
        ax.legend(loc="upper left")

        ax2 = fig.add_subplot(grid[1, 0], sharex=ax)
        ax2.axhspan(-1, 1, color=GRID, alpha=0.6, lw=0)
        ax2.axhline(0, color=MUTED, lw=0.6)
        ax2.plot(f, (data.eta - med) / data.sigma, "o-", ms=4, lw=1.0, color=INK)
        ax2.set(xlabel="Frequency (mHz)", ylabel="Residual (sigma)",
                title=f"Misfit: chi-square {quality['chi2']:.1f} for {quality['n_data']} data")  # fmt: skip
        bound = max(3.0, 1.1 * np.max(np.abs((data.eta - med) / data.sigma)))
        ax2.set_ylim(-bound, bound)

        deepest = max_depth or model.moho_bounds[1] + 500.0
        z = np.unique(np.concatenate([np.geomspace(5.0, 500.0, 60),
                                      np.arange(500.0, deepest + 1.0, 25.0)]))  # fmt: skip
        prof = result.profile(z)
        q = list(np.round(prof["q"], 2))
        vs, pri = prof["vs"] / 1e3, prof["vs_prior"] / 1e3
        axv = fig.add_subplot(grid[:, 1])
        axv.fill_betweenx(z, vs[q.index(0.05)], vs[q.index(0.95)], color=POSTERIOR, alpha=0.18,
                          lw=0, label="posterior, 5-95 %")  # fmt: skip
        axv.fill_betweenx(z, vs[q.index(0.16)], vs[q.index(0.84)], color=POSTERIOR, alpha=0.38,
                          lw=0, label="posterior, 16-84 %")  # fmt: skip
        axv.plot(vs[q.index(0.5)], z, color=POSTERIOR, lw=1.8, label="posterior median")
        for key in (0.16, 0.84):
            axv.plot(pri[q.index(key)], z, color=PRIOR, lw=1.0, ls="--",
                     label="prior, 16-84 %" if key == 0.16 else None)  # fmt: skip
        ref = model.reference
        tops = np.concatenate([[0.0], np.cumsum(ref[:-1, 0])])
        zr = np.repeat(np.append(tops, deepest * 2), 2)[1:-1]
        axv.plot(np.repeat(ref[:, 3], 2) / 1e3, zr + prof["sediment"], color=REFERENCE, lw=1.2,
                 label="reference (below sediment)")  # fmt: skip
        res = prof["resolved"]
        edges = np.flatnonzero(np.diff(np.concatenate([[0], res.astype(int), [0]])))
        for a, b in zip(edges[::2], edges[1::2]):
            axv.plot([0.02, 0.02], [z[a], z[b - 1]], transform=axv.get_yaxis_transform(),
                     color=INK, lw=4, solid_capstyle="butt",
                     label="resolved by the data" if a == edges[0] else None)  # fmt: skip
        axv.axhline(prof["reach"], color=MUTED, lw=0.8, ls=":")
        axv.text(0.98, prof["reach"], "90 % of sensitivity above", ha="right", va="bottom",
                 transform=axv.get_yaxis_transform(), fontsize=7, color=MUTED)  # fmt: skip
        axv.set_yscale("log")
        axv.set_ylim(deepest, 5.0)
        axv.set_xlim(0, 5.0)
        axv.set(xlabel="Shear velocity (km/s)", ylabel="Depth below seafloor (m)",
                title="Shear velocity")  # fmt: skip
        axv.legend(loc="upper right")

        axr = fig.add_subplot(grid[:, 2])
        res_p = result.parameter_resolution()
        names = [n for n in post.names if np.isfinite(res_p.get(n, np.nan))]
        y = np.arange(len(names))[::-1]
        values = np.array([res_p[n] for n in names])
        axr.axvline(0.6, color=MUTED, lw=0.8, ls=":")
        axr.axvline(1.0, color=MUTED, lw=0.6)
        axr.hlines(y, 0, np.minimum(values, 1.6), color=GRID, lw=1.0)
        axr.plot(np.minimum(values, 1.6), y, "o", ms=6, color=POSTERIOR)
        axr.set_yticks(y, [_short(n) for n in names])
        axr.set_xlim(0, 1.65)
        axr.set(xlabel="Posterior width / prior width",
                title="What the data constrain")  # fmt: skip
        axr.text(0.6, len(names) - 0.4, "resolved", ha="right", va="bottom", fontsize=7,
                 color=MUTED)  # fmt: skip
        axr.text(1.0, len(names) - 0.4, " prior only", ha="left", va="bottom", fontsize=7,
                 color=MUTED)  # fmt: skip
        axr.grid(axis="y", visible=False)

        label = title or (
            f"{data.station or 'station'}: water depth {data.depth:.0f} m, {quality['label']}"
            + (f" ({'; '.join(quality['reasons'])})" if quality["reasons"] else "")
        )
        fig.suptitle(label, x=0.07, ha="left", fontsize=11, color=INK)
        if path is not None:
            fig.savefig(path, dpi=160)
            plt.close(fig)
        return fig
