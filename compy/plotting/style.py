"""Matplotlib style for ComPy figures (STIX fonts, large type for slides and papers)."""

import matplotlib.pyplot as plt

_BASE = {
    "mathtext.fontset": "stix",
    "font.family": ["STIXGeneral"],
    "font.weight": "normal",
    "mathtext.default": "regular",
    "axes.grid": True,
    "axes.grid.which": "major",
    "axes.labelcolor": "black",
    "axes.labelpad": 4.0,
    "axes.labelweight": "normal",
    "axes.linewidth": 0.6,
    "axes.titlecolor": "black",
    "axes.titlepad": 5.0,
    "axes.titlesize": "large",
    "axes.titleweight": "normal",
    "axes.xmargin": 0.05,
    "axes.ymargin": 0.05,
    "grid.alpha": 1,
    "grid.color": "#b0b0b0",
    "grid.linestyle": "--",
    "grid.linewidth": 0.6,
    "xtick.color": "black",
    "xtick.direction": "out",
    "xtick.major.pad": 2.0,
    "xtick.minor.pad": 2.0,
    "xtick.minor.visible": True,
    "ytick.color": "black",
    "ytick.direction": "out",
    "ytick.major.pad": 2.0,
    "ytick.minor.pad": 2.0,
    "ytick.minor.visible": True,
    "xtick.minor.size": 10,
    "ytick.minor.size": 10,
    "ytick.minor.width": 1,
    "xtick.major.size": 10,
    "ytick.major.size": 10,
    "xtick.major.width": 2,
    "ytick.major.width": 2,
    "xtick.minor.width": 1,
    "legend.borderaxespad": 0.0,
    "legend.borderpad": 0.5,
    "legend.columnspacing": 1.5,
    "legend.edgecolor": "gray",
    "legend.facecolor": "white",
    "legend.fancybox": False,
    "legend.fontsize": "small",
    "legend.framealpha": 0.8,
    "legend.handleheight": 1.0,
    "legend.handlelength": 2.0,
    "legend.handletextpad": 0.5,
}


def set_plot_style(font_size=40, label_size="medium", tick_label_size=35):
    """Apply the ComPy figure style to Matplotlib's global ``rcParams``.

    The defaults are the paper figures' sizes; ``set_plot_style_large`` gives
    the larger type used by the spectral plots.
    """
    plt.rcParams.update(_BASE)
    plt.rcParams["font.size"] = font_size
    plt.rcParams["axes.labelsize"] = label_size
    plt.rcParams["xtick.labelsize"] = tick_label_size
    plt.rcParams["ytick.labelsize"] = tick_label_size


def set_plot_style_large():
    """ComPy style with 50 pt type everywhere (spectral and coherogram figures)."""
    set_plot_style(font_size=50, label_size=50, tick_label_size=50)
