"""Figures for ComPy: ``style``, ``spectra`` (PSD/coherence), ``inversion`` (posteriors), ``qc``.

Submodules are imported on demand so that computation never pulls in
plotting code.
"""

from compy.plotting.style import set_plot_style, set_plot_style_large

__all__ = ["set_plot_style", "set_plot_style_large"]
