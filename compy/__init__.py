"""ComPy: seafloor compliance, pressure-gauge calibration and shear-velocity inversion.

Processing chain, one module per step (all SI units):

``compy.streams``        ObsPy window helpers (split, trim, validate)
``compy.cleaning``       transient and periodic-glitch removal
``compy.preprocessing``  tilt correction and coherent-noise removal
``compy.compliance``     spectra, window selection, stacked compliance
``compy.inversion``      Vs inversion (DRAM, Metropolis, Occam), forward model
``compy.diagnostics``    convergence diagnostics and reproducible run records
``compy.plotting``       figures
``compy.numerics``       constants and shared numerical primitives

The most used functions are also available here, e.g. ``compy.correct_tilt``.
Names from ComPy 2.x (``compy.Calculate_Compliance``, ``inv_compy``,
``ffplot``...) still work and emit a ``DeprecationWarning`` with the new name.
"""

from compy._deprecation import legacy_getattr
from compy._legacy import LEGACY

__version__ = "2.5.0"

_SUBMODULES = {
    "cleaning",
    "compliance",
    "diagnostics",
    "inversion",
    "numerics",
    "plotting",
    "preprocessing",
    "streams",
}

# Public shortcuts: name -> "module:name". Imported lazily on first access.
_PUBLIC = {
    # streams
    "split_stream": "compy.streams:split_stream",
    "split_and_save_stream": "compy.streams:split_and_save_stream",
    "cut_stream_with_overlap": "compy.streams:cut_stream_with_overlap",
    "trim_streams_to_same_length": "compy.streams:trim_streams_to_same_length",
    "validate_stream": "compy.streams:validate_stream",
    # cleaning and preprocessing
    "remove_periodic_transient": "compy.cleaning:remove_periodic_transient",
    "rotate_smoothly": "compy.cleaning:rotate_smoothly",
    "correct_tilt": "compy.preprocessing:correct_tilt",
    "estimate_tilt_angles": "compy.preprocessing:estimate_tilt_angles",
    "remove_coherent_noise": "compy.preprocessing:remove_coherent_noise",
    # compliance
    "estimate_compliance": "compy.compliance:estimate_compliance",
    "stacked_compliance": "compy.compliance:stacked_compliance",
    "gravity_corrected_compliance": "compy.compliance:gravity_corrected_compliance",
    "gravitational_attraction": "compy.compliance:gravitational_attraction",
    "calculate_compliance": "compy.compliance:calculate_compliance",
    "calculate_compliance_gated": "compy.compliance:calculate_compliance_gated",
    "ComplianceResult": "compy.compliance:ComplianceResult",
    "SelectionCriteria": "compy.compliance:SelectionCriteria",
    # numerics
    "wavenumber": "compy.numerics:wavenumber",
    "compliance_uncertainty": "compy.numerics:compliance_uncertainty",
    "sliding_window": "compy.numerics:sliding_window",
    "rms": "compy.numerics:rms",
    "GRAVITY": "compy.numerics:GRAVITY",
    "GRAVITATIONAL_CONSTANT": "compy.numerics:GRAVITATIONAL_CONSTANT",
    # figures
    "set_plot_style": "compy.plotting.style:set_plot_style",
}

__all__ = ["__version__", *sorted(_SUBMODULES), *_PUBLIC]

__getattr__ = legacy_getattr(
    "compy",
    {**LEGACY["compy"], **_PUBLIC},
    quiet=set(_PUBLIC),
    submodules=_SUBMODULES,
)


def __dir__():
    return __all__
