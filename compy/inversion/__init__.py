"""Shear-velocity inversion of normalized seafloor compliance.

Submodules:

``auto``        the default inversion since ComPy 3: one call from a compliance
                curve to Vs with intervals, resolved depths and a quality label
``errors``      covariance of a stacked compliance estimate from its windows
``cli``         the ``compy-invert`` command
``dram``        sampler, forward model, data container and the ComPy 2.3-2.8
                model (``CrustalModel``), kept for earlier runs and as a check
``metropolis``  random-walk Metropolis sampler of Aminian et al. (2025),
                with parallel tempering and an L-curve scan
``occam``       smoothest-model (Occam) inversion on ln(Vs)
``minimum_structure``  minimum-structure inversion after Crawford (2004)
``literature_mcmc``    Monte Carlo inversion with the parameterizations of
                Zha & Webb (2016) and Zheng et al. (2024)
``impedance``   independent forward model: surface-impedance recursion
``propagator``  forward model: normalized compliance of a layered seafloor
``models``      starting and reference models, empirical Vp/density relations
"""

from compy.inversion.auto import (
    InversionResult,
    SeafloorModel,
    compliance_data,
    invert,
    load_result,
    model_from_crust1,
    reference_from_crust1,
)
from compy.inversion.dram import (
    ComplianceData,
    CrustalModel,
    Posterior,
    crust1_column,
    crust1_model,
    crust1_start,
    forward,
    sample,
    sensitivity_depth,
    serpentinization,
    vs_kernels,
)
from compy.inversion.impedance import compliance_by_impedance
from compy.inversion.literature_mcmc import (
    Zheng2024Model,
    ZhaWebb2016Model,
    sample_literature_model,
)
from compy.inversion.metropolis import (
    invert_compliance,
    invert_compliance_broad_bounds,
    invert_tempered,
    l_curve,
)
from compy.inversion.models import crust1_starting_model, sediment_crust_model
from compy.inversion.minimum_structure import (
    crawford_thin_layers,
    invert_minimum_structure,
)
from compy.inversion.occam import invert_occam
from compy.inversion.propagator import compliance_kernels, normalized_compliance

__all__ = [
    "InversionResult",
    "SeafloorModel",
    "compliance_data",
    "invert",
    "load_result",
    "model_from_crust1",
    "reference_from_crust1",
    "ComplianceData",
    "CrustalModel",
    "Posterior",
    "ZhaWebb2016Model",
    "Zheng2024Model",
    "compliance_by_impedance",
    "compliance_kernels",
    "crust1_column",
    "crust1_model",
    "crust1_start",
    "crawford_thin_layers",
    "crust1_starting_model",
    "forward",
    "invert_compliance",
    "invert_compliance_broad_bounds",
    "invert_minimum_structure",
    "invert_occam",
    "invert_tempered",
    "l_curve",
    "normalized_compliance",
    "sample",
    "sample_literature_model",
    "sediment_crust_model",
    "sensitivity_depth",
    "serpentinization",
    "vs_kernels",
]
