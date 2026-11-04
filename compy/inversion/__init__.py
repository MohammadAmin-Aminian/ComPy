"""Shear-velocity inversion of normalized seafloor compliance.

Submodules:

``dram``        current Bayesian inversion (DRAM sampler, CRUST1.0 models,
                full data covariance, posterior summaries, serpentinization)
``metropolis``  random-walk Metropolis sampler of Aminian et al. (2025),
                with parallel tempering and an L-curve scan
``occam``       smoothest-model (Occam) inversion
``propagator``  forward model: normalized compliance of a layered seafloor
``models``      starting and reference models, empirical Vp/density relations
"""

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
from compy.inversion.metropolis import (
    invert_compliance,
    invert_compliance_broad_bounds,
    invert_tempered,
    l_curve,
)
from compy.inversion.models import crust1_starting_model, sediment_crust_model
from compy.inversion.occam import invert_occam
from compy.inversion.propagator import compliance_kernels, normalized_compliance

__all__ = [
    "ComplianceData",
    "CrustalModel",
    "Posterior",
    "compliance_kernels",
    "crust1_column",
    "crust1_model",
    "crust1_start",
    "crust1_starting_model",
    "forward",
    "invert_compliance",
    "invert_compliance_broad_bounds",
    "invert_occam",
    "invert_tempered",
    "l_curve",
    "normalized_compliance",
    "sample",
    "sediment_crust_model",
    "sensitivity_depth",
    "serpentinization",
    "vs_kernels",
]
