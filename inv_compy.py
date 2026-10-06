#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:04:23 2023

@author: mohammadamin

"""

from compy_numerics import positive, integer, residual, log_likelihood
from compy_numerics import roughness as Roughness, wavenumber


def gravd(W, h):
    """Gravity-wave wavenumber in rad/m; scalar input retains legacy 1D output."""
    return np.atleast_1d(wavenumber(W, h))


import random
import numpy as np
import matplotlib.pyplot as plt
# %%


# %%
def _depth_profile(model, depth_grid):
    """Sample a layered model on a fixed depth grid; last row is half-space."""
    interfaces = np.cumsum(model[:-1, 0])
    return model[np.searchsorted(interfaces, depth_grid, side="right"), 3]


def _inversion(
    Data,
    f,
    depth_s,
    starting_model,
    s,
    n_layer,
    sigma_v,
    sigma_h,
    iteration,
    alpha,
    sta,
    seed,
    beta,
    return_profiles,
):
    """Symmetric random-walk Metropolis sampler with bounded model support."""
    iteration = integer(iteration, "iteration", minimum=2)
    sigma_v, sigma_h = positive(sigma_v, "sigma_v"), positive(sigma_h, "sigma_h")
    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")
    Data, f = np.asarray(Data, float), np.asarray(f, float)
    if Data.ndim != 1 or f.shape != Data.shape or np.any(f <= 0):
        raise ValueError("data and positive frequencies must be equal-length vectors")
    if s is None:
        raise ValueError("s is required: supply positive measurement uncertainties")
    residual(Data, Data, s)
    if starting_model is None:
        initial = Model_V2(1, n_layer=n_layer, sta=sta)[0][:, :, 0]
    else:
        initial = np.asarray(starting_model, float)
        if initial.ndim == 3:
            initial = initial[:, :, 0]
    initial = _validate_model(initial).copy()
    rng = np.random.default_rng(seed)
    chain = np.empty((*initial.shape, iteration))
    prediction = np.empty((iteration, len(f)))
    misfits = np.empty((1, iteration))
    likelihood = np.empty((1, iteration))
    depth_grid = np.arange(max(1, int(np.ceil(np.sum(initial[:-1, 0])))))
    prior = _depth_profile(initial, depth_grid)[None, :, None]
    profiles = np.empty((iteration, len(depth_grid), 1)) if return_profiles else None
    lower = initial[:, 3] * (0.15 if beta else 0.5)
    upper = initial[:, 3] * (1.2 if beta else 1.1)
    if beta:
        lower[-2], upper[-2] = initial[-2, 3] * np.array([0.9, 1.1])

    def evaluate(model):
        pred = calc_norm_compliance(depth_s, f, model)
        profile = _depth_profile(model, depth_grid)
        log_data = log_likelihood(Data, pred, s)
        log_posterior = log_data - 0.5 * alpha * Roughness(profile, 2)
        return pred, profile, log_data, log_posterior

    current = initial.copy()
    pred, profile, log_data, log_post = evaluate(current)
    accepted = 0
    for i in range(iteration):
        if i:
            candidate = current.copy()
            layer = rng.integers(len(initial) - 1)  # fixed elastic half-space
            if rng.integers(2) == 0:
                candidate[layer, 3] += rng.normal(0, sigma_v)
                # Hold Vp and rho fixed, as specified by Aminian et al. (2025).
            else:
                # Transfer thickness from the deepest finite layer. Total
                # finite depth is conserved; the half-space is never altered.
                if len(initial) == 2:
                    chain[:, :, i] = current
                    prediction[i] = pred
                    misfits[0, i] = np.linalg.norm(residual(Data, pred, s))
                    likelihood[0, i] = np.exp(log_data)
                    if profiles is not None:
                        profiles[i, :, 0] = profile
                    continue
                if layer == len(initial) - 2:
                    layer = rng.integers(len(initial) - 2)
                step = rng.normal(0, sigma_h)
                candidate[layer, 0] += step
                candidate[-2, 0] -= step
            valid = (
                np.all(candidate[:-1, 0] > 0)
                and np.all(candidate[:, 3] >= lower)
                and np.all(candidate[:, 3] <= upper)
                and np.all(candidate[:, 3] < candidate[:, 2])
            )
            if valid:
                new_pred, new_profile, new_log_data, new_log_post = evaluate(candidate)
                if (
                    np.log(max(rng.random(), np.finfo(float).tiny))
                    < new_log_post - log_post
                ):
                    current = candidate
                    pred, profile, log_data, log_post = (
                        new_pred,
                        new_profile,
                        new_log_data,
                        new_log_post,
                    )
                    accepted += 1
            # Out-of-support proposals are rejected, not resampled/reflected:
            # this keeps the Gaussian proposal symmetric at the boundaries.
        chain[:, :, i] = current
        prediction[i] = pred
        misfits[0, i] = np.linalg.norm(residual(Data, pred, s))
        likelihood[0, i] = np.exp(log_data)  # may underflow; never used for acceptance
        if profiles is not None:
            profiles[i, :, 0] = profile
    return (
        chain,
        profiles,
        prior,
        misfits,
        prediction,
        likelihood,
        accepted / (iteration - 1),
    )


def invert_compliace(
    Data,
    f,
    depth_s,
    starting_model=None,
    s=None,
    n_layer=3,
    sediment_thickness=80,
    n_sediment_layer=3,
    sigma_v=5,
    sigma_h=5,
    iteration=100000,
    alpha=0.25,
    sta="RR38",
    *,
    seed=0,
    return_profiles=True,
):
    """Legacy name for the v2 sampler; see README for units and priors.

    Supply a custom starting_model to control sediments. The station templates
    do not use sediment_thickness or n_sediment_layer. s must be supplied.
    Set return_profiles=False to avoid the large metre-grid profile array.
    """
    chain, profiles, _, misfits, prediction, likelihood, rate = _inversion(
        Data,
        f,
        depth_s,
        starting_model,
        s,
        n_layer,
        sigma_v,
        sigma_h,
        iteration,
        alpha,
        sta,
        seed,
        False,
        return_profiles,
    )
    return chain, profiles, misfits, prediction, likelihood, rate


def invert_compliace_beta(
    Data,
    f,
    depth_s,
    starting_model=None,
    s=None,
    n_layer=3,
    sediment_thickness=80,
    n_sediment_layer=3,
    sigma_v=5,
    sigma_h=5,
    iteration=100000,
    alpha=0.25,
    sta="RR38",
    *,
    seed=0,
    return_profiles=True,
):
    """V2 sampler with broader shallow Vs bounds; returns legacy seven-tuple."""
    return _inversion(
        Data,
        f,
        depth_s,
        starting_model,
        s,
        n_layer,
        sigma_v,
        sigma_h,
        iteration,
        alpha,
        sta,
        seed,
        True,
        return_profiles,
    )


# Correctly spelled aliases; historical names remain available.
invert_compliance = invert_compliace
invert_compliance_beta = invert_compliace_beta


def _validate_model(model):
    model = np.asarray(model, dtype=float)
    if model.ndim != 2 or model.shape[1] != 4 or len(model) < 2:
        raise ValueError(
            "model must have at least two rows and columns [thickness, rho, Vp, Vs]"
        )
    if (
        np.any(~np.isfinite(model))
        or np.any(model[:, 0] < 0)
        or np.any(model[:, 1:] <= 0)
    ):
        raise ValueError(
            "model must be finite, with nonnegative thickness and positive rho, Vp, Vs"
        )
    if np.any(model[:-1, 0] <= 0) or np.any(model[:, 2] <= model[:, 3]):
        raise ValueError(
            "finite layers must have positive thickness and Vp must exceed Vs"
        )
    return model


def Lcurve(
    Data,
    f,
    depth_s,
    starting_model=None,
    sigma_v=5,
    sigma_h=1,
    iteration=100000,
    *,
    s=None,
):
    """Explore regularization weights; returns misfit chains and alpha values.

    This diagnostic sweep is not an automatic optimal-alpha estimator.
    """
    if s is None:
        raise ValueError("s is required for the regularization sweep")
    alpha = np.logspace(-2, 1, 10)
    chains = []
    for weight in alpha:
        result = invert_compliance(
            Data,
            f,
            depth_s,
            starting_model=starting_model,
            s=s,
            sigma_v=sigma_v,
            sigma_h=sigma_h,
            iteration=iteration,
            alpha=weight,
            return_profiles=False,
        )
        chains.append(result[2])
    burnin = int(0.8 * iteration)
    plt.loglog(alpha, [np.median(chain[0, burnin:]) for chain in chains])
    plt.xlabel("Regularization weight")
    plt.ylabel("Median normalized residual norm")
    return chains, alpha


# %%
def model_exp(iteration, first_layer=200, n_layer=15, power_factor=1.15):
    """Generic exponential-thickness crustal model without sediments."""
    return Model(iteration, first_layer, n_layer, power_factor, n_sediment_layer=0)


# %%


def Model(
    iteration,
    first_layer=200,
    n_layer=10,
    power_factor=1.17,
    sediment_thickness=80,
    n_sediment_layer=1,
):
    """Generic sediment plus graded-crust starting model (SI units).

    Use Model_V2 for the station-specific CRUST templates. The last row is
    the elastic half-space, represented by zero thickness.
    """
    iteration = integer(iteration, "iteration")
    n_layer = integer(n_layer, "n_layer")
    n_sediment_layer = integer(n_sediment_layer, "n_sediment_layer", minimum=0)
    first_layer = positive(first_layer, "first_layer")
    power_factor = positive(power_factor, "power_factor")
    sediment = np.empty((0, 4))
    if n_sediment_layer:
        sediment_thickness = positive(sediment_thickness, "sediment_thickness")
        weights = np.geomspace(1, sediment_thickness, n_sediment_layer)
        vs = np.geomspace(340, 600, n_sediment_layer)
        sediment = np.column_stack(
            [sediment_thickness * weights / weights.sum(), density(vs), velp(vs), vs]
        )
    vs = np.linspace(2700, 4300, n_layer)
    crust = np.column_stack(
        [first_layer * power_factor ** np.arange(n_layer), density(vs), velp(vs), vs]
    )
    model = np.vstack([sediment, crust, [0, 3340, 8120, 4510]])
    _validate_model(model)
    chain = np.zeros((*model.shape, iteration))
    chain[:, :, 0] = model
    grid = np.arange(max(1, int(np.ceil(model[:-1, 0].sum()))))
    prior = _depth_profile(model, grid)[None, :, None]
    return chain, prior, np.zeros_like(prior)


# %%
def Model_V2(iteration, n_layer=3, sta="RR38"):
    iteration = integer(iteration, "iteration")
    if n_layer not in (3, 6, 9, 12):
        raise ValueError("n_layer must be 3, 6, 9, or 12")
    if sta not in (
        "RR28",
        "RR29",
        "RR34",
        "RR36",
        "RR38",
        "RR40",
        "RR50",
        "RR52",
    ) and n_layer not in (3, 6):
        raise ValueError("generic station template supports only 3 or 6 layers")
    dep = 0
    if sta == "RR28":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            # Thickness, Density, Vp, Vs
            starting_model[:, :, 0] = np.array(
                [
                    [260.0, 1820.0, 1750.0, 340.0],
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [260.0, 1820.0, 1750.0, 340.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [260.0, 1820.0, 1750.0, 340.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [260.0, 1820.0, 1750.0, 340.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

    elif sta == "RR29":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

    elif sta == "RR34":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7530.0, 4370.0],
                    [10000.0, 3230.0, 7530.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7530.0, 4370.0],
                    [10000.0, 3230.0, 7530.0, 4370.0],
                ]
            )

    elif sta == "RR36":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

    elif sta == "RR38":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

    elif sta == "RR40":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

    elif sta == "RR50":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )
    elif sta == "RR52":
        sediment_thickness = 0
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [1170.0, 2400.0, 4700.0, 2540.0],
                    [1170.0, 2760.0, 6300.0, 3590.0],
                    [4560.0, 2960.0, 6900.0, 3940.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [585.0, 2400.0, 4700.0, 2540.0],
                    [585.0, 2400.0, 4700.0, 2540.0],
                    [585.0, 2760.0, 6300.0, 3590.0],
                    [585.0, 2760.0, 6300.0, 3590.0],
                    [2280.0, 2960.0, 6900.0, 3940.0],
                    [2280.0, 2960.0, 6900.0, 3940.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [390.0, 2400.0, 4700.0, 2540.0],
                    [390.0, 2400.0, 4700.0, 2540.0],
                    [390.0, 2400.0, 4700.0, 2540.0],
                    [390.0, 2760.0, 6300.0, 3590.0],
                    [390.0, 2760.0, 6300.0, 3590.0],
                    [390.0, 2760.0, 6300.0, 3590.0],
                    [1520.0, 2960.0, 6900.0, 3940.0],
                    [1520.0, 2960.0, 6900.0, 3940.0],
                    [1520.0, 2960.0, 6900.0, 3940.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [292.0, 2400.0, 4700.0, 2540.0],
                    [292.0, 2400.0, 4700.0, 2540.0],
                    [292.0, 2400.0, 4700.0, 2540.0],
                    [292.0, 2400.0, 4700.0, 2540.0],
                    [292.0, 2760.0, 6300.0, 3590.0],
                    [292.0, 2760.0, 6300.0, 3590.0],
                    [292.0, 2760.0, 6300.0, 3590.0],
                    [292.0, 2760.0, 6300.0, 3590.0],
                    [1140.0, 2960.0, 6900.0, 3940.0],
                    [1140.0, 2960.0, 6900.0, 3940.0],
                    [1140.0, 2960.0, 6900.0, 3940.0],
                    [1140.0, 2960.0, 6900.0, 3940.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )
    # elif sta == "A422A":
    else:
        # sediment_thickness = 0
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 4, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [2000.0, 2400.0, 4700.0, 2000.0],
                    [2000.0, 2760.0, 6300.0, 2500.0],
                    [3000.0, 2760.0, 6300.0, 3200.0],
                    [3000.0, 2960.0, 6900.0, 3800.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 5, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [1000.0, 2400.0, 4700.0, 2000.0],
                    [1000.0, 2400.0, 4700.0, 2000.0],
                    [1000.0, 2760.0, 6300.0, 2500.0],
                    [1000.0, 2760.0, 6300.0, 2500.0],
                    [1500.0, 2960.0, 6900.0, 3200.0],
                    [1500.0, 2960.0, 6900.0, 3200.0],
                    [1500.0, 2960.0, 6900.0, 3800.0],
                    [1500.0, 2960.0, 6900.0, 3800.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

    vs0 = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    for i in range(0, len(starting_model)):
        vs0[
            0,
            int(np.sum(starting_model[:, :, 0][:, 0][0:i])) : int(
                np.sum(starting_model[:, :, 0][:, 0][0 : i + 1])
            ),
            0,
        ] = starting_model[:, 3][i][0]

    vsi = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])

    # print("Depth above half-space is "+ str(dep) + " m")
    # print("Depth above half-space is "+ str(dep - starting_model[len(starting_model)-1][0][0]) + " m")
    return (starting_model, vs0, vsi)


# %%
def velp(vs, p=0.25):
    """Vp in m/s from Vs in m/s and a stable isotropic Poisson ratio."""
    vs = np.asarray(vs, float)
    if (
        not np.isfinite(p)
        or not -1 < p < 0.5
        or np.any(~np.isfinite(vs))
        or np.any(vs <= 0)
    ):
        raise ValueError(
            "Vs must be finite and positive; Poisson ratio must lie in (-1, .5)"
        )
    return vs * np.sqrt((1 - p) / (0.5 - p))


# %%


def density(vs, p=0.25):
    """Gardner density in kg/m³ for Vp expressed in m/s.

    Empirical relation; it is not automatically applied by the v2 sampler.
    """
    return 1000 * 0.31 * velp(vs, p) ** 0.25


# %%
def misfit(d, m, l=2, s=1):
    """Lp norm of residuals normalized by positive observation uncertainty."""
    return np.linalg.norm(residual(d, m, s).ravel(), ord=l)


# %%
def liklihood(d, m, k=1, s=1):
    """Unnormalized Gaussian likelihood; use log_likelihood for MCMC."""
    return positive(k, "normalization") * np.exp(log_likelihood(d, m, s))


# %%
def liklihood_all(
    d, m, vs, vs_prior, k=1, s=1, sm=1, alpha=1, beta=1, lamda=1, order=2
):
    """Gaussian data fit with profile, roughness and damping penalties."""
    weights = np.asarray([alpha, beta, lamda], float)
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("penalty weights must be finite and nonnegative")
    model_residual = residual(vs_prior, vs, sm)
    damping = residual(np.zeros_like(m), m, s)
    penalty = (
        alpha * Roughness(vs, order)
        + beta * np.sum(model_residual**2)
        + lamda * np.sum(damping**2)
    )
    return positive(k, "normalization") * np.exp(
        log_likelihood(d, m, s) - 0.5 * penalty
    )


# %%
def liklihood_roughness(d, m, vs, k=1, s=1, alpha=1, order=2):
    """Gaussian likelihood times a nonnegative roughness prior."""
    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")
    return positive(k, "normalization") * np.exp(
        log_likelihood(d, m, s) - 0.5 * alpha * Roughness(vs, order)
    )


# %%


# %%
# stable hyperbolic tangent
def dtanh(x):
    """Stable hyperbolic tangent, including x=0 and negative arguments."""
    return np.tanh(x)


# %%


# %%


def argdtray(wd, h):

    hh = np.sqrt(abs(h))  # % magnitude of wavenumber/freq
    th = wd * hh  # % number of waves (or e-foldings) in layer in radians
    if th >= 1.5e-14:
        if h <= 0:  # % propagating wave
            c = np.cos(th)
            s = -np.sin(th) / hh
        else:  # % evenescent wave
            d = np.exp(th)
            c = 0.5 * (d + 1 / d)
            s = -0.5 * (d - 1 / d) / hh
    else:
        c = 1
        s = -wd

    return c, s


# %%


# %%
def raydep(P, om, d, ro, vp2, vs2):
    """
        %RAYDEP	Propagator matrix sol'n for P-SV waves, minor vector method
    %   u   =  horizontal velocity AT TOP OF EACH LAYER
    %   v   =  vertical velocity AT TOP OF EACH LAYER
    %   sigzx = horizontal stress AT TOP OF EACH LAYER
    %   sigzz = vertical stress AT TOP OF EACH LAYER
    %   {  Normalized compliance = -k*v/(omega*sigzz)  }
    % [v u sigzz sigzx] = raydep(p,omega,d,ro,vp2,vs2)
    %   p     = slowness (s/m) of surface wave
    %   omega = angular frequency (radians/sec) of surface wave
    %   d     = thicknesses of the model layers (meters?)
    %   rho   = density of the layer (kg/m^3) (= gm/cc * 1000)
    %   vp2   = compressional velocity squared (m/s)^2
    %   vs2   = shear velocity squared (m/s)^2

    % W.C. Crawford 1-5-89
    %   x = stress-displacement vector:
    %		( vertical velocity,  horiz velocity, vert stress, horiz stress )
    %	y = minor vector matrix
    """

    mu = ro * vs2
    n = len(d)
    ist = n - 1
    ysav = 0
    psq = P * P
    r2 = 2 * mu[ist] * P
    # % R and S are the "Wavenumbers" of compress and shear waves in botlayer
    # % RoW and SoW are divided by ang freq
    RoW = np.sqrt(psq - 1 / vp2[ist])
    SoW = np.sqrt(psq - 1 / vs2[ist])
    ym = np.zeros((ist + 1, 5))
    i = ist
    y = np.zeros((5,))
    x = np.zeros((i + 1, 4))
    y[3 - 1] = RoW
    y[4 - 1] = -SoW
    y[1 - 1] = (RoW * SoW - psq) / ro[i]
    y[2 - 1] = r2 * y[1 - 1] + P
    y[5 - 1] = ro[i] - r2 * (P + y[2 - 1])
    ym[i, :] = y
    # %*****PROPAGATE UP LAYERS*********
    while i > 0:
        i = i - 1
        ha = psq - 1 / vp2[i]
        ca, sa = argdtray(om * d[i], ha)
        hb = psq - 1 / vs2[i]
        cb, sb = argdtray(om * d[i], hb)
        hbs = hb * sb
        has = ha * sa
        r1 = 1 / ro[i]
        r2 = 2 * mu[i] * P
        b1 = r2 * y[1 - 1] - y[2 - 1]
        g3 = (y[5 - 1] + r2 * (y[2 - 1] - b1)) * r1
        g1 = b1 + P * g3
        g2 = ro[i] * y[1 - 1] - P * (g1 + b1)
        e1 = cb * g2 - hbs * y[3 - 1]
        e2 = -sb * g2 + cb * y[3 - 1]
        e3 = cb * y[4 - 1] + hbs * g3
        e4 = sb * y[4 - 1] + cb * g3
        y[3 - 1] = ca * e2 - has * e4
        y[4 - 1] = sa * e1 + ca * e3
        g3 = ca * e4 - sa * e2
        b1 = g1 - P * g3
        y[1 - 1] = (ca * e1 + has * e3 + P * (g1 + b1)) * r1
        y[2 - 1] = r2 * y[1 - 1] - b1
        y[5 - 1] = ro[i] * g3 - r2 * (y[2 - 1] - b1)
        ym[i, :] = y

    de = y[5 - 1] / np.sqrt(y[1 - 1] * y[1 - 1] + y[2 - 1] * y[2 - 1])
    ynorm = 1 / y[3 - 1]
    y[1 - 1 : 4] = np.array([0, -ynorm, 0, 0])
    # %*****PROPAGATE BACK DOWN LAYERS*********
    while i <= ist:
        x[i, 1 - 1] = (
            -ym[i, 2 - 1] * y[1 - 1] - ym[i, 3 - 1] * y[2 - 1] + ym[i, 1 - 1] * y[4 - 1]
        )
        x[i, 2 - 1] = (
            -ym[i, 4 - 1] * y[1 - 1] + ym[i, 2 - 1] * y[2 - 1] - ym[i, 1 - 1] * y[3 - 1]
        )
        x[i, 3 - 1] = (
            -ym[i, 5 - 1] * y[2 - 1] - ym[i, 2 - 1] * y[3 - 1] - ym[i, 4 - 1] * y[4 - 1]
        )
        x[i, 4 - 1] = (
            ym[i, 5 - 1] * y[1 - 1] - ym[i, 3 - 1] * y[3 - 1] + ym[i, 2 - 1] * y[4 - 1]
        )
        ls = i
        if i >= 2 - 1:
            sum = np.hypot(x[i, 0], x[i, 1])
            pbsq = 1 / vs2[i]
            if sum < 1e-4:
                break

        ha = psq - 1 / vp2[i]
        ca, sa = argdtray(om * d[i], ha)
        hb = psq - 1 / vs2[i]
        cb, sb = argdtray(om * d[i], hb)
        hbs = hb * sb
        has = ha * sa
        r2 = 2 * P * mu[i]
        e2 = r2 * y[2 - 1] - y[3 - 1]
        e3 = ro[i] * y[2 - 1] - P * e2
        e4 = r2 * y[1 - 1] - y[4 - 1]
        e1 = ro[i] * y[1 - 1] - P * e4
        e6 = ca * e2 - sa * e1
        e8 = cb * e4 - sb * e3
        y[1 - 1] = (ca * e1 - has * e2 + P * e8) / ro[i]
        y[2 - 1] = (cb * e3 - hbs * e4 + P * e6) / ro[i]
        y[3 - 1] = r2 * y[2 - 1] - e6
        y[4 - 1] = r2 * y[1 - 1] - e8
        i = i + 1
    #
    # if x(1,3) == 0
    #  error('vertical surface stress = 0 in DETRAY');
    # end
    ist = ls
    v = x[:, 1 - 1]
    u = x[:, 2 - 1]
    zz = x[:, 3 - 1]
    zx = x[:, 4 - 1]
    return v, u, zz, zx


# %%
def calc_norm_compliance(depth, freq, model):
    """calculate normalized compliance for a model and water depth

    model=[thick(m) rho(kg/m^3) vp(m/s) vs(m/s)]; last row is half-space
    """

    model = _validate_model(model)
    depth = positive(depth, "water depth")
    freq = np.asarray(freq, float)
    if (
        freq.ndim != 1
        or not len(freq)
        or np.any(~np.isfinite(freq))
        or np.any(freq <= 0)
    ):
        raise ValueError("frequencies must be a nonempty positive finite vector")
    thick = model[:, 0]
    rho = model[:, 1]
    vpsq = model[:, 2] * model[:, 2]
    vssq = model[:, 3] * model[:, 3]
    omega = 2 * np.pi * freq
    k = gravd(omega, depth)
    p = k / omega
    ncomp = np.zeros((len(p)))

    for i in np.arange((len(p))):
        v, u, sigzz, sigzx = raydep(p[i], omega[i], thick, rho, vpsq, vssq)
        ncomp[i] = -k[i] * v[1 - 1] / (omega[i] * sigzz[1 - 1])

    #   u   =  horizontal velocity AT TOP OF EACH LAYER
    #   v   =  vertical velocity AT TOP OF EACH LAYER
    #   sigzx = horizontal stress AT TOP OF EACH LAYER
    #   sigzz = vertical stress AT TOP OF EACH LAYER

    if not np.all(np.isfinite(ncomp)):
        raise ValueError("model is outside the finite real-valued propagator domain")
    return ncomp


# %%
def plot_inversion_v2(
    starting_model,
    vs,
    mis_fit,
    ncompl,
    Data,
    likelihood_data,
    freq,
    sta,
    iteration,
    s,
    sigma_v,
    sigma_h,
    n_layer,
    alpha,
    burnin=50000,
    mis_fit_trsh=1,
):

    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)
    # var_vs = np.zeros([1,vs.shape[0]])

    # for i in range(0, vs.shape[0]):
    #     var_vs[0,i] = np.var(vs[i])
    layer_number = -2
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    vs_filtered = []

    for i in range(burnin, iteration):
        if mis_fit[0][i] < mis_fit_trsh:
            # if starting_model[layer_number][3][i] < 50000:
            Vs_final = Vs_final + vs[i]
            N = N + 1
            vs_filtered.append(vs[i])
    vs_filtered = np.array(vs_filtered)

    Vs_final = Vs_final / N
    print(N)

    plt.figure(dpi=300, figsize=(15, 25))
    plt.plot(Vs_final, depth, color="black", label="Final Result", linewidth=3)

    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    dd = np.zeros([len(vs[0]), 1000])

    for i in range(0, len(vs[0])):
        dd[i] = np.histogram(vs_filtered[:, i], bins=5000, range=([0, 5000]))[0]
        dd[i] = dd[i] / np.max(dd[i])
    plt.rcParams.update({"font.size": 40})
    nn = (
        int((iteration - burnin) / 1000) + 1
    )  # I want to see just 1000 points of the data
    # nn = 1
    plt.figure(dpi=300, figsize=(35, 25))
    plt.suptitle(
        "Sigma V = "
        + str(sigma_v)
        + ", Sigma h = "
        + str(sigma_h)
        + ", Alpha = "
        + str(alpha)
        + ",N of Layer ="
        + str(n_layer)
    )
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color="green", linewidth=0.25)

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Normalized Compliance")

        # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(
        freq,
        Data,
        yerr=s,
        ecolor=("black"),
        color="black",
        fmt="none",
        linewidth=5,
        label="Measured Compliance",
        capsize=15,
    )

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
    #             label='Median of Brun-in')
    plt.plot(
        freq,
        ncompl[1],
        color="black",
        label="Start Compliance",
        linewidth=5,
        linestyle="dashed",
    )
    # plt.ylim([10e-14,10e-10])
    plt.ylim([10e-13, 10e-11])

    plt.grid(True)
    plt.legend(loc="upper left", fontsize=30)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')

    plt.subplot(122)
    plt.title(sta)
    plt.imshow(dd, aspect="auto")
    plt.colorbar()
    # plt.plot(vs[0], depth, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

    # Add labels, legend, and colorbar as in your code
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    plt.figure(dpi=300, figsize=(16, 24))
    plt.subplot(211)
    # plt.plot(likeli_hood[0, 1:iteration-1],color='Blue',label="Total Likelihood")
    # plt.plot(likelihood_Model[0, 1:iteration-1],color='red',label="Model Likelihood")
    plt.plot(
        likelihood_data[0, 1 : iteration - 1], color="green", label="Data Likelihood"
    )
    plt.xscale("log")
    plt.vlines(
        x=burnin, ymin=0, ymax=1, color="r", label="Burn-in Region", linestyles="dashed"
    )
    plt.legend(loc="upper left")

    plt.title("Likelihood")
    plt.xlabel("Iteration")
    plt.ylabel("Liklihood")

    plt.grid(True)

    plt.subplot(212)
    plt.plot(mis_fit[0, 1 : iteration - 1])
    plt.xscale("log")
    # plt.yscale('log')

    plt.vlines(
        x=burnin,
        ymin=0,
        ymax=np.max(mis_fit[0, 1 : iteration - 1]),
        color="r",
        label="Burn-in Region",
        linestyles="dashed",
        linewidth=3,
    )

    plt.title("Data Misfit")
    plt.xlabel("Iteration")
    plt.ylabel("$\\chi^2$")
    plt.grid(True)
    plt.minorticks_on()

    plt.tight_layout()

    vs_burnin = np.zeros([iteration, int(np.sum(starting_model[:, 0, 0])), 1])
    plt.legend(loc="lower left")
    plt.tight_layout()
    print(N)


# %%
def plot_inversion(
    starting_model,
    vs,
    mis_fit,
    ncompl,
    Data,
    likelihood_data,
    freq,
    sta,
    iteration,
    s,
    sigma_v,
    sigma_h,
    n_layer,
    alpha,
    burnin=50000,
    mis_fit_trsh=1,
):
    linewidth_models = 0.05
    linewidth_compliance = 0.5
    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)
    # var_vs = np.zeros([1,vs.shape[0]])

    # for i in range(0, vs.shape[0]):
    #     var_vs[0,i] = np.var(vs[i])
    layer_number = -2
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    for i in range(burnin, iteration):
        if mis_fit[0][i] < mis_fit_trsh:
            # if starting_model[layer_number][3][i] < 50000:
            Vs_final = Vs_final + vs[i]
            N = N + 1

    Vs_final = Vs_final / N
    print(N)

    plt.figure(dpi=300, figsize=(15, 25))
    plt.title(sta)

    plt.plot(Vs_final, depth, color="black", label="Final Result", linewidth=3)

    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    # plt.show()  # Don't forget to show the plot

    # plt.rcParams.update({'font.size': 30})
    # plt.figure(dpi=300, figsize=(10, 10))
    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:
    #         plt.plot(vs[i], depth, color='grey', linewidth= 1)

    plt.rcParams.update({"font.size": 40})
    nn = int((iteration - burnin) / 1000)  # I want to see just 1000 points of the data
    # nn = 1
    plt.figure(dpi=300, figsize=(35, 25))
    plt.suptitle(
        "Sigma V = "
        + str(sigma_v)
        + ", Sigma h = "
        + str(sigma_h)
        + ", Alpha = "
        + str(alpha)
        + ",N of Layer ="
        + str(n_layer)
    )
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color="green", linewidth=linewidth_compliance)

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Normalized Compliance")

        # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(
        freq,
        Data,
        yerr=s,
        ecolor=("black"),
        color="black",
        fmt="none",
        linewidth=5,
        label="Measured Compliance",
        capsize=15,
    )

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
    #             label='Median of Brun-in')
    plt.plot(
        freq,
        ncompl[1],
        color="black",
        label="Start Compliance",
        linewidth=5,
        linestyle="dashed",
    )
    # plt.ylim([10e-14,10e-10])
    # plt.ylim([10e-13,10e-11])

    plt.grid(True)
    plt.legend(loc="upper left", fontsize=30)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')

    plt.subplot(122)
    plt.title("YV." + sta)
    # plt.figure(dpi=300, figsize=(35, 25))
    for i in range(burnin, vs.shape[0]):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(vs[i], depth, color="green", linewidth=linewidth_models)
    plt.plot(
        vs[0],
        depth,
        color="black",
        label="Start Model",
        linewidth=5,
        linestyle="dashed",
    )

    # Add labels, legend, and colorbar as in your code
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    # plt.show()  # Don't forget to show the plot

    # plt.rcParams.update({'font.size': 30})
    # plt.figure(dpi=300, figsize=(10, 10))
    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:
    #         plt.plot(vs[i], depth, color='grey', linewidth= 1)

    # plt.plot(Vs_final, depth, color='black', label='Final Result',linewidth=3)
    # plt.plot(vs[0], depth, color='green', label='Start Model',linewidth=5,linestyle='dashed')

    # plt.grid(True)
    # plt.xlabel('Shear Velocity [m/s]')
    # plt.ylabel('Depth [m]')
    # plt.ylim([-int(np.sum(starting_model[:, 0, 0]))+2000, 0])

    # plt.ylim([-10000, 0])
    # plt.xlim([0,4500])
    plt.figure(dpi=300, figsize=(16, 24))

    plt.subplot(211)
    plt.suptitle("YV." + sta)

    # plt.plot(likeli_hood[0, 1:iteration-1],color='Blue',label="Total Likelihood")
    # plt.plot(likelihood_Model[0, 1:iteration-1],color='red',label="Model Likelihood")
    plt.plot(
        likelihood_data[0, 1 : iteration - 1], color="green", label="Data Likelihood"
    )
    plt.xscale("log")
    plt.vlines(
        x=burnin, ymin=0, ymax=1, color="r", label="Burn-in Region", linestyles="dashed"
    )
    plt.legend(loc="upper left")

    plt.title("Likelihood")
    plt.xlabel("Iteration")
    plt.ylabel("Liklihood")

    plt.grid(True)

    plt.subplot(212)
    plt.plot(mis_fit[0, 1 : iteration - 1])
    plt.xscale("log")
    # plt.yscale('log')

    plt.vlines(
        x=burnin,
        ymin=0,
        ymax=np.max(mis_fit[0, 1 : iteration - 1]),
        color="r",
        label="Burn-in Region",
        linestyles="dashed",
        linewidth=3,
    )

    plt.title("Data Misfit")
    plt.xlabel("Iteration")
    plt.ylabel("$\\chi^2$")
    plt.grid(True)
    plt.minorticks_on()

    plt.tight_layout()

    vs_burnin = np.zeros([iteration, int(np.sum(starting_model[:, 0, 0])), 1])
    plt.legend(loc="lower left")
    plt.tight_layout()
    print(N)


# %%
def plot_inversion_beta(
    starting_model,
    vs,
    vs0,
    mis_fit,
    ncompl,
    Data,
    likelihood_data,
    freq,
    sta,
    iteration,
    s,
    sigma_v,
    sigma_h,
    n_layer,
    alpha,
    burnin=50000,
    mis_fit_trsh=1,
):
    linewidth_models = 0.05
    linewidth_compliance = 0.5
    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)
    # var_vs = np.zeros([1,vs.shape[0]])

    # for i in range(0, vs.shape[0]):
    #     var_vs[0,i] = np.var(vs[i])
    layer_number = -2
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    for i in range(burnin, iteration):
        if mis_fit[0][i] < mis_fit_trsh:
            # if starting_model[layer_number][3][i] < 50000:
            Vs_final = Vs_final + vs[i]
            N = N + 1

    Vs_final = Vs_final / N
    print(N)

    plt.figure(dpi=300, figsize=(15, 25))
    plt.title(sta)

    plt.plot(Vs_final, depth, color="black", label="Final Result", linewidth=3)

    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.xlim([0, 5000])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    # plt.show()  # Don't forget to show the plot

    # plt.rcParams.update({'font.size': 30})
    # plt.figure(dpi=300, figsize=(10, 10))
    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:
    #         plt.plot(vs[i], depth, color='grey', linewidth= 1)

    plt.rcParams.update({"font.size": 40})

    nn = int((iteration - burnin) / 1000)  # I want to see just 1000 points of the data
    # nn = 1
    plt.figure(dpi=300, figsize=(35, 25))
    plt.suptitle(
        "Sigma V = "
        + str(sigma_v)
        + ", Sigma h = "
        + str(sigma_h)
        + ", Alpha = "
        + str(alpha)
        + ",N of Layer ="
        + str(n_layer)
    )
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            if np.mean(vs[i]) < np.mean(vs0[:, 0 : len(vs[0])]):
                plt.plot(freq, ncompl[i], color="red", linewidth=linewidth_compliance)
            else:
                plt.plot(freq, ncompl[i], color="green", linewidth=linewidth_compliance)

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Normalized Compliance")

        # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(
        freq,
        Data,
        yerr=s,
        ecolor=("black"),
        color="black",
        fmt="none",
        linewidth=5,
        label="Measured Compliance",
        capsize=15,
    )

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
    #             label='Median of Brun-in')
    plt.plot(
        freq,
        ncompl[1],
        color="black",
        label="Start Compliance",
        linewidth=5,
        linestyle="dashed",
    )
    # plt.ylim([10e-14,10e-10])
    # plt.ylim([10e-13,10e-11])

    plt.grid(True)
    plt.legend(loc="upper left", fontsize=30)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')

    plt.subplot(122)
    plt.title("YV." + sta)
    # plt.figure(dpi=300, figsize=(35, 25))
    for i in range(burnin, vs.shape[0]):
        if mis_fit[0][i] < mis_fit_trsh:
            if np.mean(vs[i]) < np.mean(vs0[:, 0 : len(vs[0])]):
                plt.plot(vs[i], depth, color="red", linewidth=linewidth_models)
            else:
                plt.plot(vs[i], depth, color="green", linewidth=linewidth_models)

    plt.plot(
        vs[0],
        depth,
        color="black",
        label="Start Model",
        linewidth=5,
        linestyle="dashed",
    )

    # Add labels, legend, and colorbar as in your code
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    # plt.show()  # Don't forget to show the plot

    # plt.rcParams.update({'font.size': 30})
    # plt.figure(dpi=300, figsize=(10, 10))
    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:
    #         plt.plot(vs[i], depth, color='grey', linewidth= 1)

    # plt.plot(Vs_final, depth, color='black', label='Final Result',linewidth=3)
    # plt.plot(vs[0], depth, color='green', label='Start Model',linewidth=5,linestyle='dashed')

    # plt.grid(True)
    # plt.xlabel('Shear Velocity [m/s]')
    # plt.ylabel('Depth [m]')
    # plt.ylim([-int(np.sum(starting_model[:, 0, 0]))+2000, 0])

    # plt.ylim([-10000, 0])
    # plt.xlim([0,4500])
    plt.figure(dpi=300, figsize=(16, 24))

    plt.subplot(211)
    plt.suptitle("YV." + sta)

    # plt.plot(likeli_hood[0, 1:iteration-1],color='Blue',label="Total Likelihood")
    # plt.plot(likelihood_Model[0, 1:iteration-1],color='red',label="Model Likelihood")
    plt.plot(
        likelihood_data[0, 1 : iteration - 1], color="green", label="Data Likelihood"
    )
    plt.xscale("log")
    plt.vlines(
        x=burnin, ymin=0, ymax=1, color="r", label="Burn-in Region", linestyles="dashed"
    )
    plt.legend(loc="upper left")

    plt.title("Likelihood")
    plt.xlabel("Iteration")
    plt.ylabel("Liklihood")

    plt.grid(True)

    plt.subplot(212)
    plt.plot(mis_fit[0, 1 : iteration - 1])
    plt.xscale("log")
    # plt.yscale('log')

    plt.vlines(
        x=burnin,
        ymin=0,
        ymax=np.max(mis_fit[0, 1 : iteration - 1]),
        color="r",
        label="Burn-in Region",
        linestyles="dashed",
        linewidth=3,
    )

    plt.title("Data Misfit")
    plt.xlabel("Iteration")
    plt.ylabel("$\\chi^2$")
    plt.grid(True)
    plt.minorticks_on()

    plt.tight_layout()

    vs_burnin = np.zeros([iteration, int(np.sum(starting_model[:, 0, 0])), 1])
    plt.legend(loc="lower left")
    plt.tight_layout()
    print(N)


# %%


def plot_inversion_density(
    vs, vs0, mis_fit, Data, s, freq, sta, burnin, ncompl, iteration, mis_fit_trsh=8
):
    start_model = start_model_plot(sta) / 50
    # start_model = start_model_plot_mean(sta)/50
    bins = 100
    jj = 0
    vs_good = []

    for ii in range(0, len(mis_fit[0])):
        if mis_fit[0][ii] < mis_fit_trsh:
            # if vs[ii][9000][0] < vs[ii][6000][0]:

            vs_good.append(vs[ii])
            jj = jj + 1

    vs_good = np.array(vs_good)

    linewidth_compliance = 0.5
    nn = int((iteration - burnin) / 1000)

    # a = []
    # b = plt.hist(vs[:,100,0],bins=100,range=([0,5000]))[1]
    # vs_good = np.array(vs_good)

    # for i in range(0,len(vs[0])):
    #     a.append(plt.hist(vs_good[:,i,0],bins=100,range=([0,5000]))[0])
    #     print(len(vs[0]) - i)

    # plt.figure(dpi=300,figsize=(20,20))
    # selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
    # selected_labels = b[selected_indices]
    # selected_labels = selected_labels.astype(int)

    # plt.imshow(a / 1000, aspect='auto', cmap='jet', norm=plt.Normalize(vmin=0, vmax=1.5))

    # # Setting custom y-axis ticks and labels
    # plt.xticks(ticks=selected_indices, labels=selected_labels)

    # plt.colorbar()
    # plt.show()

    vs_good = vs_good[:, 0:10000, 0]

    downsample_factor = vs_good.shape[1] // 100

    # Reshape the array to prepare for averaging
    # New shape will be (20000, 100, downsample_factor, 1)
    vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

    # Take the mean along the downsample_factor dimension
    vs_downsampled = vs_reshaped.mean(axis=2)

    print(vs_downsampled.shape)  # This should print (20000, 100, 1)

    # a = np.zeros([100,100])
    c = np.zeros([bins, bins])
    b = plt.hist(vs_reshaped[:, 10, 0], bins=bins, range=([0, 5000]))[1]

    for i in range(0, len(vs_downsampled[0])):
        c[i] = plt.hist(
            vs_downsampled[:, i], bins=bins, range=([0, 5000]), density=True, log=False
        )[0]
        c[i] = c[i] / np.max(c[i])
        print(i)

    from matplotlib.colors import LinearSegmentedColormap

    # Define the colors for the colormap (from white to red to black)
    colors = ["white", "red", "black"]

    # Define the transition points for the colors
    n_bins = [0, 0.4, 1]  # You can adjust these thresholds based on your data

    # Create the custom colormap
    custom_colormap = LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(n_bins, colors))
    )

    plt.figure(dpi=300, figsize=(30, 20))
    plt.suptitle("YV." + str(sta))
    plt.subplot(121)
    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color="green", linewidth=linewidth_compliance)

    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized Compliance")

    # plt.plot(freq, Data, color='black', label='Measured Compliance')
    plt.errorbar(
        freq,
        Data,
        yerr=s,
        ecolor=("black"),
        color="black",
        fmt="none",
        linewidth=5,
        label="Measured Compliance",
        capsize=15,
    )
    # plt.plot(freq, ncompl[1],color='blue', label='Starting Compliance',linewidth= 5,linestyle='dashed')
    plt.legend(loc="upper left", fontsize=30)

    plt.subplot(122)
    selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
    selected_labels = b[selected_indices]
    selected_labels = selected_labels.astype(int)

    depth = np.arange(0, 100, 1)

    ff = plt.imshow(
        c, aspect="auto", cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1)
    )
    # plt.plot(start_model,depth,color='blue',linewidth=8,linestyle='dashed',label="Starting Model ")

    # plt.imshow(a , aspect='auto', cmap='jet')

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
    #             label='Median of Brun-in')
    # plt.ylim([10e-14,10e-10])
    # plt.ylim([10e-13,10e-11])

    plt.grid(True)
    # plt.ylim([1e-11,6e-11])
    # plt.yscale('log')

    # Setting custom y-axis ticks and labels
    plt.xticks(ticks=selected_indices, labels=selected_labels)
    cbar = plt.colorbar(ff)
    # depth = np.arange(0, -10000, -1)

    # plt.plot(vs[0][0:10000], depth/100, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

    selected_indices_y = [0, 12, 24, 36, 48, 60, 72, 84]
    selected_labels_y = [0, -2000, -4000, -6000, -8000, -10000, -12000, -14000]
    plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m]")
    cbar.set_label("Probability")
    plt.legend(loc="lower left", fontsize=35)
    plt.ylim(72, 0)
    # plt.colorbar()
    plt.tight_layout()
    plt.show()


# %%


def plot_inversion_density_all(Inversion_container):

    plt.figure(dpi=300, figsize=(40, 25))

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        # vs0 = Inversion_container[ii]["Shear Velocity Starting"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        # Data = Inversion_container[ii]["compliance Measured"]
        # s = Inversion_container[ii]["uncertainty"]
        # freq = Inversion_container[ii]["compliance Frequency"]
        sta = Inversion_container[ii]["Station"]
        # burnin = Inversion_container[ii]["burnin"]
        # ncompl = Inversion_container[ii]["compliance Forward"]
        # iteration = Inversion_container[ii]["iteration"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        start_model = start_model_plot(sta) / 50

        start_model = (
            start_model_plot_mean(sta) / 50
        )  # mean of all models,sedimental models differ from rocky models

        # start_model = start_model_plot_mean(sta)/50
        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        print(vs_downsampled.shape)  # This should print (20000, 100, 1)

        # a = np.zeros([100,100])
        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            # c[i] = plt.hist(vs_downsampled[:,i],bins=bins,range=([0,5000]),density=True,log=False)[0]
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])
            print(i)

        from matplotlib.colors import LinearSegmentedColormap

        # Define the colors for the colormap (from white to red to black)
        colors = ["white", "red", "black"]
        viridis = plt.cm.get_cmap("viridis", 256)

        # Sample the colormap
        colors = [
            viridis(0.0),  # dark blue/lilac at the start of the colormap
            viridis(0.5),  # greenish in the middle of the colormap
            viridis(1.0),
        ]  # bright yellow at the end of the colormap

        # Define the transition points for the colors
        n_bins = [0, 0.4, 1]  # You can adjust these thresholds based on your data

        # Create the custom colormap
        custom_colormap = LinearSegmentedColormap.from_list(
            "custom_cmap", list(zip(n_bins, colors))
        )

        plt.subplot(2, 4, int(ii + 1))
        plt.title(str("YV.") + sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        ff = plt.imshow(
            c, aspect="auto", cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1)
        )
        # plt.plot(start_model,depth,color='blue',linewidth=8,linestyle='dashed',label="Mean Starting Model ")
        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            color="red",
            linewidth=8,
            linestyle="dashed",
            label="Median",
        )
        # plt.imshow(a , aspect='auto', cmap='jet')

        # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
        #             label='Median of Brun-in')
        # plt.ylim([10e-14,10e-10])
        # plt.ylim([10e-13,10e-11])

        plt.grid(True)
        # plt.ylim([1e-11,6e-11])
        # plt.yscale('log')

        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)
        cbar = plt.colorbar(ff)
        # depth = np.arange(0, -10000, -1)

        # plt.plot(vs[0][0:10000], depth/100, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

        selected_indices_y = [0, 12, 24, 36, 48, 60, 72, 84]
        selected_labels_y = [0, -2000, -4000, -6000, -8000, -10000, -12000, -14000]
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == 4:
            plt.ylabel("Depth [m]")
        if ii == 4 or ii == 5 or ii == 6 or ii == 7:
            plt.xlabel("Shear Velocity [m/s]")
        if ii == 3 or ii == 7:
            cbar.set_label("Probability")
        if ii == 7:
            plt.legend(loc="lower left", fontsize=30)

        if ii == 0:
            plt.text(
                0.5, 0.95, "a)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 1:
            plt.text(
                0.05,
                0.95,
                "b)",
                fontsize=40,
                fontweight="bold",
                va="top",
                color="white",
            )
        if ii == 2:
            plt.text(
                0.001,
                0.95,
                "c)",
                fontsize=40,
                fontweight="bold",
                va="top",
                color="white",
            )
        if ii == 3:
            plt.text(
                0.4, 0.95, "d)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 4:
            plt.text(
                0.9, 0.95, "e)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 5:
            plt.text(
                0.1, 0.95, "f)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 6:
            plt.text(
                0.1, 0.95, "g)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 7:
            plt.text(
                0.1, 0.95, "h)", fontsize=40, fontweight="bold", va="top", color="white"
            )

        plt.ylim(48, 0)
        # plt.ylim(72,0)
        # plt.colorbar()
    plt.tight_layout()


# %%
def final_plot(Inversion_container, *, image_dir):
    """Publication layout; image_dir contains the six named bathymetry images."""
    from pathlib import Path

    image_dir = Path(image_dir)
    from matplotlib import gridspec
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import matplotlib.image as mpimg  # Import the image module

    # Create a new figure with specified dimensions
    fig = plt.figure(figsize=(35, 30), constrained_layout=True)

    # Define the grid layout with customized row heights for the third column
    gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 2], height_ratios=[1, 1, 1])

    img1 = mpimg.imread(image_dir / "RR36.png")
    img2 = mpimg.imread(image_dir / "RR38.png")
    img3 = mpimg.imread(image_dir / "RR40.png")
    img4 = mpimg.imread(image_dir / "RR50.png")
    img5 = mpimg.imread(image_dir / "RR52.png")
    img6 = mpimg.imread(image_dir / "Rift_valley.png")

    # First plot
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title("RR36")
    ax1.imshow(img1, extent=[0, 100, 0, 100], aspect="auto")
    ax1.set_xlabel("Distance East/West [km]")
    ax1.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "a)",
        transform=ax1.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Second plot
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title("RR38")
    ax2.imshow(img2, extent=[0, 100, 0, 100], aspect="auto")
    ax2.set_xlabel("Distance East/West [km]")
    ax2.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "b)",
        transform=ax2.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Third plot
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_title("RR40")
    ax3.imshow(img3, extent=[0, 100, 0, 100], aspect="auto")
    ax3.set_xlabel("Distance East/West [km]")
    ax3.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "c)",
        transform=ax3.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Fourth plot
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.set_title("RR50")
    ax4.imshow(img4, extent=[0, 100, 0, 100], aspect="auto")
    ax4.set_xlabel("Distance East/West [km]")
    ax4.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "d)",
        transform=ax4.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Fifth plot
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.set_title("RR52")
    ax5.imshow(img5, extent=[0, 100, 0, 100], aspect="auto")
    ax5.set_xlabel("Distance East/West [km]")
    ax5.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "e)",
        transform=ax5.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Sixth plot
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.imshow(img6, extent=[-25, 25, -5, 5], aspect="auto")

    ax6.set_xlabel("Distance [km]")
    ax6.set_ylabel("Depth [km]")
    plt.text(
        0.05,
        0.95,
        "f)",
        transform=ax6.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )
    ax6.grid(False)

    # Seventh plot
    ax7 = fig.add_subplot(gs[0:3, 2])  # Spanning the top two rows of the third column

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        # vs0 = Inversion_container[ii]["Shear Velocity Starting"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        # Data = Inversion_container[ii]["compliance Measured"]
        # s = Inversion_container[ii]["uncertainty"]
        # freq = Inversion_container[ii]["compliance Frequency"]
        sta = Inversion_container[ii]["Station"]
        # burnin = Inversion_container[ii]["burnin"]
        # ncompl = Inversion_container[ii]["compliance Forward"]
        # iteration = Inversion_container[ii]["iteration"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        # start_model = start_model_plot(sta)/50
        # refrence_model = refrence_models()/50
        # start_model = start_model_plot_mean(sta)/50 #mean of all models,sedimental models differ from rocky models

        # start_model = start_model_plot_mean(sta)/50
        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        print(vs_downsampled.shape)  # This should print (20000, 100, 1)

        # a = np.zeros([100,100])
        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            # c[i] = plt.hist(vs_downsampled[:,i],bins=bins,range=([0,5000]),density=True,log=False)[0]
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])
            print(i)

        # plt.title(str("YV.")+sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        # ff = plt.imshow(c , aspect='auto', cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1))
        plt.vlines(90, 0, 100, linestyles="dashed", colors="grey", linewidth=3)
        plt.text(
            95,
            25,
            "Fresh peridotite",
            fontsize=35,
            color="grey",
            ha="center",
            rotation=90,
        )

        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            linewidth=5,
            linestyle="solid",
            label=sta,
        )
        # medians = np.median(vs_downsampled, axis=0) / 50
        # std_dev = np.std(vs_downsampled, axis=0) / 50

        # plt.errorbar(np.median(vs_downsampled,axis=0)/50, depth, xerr=np.std(vs_downsampled,axis=0)/50)
        # plt.fill_betweenx(depth, medians - std_dev, medians + std_dev, color='gray', alpha=0.5)
        plt.grid(True)

        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)

        selected_indices_y = [0, 12, 24, 36, 48, 60, 72, 84]
        selected_labels_y = [0, -2000, -4000, -6000, -8000, -10000, -12000, -14000]
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        plt.legend(loc="upper right")
        plt.ylim(48, 0)

        ax7.set_xlabel("Shear Velocity [m/s]")
        ax7.set_ylabel("Depth [m]")
        plt.text(
            0.025,
            0.98,
            "g)",
            transform=ax7.transAxes,
            fontsize=40,
            fontweight="bold",
            va="top",
        )

    # Eighth plot as inset
    ax8 = inset_axes(ax7, width="40%", height="40%", loc="lower left", borderpad=3.5)

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        # vs0 = Inversion_container[ii]["Shear Velocity Starting"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        # Data = Inversion_container[ii]["compliance Measured"]
        # s = Inversion_container[ii]["uncertainty"]
        # freq = Inversion_container[ii]["compliance Frequency"]
        sta = Inversion_container[ii]["Station"]
        # burnin = Inversion_container[ii]["burnin"]
        # ncompl = Inversion_container[ii]["compliance Forward"]
        # iteration = Inversion_container[ii]["iteration"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        # start_model = start_model_plot(sta)/50
        # refrence_model = refrence_models()/50
        # start_model = start_model_plot_mean(sta)/50 #mean of all models,sedimental models differ from rocky models

        # start_model = start_model_plot_mean(sta)/50
        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        print(vs_downsampled.shape)  # This should print (20000, 100, 1)

        # a = np.zeros([100,100])
        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            # c[i] = plt.hist(vs_downsampled[:,i],bins=bins,range=([0,5000]),density=True,log=False)[0]
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])
            print(i)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        serpentinization_percentage = (
            9000 - 2 * (np.median(vs_downsampled, axis=0))
        ) / 30
        serpentinization_percentage = np.minimum(serpentinization_percentage, 100)
        plt.plot(
            serpentinization_percentage,
            depth,
            linewidth=5,
            linestyle="solid",
            label=sta,
        )

        # medians = np.median(vs_downsampled, axis=0) / 50
        # std_dev = np.std(vs_downsampled, axis=0) / 50

        # plt.errorbar(np.median(vs_downsampled,axis=0)/50, depth, xerr=np.std(vs_downsampled,axis=0)/50)
        # plt.fill_betweenx(depth, serpentinization_percentage - std_dev, serpentinization_percentage + std_dev, color='blue', alpha=0.5)

        plt.grid(True)
        plt.xlim([0, 100])
        selected_indices_y = [0, 12, 24, 36, 48, 60, 72, 84]
        selected_labels_y = [0, -2000, -4000, -6000, -8000, -10000, -12000, -14000]
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        plt.ylim(48, 0)
        plt.legend(loc="lower right")

        ax8.set_xlabel("Serpentinization (%)")
        # ax8.set_ylabel("Depth (m)")
    plt.text(
        0.05,
        0.95,
        "h)",
        transform=ax8.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    plt.tight_layout()
    plt.show()


# %%
def plot_inversion_all(Inversion_container):
    refrence_model = refrence_models() / 50

    plt.figure(dpi=300, figsize=(15, 25))

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        # vs0 = Inversion_container[ii]["Shear Velocity Starting"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        # Data = Inversion_container[ii]["compliance Measured"]
        # s = Inversion_container[ii]["uncertainty"]
        # freq = Inversion_container[ii]["compliance Frequency"]
        sta = Inversion_container[ii]["Station"]
        # burnin = Inversion_container[ii]["burnin"]
        # ncompl = Inversion_container[ii]["compliance Forward"]
        # iteration = Inversion_container[ii]["iteration"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        start_model = start_model_plot(sta) / 50

        start_model = (
            start_model_plot_mean(sta) / 50
        )  # mean of all models,sedimental models differ from rocky models

        # start_model = start_model_plot_mean(sta)/50
        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        print(vs_downsampled.shape)  # This should print (20000, 100, 1)

        # a = np.zeros([100,100])
        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            # c[i] = plt.hist(vs_downsampled[:,i],bins=bins,range=([0,5000]),density=True,log=False)[0]
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])
            print(i)

        # plt.title(str("YV.")+sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        # plt.plot(start_model,depth,color='blue',linewidth=8,linestyle='dashed',label="Mean Starting Model ")
        plt.plot(np.median(vs_downsampled, axis=0) / 50, depth, linewidth=5, label=sta)
        # plt.imshow(a , aspect='auto', cmap='jet')

        # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
        #             label='Median of Brun-in')
        # plt.ylim([10e-14,10e-10])
        # plt.ylim([10e-13,10e-11])

        plt.grid(True)
        # plt.ylim([1e-11,6e-11])
        # plt.yscale('log')

        # Setting custom y-axis ticks and labels
    plt.xticks(ticks=selected_indices, labels=selected_labels)
    # depth = np.arange(0, -10000, -1)

    # plt.plot(vs[0][0:10000], depth/100, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

    selected_indices_y = [0, 12, 24, 36, 48, 60, 72, 84]
    selected_labels_y = [0, -2000, -4000, -6000, -8000, -10000, -12000, -14000]
    plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
    plt.ylabel("Depth [m]")
    plt.xlabel("Shear Velocity [m/s]")
    plt.plot(
        refrence_model[0],
        depth,
        color="red",
        linewidth=5,
        linestyle="dashed",
        label="SWIR 64°30'E (EW)",
    )
    plt.plot(
        refrence_model[1],
        depth,
        color="purple",
        linewidth=5,
        linestyle="dotted",
        label="SWIR 64°30'E (NS) ",
    )
    plt.plot(
        refrence_model[2],
        depth,
        color="green",
        linewidth=5,
        linestyle="dotted",
        label="SWIR 65-64°E ",
    )
    plt.legend(loc="lower left", fontsize=30)
    plt.ylim(48, 0)
    # plt.ylim(72,0)
    # plt.colorbar()
    plt.tight_layout()


# %%


def plot_inversion_density_mean_all(Inversion_container):

    plt.figure(dpi=300, figsize=(40, 25))

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        # vs0 = Inversion_container[ii]["Shear Velocity Starting"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        # Data = Inversion_container[ii]["compliance Measured"]
        # s = Inversion_container[ii]["uncertainty"]
        # freq = Inversion_container[ii]["compliance Frequency"]
        sta = Inversion_container[ii]["Station"]
        # burnin = Inversion_container[ii]["burnin"]
        # ncompl = Inversion_container[ii]["compliance Forward"]
        # iteration = Inversion_container[ii]["iteration"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        start_model = start_model_plot(sta) / 50
        refrence_model = refrence_models() / 50
        # start_model = start_model_plot_mean(sta)/50 #mean of all models,sedimental models differ from rocky models

        # start_model = start_model_plot_mean(sta)/50
        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        print(vs_downsampled.shape)  # This should print (20000, 100, 1)

        # a = np.zeros([100,100])
        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            # c[i] = plt.hist(vs_downsampled[:,i],bins=bins,range=([0,5000]),density=True,log=False)[0]
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])
            print(i)

        # Define the colors for the colormap (from white to red to black)

        # Define the transition points for the colors

        # Create the custom colormap

        plt.subplot(2, 4, int(ii + 1))
        plt.title(str("YV.") + sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        # ff = plt.imshow(c , aspect='auto', cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1))
        plt.plot(
            start_model,
            depth,
            color="blue",
            linewidth=5,
            linestyle="dashdot",
            label="Crust 1 ",
        )
        if sta == "RR40":
            plt.plot(
                refrence_model[0],
                depth,
                color="red",
                linewidth=5,
                linestyle="dashed",
                label="SWIR 64°30'E (EW)",
            )
            plt.plot(
                refrence_model[1],
                depth,
                color="purple",
                linewidth=5,
                linestyle="dotted",
                label="SWIR 64°30'E (NS) ",
            )
            plt.plot(
                refrence_model[2],
                depth,
                color="green",
                linewidth=5,
                linestyle="dotted",
                label="SWIR 65-64°E ",
            )

        if sta == "RR36" or sta == "RR38":
            plt.plot(
                refrence_model[4],
                depth,
                color="red",
                linewidth=5,
                linestyle="dashed",
                label="Atlantis Bank",
            )

        # plt.plot(refrence_models[2],depth,color='blue',linewidth=8,linestyle='dashed',label="Crust 1 ")

        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            color="black",
            linewidth=6,
            linestyle="solid",
            label="Median (This Study)",
        )
        # plt.imshow(a , aspect='auto', cmap='jet')

        # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
        #             label='Median of Brun-in')
        # plt.ylim([10e-14,10e-10])
        # plt.ylim([10e-13,10e-11])

        plt.grid(True)
        # plt.ylim([1e-11,6e-11])
        # plt.yscale('log')

        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)
        # cbar = plt.colorbar(ff)
        # depth = np.arange(0, -10000, -1)

        # plt.plot(vs[0][0:10000], depth/100, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

        selected_indices_y = [0, 12, 24, 36, 48, 60, 72, 84]
        selected_labels_y = [0, -2000, -4000, -6000, -8000, -10000, -12000, -14000]
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == 4:
            plt.ylabel("Depth [m]")
        if ii == 4 or ii == 5 or ii == 6 or ii == 7:
            plt.xlabel("Shear Velocity [m/s]")
        # if ii == 3 or ii == 7:
        # cbar.set_label('Probability')
        # if ii == 7:
        # plt.legend(loc='lower left',fontsize=30)
        plt.legend(loc="lower left", fontsize=30)
        plt.ylim(48, 0)
        # plt.ylim(72,0)
        # plt.colorbar()
    plt.tight_layout()
    # plt.savefig(file_path_save + "Inversion_Table.pdf")


# %%
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plot_inversion_serpentinization(Inversion_container):

    plt.figure(dpi=300, figsize=(45, 25))

    for ii in range(len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        sta = Inversion_container[ii]["Station"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        vs_good = []
        for i in range(len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                vs_good.append(vs[i])

        vs_good = np.array(vs_good)[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)
        vs_downsampled = vs_reshaped.mean(axis=2)

        ax_main = plt.subplot(2, len(Inversion_container), ii + 1)
        plt.title(str("YV.") + sta)
        plt.vlines(90, 0, 100, linestyles="dashed", colors="grey", linewidth=3)
        plt.text(
            95,
            25,
            "Fresh peridotite",
            fontsize=35,
            color="grey",
            ha="center",
            rotation=90,
        )

        depth = np.arange(0, 100, 1)
        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            color="black",
            linewidth=5,
            linestyle="solid",
            label="Median",
        )
        medians = np.median(vs_downsampled, axis=0) / 50
        std_dev = np.std(vs_downsampled, axis=0) / 50
        plt.fill_betweenx(
            depth, medians - std_dev, medians + std_dev, color="gray", alpha=0.5
        )
        plt.grid(True)
        # plt.legend(loc='lower left', fontsize=30)
        plt.ylim(100, 0)
        plt.xlabel("Shear Velocity [m/s]")
        plt.ylabel("Depth [m]")

        # Adding inset plot at a new position
        ax_inset = inset_axes(
            ax_main, width="30%", height="50%", loc="lower left", borderpad=1.5
        )
        serpentinization_percentage = (
            9000 - 2 * np.median(vs_downsampled, axis=0)
        ) / 30
        serpentinization_percentage = np.minimum(serpentinization_percentage, 100)
        ax_inset.plot(
            serpentinization_percentage,
            depth,
            color="blue",
            linewidth=5,
            linestyle="solid",
        )
        ax_inset.fill_betweenx(
            depth,
            serpentinization_percentage - std_dev,
            serpentinization_percentage + std_dev,
            color="blue",
            alpha=0.5,
        )
        ax_inset.set_xlabel("Serpentinization [%]")
        ax_inset.set_ylabel("Depth [m]")
        ax_inset.set_ylim(100, 0)
        ax_inset.grid(True)

    plt.tight_layout()
    plt.show()


# Assuming Inversion_container is predefined somewhere in your code with the necessary data
# plot_inversion_serpentinization(Inversion_container)


# %%
def plot_inversion_serpentinization1(Inversion_container):

    plt.figure(dpi=300, figsize=(45, 25))

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        # vs0 = Inversion_container[ii]["Shear Velocity Starting"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        # Data = Inversion_container[ii]["compliance Measured"]
        # s = Inversion_container[ii]["uncertainty"]
        # freq = Inversion_container[ii]["compliance Frequency"]
        sta = Inversion_container[ii]["Station"]
        # burnin = Inversion_container[ii]["burnin"]
        # ncompl = Inversion_container[ii]["compliance Forward"]
        # iteration = Inversion_container[ii]["iteration"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        # start_model = start_model_plot(sta)/50
        # refrence_model = refrence_models()/50
        # start_model = start_model_plot_mean(sta)/50 #mean of all models,sedimental models differ from rocky models

        # start_model = start_model_plot_mean(sta)/50
        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        print(vs_downsampled.shape)  # This should print (20000, 100, 1)

        # a = np.zeros([100,100])
        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            # c[i] = plt.hist(vs_downsampled[:,i],bins=bins,range=([0,5000]),density=True,log=False)[0]
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])
            print(i)

        # Define the colors for the colormap (from white to red to black)

        # Define the transition points for the colors
        n_bins = [0, 0.4, 1]  # You can adjust these thresholds based on your data

        # Create the custom colormap
        plt.subplot(2, len(Inversion_container), int(ii + 1))
        plt.title(str("YV.") + sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        # ff = plt.imshow(c , aspect='auto', cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1))
        plt.vlines(90, 0, 100, linestyles="dashed", colors="grey", linewidth=3)
        plt.text(
            95,
            25,
            "Fresh peridotite",
            fontsize=35,
            color="grey",
            ha="center",
            rotation=90,
        )

        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            color="black",
            linewidth=5,
            linestyle="solid",
            label="Median",
        )
        medians = np.median(vs_downsampled, axis=0) / 50
        std_dev = np.std(vs_downsampled, axis=0) / 50

        # plt.errorbar(np.median(vs_downsampled,axis=0)/50, depth, xerr=np.std(vs_downsampled,axis=0)/50)
        plt.fill_betweenx(
            depth, medians - std_dev, medians + std_dev, color="gray", alpha=0.5
        )
        plt.grid(True)

        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)
        # cbar = plt.colorbar(ff)
        # depth = np.arange(0, -10000, -1)

        # plt.plot(vs[0][0:10000], depth/100, color='black', label='Start Model',linewidth= 5 ,linestyle='dashed')

        selected_indices_y = [0, 12, 24, 36, 48, 60, 72, 84]
        selected_labels_y = [0, -2000, -4000, -6000, -8000, -10000, -12000, -14000]
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == 4:
            plt.ylabel("Depth [m]")
        if ii in range(len(Inversion_container)):
            plt.xlabel("Shear Velocity [m/s]")
        # if ii == 3 or ii == 7:
        # cbar.set_label('Probability')
        # if ii == 7:
        # plt.legend(loc='lower left',fontsize=30)
        plt.legend(loc="lower left", fontsize=30)
        plt.ylim(48, 0)
        # plt.ylim(72,0)
        # plt.colorbar()

        plt.subplot(2, len(Inversion_container), int(ii + len(Inversion_container) + 1))
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        # ff = plt.imshow(c , aspect='auto', cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1))
        serpentinization_percentage = (
            9000 - 2 * (np.median(vs_downsampled, axis=0))
        ) / 30
        serpentinization_percentage = np.minimum(serpentinization_percentage, 100)
        plt.plot(
            serpentinization_percentage,
            depth,
            color="blue",
            linewidth=5,
            linestyle="solid",
        )

        medians = np.median(vs_downsampled, axis=0) / 50
        std_dev = np.std(vs_downsampled, axis=0) / 50

        # plt.errorbar(np.median(vs_downsampled,axis=0)/50, depth, xerr=np.std(vs_downsampled,axis=0)/50)
        plt.fill_betweenx(
            depth,
            serpentinization_percentage - std_dev,
            serpentinization_percentage + std_dev,
            color="blue",
            alpha=0.5,
        )

        plt.grid(True)
        plt.xlim([0, 100])
        selected_indices_y = [0, 12, 24, 36, 48, 60, 72, 84]
        selected_labels_y = [0, -2000, -4000, -6000, -8000, -10000, -12000, -14000]
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == len(Inversion_container):
            plt.ylabel("Depth [m]")
        if ii in range(len(Inversion_container)):
            plt.xlabel("Serpentinization [%]")
        # if ii == 3 or ii == 7:
        # cbar.set_label('Probability')
        # if ii == 7:
        # plt.legend(loc='lower left',fontsize=30)
        # plt.legend(loc='lower left',fontsize=30)
        plt.ylim(48, 0)
        # plt.ylim(72,0)
        # plt.colorbar()

    plt.tight_layout()


# %%
def start_model_plot_mean(sta):
    start_model = np.zeros(100)
    if sta == "RR28" or sta == "RR29" or sta == "RR34":
        start_model[0:2] = 340
        start_model[2:7] = 2540
        start_model[7:16] = 3590
        start_model[16:44] = 3940
        start_model[44:100] = 4310

    # elif  sta == "RR29":
    #     start_model[0:1] = 340
    #     start_model[1:6] = 2540
    #     start_model[6:15] = 3590
    #     start_model[15:43] = 3940
    #     start_model[43:100] = 4310

    elif (
        sta == "RR52"
        or sta == "RR50"
        or sta == "RR40"
        or sta == "RR38"
        or sta == "RR36"
    ):
        start_model[0:5] = 2540
        start_model[5:14] = 3590
        start_model[14:42] = 3940
        start_model[42:100] = 4310

    return start_model


# %%
def refrence_models(p=0.25):
    # Possion Ratio
    ref_models = np.zeros([5, 100])

    # SWIR-64-EW_profile
    ref_models[0][0:2] = 3500
    ref_models[0][2:3] = 4070
    ref_models[0][3:5] = 4730
    ref_models[0][5:7] = 5200
    ref_models[0][7:10] = 5900
    ref_models[0][10:13] = 6200
    ref_models[0][13:15] = 6600
    ref_models[0][15:17] = 6900
    ref_models[0][17:20] = 7120
    ref_models[0][20:23] = 7230
    ref_models[0][23:39] = 7600
    ref_models[0][39:100] = 7780

    # SWIR-64-NS_profile.txt
    ref_models[1][0:1] = 2938
    ref_models[1][1:4] = 4214
    ref_models[1][4:7] = 4994
    ref_models[1][7:10] = 5881
    ref_models[1][10:13] = 6414
    ref_models[1][13:19] = 7018
    ref_models[1][19:23] = 7374
    ref_models[1][23:28] = 7589
    ref_models[1][28:35] = 7770
    ref_models[1][35:100] = 7846

    # SWIR65-66E.txt
    ref_models[2][0:12] = 3680
    ref_models[2][12:34] = 6540
    ref_models[2][34:100] = 6910

    # BestFit
    ref_models[3][0:1] = 3500
    ref_models[3][1:2] = 3730
    ref_models[3][2:3] = 3910
    ref_models[3][3:5] = 4060
    ref_models[3][5:7] = 4230
    ref_models[3][7:8] = 4420
    ref_models[3][8:19] = 4630
    ref_models[3][9:10] = 4870
    ref_models[3][10:11] = 5080
    ref_models[3][11:12] = 5280
    ref_models[3][12:13] = 5460
    ref_models[3][13:14] = 5610

    ref_models[3][14:15] = 5760
    ref_models[3][15:16] = 5910
    ref_models[3][16:17] = 6060
    ref_models[3][17:18] = 6210
    ref_models[3][18:20] = 6360
    ref_models[3][20:21] = 6520
    ref_models[3][21:22] = 6670
    ref_models[3][22:23] = 6850

    ref_models[3][23:24] = 7030
    ref_models[3][24:25] = 7200
    ref_models[3][25:27] = 7370
    ref_models[3][27:28] = 7510
    ref_models[3][28:35] = 7640
    ref_models[3][35:100] = 8000

    # Atlantis Bank M uller 1997
    ref_models[4][0:1] = 5900
    ref_models[4][1:4] = 6000
    ref_models[4][4:7] = 6200
    ref_models[4][7:10] = 6300
    ref_models[4][10:13] = 6414
    ref_models[4][13:19] = 6500
    ref_models[4][19:23] = 6600
    ref_models[4][23:28] = 6789
    ref_models[4][28:30] = 6800
    ref_models[4][30:100] = 8000

    return ref_models / np.sqrt((1 - p) / (0.5 - p))


# %%
def start_model_plot(sta):
    start_model = np.zeros(100)
    if sta == "RR28":
        start_model[0:1] = 340
        start_model[1:6] = 2700
        start_model[6:16] = 3700
        start_model[16:45] = 4050
        start_model[45:99] = 4510

    elif sta == "RR29":
        start_model[0:1] = 340
        start_model[1:5] = 2700
        start_model[5:14] = 3700
        start_model[14:44] = 4050
        start_model[44:99] = 4510

    elif sta == "RR34":
        start_model[0:2] = 340
        start_model[2:7] = 2700
        start_model[7:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4500

    elif sta == "RR36":
        start_model[0:4] = 2700
        start_model[4:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4370

    elif sta == "RR38":
        start_model[0:4] = 2700
        start_model[4:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4190

    elif sta == "RR40":
        start_model[0:4] = 2700
        start_model[4:13] = 3700
        start_model[13:42] = 4050
        start_model[42:100] = 4360

    elif sta == "RR50":
        start_model[0:7] = 2540
        start_model[7:14] = 3590
        start_model[14:42] = 3940
        start_model[42:100] = 4190

    elif sta == "RR52":
        start_model[0:7] = 2660
        start_model[7:14] = 3670
        start_model[14:42] = 4022
        start_model[42:100] = 4330

    else:
        # elif  sta == "A422A":
        start_model[0:13] = 2200
        start_model[13:25] = 2700
        start_model[25:42] = 3500
        start_model[42:60] = 3800
        start_model[60:100] = 4310
    return start_model


# %%
def plot_hist(starting_model, burnin, mis_fit, mis_fit_trsh):

    starting_model_opt = []
    for i in range(burnin, len(starting_model[0][0][:]) - 1):
        if mis_fit[0][i] < mis_fit_trsh:
            starting_model_opt.append(starting_model[:, :, i])
    starting_model_opt = np.array(starting_model_opt)
    plt.rcParams.update({"font.size": 40})

    plt.figure(dpi=300, figsize=(35, 70))
    plt.suptitle("Shear Velocity Distribution of Layers", y=0.99)
    for i in range(0, (len(starting_model_opt[0]) - 2)):
        plt.subplot((len(starting_model_opt[0]) - 1), 4, i + 1)

        # plt.hist(starting_model_opt[i-1,3,burnin:-1], 100, density=True, facecolor='k', alpha=1)
        plt.hist(
            starting_model_opt[:, i, 3],
            bins=20,
            density=True,
            histtype="barstacked",
            facecolor="k",
            alpha=1,
            stacked=True,
            label="Layer " + str(i + 1),
        )
        plt.xlim([0, 4500])
        plt.legend(loc="upper left")
        # plt.title('Layer ' + str(i+1))
        plt.tight_layout()

    plt.figure(dpi=300, figsize=(35, 70))
    plt.suptitle("Thickness Distribution of Layers", y=0.99)
    for i in range(0, (len(starting_model_opt[0]) - 2)):
        plt.subplot((len(starting_model_opt[0]) - 1), 4, i + 1)

        # plt.hist(starting_model_opt[i-1,3,burnin:-1], 100, density=True, facecolor='k', alpha=1)
        plt.hist(
            starting_model_opt[:, i, 0],
            bins=20,
            histtype="barstacked",
            facecolor="k",
            alpha=1,
            stacked=True,
            label="Layer " + str(i + 1),
        )
        plt.xlim([0, 2000])
        plt.legend(loc="upper left")
        # plt.title('Layer ' + str(i+1))
        plt.tight_layout()


# %%
def plot_hist2d(starting_model, burnin, mis_fit, mis_fit_trsh):

    starting_model_opt = []
    for i in range(burnin, len(starting_model[0][0][:]) - 1):
        if mis_fit[0][i] < mis_fit_trsh:
            starting_model_opt.append(starting_model[:, :, i])
    starting_model_opt = np.array(starting_model_opt)

    plt.rcParams.update({"font.size": 40})
    plt.figure(dpi=300, figsize=(50, 100))
    plt.suptitle("Shear Velocity and Thickness Distribution of Layers", y=0.99)
    for i in range(0, (len(starting_model_opt[0]) - 2)):
        plt.subplot((len(starting_model_opt[0]) - 1), 4, i + 1)

        # plt.hist(starting_model_opt[i-1,3,burnin:-1], 100, density=True, facecolor='k', alpha=1)
        plt.hist2d(
            starting_model_opt[:, i, 0],
            starting_model_opt[:, i, 3],
            density=True,
            bins=10,
            cmap="jet",
        )
        # bins = 50, density = True, histtype = "barstacked",facecolor='k', alpha=1,stacked = True,label='Layer ' + str(i+1))
        plt.xlabel("Thickness[m]")
        plt.ylabel("Vs [m/s]")
        # plt.legend(loc='upper left')
        # plt.title('Layer ' + str(i+1))
        plt.colorbar()
        plt.tight_layout()


# %%
def autocorreletion(starting_model, N):

    cc = np.zeros([len(starting_model), len(starting_model[0, 3, :])])

    for i in range(0, len(cc) - 1):
        s = starting_model[i, 3, :]
        dd = np.correlate(s - np.mean(s), s - np.mean(s), "full") / np.sum(
            (s - np.mean(s)) ** 2
        )
        # Auto-correlations.
        cc[i] = dd[N - 1 :]

    # Estimate of the effective sample size (Gelman et al., 2013).
    Neff = np.zeros([len(starting_model)])
    for i in range(0, len(cc)):
        for j in range(N - 1):
            if cc[i][j] + cc[i][j + 1] > 0.0:
                Neff[i] += cc[i][j]

    Neff = N / (1.0 + 2.0 * Neff)
    # for i in range(0,len(cc)):
    #     print('Effective Sample Size (parameter 1): %f' % Neff[i])

    # Plot autocorrelation function.
    # plt.figure(dpi = 300,figsize=(40,40))
    # plt.rcParams.update({'font.size': 30})
    # for i in range(0,len(cc)-1):
    #     plt.plot(cc[i][0:N],'k',linewidth=3,label='Neff = : %f' % Neff[i])
    #     plt.xlabel('Iteration',labelpad=15)
    #     plt.xlim([0,N])
    #     plt.title('Layer ' + str(i+1))
    #     plt.legend(loc ="upper right",fontsize = 15)
    #     plt.grid()
    #     plt.tight_layout()

    plt.figure(dpi=300, figsize=(30, 40))
    plt.rcParams.update({"font.size": 30})
    for i in range(0, len(cc) - 1):
        plt.subplot(len(cc) - 1, 3, i + 1)
        plt.plot(cc[i][0:N], "k", linewidth=3, label="Neff = : %f" % Neff[i])
        plt.xlabel("Iteration", labelpad=15)
        plt.xlim([0, N])
        plt.title("Layer " + str(i + 1))
        plt.legend(loc="upper right", fontsize=15)
        plt.grid()
        plt.tight_layout()

    plt.figure(dpi=300, figsize=(30, 40))
    plt.rcParams.update({"font.size": 30})
    plt.suptitle("Vs")
    for i in range(0, len(cc) - 1):
        plt.subplot(len(cc) - 1, 3, i + 1)
        plt.plot(starting_model[i, 3, :], "k", linewidth=2)
        plt.xlabel("Iteration", labelpad=15)
        plt.xlim([0, N])
        plt.title("Layer " + str(i + 1))
        plt.grid()
        plt.tight_layout()

    plt.figure(dpi=300, figsize=(30, 40))
    plt.rcParams.update({"font.size": 30})
    plt.suptitle("Thickness")
    for i in range(0, len(cc) - 1):
        plt.subplot(len(cc) - 1, 3, i + 1)
        plt.plot(starting_model[i, 0, :], "k", linewidth=2)
        plt.xlabel("Iteration", labelpad=15)
        plt.xlim([0, N])
        plt.title("Layer " + str(i + 1))
        plt.grid()
        plt.tight_layout()

    return Neff
