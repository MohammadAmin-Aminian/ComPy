"""Random-walk Metropolis inversion of normalized compliance (Aminian et al. 2025).

``invert_compliance`` is the sampler used in the GJI paper;
``invert_tempered`` wraps it in parallel tempering and ``l_curve`` scans
the roughness weight. Likelihoods and misfits are the ones the sampler
uses. For the newer DRAM sampler see ``compy.inversion.dram``.
"""

import matplotlib.pyplot as plt
import numpy as np

from compy._deprecation import renamed

from compy.inversion.models import (
    _depth_profile,
    _validate_model,
    crust1_starting_model,
    rhumrum_starting_model,
)
from compy.inversion.propagator import normalized_compliance
from compy.numerics import integer, log_likelihood, positive, residual, roughness


def _step_sizes(value, name):
    """Validate a proposal scale: a positive scalar or a 1D positive array."""
    if np.ndim(value) == 0:
        return positive(value, name)
    value = np.asarray(value, dtype=float)
    if (
        value.ndim != 1
        or value.size == 0
        or np.any(~np.isfinite(value))
        or np.any(value <= 0)
    ):
        raise ValueError(f"{name} must be finite and positive")
    return value


def _explicit_bounds(vs_bounds, initial):
    """Uniform Vs prior bounds for every row; the half-space is pinned."""
    try:
        lo, hi = vs_bounds
    except (TypeError, ValueError) as exc:
        raise ValueError("vs_bounds must be a (lower, upper) pair") from exc
    n = len(initial) - 1
    lo = (
        np.broadcast_to(np.asarray(lo, float), (n,))
        if np.ndim(lo) == 0
        else np.asarray(lo, float)
    )
    hi = (
        np.broadcast_to(np.asarray(hi, float), (n,))
        if np.ndim(hi) == 0
        else np.asarray(hi, float)
    )
    if lo.shape != (n,) or hi.shape != (n,):
        raise ValueError(f"vs_bounds entries must be scalars or have {n} values")
    if (
        np.any(~np.isfinite(lo))
        or np.any(~np.isfinite(hi))
        or np.any(lo <= 0)
        or np.any(hi <= lo)
    ):
        raise ValueError("vs_bounds must be finite with 0 < lower < upper")
    half = initial[-1, 3]
    return np.append(lo, half), np.append(hi, half)


def _broadcast_steps(value, n, name):
    """Return a fresh, writable per-layer copy of a validated proposal scale."""
    if np.ndim(value) == 0:
        return np.full(n, value, dtype=float)
    if len(value) != n:
        raise ValueError(f"{name} must be a scalar or have {n} entries")
    return np.array(value, dtype=float)


def _metropolis(
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
    *,
    adapt_until=0,
    target_acceptance=0.3,
    initial_state=None,
    diagnostics=None,
    vs_bounds=None,
    vary_thickness=True,
    alpha_bounds=None,
    interface_move="deepest",
    data_covariance=None,
    elastic_from_vs=None,
):
    """Symmetric random-walk Metropolis sampler with bounded model support.

    Optional, keyword-only extensions (defaults reproduce the v2 sampler
    exactly, draw for draw):

    sigma_v, sigma_h
        Scalars, or per-layer arrays: ``sigma_v`` has one entry per finite
        layer; ``sigma_h`` one entry per finite layer except the deepest
        (which absorbs every thickness transfer).
    adapt_until
        Number of initial iterations during which each move's step size is
        tuned by Robbins-Monro on the log scale towards ``target_acceptance``.
        Step sizes are frozen afterwards, so iterations at or beyond
        ``adapt_until`` form an ordinary Metropolis chain. Discard at least
        ``adapt_until`` iterations as burn-in.
    initial_state
        Model at which the chain starts. Prior bounds, Vp, density, the
        half-space and the conserved total finite thickness are still defined
        by ``starting_model``; the state must lie inside that support. Use
        different states per chain for meaningful convergence diagnostics.
    diagnostics
        Optional dict, filled with the final step sizes and the per-move
        acceptance rates after adaptation.
    vs_bounds
        Optional ``(lower, upper)`` uniform Vs prior in m/s, as scalars or
        arrays with one entry per finite layer (the half-space is fixed).
        Replaces the bounds otherwise derived from ``starting_model``, so the
        prior no longer depends on where a chain starts.
    vary_thickness
        If False, only Vs is sampled and the layer geometry stays fixed: the
        many-thin-layer "minimum structure" parameterization of Crawford et
        al. (1991), regularized through ``alpha``.
    alpha_bounds
        Optional ``(alpha_min, alpha_max)``. The smoothing weight then becomes
        an unknown (hierarchical Bayes, Malinverno & Briggs 2004) with a
        Jeffreys prior on that interval, starting at ``alpha``. The roughness
        term is a Gaussian prior on Vs jumps between the k = n_finite - 1
        interfaces, so alpha is Gibbs-sampled each iteration from its
        conditional Gamma(k/2, rate = roughness/2), truncated to the bounds.
        The data, not a hand-picked value, then set how smooth models are.
        The trace is returned in ``diagnostics["alpha"]``.
    interface_move
        ``"deepest"`` (default, v2 behaviour) transfers thickness between a
        layer and the deepest finite layer, which shifts every interface in
        between. ``"adjacent"`` moves a single interface by exchanging
        thickness between the two layers that share it, leaving all other
        interfaces in place; free-interface chains mix much faster with it.
    data_covariance
        Optional full data covariance (N x N, 1/Pa**2). Use it for errors that
        are correlated between frequencies, e.g. a pressure-gain uncertainty
        common to all bins: C = diag(sigma**2) + (g * eta)(g * eta)^T. The
        Gaussian likelihood then uses the Cholesky-whitened residual and the
        returned misfit is its norm; ``s`` is ignored and may be None.
    elastic_from_vs
        Optional callable ``vs -> (vp, rho)`` (m/s, kg/m**3), e.g.
        ``brocher_elastic``. Vp and density of every finite layer then follow
        Vs instead of staying fixed, which removes the stiffness floor and
        the weak or reversed shallow Vs sensitivity of a fixed-Vp column.
        The half-space keeps its given properties.
    """
    iteration = integer(iteration, "iteration", minimum=2)
    sigma_v, sigma_h = _step_sizes(sigma_v, "sigma_v"), _step_sizes(sigma_h, "sigma_h")
    adapt_until = integer(adapt_until, "adapt_until", minimum=0)
    if interface_move not in ("deepest", "adjacent"):
        raise ValueError("interface_move must be 'deepest' or 'adjacent'")
    if not 0 < target_acceptance < 1:
        raise ValueError("target_acceptance must lie strictly between 0 and 1")
    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")
    Data, f = np.asarray(Data, float), np.asarray(f, float)
    if Data.ndim != 1 or f.shape != Data.shape or np.any(f <= 0):
        raise ValueError("data and positive frequencies must be equal-length vectors")
    if data_covariance is not None:
        from scipy.linalg import solve_triangular

        cov = np.asarray(data_covariance, float)
        if (
            cov.shape != (len(Data), len(Data))
            or not np.all(np.isfinite(cov))
            or not np.allclose(cov, cov.T, rtol=1e-10, atol=0)
        ):
            raise ValueError("data_covariance must be a finite symmetric N x N matrix")
        try:
            chol = np.linalg.cholesky(cov)
        except np.linalg.LinAlgError as exc:
            raise ValueError("data_covariance must be positive definite") from exc

        def whitened(pred):
            return solve_triangular(chol, Data - pred, lower=True)

        def data_log_likelihood(pred):
            r = whitened(pred)
            return -0.5 * float(r @ r)

        def data_misfit(pred):
            return float(np.linalg.norm(whitened(pred)))
    else:
        if s is None:
            raise ValueError("s is required: supply positive measurement uncertainties")
        residual(Data, Data, s)

        def data_log_likelihood(pred):
            return log_likelihood(Data, pred, s)

        def data_misfit(pred):
            return np.linalg.norm(residual(Data, pred, s))

    if starting_model is None:
        initial = rhumrum_starting_model(1, n_layer=n_layer, sta=sta)[0][:, :, 0]
    elif isinstance(starting_model, dict):  # a compy.inversion.crust1_column
        initial = crust1_starting_model(starting_model, sublayers=max(1, n_layer // 3))
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
    n_finite = len(initial) - 1
    if vs_bounds is not None:
        lower, upper = _explicit_bounds(vs_bounds, initial)
    # Adaptation multiplies by exp(.), so scales stay positive and untouched
    # scales are passed to the RNG exactly as supplied.
    step_size = {
        "vs": _broadcast_steps(sigma_v, n_finite, "sigma_v"),
        "h": _broadcast_steps(sigma_h, max(n_finite - 1, 1), "sigma_h"),
    }
    tuned = {key: np.zeros(len(value), int) for key, value in step_size.items()}
    proposed = {key: np.zeros(len(value), int) for key, value in step_size.items()}
    kept = {key: np.zeros(len(value), int) for key, value in step_size.items()}

    def in_support(model):
        return bool(
            np.all(model[:-1, 0] > 0)
            and np.all(model[:, 3] >= lower)
            and np.all(model[:, 3] <= upper)
            and np.all(model[:, 3] < model[:, 2])
        )

    hierarchical = alpha_bounds is not None
    if hierarchical:
        from scipy.stats import gamma as gamma_dist

        try:
            a_lo, a_hi = (positive(v, "alpha_bounds") for v in alpha_bounds)
        except (TypeError, ValueError) as exc:
            raise ValueError("alpha_bounds must be two positive values") from exc
        if not a_lo < a_hi or not a_lo <= alpha <= a_hi:
            raise ValueError(
                "alpha_bounds must satisfy lower < upper and contain alpha"
            )
        n_jumps = max(n_finite - 1, 1)
        alpha_trace = np.empty(iteration)
    state = {"alpha": float(alpha)}

    def evaluate(model):
        pred = normalized_compliance(depth_s, f, model)
        profile = _depth_profile(model, depth_grid)
        log_data = data_log_likelihood(pred)
        log_posterior = log_data - 0.5 * state["alpha"] * roughness(profile, 2)
        return pred, profile, log_data, log_posterior

    def gibbs_alpha(profile):
        """Draw alpha from Gamma(k/2, rate=R/2) truncated to alpha_bounds."""
        rate = 0.5 * roughness(profile, 2)
        shape = 0.5 * n_jumps
        if rate <= np.finfo(float).tiny:  # flat profile: alpha^(k/2-1) on bounds
            lo, hi = a_lo**shape, a_hi**shape
            return (lo + rng.random() * (hi - lo)) ** (1 / shape)
        dist = gamma_dist(shape, scale=1 / rate)
        c_lo, c_hi = dist.cdf(a_lo), dist.cdf(a_hi)
        if c_hi - c_lo < 1e-12:  # all mass beyond one bound
            return a_lo if c_lo >= 1 - 1e-12 else a_hi
        return float(np.clip(dist.ppf(c_lo + rng.random() * (c_hi - c_lo)), a_lo, a_hi))

    def couple(model):
        if elastic_from_vs is not None:
            vp, rho = elastic_from_vs(model[:-1, 3])
            model[:-1, 2], model[:-1, 1] = vp, rho
        return model

    if elastic_from_vs is not None:
        initial = couple(initial)
        if not in_support(initial):
            raise ValueError(
                "elastic_from_vs must give Vp above Vs for the starting model"
            )
    if initial_state is None:
        current = initial.copy()
    else:
        current = couple(_validate_model(initial_state).copy())
        if (
            current.shape != initial.shape
            or not np.array_equal(current[:, 1:3], initial[:, 1:3])
            and elastic_from_vs is None
            or not np.array_equal(current[-1, 1:], initial[-1, 1:])
            or current[-1, 3] != initial[-1, 3]
            or not np.isclose(current[:-1, 0].sum(), initial[:-1, 0].sum())
            or not in_support(current)
        ):
            raise ValueError(
                "initial_state must match starting_model in shape, density, Vp, "
                "half-space Vs and total finite thickness, and lie within the "
                "Vs prior bounds"
            )
    pred, profile, log_data, log_post = evaluate(current)
    accepted = 0
    for i in range(iteration):
        if i:
            candidate = current.copy()
            layer = rng.integers(len(initial) - 1)  # fixed elastic half-space
            if not vary_thickness or rng.integers(2) == 0:
                move = "vs"
                candidate[layer, 3] += rng.normal(0, step_size[move][layer])
                if elastic_from_vs is not None and candidate[layer, 3] > 0:
                    couple(candidate)
                # Hold Vp and rho fixed, as specified by Aminian et al. (2025).
            else:
                # Transfer thickness from the deepest finite layer. Total
                # finite depth is conserved; the half-space is never altered.
                if len(initial) == 2:
                    chain[:, :, i] = current
                    prediction[i] = pred
                    misfits[0, i] = data_misfit(pred)
                    likelihood[0, i] = np.exp(log_data)
                    if profiles is not None:
                        profiles[i, :, 0] = profile
                    if hierarchical:
                        alpha_trace[i] = state["alpha"]
                    continue
                if layer == len(initial) - 2:
                    layer = rng.integers(len(initial) - 2)
                move = "h"
                step = rng.normal(0, step_size[move][layer])
                candidate[layer, 0] += step
                # Partner layer: the deepest finite layer (v2) or the layer
                # directly below, so only interface ``layer`` moves.
                candidate[layer + 1 if interface_move == "adjacent" else -2, 0] -= step
            took = False
            if in_support(candidate):
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
                    took = True
            # Out-of-support proposals are rejected, not resampled/reflected:
            # this keeps the Gaussian proposal symmetric at the boundaries.
            if i < adapt_until:
                tuned[move][layer] += 1
                gain = tuned[move][layer] ** -0.6  # diminishing Robbins-Monro gain
                step_size[move][layer] *= np.exp(gain * (took - target_acceptance))
            else:  # frozen kernel: record per-move acceptance only
                proposed[move][layer] += 1
                kept[move][layer] += took
            if hierarchical:  # Gibbs step for the smoothing weight
                state["alpha"] = gibbs_alpha(profile)
                log_post = log_data - 0.5 * state["alpha"] * roughness(profile, 2)
        if hierarchical:
            alpha_trace[i] = state["alpha"]
        chain[:, :, i] = current
        prediction[i] = pred
        misfits[0, i] = data_misfit(pred)
        likelihood[0, i] = np.exp(log_data)  # may underflow; never used for acceptance
        if profiles is not None:
            profiles[i, :, 0] = profile
    if diagnostics is not None:
        with np.errstate(invalid="ignore", divide="ignore"):
            diagnostics.update(
                sigma_v=step_size["vs"].copy(),
                sigma_h=step_size["h"].copy(),
                acceptance_vs=kept["vs"] / proposed["vs"],
                acceptance_h=kept["h"] / proposed["h"],
                adapt_until=adapt_until,
            )
            if hierarchical:
                diagnostics["alpha"] = alpha_trace
    return (
        chain,
        profiles,
        prior,
        misfits,
        prediction,
        likelihood,
        accepted / (iteration - 1),
    )


def invert_compliance(
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
    adapt_until=0,
    target_acceptance=0.3,
    initial_state=None,
    diagnostics=None,
    vs_bounds=None,
    vary_thickness=True,
    alpha_bounds=None,
    interface_move="deepest",
    data_covariance=None,
    elastic_from_vs=None,
):
    """Legacy name for the v2 sampler; see README for units and priors.

    starting_model is a layered model or a CRUST1.0 column
    (``compy.inversion.crust1_column``), which ``crust1_starting_model``
    turns into a model with n_layer // 3 sublayers per crustal layer. Without
    it, the deprecated ``rhumrum_starting_model`` template of ``sta`` is used. s must be supplied.
    Set return_profiles=False to avoid the large metre-grid profile array.
    adapt_until, target_acceptance, initial_state and diagnostics are optional
    sampler-tuning keywords documented in ``_inversion``.
    """
    chain, profiles, _, misfits, prediction, likelihood, rate = _metropolis(
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
        adapt_until=adapt_until,
        target_acceptance=target_acceptance,
        initial_state=initial_state,
        diagnostics=diagnostics,
        vs_bounds=vs_bounds,
        vary_thickness=vary_thickness,
        alpha_bounds=alpha_bounds,
        interface_move=interface_move,
        data_covariance=data_covariance,
        elastic_from_vs=elastic_from_vs,
    )
    return chain, profiles, misfits, prediction, likelihood, rate


def invert_compliance_broad_bounds(
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
    adapt_until=0,
    target_acceptance=0.3,
    initial_state=None,
    diagnostics=None,
    vs_bounds=None,
    vary_thickness=True,
    alpha_bounds=None,
    interface_move="deepest",
    data_covariance=None,
    elastic_from_vs=None,
):
    """Metropolis sampler with broader shallow Vs bounds; returns the seven-tuple of ``invert_compliance``.

    adapt_until, target_acceptance, initial_state and diagnostics are optional
    sampler-tuning keywords documented in ``_inversion``.
    """
    return _metropolis(
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
        adapt_until=adapt_until,
        target_acceptance=target_acceptance,
        initial_state=initial_state,
        diagnostics=diagnostics,
        vs_bounds=vs_bounds,
        vary_thickness=vary_thickness,
        alpha_bounds=alpha_bounds,
        interface_move=interface_move,
        data_covariance=data_covariance,
        elastic_from_vs=elastic_from_vs,
    )


def l_curve(
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


def misfit(d, m, l=2, s=1):
    """Lp norm of residuals normalized by positive observation uncertainty."""
    return np.linalg.norm(residual(d, m, s).ravel(), ord=l)


def likelihood(d, m, k=1, s=1):
    """Unnormalized Gaussian likelihood; use log_likelihood for MCMC."""
    return positive(k, "normalization") * np.exp(log_likelihood(d, m, s))


def regularized_likelihood(
    d, m, vs, vs_prior, k=1, s=1, sm=1, alpha=1, beta=1, lamda=1, order=2
):
    """Gaussian data fit with profile, roughness and damping penalties."""
    weights = np.asarray([alpha, beta, lamda], float)
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("penalty weights must be finite and nonnegative")
    model_residual = residual(vs_prior, vs, sm)
    damping = residual(np.zeros_like(m), m, s)
    penalty = (
        alpha * roughness(vs, order)
        + beta * np.sum(model_residual**2)
        + lamda * np.sum(damping**2)
    )
    return positive(k, "normalization") * np.exp(
        log_likelihood(d, m, s) - 0.5 * penalty
    )


def roughness_likelihood(d, m, vs, k=1, s=1, alpha=1, order=2):
    """Gaussian likelihood times a nonnegative roughness prior."""
    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")
    return positive(k, "normalization") * np.exp(
        log_likelihood(d, m, s) - 0.5 * alpha * roughness(vs, order)
    )


def _tempered_segment(job):
    """Advance one replica of ``invert_tempered`` (module level so it pickles)."""
    data, f, depth_s, start, s, cov, t, state, alpha, sv, sh, n, seed, adapt, kw = job
    info = {}
    chain, _, _, misfits, pred, _, rate = _metropolis(
        data,
        f,
        depth_s,
        start,
        None if s is None else s * np.sqrt(t),
        3,
        sv,
        sh,
        n + 1,
        alpha,
        "RR38",
        seed,
        True,
        False,
        adapt_until=adapt,
        initial_state=state,
        diagnostics=info,
        data_covariance=None if cov is None else cov * t,
        **kw,
    )
    alphas = info["alpha"][1:] if "alpha" in info else np.full(n, float(alpha))
    # Sample 0 repeats the incoming state; chi-square is reported at T = 1.
    return (
        chain[:, :, 1:],
        pred[1:],
        misfits[0, 1:] ** 2 * t,
        alphas,
        info["sigma_v"],
        info["sigma_h"],
        rate,
    )


def invert_tempered(
    Data,
    f,
    depth_s,
    starting_model,
    initial_states,
    *,
    s=None,
    data_covariance=None,
    temperatures=(1.0, 2.0, 4.0, 8.0, 16.0),
    n_ensembles=4,
    adapt_until=8000,
    segment=200,
    rounds=200,
    burn_rounds=30,
    sigma_v=100.0,
    sigma_h=50.0,
    alpha=0.25,
    seed=0,
    processes=None,
    **sampler_kw,
):
    """Parallel tempering (replica exchange) around the Metropolis sampler.

    Shallow compliance posteriors can be multimodal: e.g. slow-over-fast and
    fast-over-slow layering below the resolution depth fit equally well, and
    single random-walk chains rarely cross between them (high R-hat). Each of
    ``n_ensembles`` independent ensembles runs one replica per temperature on
    likelihood**(1/T) (for a Gaussian likelihood: uncertainties times sqrt(T),
    covariance times T) with the prior unchanged, and swaps adjacent replicas
    after every ``segment`` iterations with the replica-exchange acceptance
    (Geyer 1991; Sambridge 2014). The T = 1 replicas then sample the posterior,
    and the independent ensembles give the R-hat/ESS diagnostics.

    Each replica is first tuned for ``adapt_until`` iterations at its own
    temperature; step sizes then stay frozen, so every later segment is an
    ordinary Metropolis kernel. ``initial_states`` is a list of models;
    ensemble e at temperature slot j starts from entry (e + j) mod len.
    ``sampler_kw`` is passed to the sampler (vs_bounds, vary_thickness,
    alpha_bounds, elastic_from_vs, ...); the first ``burn_rounds`` rounds after
    tuning are discarded. ``processes`` sets the worker pool (1: serial).

    Returns a dict: ``chain`` (ensembles, rows, 4, draws) of the T = 1
    replicas, ``predictions``, ``chi2`` (at T = 1), ``alpha``,
    ``swap_acceptance`` per adjacent pair, ``temperatures`` and the tuned
    ``sigma_v``/``sigma_h`` per temperature slot.
    """
    temps = np.asarray(temperatures, float)
    if (
        temps.ndim != 1
        or temps.size < 1
        or temps[0] != 1.0
        or np.any(np.diff(temps) <= 0)
    ):
        raise ValueError("temperatures must increase strictly from 1.0")
    if (s is None) == (data_covariance is None):
        raise ValueError("give exactly one of s or data_covariance")
    for name, value in (
        ("n_ensembles", n_ensembles),
        ("segment", segment),
        ("rounds", rounds),
    ):
        integer(value, name)
    if not 0 <= burn_rounds < rounds:
        raise ValueError("burn_rounds must lie in [0, rounds)")
    states = [_validate_model(m).copy() for m in initial_states]
    if not states:
        raise ValueError("initial_states must contain at least one model")
    data, f = np.asarray(Data, float), np.asarray(f, float)
    s = None if s is None else np.broadcast_to(np.asarray(s, float), data.shape).copy()
    cov = None if data_covariance is None else np.asarray(data_covariance, float)
    start = _validate_model(starting_model)
    reps = [(e, j) for e in range(n_ensembles) for j in range(len(temps))]
    rng = np.random.default_rng(seed)

    def job(r, state, alpha_r, sv, sh, n, adapt):
        return (
            data,
            f,
            depth_s,
            start,
            s,
            cov,
            temps[r[1]],
            state,
            alpha_r,
            sv,
            sh,
            n,
            int(rng.integers(2**31)),
            adapt,
            sampler_kw,
        )

    if processes == 1:
        pool, run = None, lambda jobs: [_tempered_segment(j) for j in jobs]
    else:
        from multiprocessing import Pool

        pool = Pool(processes)
        run = lambda jobs: pool.map(_tempered_segment, jobs)  # noqa: E731
    try:
        out = run(
            [
                job(
                    r,
                    states[(r[0] + r[1]) % len(states)],
                    alpha,
                    sigma_v,
                    sigma_h,
                    max(adapt_until, 1),
                    adapt_until,
                )
                for r in reps
            ]
        )
        cur = {r: o[0][:, :, -1].copy() for r, o in zip(reps, out)}
        alph = {r: float(o[3][-1]) for r, o in zip(reps, out)}
        chi2 = {r: float(o[2][-1]) for r, o in zip(reps, out)}
        sv = {r: o[4] for r, o in zip(reps, out)}  # tied to the temperature slot
        sh = {r: o[5] for r, o in zip(reps, out)}
        keep = [dict(chain=[], pred=[], chi2=[], alpha=[]) for _ in range(n_ensembles)]
        tried, swapped = np.zeros(len(temps) - 1), np.zeros(len(temps) - 1)
        for rnd in range(rounds):
            out = run([job(r, cur[r], alph[r], sv[r], sh[r], segment, 0) for r in reps])
            for r, o in zip(reps, out):
                cur[r], alph[r], chi2[r] = (
                    o[0][:, :, -1].copy(),
                    float(o[3][-1]),
                    float(o[2][-1]),
                )
                if r[1] == 0 and rnd >= burn_rounds:
                    k = keep[r[0]]
                    k["chain"].append(o[0])
                    k["pred"].append(o[1])
                    k["chi2"].append(o[2])
                    k["alpha"].append(o[3])
            for e in range(n_ensembles):  # alternate even and odd adjacent pairs
                for j in range(rnd % 2, len(temps) - 1, 2):
                    a, b = (e, j), (e, j + 1)
                    log_acc = (chi2[a] - chi2[b]) * (
                        0.5 / temps[j] - 0.5 / temps[j + 1]
                    )
                    tried[j] += 1
                    if np.log(max(rng.random(), np.finfo(float).tiny)) < log_acc:
                        swapped[j] += 1
                        for d in (cur, alph, chi2):
                            d[a], d[b] = d[b], d[a]
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    with np.errstate(invalid="ignore"):
        swap_rate = swapped / tried
    return dict(
        chain=np.stack([np.concatenate(k["chain"], axis=2) for k in keep]),
        predictions=np.stack([np.concatenate(k["pred"]) for k in keep]),
        chi2=np.stack([np.concatenate(k["chi2"]) for k in keep]),
        alpha=np.stack([np.concatenate(k["alpha"]) for k in keep]),
        swap_acceptance=swap_rate,
        temperatures=temps,
        sigma_v=[sv[(0, j)] for j in range(len(temps))],
        sigma_h=[sh[(0, j)] for j in range(len(temps))],
    )


__getattr__ = renamed(
    __name__, {"invert_compliance_beta": "invert_compliance_broad_bounds"}
)
