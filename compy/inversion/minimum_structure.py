"""Minimum-structure Vs inversion after Crawford (2004, GJI 157, section 3.3).

Separate from ``compy.inversion.occam``, which smooths ln(Vs) and linearizes in
ln(compliance). Here, as in Crawford's equations 10 to 14, the model is Vs itself and
the data are compliance:

    m = ((W J)^T W J + mu D^T D)^-1 (W J)^T W (d - F(m0) + J m0)

where D holds second differences of Vs and W = diag(1 / sigma). Each iteration scans
mu and keeps the best fit until chi2 reaches the target, then the smoothest model
with chi2 at or below it (Constable et al. 1987). Layer thicknesses, Vp, density and
the half-space stay fixed.

Crawford gives no bounds. Here each solve is a bounded least-squares problem, with Vs
between ``vs_min`` and Vp / sqrt(2) (Poisson ratio >= 0), and by default every layer
stays below the half-space Vs, so that Vs increases across the Moho and does not drop
just above it. Clipping an unbounded solution instead stalled on real data.

When a straight line in layer number already fits, the roughness is zero and chi2
stays below the target: the data require no more structure than that.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import lsq_linear

from compy.inversion.impedance import ImpedanceForward
from compy.numerics import GRAVITY, integer, positive

MOHO_JUMP_MINIMUM = 1.0  # m/s between the last layer and the half-space


def crawford_thin_layers(reference, top=50.0, growth=1.1, vp_vs_start=1.87):
    """Thin-layer model on a reference, with Vs starting from Vp / ``vp_vs_start``.

    Crawford (2004) uses a top layer of 0.05 km and a growth factor of 1.1, and
    Crawford (1994) a starting Vp/Vs of 1.87. Layers stop at the base of the reference,
    so a CRUST1.0 column gives the Moho as the last interface.

    Args:
        reference (array): Rows ``[thickness, rho, vp, vs]``; the last is the half-space.
        top (float): Thickness of the top layer [m].
        growth (float): Ratio between successive layers.
        vp_vs_start (float): Vp/Vs of the starting finite layers.

    Returns:
        tuple: Model with the same columns, and the layer edges below the seafloor [m].
    """
    reference = np.asarray(reference, float)
    positive(top, "top")
    if growth < 1 or len(reference) < 2:
        raise ValueError("need growth >= 1 and at least one finite layer")
    base = reference[:-1, 0].sum()
    edges, t = [0.0], float(top)
    while edges[-1] + t < base:
        edges.append(edges[-1] + t)
        t *= growth
    if len(edges) > 1 and base - edges[-1] < 0.5 * (edges[-1] - edges[-2]):
        edges.pop()  # no sliver above the last interface
    edges = np.array([*edges, base])
    tops = np.r_[0.0, np.cumsum(reference[:-1, 0])]
    rows = reference[
        np.searchsorted(tops, 0.5 * (edges[:-1] + edges[1:]), side="right") - 1
    ]
    model = np.vstack(
        [
            np.column_stack(
                [np.diff(edges), rows[:, 1], rows[:, 2], rows[:, 2] / vp_vs_start]
            ),
            reference[-1],
        ]
    )
    return model, edges


@dataclass
class MinimumStructureResult:
    """Outcome of ``invert_minimum_structure``.

    Attributes:
        model (ndarray): Rows ``[thickness, rho, vp, vs]`` with the inverted Vs.
        prediction (ndarray): Compliance of the model [1/Pa].
        chi2 (float): Sum of squared normalized residuals.
        target_chi2 (float): The misfit asked for.
        reached_target (bool): False means no model of this form fits the data.
        roughness (float): Sum of squared second differences of Vs [(m/s)**2].
        mu (float): Lagrange multiplier of the last accepted step.
        iterations (int): Accepted steps.
        history (list): chi2, roughness, mu and step length per accepted step.
    """

    model: np.ndarray
    prediction: np.ndarray
    chi2: float
    target_chi2: float
    reached_target: bool
    roughness: float
    mu: float
    iterations: int
    history: list


def invert_minimum_structure(
    eta,
    freqs,
    depth,
    model,
    sigma,
    *,
    start=None,
    target_chi2=None,
    vs_min=50.0,
    moho_jump=True,
    max_iter=60,
    gravity=GRAVITY,
):
    """Smoothest Vs profile that fits compliance, after Crawford (2004).

    Args:
        eta (array): Measured normalized compliance [1/Pa].
        freqs (array): Frequencies [Hz].
        depth (float): Water depth [m].
        model (array): Rows ``[thickness, rho, vp, vs]`` (see ``crawford_thin_layers``);
            only the Vs of the finite layers changes.
        sigma (array): Standard error per frequency [1/Pa].
        start (array or None): Starting Vs of the finite layers; default the model's.
        target_chi2 (float or None): Target misfit; default the number of data.
        vs_min (float): Lower bound on Vs [m/s].
        moho_jump (bool): Keep the Vs of every finite layer below the half-space Vs.
        max_iter (int): Maximum iterations.
        gravity (float): Gravitational acceleration [m/s**2].

    Returns:
        MinimumStructureResult: The model and its fit.
    """
    eta, freqs, sigma = (np.asarray(v, float) for v in (eta, freqs, sigma))
    model = np.array(model, float)
    if eta.ndim != 1 or eta.shape != freqs.shape or eta.shape != sigma.shape:
        raise ValueError("eta, freqs and sigma must be vectors of equal length")
    if np.any(~np.isfinite(eta)) or np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("eta must be finite and sigma finite and positive")
    if model.ndim != 2 or model.shape[1] != 4 or len(model) < 4:
        raise ValueError(
            "model needs columns [thickness, rho, vp, vs] and >= 3 finite layers"
        )
    n = len(model) - 1
    target = len(eta) if target_chi2 is None else positive(target_chi2, "target_chi2")
    lower = np.full(n, positive(vs_min, "vs_min"))
    upper = model[:-1, 2] / np.sqrt(2)
    if moho_jump:
        upper = np.minimum(upper, model[-1, 3] - MOHO_JUMP_MINIMUM)
    if np.any(upper <= lower):
        raise ValueError("no Vs satisfies the bounds in at least one layer")
    forward = ImpedanceForward(depth, freqs, gravity)
    second = np.diff(np.eye(n), n=2, axis=0)
    zeros = np.zeros(n - 2)

    def predict(vs):
        thickness, rho, vp = model[:, 0].copy(), model[:, 1].copy(), model[:, 2].copy()
        return forward(thickness, rho, vp, np.r_[vs, model[-1, 3]])

    def misfit(prediction):
        return float(np.sum(((eta - prediction) / sigma) ** 2))

    def roughness(vs):
        return float(np.sum((second @ vs) ** 2))

    vs = model[:-1, 3].copy() if start is None else np.array(start, float)
    if vs.shape != (n,):
        raise ValueError("start needs one Vs per finite layer")
    vs = np.clip(vs, lower, upper)
    prediction = predict(vs)
    chi2 = misfit(prediction)
    history = [dict(chi2=chi2, roughness=roughness(vs))]
    mu_now = float("nan")
    for _ in range(integer(max_iter, "max_iter")):
        jacobian = np.empty((len(eta), n))
        for j in range(n):
            trial = vs.copy()
            trial[j] += 1.0
            jacobian[:, j] = predict(trial) - prediction
        weighted = jacobian / sigma[:, None]
        rhs = np.r_[(eta - prediction + jacobian @ vs) / sigma, zeros]
        scale = np.sum(weighted**2) / np.sum(second**2)
        trials = []
        for mu in scale * np.logspace(-8, 6, 29):
            design = np.vstack([weighted, np.sqrt(mu) * second])
            candidate = lsq_linear(design, rhs, bounds=(lower, upper), max_iter=200).x
            p = predict(candidate)
            if np.all(np.isfinite(p)):
                trials.append((mu, misfit(p), candidate, p))
        if not trials:
            break
        fitting = [t for t in trials if t[1] <= target]
        mu, new_chi2, candidate, p = (
            min(fitting, key=lambda t: roughness(t[2]))
            if fitting
            else min(trials, key=lambda t: t[1])
        )
        reached, old_roughness = chi2 <= target, roughness(vs)

        def admissible(c, v):
            return (
                (c <= target and roughness(v) < old_roughness) if reached else c < chi2
            )

        step = 1.0
        while not admissible(new_chi2, candidate) and step > 1 / 64:
            step /= 2
            candidate = vs + step * (candidate - vs)
            p = predict(candidate)
            new_chi2 = misfit(p)
        if not admissible(new_chi2, candidate):
            break
        vs, prediction, chi2, mu_now = candidate, p, new_chi2, float(mu)
        history.append(dict(chi2=chi2, roughness=roughness(vs), mu=mu_now, step=step))
        if reached and old_roughness - roughness(vs) <= 1e-3 * old_roughness:
            break
    out = model.copy()
    out[:-1, 3] = vs
    return MinimumStructureResult(
        model=out,
        prediction=prediction,
        chi2=chi2,
        target_chi2=float(target),
        reached_target=chi2 <= target,
        roughness=roughness(vs),
        mu=mu_now,
        iterations=len(history) - 1,
        history=history,
    )
