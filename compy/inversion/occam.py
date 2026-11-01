"""Occam inversion for the smoothest Vs profile that fits the data (Constable et al. 1987)."""

import numpy as np

from compy.inversion.models import _validate_model
from compy.inversion.propagator import normalized_compliance
from compy.numerics import integer, positive, residual


def invert_occam(
    Data,
    f,
    depth_s,
    model,
    s,
    *,
    target_chi2=None,
    order=2,
    max_iter=40,
    vs_min=50.0,
):
    """Occam's inversion for a smooth Vs profile (Constable et al. 1987).

    The linearized minimum-structure method used for compliance by Crawford
    et al. (1991, 1998). Finite-layer Vs values are inverted; thickness, Vp,
    density and the half-space stay as given in ``model``. The problem is
    linearized in ln(Vs) and ln(compliance), since compliance scales roughly
    as 1/Vs**2. Each iteration builds a finite-difference Jacobian, scans the
    Lagrange multiplier and keeps the smoothest model whose chi-square is at
    most ``target_chi2`` (default: the number of data); while the target is
    out of reach it keeps the best-fitting one. A halving line search stops
    a step from degrading the fit. Roughness is the ``order``-th difference of
    ln(Vs) between layers. Chi-square is computed on the compliance itself,
    ``sum(((Data - prediction) / s) ** 2)``.

    Returns a dict with ``model``, ``chi2``, ``mu``, ``roughness``,
    ``iterations`` and a per-iteration ``history``.
    """
    model = _validate_model(model).copy()
    Data, f = np.asarray(Data, float), np.asarray(f, float)
    s = np.broadcast_to(np.asarray(s, float), Data.shape)
    residual(Data, Data, s)
    if np.any(Data <= 0):
        raise ValueError("Occam inversion works on positive normalized compliance")
    n = len(model) - 1
    order = integer(order, "order")
    if n <= order:
        raise ValueError("need more finite layers than the roughness order")
    target = len(Data) if target_chi2 is None else positive(target_chi2, "target_chi2")
    lo, hi = np.log(positive(vs_min, "vs_min")), np.log(0.98 * model[:-1, 2])
    D = np.diff(np.eye(n), n=order, axis=0)
    DtD = D.T @ D
    log_data, w = np.log(Data), Data / s  # 1 / relative uncertainty

    def forward(x):
        trial = model.copy()
        trial[:-1, 3] = np.exp(x)
        return normalized_compliance(depth_s, f, trial)

    def chi2(pred):
        return float(np.sum(((Data - pred) / s) ** 2))

    def rough(x):
        return float(np.sum((D @ x) ** 2))

    x = np.clip(np.log(model[:-1, 3]), lo, hi)
    pred = forward(x)
    c_now = chi2(pred)
    history = [dict(chi2=c_now, mu=np.nan, roughness=rough(x))]
    mu_now = np.nan
    for iteration in range(1, integer(max_iter, "max_iter") + 1):
        log_pred = np.log(pred)
        J = np.empty((len(Data), n))
        for j in range(n):
            step = x.copy()
            step[j] += 1e-3
            J[:, j] = (np.log(forward(step)) - log_pred) / 1e-3
        WJ = J * w[:, None]
        JtJ = WJ.T @ WJ
        rhs = WJ.T @ (w * (log_data - log_pred + J @ x))
        scale = np.trace(JtJ) / np.trace(DtD)
        trials = []
        for mu in scale * np.logspace(-8, 6, 57):
            cand = np.clip(np.linalg.solve(mu * DtD + JtJ, rhs), lo, hi)
            try:
                cand_pred = forward(cand)
            except ValueError:
                continue
            trials.append((mu, chi2(cand_pred), cand, cand_pred))
        feasible = [t for t in trials if t[1] <= target]
        # Smoothest model that fits (evaluated with the full forward model),
        # or the best fit while the target is still out of reach.
        mu, c2, cand, cand_pred = (
            min(feasible, key=lambda t: rough(t[2]))
            if feasible
            else min(trials, key=lambda t: t[1])
        )
        # Damped step. Before the target is reached a step must improve the
        # fit; afterwards it must stay feasible AND be smoother, so phase two
        # decreases roughness monotonically (Constable et al. 1987).
        reached = c_now <= target
        r_old = rough(x)

        def admissible(c, xx):
            return (c <= target and rough(xx) < r_old) if reached else c < c_now

        t = 1.0
        while not admissible(c2, cand) and t > 1 / 64:
            t /= 2
            trial_x = x + t * (cand - x)
            trial_pred = forward(trial_x)
            c2, cand, cand_pred = chi2(trial_pred), trial_x, trial_pred
        if not admissible(c2, cand):  # no admissible step: converged
            break
        r_new = rough(cand)
        x, pred, c_now, mu_now = cand, cand_pred, c2, mu
        history.append(dict(chi2=c_now, mu=float(mu), roughness=r_new, step=t))
        if reached and r_old - r_new <= 1e-3 * r_old:
            break
    final = model.copy()
    final[:-1, 3] = np.exp(x)
    return dict(
        model=final,
        chi2=c_now,
        mu=float(mu_now),
        roughness=rough(x),
        iterations=len(history) - 1,
        history=history,
    )
