"""Monte Carlo compliance inversion with parameterizations published by other groups.

A second inversion, separate from ``compy.inversion.dram``: it shares neither that
module's model nor its sampler, and uses the forward model of
``compy.inversion.impedance``. Use it to test whether a result depends on the
parameterization.

Two published setups:

``ZhaWebb2016Model``
    Zha & Webb (2016, JGR), as described in Zha's PhD thesis (2015, section 4.3): a
    sediment layer with a fixed Vs and Vp profile and variable thickness, crustal
    layers from 100 m doubling down to a fixed Moho, a mantle half-space, Vp/Vs of
    1.85 in crust and mantle, and a density profile fixed by the starting Vp.
``Zheng2024Model``
    Zheng et al. (2024, JGR, doi:10.1029/2024JB029862, section 3.3, compliance only):
    a linear-gradient sediment layer, a monotonic crust of 4 cubic B-splines to 5 km,
    Vp from Castagna et al. (1985) and Brocher (2005), and constant densities.

Both reject any model without a positive Vs jump at the Moho, the hard prior of
Bayesian Monte Carlo crustal inversions (Shen et al. 2013), with mantle Vs up to
4.9 km/s.

Differences from the papers are in each class docstring.

The sampler follows the papers' Metropolis rule on exp(-(chi2 + penalty) / 2) with 10
chains started from prior draws. The proposal covariance adapts during burn-in only
(Haario et al. 2001). The misfit filter of Shen et al. (2013) is then applied, as in
both papers.
"""

import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from compy.diagnostics import rank_rhat_ess
from compy.inversion.impedance import ImpedanceForward
from compy.numerics import GRAVITY, integer, positive

MANTLE_VS_BOUNDS = (3800.0, 4900.0)  # m/s; 4.9 km/s is the limit of Shen et al. (2013)


def hamilton_sediment_vs(z):
    """Vs [m/s] of silt-clay sediment at depth z [m] (Hamilton 1979)."""
    z = np.asarray(z, float)
    return np.where(
        z < 36, 116 + 4.65 * z, np.where(z < 120, 237 + 1.28 * z, 322 + 0.58 * z)
    )


def hamilton_sediment_vp(z):
    """Vp [m/s] of sediment at depth z [m] (Hamilton 1979)."""
    d = np.asarray(z, float) / 1000.0
    return 1000.0 * (1.511 + 1.304 * d - 0.741 * d**2 + 0.257 * d**3)


def brocher_crust_vp(vs):
    """Vp [m/s] from Vs [m/s] (Brocher 2005 regression)."""
    s = np.asarray(vs, float) / 1000.0
    return 1000.0 * (
        0.9409 + 2.0947 * s - 0.8206 * s**2 + 0.2683 * s**3 - 0.0251 * s**4
    )


def castagna_sediment_vp(vs):
    """Vp [m/s] from Vs [m/s] (Castagna et al. 1985: Vp = 1.16 Vs + 1.36 km/s)."""
    return 1160.0 * (np.asarray(vs, float) / 1000.0) + 1360.0


def positive_moho_jump(vs_last_crust, vs_mantle):
    """True if Vs increases across the Moho."""
    return vs_mantle > vs_last_crust


class ZhaWebb2016Model:
    """Parameterization of Zha & Webb (2016).

    Parameters: sediment thickness, one Vs per crustal layer, mantle Vs.

    Differences from the paper:

    * sediment Vs is Hamilton (1979), the alternative Zha tests, not Ruan et al. (2014);
    * sediment thickness is continuous with a partial last 50 m sublayer, not on a
      50 m grid, because the grid makes chains stick in neighbouring cells;
    * mantle Vs is 3800-4900 m/s with a positive Moho jump, not 3800-4200 m/s, which
      suits a hot back-arc;
    * the paper does not state the crustal Vs bounds or the roughness weight, so
      100-4500 m/s and 1 per (km/s)**2 are defaults here;
    * the Moho and starting Vp come from the caller; the paper used local seismic models.

    Args:
        moho_depth (float): Moho depth below the seafloor [m].
        vp_tops (array): Tops of the starting Vp layers [m].
        vp_values (array): Starting Vp of those layers [m/s].
        vp_mantle (float): Starting mantle Vp [m/s].
        roughness_weight (float): Weight on squared second differences of crustal Vs
            [per (km/s)**2].
    """

    label = "zha_webb_2016"
    sediment_max, vp_vs = 2000.0, 1.85
    vs_crust_bounds = (100.0, 4500.0)

    def __init__(self, moho_depth, vp_tops, vp_values, vp_mantle, roughness_weight=1.0):
        self.moho_depth = positive(moho_depth, "moho_depth")
        self.vp_tops, self.vp_values = (
            np.asarray(vp_tops, float),
            np.asarray(vp_values, float),
        )
        self.rho_mantle = 1850.0 + 0.165 * positive(vp_mantle, "vp_mantle")
        self.roughness_weight = float(roughness_weight)
        fixed = []
        limit = self.moho_depth - self.sediment_max - 500.0
        while sum(fixed) + 100.0 * 2 ** len(fixed) <= limit:
            fixed.append(100.0 * 2 ** len(fixed))
        if not fixed:
            raise ValueError("Moho too shallow for this parameterization")
        self.fixed = np.array(fixed)
        self.n_crust = len(fixed) + 1
        self.names = [
            "sediment_thickness",
            *[f"vs_crust_{i + 1}" for i in range(self.n_crust)],
            "vs_mantle",
        ]
        self.lower = np.r_[
            0.0, np.full(self.n_crust, self.vs_crust_bounds[0]), MANTLE_VS_BOUNDS[0]
        ]
        self.upper = np.r_[
            self.sediment_max,
            np.full(self.n_crust, self.vs_crust_bounds[1]),
            MANTLE_VS_BOUNDS[1],
        ]
        self.index_last_crust = self.n_crust

    def initial(self, rng):
        """A random draw from the bounds; the sampler redraws until it is valid."""
        return rng.uniform(self.lower, self.upper)

    def is_valid(self, theta):
        """True if ``theta`` is inside the bounds and has a positive Moho jump."""
        return bool(
            np.all(theta >= self.lower)
            and np.all(theta <= self.upper)
            and positive_moho_jump(theta[self.n_crust], theta[-1])
        )

    def layers(self, theta):
        """Arrays ``thickness, rho, vp, vs``; the last row is the half-space."""
        h = float(theta[0])
        n = int(math.ceil(h / 50.0))
        sub = np.diff(np.r_[np.arange(n) * 50.0, h]) if n else np.empty(0)
        z = np.arange(n) * 50.0 + 0.5 * sub
        edges = np.r_[h, h + np.cumsum(self.fixed), self.moho_depth]
        mid = 0.5 * (edges[:-1] + edges[1:])
        vs = theta[1 : 1 + self.n_crust]
        start_vp = self.vp_values[np.searchsorted(self.vp_tops, mid, side="right") - 1]
        return (
            np.r_[sub, np.diff(edges), 0.0],
            np.r_[np.full(n, 2000.0), 1850.0 + 0.165 * start_vp, self.rho_mantle],
            np.r_[hamilton_sediment_vp(z), self.vp_vs * vs, self.vp_vs * theta[-1]],
            np.r_[hamilton_sediment_vs(z), vs, theta[-1]],
        )

    def penalty(self, theta):
        """Roughness term alpha R of Zha's likelihood."""
        return self.roughness_weight * float(
            np.sum(np.diff(theta[1 : 1 + self.n_crust] / 1000.0, n=2) ** 2)
        )

    def retained(self, chi2, n_data):
        """Zha (2015): keep chi2 / N below max(2 min, min + 0.25)."""
        norm = np.asarray(chi2, float) / n_data
        return norm < max(2 * norm.min(), norm.min() + 0.25)


class Zheng2024Model:
    """Parameterization of Zheng et al. (2024), compliance-only case.

    Parameters: sediment thickness, Vs at the sediment top and base, 4 B-spline
    coefficients of crustal Vs, mantle Vs.

    Differences from the paper:

    * the paper's model stops at 5 km; here Vs stays at the last coefficient from 5 km
      to the given Moho, above a mantle half-space with a positive jump, of density
      3000 kg/m**3 (the paper's mantle density);
    * sediment thickness is uniform in 0-2000 m, because there is no GlobSed or
      reflection prior here;
    * the sampler is Metropolis, not No-U-Turn.

    Args:
        moho_depth (float): Moho depth below the seafloor [m].
    """

    label = "zheng_2024"
    spline_bottom, sediment_max, n_sediment, n_crust = 5000.0, 2000.0, 8, 24
    vs_max = 4500.0
    vs_sediment_max = (4000.0 - 1360.0) / 1.16  # Vp <= 4 km/s in the sediment
    vs_sediment_min = 1360.0 / (15.0 - 1.16)  # Vp/Vs <= 15 in the sediment

    def __init__(self, moho_depth):
        self.moho_depth = positive(moho_depth, "moho_depth")
        if self.moho_depth < self.spline_bottom:
            raise ValueError("this parameterization needs a Moho below 5 km")
        self.names = [
            "sediment_thickness",
            "vs_sediment_top",
            "vs_sediment_base",
            "spline_1",
            "spline_2",
            "spline_3",
            "spline_4",
            "vs_mantle",
        ]
        self.lower = np.r_[
            0.0,
            self.vs_sediment_min,
            self.vs_sediment_min,
            0.0,
            0.0,
            0.0,
            0.0,
            MANTLE_VS_BOUNDS[0],
        ]
        self.upper = np.array(
            [
                self.sediment_max,
                self.vs_sediment_max,
                self.vs_sediment_max,
                *[self.vs_max] * 4,
                MANTLE_VS_BOUNDS[1],
            ]
        )
        self.index_last_crust = 6

    def initial(self, rng):
        """A random draw from the bounds, ordered so that Vs increases with depth."""
        for _ in range(100000):
            t = rng.uniform(self.lower, self.upper)
            t[1:3] = np.sort(t[1:3])
            t[3:7] = np.sort(rng.uniform(t[2], self.vs_max, 4))
            t[7] = rng.uniform(max(t[6], MANTLE_VS_BOUNDS[0]), MANTLE_VS_BOUNDS[1])
            if self.is_valid(t):
                return t
        raise RuntimeError("no valid starting model")

    def is_valid(self, theta):
        """Bounds, Vs increasing from the sediment top to the mantle, Vp/Vs >= 1.7 in the sediment."""
        if np.any(theta < self.lower) or np.any(theta > self.upper):
            return False
        if np.any(np.diff(theta[1:7]) < 0) or not positive_moho_jump(
            theta[6], theta[7]
        ):
            return False
        return bool(
            castagna_sediment_vp(theta[1]) / theta[1] >= 1.70
            and castagna_sediment_vp(theta[2]) / theta[2] >= 1.70
            and theta[6] > 0
        )

    def layers(self, theta):
        """Arrays ``thickness, rho, vp, vs``; the last row is the half-space."""
        h = float(theta[0])
        edges = np.linspace(h, self.spline_bottom, self.n_crust + 1)
        x = (0.5 * (edges[:-1] + edges[1:]) - h) / (self.spline_bottom - h)
        vs = (
            theta[3] * (1 - x) ** 3
            + 3 * theta[4] * x * (1 - x) ** 2
            + 3 * theta[5] * x**2 * (1 - x)
            + theta[6] * x**3
        )
        vs = np.maximum(vs, 1.0)
        vp = np.maximum(brocher_crust_vp(vs), 1.5 * vs)
        thick, rho, p, s = np.diff(edges), np.full(self.n_crust, 2800.0), vp, vs
        last = float(theta[6])
        if self.moho_depth > self.spline_bottom:
            thick = np.r_[thick, self.moho_depth - self.spline_bottom]
            rho = np.r_[rho, 2800.0]
            p = np.r_[p, max(float(brocher_crust_vp(last)), 1.5 * last)]
            s = np.r_[s, last]
        if h > 1.0:
            e = np.linspace(0.0, h, self.n_sediment + 1)
            vs_sed = theta[1] + (0.5 * (e[:-1] + e[1:]) / h) * (theta[2] - theta[1])
            thick = np.r_[np.diff(e), thick]
            rho = np.r_[np.full(self.n_sediment, 2000.0), rho]
            p = np.r_[castagna_sediment_vp(vs_sed), p]
            s = np.r_[vs_sed, s]
        vm = float(theta[7])
        return (
            np.r_[thick, 0.0],
            np.r_[rho, 3000.0],
            np.r_[p, max(float(brocher_crust_vp(vm)), 1.5 * vm)],
            np.r_[s, vm],
        )

    def penalty(self, theta):
        """This parameterization has no roughness term."""
        return 0.0

    def retained(self, chi2, n_data):
        """Zheng et al. (2024): keep rms below 2 min, or below min + 0.5 if min < 0.5."""
        rms = np.sqrt(np.asarray(chi2, float) / n_data)
        low = rms.min()
        return rms < (2 * low if low >= 0.5 else low + 0.5)


def vs_at_depths(model, theta, z):
    """Vs [m/s] of one parameter vector at depths ``z`` [m] below the seafloor."""
    thickness, _, _, vs = model.layers(theta)
    tops = np.r_[0.0, np.cumsum(thickness[:-1])]
    return vs[np.searchsorted(tops, np.asarray(z, float), side="right") - 1]


def monthly_standard_deviation(curves):
    """Data error from monthly compliance curves, shape (months, frequencies).

    Zheng et al. (2024) use the standard deviation of their daily compliance. This
    is the same statistic on monthly stacks, so it is larger than a standard error.
    """
    curves = np.atleast_2d(np.asarray(curves, float))
    if curves.shape[0] < 3:
        raise ValueError("need at least 3 monthly curves")
    return curves.std(axis=0, ddof=1)


def _run_chain(job):
    """One adaptive-Metropolis chain; the proposal is frozen after burn-in."""
    model, depth, freqs, eta, sigma, gravity, seed, n_iter, n_burn, thin, prior_only = (
        job
    )
    rng = np.random.default_rng(seed)
    forward = ImpedanceForward(depth, freqs, gravity)

    def evaluate(theta):
        if not model.is_valid(theta):
            return None
        penalty = model.penalty(theta)
        if prior_only:
            return 0.0, penalty, np.zeros(len(eta))
        prediction = forward(*model.layers(theta))
        if not np.all(np.isfinite(prediction)):
            return None
        return float(np.sum(((eta - prediction) / sigma) ** 2)), penalty, prediction

    x, state = model.initial(rng), None
    for _ in range(100000):
        state = evaluate(x)
        if state is not None:
            break
        x = model.initial(rng)
    if state is None:
        raise RuntimeError("no valid starting model")
    d, width = len(x), model.upper - model.lower
    chol = np.linalg.cholesky(np.diag((0.05 * width) ** 2))
    history = np.empty((n_burn, d))
    log_scale, recent = 0.0, [0, 0]
    keep, misfit, penalties, predictions, accepted = [], [], [], [], 0
    for it in range(n_iter):
        y = x + chol @ rng.standard_normal(d)
        new = evaluate(y)
        moved = False
        if new is not None:
            delta = (state[0] + state[1]) - (new[0] + new[1])
            moved = delta >= 0 or rng.random() < math.exp(delta / 2)
        if moved:
            x, state = y, new
        if it < n_burn:
            history[it] = x
            recent[0] += moved
            recent[1] += 1
            if it >= 500 and it % 100 == 0:
                log_scale += (
                    (recent[0] / recent[1] - 0.234) * 2 * (1 + (it - 500) / 100) ** -0.6
                )
                recent = [0, 0]
                cov = np.cov(history[it // 2 : it + 1], rowvar=False) + np.diag(
                    (1e-4 * width) ** 2
                )
                chol = np.linalg.cholesky(math.exp(2 * log_scale) * 2.38**2 / d * cov)
        else:
            accepted += moved
            if (it - n_burn) % thin == 0:
                keep.append(x.copy())
                misfit.append(state[0])
                penalties.append(state[1])
                predictions.append(state[2])
    return dict(
        theta=np.array(keep),
        chi2=np.array(misfit),
        penalty=np.array(penalties),
        prediction=np.array(predictions),
        acceptance=accepted / (n_iter - n_burn),
    )


@dataclass
class LiteratureMCMCResult:
    """Retained draws of ``sample_literature_model``; arrays are shaped (chains, draws, ...)."""

    model: object
    depth: float
    freqs: np.ndarray
    eta: np.ndarray
    sigma: np.ndarray
    prior_only: bool
    theta: np.ndarray
    chi2: np.ndarray
    penalty: np.ndarray
    prediction: np.ndarray
    acceptance: np.ndarray

    @property
    def retained_mask(self):
        """Pooled draws that pass the misfit filter (all of them in a prior-only run)."""
        pooled = self.chi2.ravel()
        if self.prior_only:
            return np.ones(len(pooled), bool)
        return self.model.retained(pooled, len(self.eta))

    @property
    def retained_theta(self):
        """Pooled parameter draws that pass the misfit filter."""
        return self.theta.reshape(-1, self.theta.shape[-1])[self.retained_mask]

    def diagnostics(self):
        """Rank-normalized split R-hat and bulk ESS per parameter, over all draws."""
        out = {}
        for i, name in enumerate(self.model.names):
            r, ess = rank_rhat_ess(self.theta[:, :, i])
            out[name] = dict(rhat=float(r), ess_bulk=float(ess))
        return out

    def vs_draws(self, z, max_draws=8000, seed=0):
        """Vs [m/s] of retained draws at depths ``z``, shape (draws, len(z))."""
        theta = self.retained_theta
        if len(theta) > max_draws:
            pick = np.random.default_rng(seed).choice(
                len(theta), max_draws, replace=False
            )
            theta = theta[pick]
        return np.array([vs_at_depths(self.model, t, z) for t in theta])

    def vs_quantiles(self, z, q=(0.05, 0.16, 0.5, 0.84, 0.95), max_draws=8000):
        """Pointwise Vs quantiles at depths ``z``, shape (len(q), len(z))."""
        return np.quantile(self.vs_draws(z, max_draws), q, axis=0)

    def summary(self):
        """JSON-ready summary: parameter quantiles, diagnostics, fit and the Moho jump."""
        kept = self.retained_theta
        diag = self.diagnostics()
        q = np.quantile(kept, [0.05, 0.16, 0.5, 0.84, 0.95], axis=0)
        chi2 = self.chi2.ravel()
        jump = kept[:, -1] - kept[:, self.model.index_last_crust]
        return dict(
            setup=self.model.label,
            prior_only=self.prior_only,
            n_data=len(self.eta),
            depth=self.depth,
            moho_depth=self.model.moho_depth,
            names=list(self.model.names),
            parameters={
                name: dict(
                    zip(("p05", "p16", "median", "p84", "p95"), map(float, q[:, i])),
                    **diag[name],
                )
                for i, name in enumerate(self.model.names)
            },
            chi2_min=float(chi2.min()),
            chi2_median_retained=float(np.median(chi2[self.retained_mask])),
            retained=int(self.retained_mask.sum()),
            total=len(chi2),
            max_rhat=max(v["rhat"] for v in diag.values()),
            min_ess=min(v["ess_bulk"] for v in diag.values()),
            acceptance=self.acceptance.tolist(),
            moho_jump=dict(median=float(np.median(jump)), minimum=float(jump.min())),
        )

    def save(self, path):
        """Write ``<path>.npz`` (retained draws and data) and ``<path>.json`` (summary)."""
        path = Path(path)
        pooled = self.prediction.reshape(-1, len(self.eta))[self.retained_mask]
        np.savez_compressed(
            path.with_suffix(".npz"),
            theta=self.retained_theta,
            chi2=self.chi2.ravel()[self.retained_mask],
            names=np.array(self.model.names),
            freqs=self.freqs,
            eta=self.eta,
            sigma=self.sigma,
            depth=self.depth,
            prediction_median=np.median(pooled, axis=0),
        )
        path.with_suffix(".json").write_text(json.dumps(self.summary(), indent=1))


def sample_literature_model(
    model,
    depth,
    freqs,
    eta,
    sigma,
    *,
    n_chains=10,
    n_iter=400000,
    n_burn=200000,
    thin=40,
    seed=11,
    processes=None,
    prior_only=False,
    gravity=GRAVITY,
):
    """Run the Metropolis chains of a published parameterization.

    The defaults are far above the 10 x 10,000 of Zha and 10 x 20,000 of Zheng,
    because synthetic tests had not converged at 50,000 iterations. Check
    ``result.summary()["max_rhat"]`` before using a result.

    Args:
        model (ZhaWebb2016Model or Zheng2024Model): The parameterization.
        depth (float): Water depth [m].
        freqs (array): Frequencies [Hz].
        eta (array): Measured normalized compliance [1/Pa].
        sigma (array): Standard error per frequency [1/Pa].
        n_chains, n_iter, n_burn, thin (int): Chain settings.
        seed (int): Seed of all chains; a result repeats for a given seed.
        processes (int or None): Worker processes; 1 runs serially. Workers are
            spawned, so call from a script guarded by ``if __name__ == "__main__":``.
        prior_only (bool): Sample the prior without the data.
        gravity (float): Gravitational acceleration [m/s**2].

    Returns:
        LiteratureMCMCResult: Draws, fit and diagnostics.
    """
    freqs, eta, sigma = (np.asarray(v, float) for v in (freqs, eta, sigma))
    if freqs.ndim != 1 or eta.shape != freqs.shape or sigma.shape != freqs.shape:
        raise ValueError("freqs, eta and sigma must be vectors of equal length")
    if np.any(~np.isfinite(eta)) or np.any(eta <= 0):
        raise ValueError("eta must be finite and positive")
    if np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("sigma must be finite and positive")
    n_chains, n_iter, n_burn, thin = (
        integer(n_chains, "n_chains", 2),
        integer(n_iter, "n_iter"),
        integer(n_burn, "n_burn"),
        integer(thin, "thin"),
    )
    if n_iter <= n_burn:
        raise ValueError("n_iter must exceed n_burn")
    seeds = np.random.SeedSequence(seed).generate_state(n_chains)
    jobs = [
        (
            model,
            depth,
            freqs,
            eta,
            sigma,
            gravity,
            int(s),
            n_iter,
            n_burn,
            thin,
            prior_only,
        )
        for s in seeds
    ]
    if processes == 1:
        chains = [_run_chain(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=processes) as pool:
            chains = list(pool.map(_run_chain, jobs))
    stack = {
        key: np.stack([c[key] for c in chains])
        for key in ("theta", "chi2", "penalty", "prediction")
    }
    return LiteratureMCMCResult(
        model=model,
        depth=float(depth),
        freqs=freqs,
        eta=eta,
        sigma=sigma,
        prior_only=prior_only,
        acceptance=np.array([c["acceptance"] for c in chains]),
        **stack,
    )
