"""Deep-crust synthetic compliance inversion following the compliance literature.

Usage (repository root; about 25 minutes on two cores):
    python _Example/deep_crust_inversion.py --output deep-crust
    python _Example/plot_deep_posterior.py deep-crust

Set-up, with the paper behind each choice:
* Water depth 4200 m, 16 frequencies 3.5-18.5 mHz, below the infragravity
  cut-off sqrt(g / 2 pi H) = 19.3 mHz (Crawford et al. 1991).
* 3 % Gaussian errors, a good-coherence station (Crawford 2004: 1-15 %).
* True crust: sediment, layers 2A/2B, a low-velocity zone (the target of
  Crawford et al. 1999), gabbro, mantle half-space.
* Inversion on a different, minimum-structure grid (Crawford et al. 1991):
  14 layers whose thickness grows by 1.25 to 13 km, Vs only, Vp and density
  fixed from a refraction-type reference as in Crawford's studies and
  Aminian et al. (2025). Using a grid unlike the true layering avoids the
  "inverse crime".
* Uniform Vs prior 200-4800 m/s with Vs < Vp, plus a Gaussian prior on Vs
  jumps whose scale is learned from the data (hierarchical Bayes, Malinverno
  & Briggs 2004) instead of hand-picked or L-curve chosen.
* Four chains from deliberately different starts, adaptive step sizes during
  burn-in, rank-normalized split R-hat and bulk ESS (Vehtari et al. 2021).
* Cross-checks: Occam's inversion (Constable et al. 1987), numerical Vs
  kernels (Doran & Laske 2020) and a prior-only run with the data switched
  off, so the posterior can be compared with what the prior alone allows.
"""

from __future__ import annotations

import argparse
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import inv_compy as inv
from compy_diagnostics import rank_rhat_ess

H = 4200.0
F = np.geomspace(0.0035, 0.0185, 16)
NOISE, NOISE_SEED = 0.03, 2026
TRUTH = np.array(
    [
        [300.0, 1900.0, 2000.0, 800.0],
        [900.0, 2500.0, 4200.0, 2100.0],
        [1300.0, 2800.0, 6200.0, 3500.0],
        [1500.0, 2800.0, 6000.0, 2700.0],
        [2500.0, 3000.0, 7000.0, 3900.0],
        [0.0, 3300.0, 8000.0, 4400.0],
    ]
)
N_LAYER, FIRST, GROWTH = 14, 150.0, 1.25
THICK = FIRST * GROWTH ** np.arange(N_LAYER)
EDGES = np.concatenate([[0.0], np.cumsum(THICK)])
MID = 0.5 * (EDGES[:-1] + EDGES[1:])
VS_BOUNDS = (200.0, 4800.0)
SIGMA_RANGE = (100.0, 5000.0)  # allowed Vs-jump scale, m/s
STEP = inv.Roughness(np.r_[np.zeros(100), np.ones(100)], 2)  # roughness of a unit step
ALPHA_BOUNDS = (1 / (STEP * SIGMA_RANGE[1] ** 2), 1 / (STEP * SIGMA_RANGE[0] ** 2))
ALPHA0 = 1 / (STEP * 600.0**2)


def sample_truth(column, z):
    return TRUTH[np.searchsorted(np.cumsum(TRUTH[:-1, 0]), z, side="right"), column]


def grid_model(vs):
    rows = np.column_stack([THICK, sample_truth(1, MID), sample_truth(2, MID), vs])
    return np.vstack([rows, TRUTH[-1]])


VP = sample_truth(2, MID)
START_VS = np.minimum(1200 + 2900 * (1 - np.exp(-MID / 3000)), 0.9 * VP)
STARTS = {
    "slow": np.clip(400 + 0.25 * VP, 250, 0.9 * VP),
    "fast": np.clip(0.58 * VP, 250, 4700),
    "constant": np.minimum(2500.0, 0.9 * VP),
    "generic": START_VS,
}
CLEAN = inv.calc_norm_compliance(H, F, TRUTH)
UNC = NOISE * np.abs(CLEAN)
OBS = CLEAN + np.random.default_rng(NOISE_SEED).normal(0, UNC)


def run_chain(args):
    k, iterations, adapt, prior_only = args
    name = list(STARTS)[k % len(STARTS)]
    info = {}
    out = inv.invert_compliance_beta(
        OBS,
        F,
        H,
        starting_model=grid_model(START_VS),
        s=UNC * (1e9 if prior_only else 1.0),
        iteration=iterations,
        sigma_v=100,
        alpha=ALPHA0,
        seed=(900 if prior_only else 200) + k,
        return_profiles=False,
        adapt_until=adapt,
        initial_state=grid_model(STARTS[name]),
        vs_bounds=VS_BOUNDS,
        vary_thickness=False,
        alpha_bounds=ALPHA_BOUNDS,
        diagnostics=info,
    )
    return (
        out[0][:, 3, :].astype(np.float32),
        out[4].astype(np.float32),
        out[3][0],
        out[-1],
        info["alpha"].astype(np.float32),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=Path("deep-crust"))
    parser.add_argument("--iterations", type=int, default=160_000)
    parser.add_argument("--adapt", type=int, default=20_000)
    parser.add_argument("--burnin", type=int, default=40_000)
    parser.add_argument("--prior-iterations", type=int, default=100_000)
    args = parser.parse_args()
    if not args.adapt <= args.burnin < args.iterations - 8:
        raise ValueError("need adapt <= burnin < iterations - 8")
    args.output.mkdir(parents=True, exist_ok=False)
    t0 = time.time()
    with Pool(2) as pool:
        res = pool.map(
            run_chain, [(k, args.iterations, args.adapt, False) for k in range(4)]
        )
        prior = pool.map(
            run_chain, [(k, args.prior_iterations, args.adapt, True) for k in range(2)]
        )
    vs = np.stack([r[0] for r in res])
    misfit = np.stack([r[2] for r in res])
    alphas = np.stack([r[4] for r in res])
    B = args.burnin
    occam = inv.invert_occam(OBS, F, H, grid_model(START_VS), UNC)
    z_k, kernels = inv.compliance_kernels(H, F, TRUTH, dz=100, max_depth=20000)
    layers = {}
    for j in range(N_LAYER):
        post = vs[:, j, B:].astype(float)
        r, ess = rank_rhat_ess(post)
        q = np.quantile(post, [0.05, 0.5, 0.95])
        layers[f"layer_{j + 1:02d}"] = dict(
            top_m=float(EDGES[j]),
            bottom_m=float(EDGES[j + 1]),
            rhat=r,
            ess_bulk=ess,
            p05=float(q[0]),
            median=float(q[1]),
            p95=float(q[2]),
            truth_mid=float(sample_truth(3, MID[j])),
            occam=float(occam["model"][j, 3]),
        )
    sigma = np.sqrt(1 / (STEP * alphas[:, B:].astype(float)))
    rs, ess_s = rank_rhat_ess(np.log(sigma))
    summary = dict(
        settings=dict(
            iterations=args.iterations,
            adapt_until=args.adapt,
            burnin=B,
            chains=4,
            starts=list(STARTS),
            water_depth_m=H,
            freq_hz=F.tolist(),
            noise=NOISE,
            vs_bounds=VS_BOUNDS,
            sigma_range=SIGMA_RANGE,
            n_layer=N_LAYER,
        ),
        layers=layers,
        sigma_jump=dict(
            p05=float(np.quantile(sigma, 0.05)),
            median=float(np.median(sigma)),
            p95=float(np.quantile(sigma, 0.95)),
            rhat=rs,
            ess_bulk=ess_s,
        ),
        chi2_truth=float(np.sum(((OBS - CLEAN) / UNC) ** 2)),
        chi2_posterior_median=float(np.median(misfit[:, B:] ** 2)),
        occam=dict(chi2=occam["chi2"], iterations=occam["iterations"]),
        acceptance=[r[3] for r in res],
        max_rhat=max(v["rhat"] for v in layers.values()),
        min_ess=min(v["ess_bulk"] for v in layers.values()),
        runtime_s=time.time() - t0,
    )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=1))
    np.savez_compressed(
        args.output / "chains.npz",
        vs=vs,
        preds=np.stack([r[1] for r in res]),
        misfit=misfit,
        alpha=alphas,
        occam=occam["model"],
        obs=OBS,
        clean=CLEAN,
        unc=UNC,
        truth=TRUTH,
        edges=EDGES,
        freq=F,
        water_depth=H,
        kernel_z=z_k,
        kernels=kernels,
        prior_vs=np.stack([p[0][:, args.adapt :] for p in prior]),
    )
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "max_rhat",
                    "min_ess",
                    "sigma_jump",
                    "chi2_truth",
                    "chi2_posterior_median",
                )
            },
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
