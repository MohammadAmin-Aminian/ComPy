"""Converged 4-chain synthetic inversion with the tuned ComPy sampler.

Usage (repository root):
    mkdir converged && python _Example/converged_inversion.py converged
    python _Example/plot_posterior_heatmap.py converged

Same synthetic problem as _Example/reproducible_inversion.py (truth, 12
frequencies, 3% noise with seed 2026, biased starting model Vs1 = 1800 m/s that
defines the prior). Differences: per-layer adaptive step sizes during burn-in,
dispersed chain starts, longer chains, rank-normalized R-hat and bulk ESS.
"""
import sys
import json
import time
from multiprocessing import Pool

import numpy as np
from scipy.stats import norm, rankdata

import inv_compy as inv

ITER, ADAPT, BURN = 60_000, 10_000, 15_000
OUT = sys.argv[1]

TRUTH = np.array([[100.0, 2200.0, 3000.0, 1500.0],
                  [1000.0, 2800.0, 6000.0, 3500.0],
                  [0.0, 3300.0, 8000.0, 4500.0]])
START = TRUTH.copy()
START[0, 3] = 1800.0                      # biased prior-defining model, as in the repo
F = np.linspace(0.005, 0.02, 12)
DEPTH = 4000.0
CLEAN = inv.calc_norm_compliance(DEPTH, F, TRUTH)
UNC = 0.03 * np.abs(CLEAN)
OBS = CLEAN + np.random.default_rng(2026).normal(0, UNC)

# Over-dispersed starts spanning the prior support (Vs1 in [270, 2160],
# Vs2 in [3150, 3850]); total finite thickness 1100 m is conserved.
STARTS = [(400.0, 3200.0, 60.0), (900.0, 3800.0, 180.0),
          (1400.0, 3300.0, 130.0), (2100.0, 3700.0, 80.0)]


def run_chain(k):
    vs1, vs2, h1 = STARTS[k]
    state = START.copy()
    state[0, 3], state[1, 3], state[0, 0], state[1, 0] = vs1, vs2, h1, 1100.0 - h1
    info = {}
    out = inv.invert_compliance_beta(
        OBS, F, DEPTH, starting_model=START, s=UNC, iteration=ITER,
        sigma_v=35, sigma_h=8, alpha=0, seed=100 + k, return_profiles=False,
        adapt_until=ADAPT, target_acceptance=0.3, initial_state=state,
        diagnostics=info,
    )
    return out[0], out[4], out[3][0], out[-1], info


def rank_rhat_ess(x):
    """Rank-normalized split R-hat (max of bulk/folded) and bulk ESS (Vehtari 2021)."""
    m, n = x.shape
    half = n // 2
    s = np.concatenate([x[:, :half], x[:, half:2 * half]])

    def z(a):
        r = rankdata(a, method="average").reshape(a.shape)
        return norm.ppf((r - 0.375) / (a.size + 0.25))

    def rhat(a):
        nn = a.shape[1]
        w = a.var(axis=1, ddof=1).mean()
        b = nn * a.mean(axis=1).var(ddof=1)
        return np.sqrt(((nn - 1) / nn * w + b / nn) / w)

    zs = z(s)
    folded = z(np.abs(s - np.median(s)))
    r = max(rhat(zs), rhat(folded))
    # Bulk ESS with Geyer's initial monotone sequence.
    mm, nn = zs.shape
    c = zs - zs.mean(axis=1, keepdims=True)
    fft = np.fft.rfft(c, n=2 * nn)
    acov = np.fft.irfft(fft * np.conj(fft))[:, :nn] / nn
    w = acov[:, 0].mean() * nn / (nn - 1)
    var_plus = w * (nn - 1) / nn + zs.mean(axis=1).var(ddof=1)
    rho = 1 - (w - acov.mean(axis=0)) / var_plus
    rho[0] = 1.0
    pairs = rho[:-1:2] + rho[1::2]
    keep = np.argmax(pairs < 0) if np.any(pairs < 0) else len(pairs)
    pairs = np.minimum.accumulate(pairs[:keep])
    tau = -1 + 2 * pairs.sum()
    return float(r), float(mm * nn / tau)


if __name__ == "__main__":
    t0 = time.time()
    with Pool(2) as pool:
        res = pool.map(run_chain, range(4))
    runtime = time.time() - t0
    models = np.stack([r[0] for r in res])        # chains x layers x 4 x iters
    preds = np.stack([r[1] for r in res])         # chains x iters x freqs
    misfit = np.stack([r[2] for r in res])        # chains x iters
    acc = np.array([r[3] for r in res])
    infos = [r[4] for r in res]

    params = {
        "Vs1": models[:, 0, 3], "Vs2": models[:, 1, 3],
        "h1": models[:, 0, 0], "misfit": misfit,
    }
    truth = {"Vs1": 1500, "Vs2": 3500, "h1": 100}
    table = {}
    for k, v in params.items():
        post = v[:, BURN:]
        q = np.quantile(post, [0.05, 0.16, 0.5, 0.84, 0.95])
        r, ess = rank_rhat_ess(post)
        table[k] = dict(mean=float(post.mean()), p05=q[0], p16=q[1], median=q[2],
                        p84=q[3], p95=q[4], rhat=r, ess_bulk=ess)
        if k in truth:
            table[k]["truth"] = truth[k]
            table[k]["truth_in_90pct"] = bool(q[0] <= truth[k] <= q[4])
    summary = dict(
        settings=dict(iterations=ITER, adapt_until=ADAPT, burnin=BURN, chains=4,
                      starts=STARTS, seeds=[100, 101, 102, 103], noise_seed=2026),
        parameters={k: {kk: (float(vv) if isinstance(vv, (float, np.floating)) else vv)
                        for kk, vv in d.items()} for k, d in table.items()},
        tuned_sigma_v=[i["sigma_v"].tolist() for i in infos],
        tuned_sigma_h=[i["sigma_h"].tolist() for i in infos],
        acceptance_vs_after_adapt=[i["acceptance_vs"].tolist() for i in infos],
        acceptance_h_after_adapt=[i["acceptance_h"].tolist() for i in infos],
        overall_acceptance=acc.tolist(),
        chi2_true_model=float(np.sum(((OBS - CLEAN) / UNC) ** 2)),
        runtime_s=runtime,
    )
    json.dump(summary, open(f"{OUT}/summary.json", "w"), indent=1)
    np.savez_compressed(f"{OUT}/chains.npz", models=models, preds=preds, misfit=misfit,
                        truth=TRUTH, start=START, f=F, obs=OBS, clean=CLEAN, unc=UNC,
                        starts=np.array(STARTS), burnin=BURN, adapt=ADAPT)
    print(json.dumps(summary, indent=1))
