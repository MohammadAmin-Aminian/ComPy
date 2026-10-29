"""Bayesian Vs inversion of a compliance estimate with compy_inversion (ComPy 2.3).

Usage (repository root):
    python _Example/4_Vs_inversion.py                      # offline synthetic
    python _Example/4_Vs_inversion.py RR52_compliance_compy.npz \
        --lat -20.4723 --lon 68.1094 --crust1 path/to/crust1.0

The input file is the output of compy_processing.stacked_compliance saved as
npz (keys freqs, compliance, uncertainty, band, depth). Without CRUST1.0 the
Moho is put 7 km below the seafloor. About 15 s on 8 cores with Numba.
"""

import argparse

import numpy as np

import compy_inversion as ci


def synthetic():
    depth, f = 2918.0, np.arange(6, 15) * 1.220703125e-3
    model = ci.CrustalModel(moho_depth=6900.0)
    truth = np.array([150.0, 1800, 2200, 2900, 3500, 3800, 3950, 4300])
    clean = ci.forward(depth, f, model.layers(truth))
    noise = ci.ComplianceData.build(
        f, clean, depth, random=0.003 * clean, gain_rel=0.05, floor_rel=0.01
    )
    obs = clean + np.linalg.cholesky(noise.covariance) @ np.random.default_rng(
        1
    ).standard_normal(len(f))
    data = ci.ComplianceData.build(f, obs, depth, random=0.003 * obs, gain_rel=0.05,
                                   floor_rel=0.01, station="synthetic")  # fmt: skip
    print("true parameters:", dict(zip(model.names, truth)))
    return data, model


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("npz", nargs="?")
    parser.add_argument("--lat", type=float)
    parser.add_argument("--lon", type=float)
    parser.add_argument("--crust1")
    parser.add_argument("--output", default="vs_inversion")
    args = parser.parse_args()
    if args.npz is None:
        data, model = synthetic()
    else:
        data = ci.ComplianceData.from_npz(args.npz, gain_rel=0.05, floor_rel=0.01)
        moho = 7000.0
        if args.crust1:
            moho = ci.crust1_column(args.lat, args.lon, args.crust1)["moho_depth"]
        model = ci.CrustalModel(moho_depth=moho)
    post = ci.sample(
        data, model, n_chains=8, n_iter=60000, n_burn=20000, thin=10, seed=0
    )
    s = post.summary()
    print(f"chi2 median {s['chi2']['median']:.1f} for N={s['n_data']}, max R-hat {s['max_rhat']:.3f}, "
          f"min ESS {s['min_ess']:.0f}")  # fmt: skip
    for name, p in s["parameters"].items():
        print(f"  {name:20s} {p['median']:7.0f}  (90 %: {p['p05']:.0f}-{p['p95']:.0f})")
    z = np.array([0.1, 0.5, 1.0, 2.0, 3.0]) * 1e3
    q = post.vs_quantiles(z, q=(0.05, 0.5, 0.95))
    for depth, lo, med, hi in zip(z, *q):
        print(f"  Vs at {depth / 1e3:.1f} km: {med:.0f} m/s (90 %: {lo:.0f}-{hi:.0f})")
    post.save(args.output)
    print(f"saved {args.output}.npz and {args.output}.json")


if __name__ == "__main__":
    main()
