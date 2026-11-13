"""E3: simulation-based calibration and misspecification tests of ``compy.inversion.auto``.

For a station template (its frequencies, water depth, relative covariance and
CRUST1.0 reference), repeat: draw a true model, compute its compliance, add noise
with the covariance, invert, and record where the truth falls in the posterior.

``--truth prior``     truth drawn from the prior of the model (calibration proper:
                      ranks must be uniform, intervals must have nominal coverage).
``--truth gradient``  truth from a different generator (smooth gradients on 60 thin
                      layers, depth-dependent Vp/Vs, another density rule, crustal
                      thickness and sediment of its own): tests the layered model,
                      the Vp rule and the prior against an Earth they do not describe.
``--tilt X``          adds an unmodelled tilt of X (relative, edge to edge) to the data.
``--gain-error X``    multiplies the data by a gain error drawn with sigma X instead
                      of the prior sigma.

Run: PYTHONPATH=<worktree> python e03_sbc.py OUT --station YV.RR52 --n 200 --workers 12
"""

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

DEPTHS = np.array([50.0, 150.0, 300.0, 700.0, 1500.0, 3000.0, 5000.0, 7000.0])
SAMPLER = dict(n_chains=4, n_iter=30000, n_burn=10000, thin=10)


def gradient_truth(rng, moho_ref):
    """A smooth oceanic model that the inversion's layers cannot represent exactly."""
    h_sed = rng.choice([0.0, rng.uniform(5, 400)])
    v_top = rng.uniform(900, 2600)  # layer 2A at the seafloor
    z2a = rng.uniform(200, 800)  # base of the steep gradient
    v2 = rng.uniform(2700, 3500)
    z3 = z2a + rng.uniform(800, 2000)  # top of layer 3
    v3 = rng.uniform(3500, 3950)
    moho = np.clip(rng.normal(moho_ref, 0.2 * moho_ref), 3000, 1.5 * moho_ref)
    v_bot = v3 + rng.uniform(0, 250)
    vm = rng.uniform(max(v_bot + 100, 4200), 4800)
    edges = np.unique(
        np.concatenate(
            [
                np.linspace(0, z2a, 16),
                np.linspace(z2a, z3, 16),
                np.linspace(z3, moho, 28),
            ]
        )
    )
    mid = 0.5 * (edges[:-1] + edges[1:])
    vs = np.interp(mid, [0, z2a, z3, moho], [v_top, v2, v3, v_bot])
    ratio = np.interp(
        mid,
        [0, z2a, z3],
        [rng.uniform(1.9, 2.5), rng.uniform(1.8, 2.0), rng.uniform(1.76, 1.9)],
    )
    vp = ratio * vs
    rho = 1850.0 + 0.165 * vp  # Christensen & Shaw (1970), not the inversion's rule
    rows = [np.column_stack([np.diff(edges), rho, vp, vs])]
    if h_sed > 0:
        e = np.linspace(0, h_sed, 13)
        m = 0.5 * (e[:-1] + e[1:])
        from compy.inversion import dram

        svs = rng.uniform(0.7, 1.5) * dram.hamilton_vs(m)
        svp = dram.hamilton_vp(m) * rng.uniform(0.95, 1.05)
        rows.insert(
            0,
            np.column_stack(
                [np.diff(e), np.full(12, rng.uniform(1600, 2000)), svp, svs]
            ),
        )
    rows.append(np.array([[0.0, 3300.0, 1.78 * vm, vm]]))
    return np.vstack(rows), dict(
        sediment_thickness=float(h_sed), moho_depth=float(moho)
    )


def vs_of(layers, z):
    tops = np.concatenate([[0.0], np.cumsum(layers[:-1, 0])])
    return layers[np.searchsorted(tops, z, side="right") - 1, 3]


def replicate(job):
    (i, template, truth, tilt, gain_error, seed) = job
    from compy.inversion import auto, dram
    from compy.inversion.errors import log_frequency_ramp

    t0 = time.time()
    freqs, depth, cov_rel = (
        np.array(template[k]) for k in ("freqs", "depth", "cov_rel")
    )
    depth = float(depth)
    model = auto.SeafloorModel(
        np.array(template["reference"]), vs_sigma=template["vs_sigma"]
    )
    sigma = template["gain_sigma"]
    rng = np.random.default_rng([seed, i])
    forward = dram.Forward(depth, freqs)
    attraction = dram.attraction_compliance(freqs, depth)
    for _ in range(1000):
        if truth == "prior":
            theta = model.prior_draw(rng)
            layers, extra = model.layers(theta), {}
        else:
            layers, extra = gradient_truth(rng, model.moho_depth)
            theta = None
        try:
            eta = forward(layers)
            break
        except ValueError:
            continue
    g = rng.normal(1.0, sigma if gain_error is None else gain_error)
    clean = attraction + g * (eta - attraction)
    cov = cov_rel * np.outer(clean, clean)
    d = clean + np.linalg.cholesky(cov) @ rng.standard_normal(len(freqs))
    if tilt:
        d = d * (1 + tilt * log_frequency_ramp(freqs))
    data = dram.ComplianceData(
        freqs, d, depth, cov_rel * np.outer(d, d), station=f"sbc{i}"
    )
    res = auto.invert(
        data,
        model,
        gain_sigma=sigma,
        seed=int(rng.integers(2**31)),
        processes=1,
        max_rounds=1,
        n_prior=1000,
        **SAMPLER,
    )
    post = res.posterior
    flat = post.theta.reshape(-1, post.theta.shape[-1])
    sub = flat[np.linspace(0, len(flat) - 1, 400).astype(int)]
    s = post.summary()
    out = dict(
        i=i,
        seconds=time.time() - t0,
        max_rhat=s["max_rhat"],
        min_ess=s["min_ess"],
        chi2=s["chi2"]["median"],
        n_data=len(freqs),
        gain_true=float(g),
        extra=extra,
    )
    vs_true = vs_of(layers, DEPTHS)
    vs_post = np.array([model.vs_at(t[: len(model.names)], DEPTHS) for t in sub])
    vs_prior = np.array(
        [model.vs_at(t[: len(model.names)], DEPTHS) for t in res.prior_theta[:400]]
    )
    q = np.quantile(vs_post, [0.05, 0.25, 0.5, 0.75, 0.95], axis=0)
    w = np.log(np.quantile(vs_post, 0.84, 0) / np.quantile(vs_post, 0.16, 0))
    wp = np.log(np.quantile(vs_prior, 0.84, 0) / np.quantile(vs_prior, 0.16, 0))
    out.update(
        vs_true=vs_true.tolist(),
        vs_median=q[2].tolist(),
        vs_in50=((vs_true >= q[1]) & (vs_true <= q[3])).tolist(),
        vs_in90=((vs_true >= q[0]) & (vs_true <= q[4])).tolist(),
        vs_rank=(vs_post < vs_true).mean(0).tolist(),
        width_ratio=(w / wp).tolist(),
        reach=float(res.reach(0.9)),
    )
    if theta is not None:
        full = np.append(theta, g)
        out["rank"] = (sub < full).mean(0).tolist()
        pq = np.quantile(flat, [0.05, 0.25, 0.75, 0.95], axis=0)
        out["in50"] = ((full >= pq[1]) & (full <= pq[2])).tolist()
        out["in90"] = ((full >= pq[0]) & (full <= pq[3])).tolist()
        out["names"] = post.names
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out")
    parser.add_argument("--station", default="YV.RR52")
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--truth", default="prior", choices=("prior", "gradient"))
    parser.add_argument("--tilt", type=float, default=0.0)
    parser.add_argument("--gain-error", type=float, default=None)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--name", default=None)
    args = parser.parse_args()
    import archive

    from compy.inversion import auto

    data, column, info = archive.station_inputs(args.station)
    template = dict(
        freqs=data.freqs.tolist(),
        depth=data.depth,
        cov_rel=(data.covariance / np.outer(data.eta, data.eta)).tolist(),
        reference=auto.reference_from_crust1(column).tolist(),
        vs_sigma=auto.model_from_crust1(column).vs_sigma,
        gain_sigma=auto.GAIN_SIGMA[info["gain_source"]],
    )
    jobs = [
        (i, template, args.truth, args.tilt, args.gain_error, args.seed)
        for i in range(args.n)
    ]
    name = args.name or f"{args.station}_{args.truth}"
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    results = []
    t0 = time.time()
    with mp.get_context("spawn").Pool(args.workers) as pool:
        for r in pool.imap_unordered(replicate, jobs):
            results.append(r)
            if len(results) % 20 == 0:
                print(
                    f"{name}: {len(results)}/{args.n} after {time.time() - t0:.0f} s",
                    flush=True,
                )
                (out / f"{name}.json").write_text(
                    json.dumps(
                        dict(template=template, args=vars(args), results=results)
                    )
                )
    (out / f"{name}.json").write_text(
        json.dumps(dict(template=template, args=vars(args), results=results))
    )
    print(f"{name}: done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
