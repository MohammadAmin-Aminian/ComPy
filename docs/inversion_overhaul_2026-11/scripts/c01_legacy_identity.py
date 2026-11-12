"""C1: the ComPy 2.x stage E variant ``final`` still gives the archived draws.

Rebuilds the RR52 run of 12_inversion_all_stations (ComPy 2.8.0, variant ``final``)
from its saved data and compares every draw. The sampler changes of ComPy 3
(model log prior, Gaussian gain prior, new covariance terms) must not change a
2.x run.

Run: PYTHONPATH=<worktree> python c01_legacy_identity.py [--processes 2]
"""

import argparse
import json
from pathlib import Path

import numpy as np

import compy.inversion.dram as ci

OLD = Path.home() / "Desktop/ComPy/12_inversion_all_stations/YV/results/RR52/inversion"
CRUST1 = Path.home() / "Desktop/BRUIT-FM/compliance_yv_2026-10-28/literature/crust1"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=2)
    args = parser.parse_args()
    d = np.load(OLD / "RR52_final.npz")
    record = json.loads((OLD / "RR52_final.json").read_text())
    data = ci.ComplianceData(d["freqs"], d["eta"], float(d["depth"]), d["covariance"], station="RR52")
    column = ci.crust1_column(-20.4723, 68.1094, CRUST1)
    model = ci.CrustalModel(moho_depth=column["moho_depth"], sediment_max=1000.0,
                            vp_vs_range=(1.7, 2.4), require_moho_jump=True)
    start = ci.crust1_start(model, column)
    post = ci.sample(data, model, seed=0, processes=args.processes, start=start, start_spread=0.1,
                     n_chains=8, n_iter=80000, n_burn=30000, thin=10, gain_bounds=(0.95, 1.05),
                     temperatures=1.6 ** np.arange(8))
    same = np.array_equal(post.theta, d["theta"])
    diff = float(np.max(np.abs(post.theta - d["theta"]))) if post.theta.shape == d["theta"].shape else None
    out = dict(identical=bool(same), max_abs_difference=diff, shape=list(post.theta.shape),
               archived_version="2.8.0", archived_chi2=record["chi2"]["median"],
               chi2=post.summary()["chi2"]["median"])
    print(json.dumps(out, indent=1))
    target = Path.home() / "Desktop/ComPy/14_inversion_overhaul/c01_legacy_identity.json"
    target.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
