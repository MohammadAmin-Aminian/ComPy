"""E4: summary of the calibration runs of ``e03_sbc.py``.

For every JSON in the folder: coverage of the 50 % and 90 % intervals per parameter
(prior truth) and of Vs at fixed depths (all truths), with the binomial 95 % range
expected for the number of replicates; a chi-square test of rank uniformity
(10 bins); the median error of the Vs median; coverage separately for depths the
run counts as resolved (layer width ratio <= 0.6 and above the reach) and not.

Run: python e04_sbc_summary.py <folder> [--figure out.png]
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import binom, chi2

DEPTHS = np.array([50.0, 150.0, 300.0, 700.0, 1500.0, 3000.0, 5000.0, 7000.0])


def expected(n, p):
    lo, hi = binom.ppf([0.025, 0.975], n, p)
    return lo / n, hi / n


def uniformity(ranks, bins=10):
    counts, _ = np.histogram(ranks, np.linspace(0, 1, bins + 1))
    e = len(ranks) / bins
    stat = float(np.sum((counts - e) ** 2 / e))
    return float(chi2.sf(stat, bins - 1))


def summarize(path):
    run = json.loads(Path(path).read_text())
    res = run["results"]
    n = len(res)
    out = dict(name=Path(path).stem, n=n, truth=run["args"]["truth"], tilt=run["args"]["tilt"],
               gain_error=run["args"]["gain_error"],
               converged=float(np.mean([r["max_rhat"] <= 1.05 for r in res])),
               chi2_over_n=float(np.median([r["chi2"] / r["n_data"] for r in res])),
               expected50=expected(n, 0.5), expected90=expected(n, 0.9))
    if "rank" in res[0]:
        names = res[0]["names"]
        rank = np.array([r["rank"] for r in res])
        in50 = np.array([r["in50"] for r in res])
        in90 = np.array([r["in90"] for r in res])
        out["parameters"] = {nm: dict(cover50=float(in50[:, i].mean()), cover90=float(in90[:, i].mean()),
                                      p_uniform=uniformity(rank[:, i]))
                             for i, nm in enumerate(names)}
    v50 = np.array([r["vs_in50"] for r in res], float)
    v90 = np.array([r["vs_in90"] for r in res], float)
    err = np.array([np.array(r["vs_median"]) / np.array(r["vs_true"]) - 1 for r in res])
    ratio = np.array([r["width_ratio"] for r in res])
    reach = np.array([r["reach"] for r in res])
    resolved = (ratio <= 0.6) & (DEPTHS[None, :] <= reach[:, None])
    out["vs"] = {f"{int(z)}": dict(cover50=float(v50[:, k].mean()), cover90=float(v90[:, k].mean()),
                                   bias=float(np.median(err[:, k])),
                                   abs_error=float(np.median(np.abs(err[:, k]))),
                                   pointwise_ratio=float(np.median(ratio[:, k])))
                 for k, z in enumerate(DEPTHS)}
    out["vs_all"] = dict(cover50=float(v50.mean()), cover90=float(v90.mean()))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder")
    args = parser.parse_args()
    folder = Path(args.folder)
    summaries = [summarize(p) for p in sorted(folder.glob("*.json")) if p.name != "summary.json"]
    (folder / "summary.json").write_text(json.dumps(summaries, indent=1))
    for s in summaries:
        print(f"\n== {s['name']}: n={s['n']}, truth={s['truth']}, tilt={s['tilt']}, gain error={s['gain_error']}, "
              f"converged {100 * s['converged']:.0f} %, median chi2/N {s['chi2_over_n']:.2f}")
        print(f"   expected coverage (95 % range): 50 %: {s['expected50'][0]:.2f}-{s['expected50'][1]:.2f}, "
              f"90 %: {s['expected90'][0]:.2f}-{s['expected90'][1]:.2f}")
        for name, p in s.get("parameters", {}).items():
            flag = "" if (s["expected50"][0] <= p["cover50"] <= s["expected50"][1]
                          and s["expected90"][0] <= p["cover90"] <= s["expected90"][1]) else "  <-- outside"
            print(f"   {name:20s} 50 %: {p['cover50']:.2f}  90 %: {p['cover90']:.2f}  rank uniformity p={p['p_uniform']:.2f}{flag}")
        print("   Vs at depth        50 %   90 %   median bias  median |error|")
        for z, v in s["vs"].items():
            print(f"   {z:>6s} m          {v['cover50']:.2f}   {v['cover90']:.2f}   {100 * v['bias']:+6.1f} %     {100 * v['abs_error']:5.1f} %")


if __name__ == "__main__":
    main()
