"""Stage J: tables for the gain-sampled 'final' runs.

results/tables/final.json and final.md: per station, the fit of the final
run (H1) and of its H2 bracket, the sampled gain factors and the implied
absolute gains next to the Rayleigh-wave gains, the crustal Vp/Vs, the
sediment and Vs at fixed depths, and the after/stable compliance ratio of
Stage G.
"""

import json
import sys

import numpy as np

from a_invert import EPOCH, ROOT
from i_paper_figures import STATIONS, load

DEPTHS = (0.5, 1.0, 2.0, 3.0)


def q(x):
    a = np.quantile(x, [0.05, 0.5, 0.95])
    return dict(p05=float(a[0]), median=float(a[1]), p95=float(a[2]))


def main():
    rows = []
    for s in STATIONS:
        res = ROOT / "results" / s
        final = json.loads((res / f"{s}_final.json").read_text())
        row = dict(station=s, n=final["n_data"], chi2=final["chi2"]["median"], max_rhat=final["max_rhat"],
                   min_ess=final["min_ess"], groups=final["groups"], gains_used=final["input"]["gains"])  # fmt: skip
        for variant in ("final_h2", "final_stable", "gain", "main"):
            path = res / f"{s}_{variant}.json"
            if path.exists():
                j = json.loads(path.read_text())
                row[f"chi2_{variant}"] = j["chi2"]["median"]
                row[f"n_{variant}"] = j["n_data"]
        post, _ = load(s)
        flat = post.theta.reshape(-1, post.theta.shape[-1])
        names = post.names
        row["vp_vs"] = q(flat[:, names.index("vp_vs")])
        row["sediment"] = q(flat[:, names.index("sediment_thickness")])
        row["gain_factor"], row["gain_abs"] = {}, {}
        for g, used in zip(final["groups"], final["input"]["gains"]):
            key = "gain" if len(final["groups"]) == 1 else f"gain_{g}"
            row["gain_factor"][g] = q(flat[:, names.index(key)])
            row["gain_abs"][g] = q(used * flat[:, names.index(key)])
        if len(final["groups"]) == 2:
            i, j = names.index("gain_stable"), names.index("gain_after")
            ratio = final["input"]["gains"][1] * flat[:, j] / (final["input"]["gains"][0] * flat[:, i])
            row["gain_ratio_after_stable"] = q(ratio)
            row["rayleigh_ratio"] = final["input"]["gains"][1] / final["input"]["gains"][0]
        vs = post.vs_draws(np.array(DEPTHS) * 1e3, max_draws=6000)
        row["vs"] = {f"{d:g}": q(vs[:, k]) for k, d in enumerate(DEPTHS)}
        if s in EPOCH:
            g = json.loads((ROOT / "compliance_after" / s / f"{s}_summary.json").read_text())
            row["after_over_stable"] = [min(g["ratio_after_over_stable"]), max(g["ratio_after_over_stable"])]
            row["after_windows"] = g["n_windows"]
        rows.append(row)
    (ROOT / "results" / "tables" / "final.json").write_text(json.dumps(rows, indent=1))
    head = ("| Station | χ²/N (H1) | χ²/N (H2) | R-hat | Vp/Vs | Gain used → sampled | Sediment [m] | "
            + " | ".join(f"Vs {d:g} km" for d in DEPTHS) + " |")  # fmt: skip
    lines = [head, "|" + "---|" * (7 + len(DEPTHS))]
    for r in rows:
        gains = "; ".join(f"{g}: {u:.3f} → {r['gain_abs'][g]['median']:.3f}"
                          for g, u in zip(r["groups"], r["gains_used"]))  # fmt: skip
        h2 = f"{r['chi2_final_h2']:.1f}/{r['n_final_h2']}" if "chi2_final_h2" in r else "–"
        vs = " | ".join(f"{v['median']:.0f} ({v['p05']:.0f}–{v['p95']:.0f})" for v in r["vs"].values())
        lines.append(f"| {r['station']} | {r['chi2']:.1f}/{r['n']} | {h2} | {r['max_rhat']:.3f} | "
                     f"{r['vp_vs']['median']:.2f} ({r['vp_vs']['p05']:.2f}–{r['vp_vs']['p95']:.2f}) | {gains} | "
                     f"{r['sediment']['median']:.0f} ({r['sediment']['p05']:.0f}–{r['sediment']['p95']:.0f}) | {vs} |")  # fmt: skip
    (ROOT / "results" / "tables" / "final.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    for r in rows:
        if "gain_ratio_after_stable" in r:
            g = r["gain_ratio_after_stable"]
            print(f"{r['station']}: after/stable gain sampled {g['median']:.3f} ({g['p05']:.3f}–{g['p95']:.3f}) "
                  f"vs Rayleigh {r['rayleigh_ratio']:.3f}; compliance after/stable {r['after_over_stable']}")  # fmt: skip


if __name__ == "__main__":
    sys.exit(main())
