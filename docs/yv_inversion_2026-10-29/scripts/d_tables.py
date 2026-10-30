"""Stage D: variant comparison table from the run summaries.

Writes results/tables/variants.json and variants.md: per station and variant,
the chi-square, convergence, sediment thickness, the crustal Vp/Vs (vpvs) and
the median Vs at fixed depths.
"""

import json
import sys

import numpy as np

from c_figures import STATIONS, VARIANT_LABELS, load
from a_invert import ROOT

DEPTHS = (0.25, 0.5, 1.0, 2.0, 3.0)


def main():
    rows = []
    for station in STATIONS:
        for variant in VARIANT_LABELS:
            path = ROOT / "results" / station / f"{station}_{variant}.json"
            if not path.exists():
                continue
            post, _, s = load(station, variant)
            vs = np.median(
                post.vs_draws(np.array(DEPTHS) * 1e3, max_draws=4000), axis=0
            )
            p = s["parameters"]
            rows.append(dict(
                station=station, variant=variant, chi2=s["chi2"]["median"], n=s["n_data"],
                max_rhat=s["max_rhat"], min_ess=s["min_ess"],
                sediment=p["sediment_thickness"]["median"],
                vp_vs=p["vp_vs"]["median"] if "vp_vs" in p else None,
                vp_vs_90=[p["vp_vs"]["p05"], p["vp_vs"]["p95"]] if "vp_vs" in p else None,
                vs=dict(zip(map(str, DEPTHS), vs.round().tolist())),
                converged=s["max_rhat"] < 1.05 and s["min_ess"] > 100,
            ))  # fmt: skip
    (ROOT / "results" / "tables" / "variants.json").write_text(
        json.dumps(rows, indent=1)
    )
    head = (
        "| Station | Variant | χ²/N | R-hat | Sed. [m] | "
        + " | ".join(f"Vs {d:g} km" for d in DEPTHS)
        + " |"
    )
    lines = [head, "|" + "---|" * (5 + len(DEPTHS))]
    for r in rows:
        flag = "" if r["converged"] else " ✗"
        vpvs = f" (Vp/Vs {r['vp_vs']:.2f})" if r["vp_vs"] else ""
        lines.append(f"| {r['station']} | {VARIANT_LABELS[r['variant']]}{vpvs} | {r['chi2']:.1f}/{r['n']} | "
                     f"{r['max_rhat']:.2f}{flag} | {r['sediment']:.0f} | "
                     + " | ".join(f"{v:.0f}" for v in r["vs"].values()) + " |")  # fmt: skip
    (ROOT / "results" / "tables" / "variants.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    main_rows = json.loads(
        (ROOT / "results" / "tables" / "inversion_summary.json").read_text()
    )
    out = ["| Station | Record | Band [mHz] | χ²/N | R-hat | ESS | Sediment [m] | Vs 0.5 km | Vs 1 km | Vs 2 km "
           "| 50/90 % sensitivity [km] |", "|" + "---|" * 11]  # fmt: skip
    for r in main_rows:
        out.append(f"| {r['station']} | {'stable epoch' if r['epoch'] else 'whole'} | {r['band_mhz']} | "
                   f"{r['chi2_median']:.1f}/{r['n_data']} | {r['max_rhat']:.3f} | {r['min_ess']} | {r['sediment_m']} | "
                   f"{r['vs_0.5km']} | {r['vs_1km']} | {r['vs_2km']} | {r['sens50_km']:.1f}/{r['sens90_km']:.1f} |")  # fmt: skip
    (ROOT / "results" / "tables" / "main.md").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    sys.exit(main())
