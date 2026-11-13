"""Collect the results of ``r02_batch.py`` into one table.

Run: python r03_collect.py OUTDIR [--tag auto]
Writes OUTDIR/tables/results_<tag>.csv: one row per station with the quality label,
fit, convergence, key parameters with their priors, the resolved depth range, Vs at
fixed depths (median, 16 %, 84 %, prior median, resolved flag) and, where the 2.8
all-station run has the station, its Vs at the same depths.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

DEPTHS = (100, 300, 500, 1000, 2000, 3000, 5000)
OLD = (
    Path.home() / "Desktop/ComPy/12_inversion_all_stations/tables/inversion_results.csv"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out")
    parser.add_argument("--tag", default="auto")
    args = parser.parse_args()
    out = Path(args.out)
    old = {}
    if OLD.exists():
        with open(OLD) as fh:
            old = {r["sid"]: r for r in csv.DictReader(fh)}
    rows = []
    for path in sorted(out.glob(f"*/*/*_{args.tag}.json")):
        s = json.loads(path.read_text())
        info, q = s["info"], s["quality"]
        sid = info["sid"]
        with np.load(path.with_suffix(".npz")) as d:
            chi2_min = float(-2 * d["loglik"].max())
        prof = np.genfromtxt(
            path.with_name(path.stem + "_profile.csv"), delimiter=",", names=True
        )
        p = s["parameters"]
        row = dict(
            sid=sid,
            tier=info["tier"],
            depth_m=info["depth"],
            latitude=info["latitude"],
            longitude=info["longitude"],
            gain_source=info["gain_source"],
            response=info["response"],
            tilt_rel=info["tilt_rel"],
            n_windows=info["n_windows"],
            n_data=info["n_data"],
            f_lo=info["band"][0],
            f_hi=info["band"][1],
            quality=q["label"],
            reasons=" | ".join(q["reasons"]),
            chi2_median=q["chi2"],
            chi2_min=chi2_min,
            chi2_limit=q["chi2_limit"],
            max_rhat=q["max_rhat"],
            min_ess=q["min_ess"],
            rounds=s["rounds"],
            seconds=info.get("seconds"),
            gain=q["gain"],
            gain_pull=q["gain_pull"],
            reach_m=s["reach_m"],
            resolved_from_m=s["resolved_from_m"],
            resolved_to_m=s["resolved_to_m"],
            crust1_moho=s["model"]["moho_depth"],
            crust1_sediment=info["crust1"]["sediment_thickness"],
        )
        for name in (
            "sediment_thickness",
            "sediment_vs_scale",
            "sediment_modulus",
            "moho_depth",
            "vs_mantle",
            "vp_upper",
            "vp_lower",
        ):
            for key in ("p16", "median", "p84"):
                row[f"{name}_{key}"] = p[name][key]
            row[f"{name}_prior"] = s["prior"][name]["median"]
        for z in DEPTHS:
            i = int(np.argmin(np.abs(prof["depth_m"] - z)))
            row[f"vs_{z}"] = prof["vs_p50"][i]
            row[f"vs_{z}_p16"] = prof["vs_p16"][i]
            row[f"vs_{z}_p84"] = prof["vs_p84"][i]
            row[f"prior_{z}"] = prof["prior_p50"][i]
            row[f"resolved_{z}"] = int(prof["resolved"][i])
            if sid in old and f"vs_{z}" in old[sid]:
                row[f"old_vs_{z}"] = old[sid][f"vs_{z}"]
        if sid in old:
            row.update(
                old_quality=old[sid]["quality"],
                old_chi2_over_n=old[sid]["chi2_over_n"],
                old_sediment=old[sid]["sediment_m"],
            )
        rows.append(row)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    (out / "tables").mkdir(exist_ok=True)
    with open(out / "tables" / f"results_{args.tag}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} stations -> {out / 'tables' / f'results_{args.tag}.csv'}")


if __name__ == "__main__":
    main()
