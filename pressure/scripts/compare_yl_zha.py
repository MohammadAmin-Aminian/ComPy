"""Compare saved YL gains with Table 1; never fit or write masters."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from uuid import uuid4

import numpy as np

SOURCE = "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688"
# Table 1: gain coefficient, standard deviation, retained events, depth (m).
PUBLISHED = {
    "A02W": (1.50, 0.01, 19, 2015),
    "A06W": (1.53, 0.03, 19, 2137),
    "A08W": (1.43, 0.02, 17, 2290),
    "A10W": (1.82, 0.01, 17, 1879),
    "A12W": (1.67, 0.03, 17, 2574),
    "A14W": (1.97, 0.10, 9, 884),
    "B01W": (0.19, 0.02, 18, 1723),
    "B04W": (1.41, 0.06, 19, 2197),
    "B09W": (1.32, 0.01, 19, 2626),
    "B11W": (1.90, 0.09, 17, 2060),
    "B14W": (1.62, 0.04, 19, 2219),
    "C01W": (1.85, 0.08, 17, 1655),
    "C03W": (1.55, 0.02, 19, 2084),
    "C05W": (1.77, 0.01, 18, 2289),
    "C07W": (1.99, 0.01, 11, 2648),
    "C08W": (1.78, 0.02, 19, 2785),
    "C09W": (1.84, 0.02, 19, 2537),
    "C10W": (1.51, 0.08, 7, 2162),
    "C11W": (1.83, 0.03, 16, 2214),
    "C16W": (1.79, 0.03, 19, 2213),
    "C17W": (1.43, 0.01, 19, 2822),
    "F02W": (1.47, 0.04, 17, 2602),
    "N02W": (1.73, 0.03, 19, 2717),
    "N03W": (1.80, 0.02, 19, 2304),
    "S01W": (1.81, 0.10, 16, 1980),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    workflow_path = Path("results/NETWORKS/YL_ALL_NO_TILT/_workflow/rayleigh.json")
    workflow = json.loads(workflow_path.read_text())
    if workflow.get("execution_status") != "complete":
        raise ValueError("Wait for the all-station YL run to finish")
    protected = {workflow_path, *Path("exports").glob("master_inventory.*")}
    before = {str(p): sha(p) for p in protected}
    rows, details = [], []
    for entry in sorted(workflow["entries"].values(), key=lambda e: e["seed_id"]):
        path = Path(entry["update"]["rayleigh_result_file"])
        before[str(path)] = sha(path)
        result = json.loads(path.read_text())
        station = entry["seed_id"].split(".")[1]
        published, std, count, depth = PUBLISHED[station]
        g = result.get("median_pressure_multiplier")
        coefficient = 1 / float(g) if g else None
        events = [e for e in result.get("events", []) if e["status"] == "candidate"]
        large = [e for e in events if e["magnitude"] >= 7]
        large_g = (
            float(np.median([e["pressure_multiplier"] for e in large]))
            if len(large) >= 3
            else None
        )
        row = dict(
            station=station,
            status=entry["status"],
            pressure_multiplier=g,
            reciprocal_comparison_coefficient=coefficient,
            published_coefficient=published,
            published_standard_deviation=std,
            relative_difference_percent=100 * (coefficient / published - 1)
            if coefficient
            else None,
            candidate_events=len(events),
            published_events=count,
            candidate_M_ge_7_events=len(large),
            saved_M_ge_7_subset_reciprocal=1 / large_g if large_g else None,
            depth_m=result.get("settings", {}).get("depth_m"),
            published_depth_m=depth,
        )
        rows.append(row)
        reasons = Counter(
            e.get("reason", "")
            for e in result.get("events", [])
            if e.get("reason") and e["status"] != "excluded_cluster"
        )
        details.append(
            dict(
                station=station,
                source_result=str(path),
                rejected_reasons=dict(reasons),
                candidate_pressure_multipliers=[
                    dict(
                        event=e["id"],
                        origin_utc=e["origin_utc"],
                        magnitude=e["magnitude"],
                        gain=e["pressure_multiplier"],
                    )
                    for e in events
                ],
            )
        )
    comparable = [r for r in rows if r["reciprocal_comparison_coefficient"] is not None]
    peers = [r for r in comparable if r["station"] != "C07W"]
    statistics = dict(
        stations=len(rows),
        numerical_comparisons=len(comparable),
        unresolved=[
            r["station"] for r in rows if r["reciprocal_comparison_coefficient"] is None
        ],
        peer_count_excluding_C07W=len(peers),
        peer_difference_range_percent=[
            min(r["relative_difference_percent"] for r in peers),
            max(r["relative_difference_percent"] for r in peers),
        ],
        peer_median_difference_percent=float(
            np.median([r["relative_difference_percent"] for r in peers])
        ),
    )
    assert all(sha(Path(p)) == value for p, value in before.items()), (
        "Source/master changed"
    )
    audit = dict(
        scope="Saved gains: reciprocal comparison only, not calibration acceptance",
        published_source=SOURCE,
        workflow=str(workflow_path),
        options=workflow["options"],
        statistics=statistics,
        stations=rows,
        event_details=details,
        protected_sha256=before,
        protected_inputs_preserved=True,
        correction_applied=False,
        caveats=[
            "Check nominal response versions before absolute comparisons",
            "Inverting coefficients does not equate different fit objectives",
            "Current band/model/windows and event selection differ from the paper",
            "M>=7 subset uses already accepted events; the paper's event set differs",
        ],
    )
    out = Path("results/VALIDATION/YL_Zha_comparison") / f"run_{uuid4().hex[:12]}"
    out.mkdir(parents=True, exist_ok=False)
    for root, stem in [(out, "comparison"), (Path("docs"), "yl_all_zha_comparison")]:
        (root / f"{stem}.json").write_text(json.dumps(audit, indent=2) + "\n")
        with (root / f"{stem}.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps(statistics, indent=2))
    for r in rows:
        coefficient, difference = (
            r["reciprocal_comparison_coefficient"],
            r["relative_difference_percent"],
        )
        print(
            r["station"],
            f"{coefficient:.3f}" if coefficient else "unresolved",
            r["published_coefficient"],
            f"{difference:+.1f}%" if difference is not None else "",
            f"events {r['candidate_events']}/{r['published_events']}",
        )
    print("Output:", out)


if __name__ == "__main__":
    main()
