"""Bounded saved-metadata audit; no waveform reads, fitting or source writes."""

import hashlib
import json
from collections import Counter
from pathlib import Path
from uuid import uuid4

import numpy as np
from obspy import UTCDateTime, read_inventory


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analog(stage, frequencies, reflect=False):
    """Explicit FDSN forward PAZ, including gain; reflection is diagnostic only."""
    assert stage.pz_transfer_function_type == "LAPLACE (RADIANS/SECOND)"
    s = 2j * np.pi * frequencies
    h = np.full(s.shape, stage.normalization_factor * stage.stage_gain, complex)
    for zero in stage.zeros:
        h *= s - complex(zero)
    for pole in stage.poles:
        p = complex(pole)
        if reflect and p.real > 0:
            p = complex(-p.real, p.imag)
        h /= s - p
    return h


def main():
    workflow_path = Path("results/NETWORKS/TEST/_workflow/rayleigh.json")
    workflow = json.loads(workflow_path.read_text())
    summaries = []
    protected = {workflow_path, *Path("exports").glob("master_inventory.*")}
    for entry in workflow["entries"].values():
        files = [Path(p) for p in entry.get("files", [])]
        protected.update(p for p in files if p.is_file())
        path = next((p for p in files if p.name == "comparison.json"), None)
        if path:
            summaries.append((path, json.loads(path.read_text())))
    snapshots = {}
    for path, summary in summaries:
        station = summary["station"].split(".")[-1]
        if station not in {"RR03", "RR22", "RR52"}:
            continue
        key = next(k for k in summary["epoch_keys"] if k["Seed_id"].endswith("BDH"))
        root = path.parents[2].parent
        inv_path = next(root.glob("*BDH*/inventory.xml"))
        inventory = read_inventory(str(inv_path))
        channels = [
            c
            for n in inventory
            for st in n
            for c in st
            if c.code == "BDH"
            and str(c.start_date or "") == key["channel_start_time"]
            and str(c.end_date or "") == key["channel_end_time"]
        ]
        assert len(channels) == 1, "Snapshot epoch ambiguity"
        snapshots[station] = (inv_path, channels[0].response)
        protected.add(inv_path)
    before = {str(p): sha(p) for p in sorted(protected)}
    audit = {
        "scope": "Saved response snapshots and candidate timing metadata only",
        "handoff_baseline_commit": "833cc05",
        "convention_source": "https://docs.fdsn.org/projects/stationxml/en/latest/response.html",
        "workflow_execution_status": workflow["execution_status"],
        "workflow_status_counts": workflow["status_counts"],
        "responses": [],
        "windows": [],
        "limitations": [
            "Stage PAZ agreement verifies interpretation, not physical calibration",
            "Reflected-pole ratios are algebra only; no corrected response is created",
            "First-arrival instants exclude neither body-wave duration nor coda",
            "Group-arrival ranges are predictions, not observed packet durations",
            "No new boundaries, event decisions, gains or scores are assigned",
        ],
    }
    frequencies = np.array([0.03, 0.05, 0.07])
    for station, (inv_path, response) in sorted(snapshots.items()):
        for stage in response.response_stages:
            if not hasattr(stage, "poles"):
                continue
            direct = analog(stage, frequencies)
            evaluated = response.get_evalresp_response_for_frequencies(
                frequencies,
                output="DEF",
                start_stage=stage.stage_sequence_number,
                end_stage=stage.stage_sequence_number,
            )
            relative = float(np.max(np.abs(evaluated / direct - 1)))
            assert relative < 1e-10, "Independent PAZ/evalresp mismatch"
            ratio = direct / analog(stage, frequencies, reflect=True)
            audit["responses"].append(
                {
                    "station": station,
                    "inventory": str(inv_path),
                    "stage": stage.stage_sequence_number,
                    "input_units": stage.input_units,
                    "output_units": stage.output_units,
                    "frequency_hz": frequencies.tolist(),
                    "poles": [[complex(p).real, complex(p).imag] for p in stage.poles],
                    "direct_phase_deg": np.angle(direct, deg=True).tolist(),
                    "max_complex_relative_evalresp_error": relative,
                    "original_over_reflected_amplitude": np.abs(ratio).tolist(),
                    "original_over_reflected_phase_deg": np.angle(
                        ratio, deg=True
                    ).tolist(),
                }
            )
    for path, summary in summaries:
        station = summary["station"].split(".")[-1]
        for event in summary.get("events", []):
            if event["status"] != "candidate":
                continue
            assert summary["settings"]["segment_seconds"] == 200
            assert summary["settings"]["window_seconds"] == 1800
            start, end = map(UTCDateTime, (event["start_utc"], event["end_utc"]))
            arrivals = {
                kind: UTCDateTime(event[key])
                for kind, key in (
                    ("P", "predicted_first_p_utc"),
                    ("S", "predicted_first_s_utc"),
                )
                if event.get(key)
            }
            inside = {k: t for k, t in arrivals.items() if start <= t <= end}
            cuts = sorted(
                {
                    start.timestamp,
                    end.timestamp,
                    *(t.timestamp for t in inside.values()),
                }
            )
            longest = max(b - a for a, b in zip(cuts, cuts[1:]))
            last = max(arrivals.values()) if arrivals else None
            tail = max(0.0, end - max(start, last)) if last else None
            group = event.get("group_arrival_range_seconds")
            audit["windows"].append(
                {
                    "station": station,
                    "event": event["id"],
                    "result_file": str(path),
                    "distance_km": event["distance_km"],
                    "start_utc": str(start),
                    "end_utc": str(end),
                    "body_arrivals_inside": list(inside),
                    "longest_first_arrival_free_span_seconds": longest,
                    "can_avoid_first_arrival_instants_with_1600s": longest >= 1600,
                    "seconds_after_last_predicted_first_body_arrival": tail,
                    "can_keep_1600s_after_last_first_body_arrival": tail >= 1600
                    if tail is not None
                    else None,
                    "predicted_group_range_width_seconds": group[1] - group[0]
                    if group
                    else None,
                }
            )
    flagged = [w for w in audit["windows"] if w["body_arrivals_inside"]]
    audit["counts"] = {
        "candidate_windows": len(audit["windows"]),
        "with_first_body_arrival_inside": len(flagged),
        "flagged_without_1600s_first_arrival_free_span": sum(
            not w["can_avoid_first_arrival_instants_with_1600s"] for w in flagged
        ),
        "flagged_without_1600s_after_last_first_body_arrival": sum(
            not w["can_keep_1600s_after_last_first_body_arrival"] for w in flagged
        ),
        "candidate_windows_by_event": dict(
            Counter(w["event"] for w in audit["windows"])
        ),
    }
    assert before == {str(p): sha(p) for p in sorted(protected)}, (
        "Protected product changed"
    )
    audit["protected_sha256"] = before
    audit["protected_products_preserved"] = True
    out = Path("results/VALIDATION/YV_rayleigh") / f"conventions_{uuid4().hex[:12]}"
    out.mkdir(exist_ok=False)
    (out / "validation.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"output": str(out), "counts": audit["counts"]}, indent=2))


if __name__ == "__main__":
    main()
