"""Saved-spectrum pilot of fixed YV response hypotheses; development only.

Approved by the owner on 21 October 2026 after
``docs/yv_response_confirmation.md``. Each station gets one frozen hypothesis
from the existing diagnosis helper (no divisor search, no new variant):

* RR22, RR56: H2, Table A1 conjugate pairs restored (pressure as H1);
* RR17, RR18: H4, pressure pole moved to -p/10, labelled post hoc;
* controls: RR52 identity (must reproduce g = 0.6787278140712512 exactly or
  the pilot stops), RR16 identity under H2 (trial g), RR10 negative H4.

Only saved event spectra of the completed pole-only trial are read. Before any
hypothesis is applied, H0 and H1 must replay the trial's nominal and pole
per-event outcomes. Coherence, phase, overlap, event and fit rules are the
unchanged ``evaluate``/``absolute`` path. Every input is hashed before and
after; output goes to a new exclusive folder. No response, XML, master, saved
result, default, score or rating is changed, and no gain is applied anywhere.
Events here were inspected during the diagnosis: results are development
evidence, not validation.
"""

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from compy_pressure.calibration_experiment import apply_response_hypothesis, evaluate
from scripts.diagnose_yv_rayleigh import (
    alternative_stages,
    laplace_stages,
    pz_value,
    sha256,
    snapshot_channel,
    stage_ratio,
    wrap,
)

RR52_G = 0.6787278140712512
PLAN = {
    "RR22": {"hypotheses": ["H2"], "role": "pilot", "label": "supported metadata"},
    "RR56": {"hypotheses": ["H2"], "role": "pilot", "label": "supported metadata"},
    "RR17": {"hypotheses": ["H4"], "role": "pilot", "label": "post hoc"},
    "RR18": {"hypotheses": ["H4"], "role": "pilot", "label": "post hoc"},
    "RR52": {
        "hypotheses": ["H2", "H4"],
        "role": "identity_control",
        "expected_g": RR52_G,
        "exact": True,
    },
    "RR16": {
        "hypotheses": ["H2"],
        "role": "identity_control",
        "expected_g": 1.078308697242057,
        "exact": False,
    },
    "RR10": {"hypotheses": ["H4"], "role": "negative_control"},
}
TRIAL_VARIANTS = {"H0": "nominal_absolute", "H1": "pole_absolute"}
REPLAY_ATOL = 1e-10


class Inputs:
    """Hash every consumed file on first use and recheck before writing."""

    def __init__(self):
        self.hashes = {}

    def bind(self, path, expected=None):
        path = str(Path(path).resolve())
        if path not in self.hashes:
            self.hashes[path] = sha256(path)
        if expected is not None and self.hashes[path] != expected:
            raise ValueError(f"Input hash differs from recorded value: {path}")
        return path

    def json(self, path, expected=None):
        return json.loads(Path(self.bind(path, expected)).read_text(encoding="utf-8"))

    def verify(self):
        for path, digest in self.hashes.items():
            if sha256(path) != digest:
                raise ValueError(f"Input changed during the pilot: {path}")


def same_outcome(mine, saved):
    if mine["status"] != saved["status"]:
        return False
    a, b = mine.get("pressure_multiplier"), saved.get("pressure_multiplier")
    if a is None or b is None:
        return a is b
    return bool(np.isclose(a, b, rtol=0, atol=REPLAY_ATOL))


def coherent_bins(spectra, threshold):
    """Bins of events passing evaluate's own coherence rule."""
    c = np.asarray(spectra["coherence_squared"])
    coherent = c >= threshold
    if np.median(c) < threshold or np.mean(c) < threshold or coherent.sum() < 3:
        return np.zeros(c.shape, dtype=bool)
    return coherent


def ratios(frequency, stages, name, hypotheses):
    alt_p, alt_z, polarity = hypotheses[name]
    pressure, vertical = stages
    rp = stage_ratio(frequency, pressure, alt_p) * np.exp(1j * np.deg2rad(polarity))
    return rp, stage_ratio(frequency, vertical, alt_z)


def a0_ratio(channel, alternative):
    """Product over analogue stages of A0(alternative roots) / A0 recorded.

    The hypotheses keep the recorded A0. If A0 were instead recomputed from
    the alternative roots, that channel's response would scale by this ratio.
    """
    stages = [
        s
        for s in channel.response.response_stages
        if (getattr(s, "poles", None) or getattr(s, "zeros", None))
        and getattr(s, "pz_transfer_function_type", "") == "LAPLACE (RADIANS/SECOND)"
    ]
    ratio = 1.0
    for stage, (poles, zeros) in zip(stages, alternative, strict=True):
        value = pz_value(np.array([stage.normalization_frequency]), poles, zeros)[0]
        ratio *= 1.0 / abs(value) / float(stage.normalization_factor)
    return ratio


def summarize(events, name, minimum, threshold):
    outcomes = [e["variants"][name] for e in events]
    fitted = [o["pressure_multiplier"] for o in outcomes if o["status"] == "candidate"]
    top = [
        o["top_bin_removed_gain"]
        for o in outcomes
        if o["status"] == "candidate" and "top_bin_removed_gain" in o
    ]
    phases = np.concatenate([e["phase_in_coherent_bins"][name] for e in events] or [[]])
    statuses = {}
    for outcome in outcomes:
        statuses[outcome["status"]] = statuses.get(outcome["status"], 0) + 1
    gain = float(np.median(fitted)) if len(fitted) >= minimum else None
    return {
        "candidate_events": len(fitted),
        "development_g": gain,
        "diagnostic_median": float(np.median(fitted)) if fitted else None,
        "candidate_multipliers": fitted,
        "median_top_bin_removed_g": float(np.median(top)) if top else None,
        "status_counts": statuses,
        "coherent_bins": int(phases.size),
        "phase_rms_deg": float(np.sqrt(np.mean(phases**2))) if phases.size else None,
        "status": "development_estimate"
        if gain is not None
        else "insufficient_candidates",
        "coherence_threshold": threshold,
    }


def station_pilot(station, plan, workflow, trial_root, inputs):
    entries = [
        e for e in workflow["entries"].values() if e["seed_id"].split(".")[1] == station
    ]
    if len(entries) != 1 or not entries[0]["seed_id"].endswith(".BDH"):
        raise ValueError(f"Expected one pressure workflow entry for {station}")
    entry = entries[0]
    comparison_path = inputs.bind(entry["update"]["rayleigh_result_file"])
    comparison = inputs.json(comparison_path)
    trial = inputs.json(Path(trial_root) / f"YV_{station}" / "trial.json")
    if (
        trial.get("network") != "YV"
        or trial.get("station") != station
        or str(Path(trial.get("baseline_result", "")).resolve()) != comparison_path
    ):
        raise ValueError(f"Trial baseline binding mismatch: {station}")
    settings = comparison["settings"]
    if (
        settings.get("pressure_sign") != -1
        or settings.get("tilt_method") != "none"
        or list(settings.get("band_hz", [])) != [0.03, 0.07]
        or settings.get("coherence_squared_threshold") != 0.95
        or settings.get("phase_limit_deg") != 20.0
    ):
        raise ValueError(f"{station}: not the frozen no-tilt baseline settings")
    for path, digest in comparison["source_inventory_sha256"].items():
        inputs.bind(path, digest)
    epochs = {e["Seed_id"]: e for e in comparison["epoch_keys"]}
    channels = {}
    for path, digest in comparison["response_snapshot_sha256"].items():
        seed_id = Path(path).parent.name.split("__")[0]
        if seed_id not in epochs or seed_id in channels:
            raise ValueError(f"Unbound or duplicate response snapshot: {seed_id}")
        channels[seed_id] = snapshot_channel(
            inputs.bind(path, digest), seed_id, epochs[seed_id]
        )
    vertical_ids = {
        e["tilt"]["vertical_component"] for e in comparison["events"] if "tilt" in e
    }
    pressure_id = entry["seed_id"]
    if len(vertical_ids) != 1 or next(iter(vertical_ids)) not in channels:
        raise ValueError(f"{station}: vertical partner is not unique or unbound")
    vertical_id = vertical_ids.pop()
    stages = (
        laplace_stages(channels[pressure_id]),
        laplace_stages(channels[vertical_id]),
    )
    hypotheses = alternative_stages(*stages)
    names = ["H0", "H1", *plan["hypotheses"]]
    trial_events = {e["id"]: e for e in trial["events"]}
    threshold = settings["coherence_squared_threshold"]
    events = []
    for original in comparison["events"]:
        path = original.get("files", {}).get("spectra")
        if not path:
            continue
        with np.load(inputs.bind(path), allow_pickle=False) as data:
            saved = {k: data[k].copy() for k in data.files}
        f = saved["frequency_hz"]
        overlap = bool(original.get("other_predicted_arrival_ids"))
        record = {"id": original["id"], "baseline_status": original["status"]}
        record["variants"], record["phase_in_coherent_bins"] = {}, {}
        for name in names:
            rp, rz = ratios(f, stages, name, hypotheses)
            spectra = apply_response_hypothesis(saved, rp, rz)
            outcome = evaluate(spectra, settings, overlap, np.ones(f.shape), "absolute")
            record["variants"][name] = outcome
            bins = coherent_bins(spectra, threshold)
            record["phase_in_coherent_bins"][name] = wrap(spectra["phase_deg"][bins])
        trial_event = trial_events.get(original["id"], {})
        for name, variant in TRIAL_VARIANTS.items():
            saved_outcome = trial_event.get("variants", {}).get(variant)
            if saved_outcome is None or not same_outcome(
                record["variants"][name], saved_outcome
            ):
                raise ValueError(
                    f"{station} {original['id']}: {name} does not replay {variant}"
                )
        events.append(record)
    if len(events) != sum("variants" in e for e in trial["events"]):
        raise ValueError(f"{station}: saved event set differs from the trial")
    minimum = settings["minimum_candidate_events"]
    summaries = {name: summarize(events, name, minimum, threshold) for name in names}
    for name, variant in TRIAL_VARIANTS.items():
        expected = trial["summaries"][variant]["median_pressure_multiplier"]
        actual = summaries[name]["development_g"]
        if (actual is None) != (expected is None) or (
            actual is not None
            and not np.isclose(actual, expected, rtol=0, atol=REPLAY_ATOL)
        ):
            raise ValueError(
                f"{station}: {name} station gain does not replay the trial"
            )
    a0_effect = {}
    for name in names:
        alt_p, alt_z, _ = hypotheses[name]
        cp = a0_ratio(channels[pressure_id], alt_p)
        cz = a0_ratio(channels[vertical_id], alt_z)
        # P and A would scale by 1/cp and 1/cz, so g would scale by cp/cz.
        a0_effect[name] = {
            "pressure_a0_ratio": cp,
            "vertical_a0_ratio": cz,
            "g_factor": cp / cz,
        }
    for record in events:
        record["phase_in_coherent_bins"] = {
            k: [float(x) for x in v]
            for k, v in record["phase_in_coherent_bins"].items()
        }
    return {
        "station": station,
        "role": plan["role"],
        "label": plan.get("label", plan["role"]),
        "hypotheses": plan["hypotheses"],
        "pressure_seed_id": pressure_id,
        "vertical_seed_id": vertical_id,
        "comparison": comparison_path,
        "summaries": summaries,
        "if_a0_recomputed_from_hypothesis_roots": a0_effect,
        "events": events,
    }


def control_checks(stations, plan):
    checks = []
    for station, report in stations.items():
        expected_plan = plan.get(station, {})
        role = report["role"]
        for name in report["hypotheses"]:
            g = report["summaries"][name]["development_g"]
            if role == "identity_control":
                expected = expected_plan["expected_g"]
                passed = g is not None and (
                    g == expected
                    if expected_plan["exact"]
                    else bool(np.isclose(g, expected, rtol=0, atol=REPLAY_ATOL))
                )
                checks.append(
                    {
                        "station": station,
                        "hypothesis": name,
                        "check": "identity_g",
                        "expected": expected,
                        "actual": g,
                        "passed": passed,
                    }
                )
            elif role == "negative_control":
                h1 = report["summaries"]["H1"]["phase_rms_deg"]
                mine = report["summaries"][name]["phase_rms_deg"]
                passed = h1 is not None and mine is not None and mine > h1
                checks.append(
                    {
                        "station": station,
                        "hypothesis": name,
                        "check": "phase_worsens",
                        "H1_phase_rms_deg": h1,
                        "actual": mine,
                        "passed": passed,
                    }
                )
    return checks


def summary_rows(stations):
    rows = []
    for station, report in stations.items():
        for name, summary in report["summaries"].items():
            rows.append(
                {
                    "station": station,
                    "role": report["role"],
                    "label": report["label"] if name in report["hypotheses"] else "",
                    "variant": name,
                    "candidate_events": summary["candidate_events"],
                    "development_g": summary["development_g"],
                    "diagnostic_median": summary["diagnostic_median"],
                    "median_top_bin_removed_g": summary["median_top_bin_removed_g"],
                    "phase_rms_deg": summary["phase_rms_deg"],
                    "coherent_bins": summary["coherent_bins"],
                    "status_counts": ";".join(
                        f"{k}={v}" for k, v in sorted(summary["status_counts"].items())
                    ),
                    "g_factor_if_a0_recomputed": report[
                        "if_a0_recomputed_from_hypothesis_roots"
                    ][name]["g_factor"],
                }
            )
    return rows


def run(args, plan=None):
    plan = PLAN if plan is None else plan
    out = Path(args.out)
    if out.exists():
        raise FileExistsError(f"Output folder already exists: {out}")
    inputs = Inputs()
    workflow = inputs.json(args.workflow, args.workflow_sha256)
    inputs.bind(Path(args.trial_root) / "experiment.json", args.trial_sha256)
    stations = {}
    # The RR52 identity control runs first so a failure stops everything.
    for station in sorted(plan, key=lambda name: name != "RR52"):
        stations[station] = station_pilot(
            station, plan[station], workflow, args.trial_root, inputs
        )
        if station == "RR52":
            check = control_checks({station: stations[station]}, plan)
            if not all(c["passed"] for c in check):
                raise RuntimeError(
                    f"RR52 identity control failed; pilot stopped: {check}"
                )
    checks = control_checks(stations, plan)
    inputs.verify()
    report = {
        "tool": "scripts/run_response_pilot.py",
        "code_sha256": sha256(Path(__file__)),
        "experiment_code_sha256": sha256(
            Path(__file__).resolve().parents[1]
            / "../compy_pressure/calibration_experiment.py"
        ),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "development_only_inspected_events",
        "plan": plan,
        "control_checks": checks,
        "input_sha256": inputs.hashes,
        "stations": stations,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir()
    with open(out / "pilot.json", "x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    rows = summary_rows(stations)
    with open(out / "stations.csv", "x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--workflow", required=True)
    parser.add_argument("--workflow-sha256", required=True)
    parser.add_argument("--trial-root", required=True)
    parser.add_argument("--trial-sha256", required=True, help="experiment.json SHA-256")
    parser.add_argument("--out", required=True, help="new folder; must not exist")
    report = run(parser.parse_args())
    for check in report["control_checks"]:
        print(check)


if __name__ == "__main__":
    main()
