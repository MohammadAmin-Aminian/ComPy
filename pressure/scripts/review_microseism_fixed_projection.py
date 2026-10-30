"""Read-only test of saved quiet-trained weights on disjoint saved P/H/H/Z data.

Not a new calibration, geographic rotation or generally signal-preserving
correction. The primary hypothesis is chosen by date before target amplitudes
are examined. No raw waveform or new model/fit is computed.
"""

import argparse
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np
from obspy import UTCDateTime

from compy_pressure.godin_pilot import write_csv
from compy_pressure.microseism_fixed_projection import (
    RULE_ID,
    ideal_rayleigh_distortion,
    prior_projections,
    projection_products,
)
from compy_pressure.ppsd import atomic_json
from scripts.combine_microseism_weeks import InputGuard

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

ROOT = Path(__file__).resolve().parents[1]


def finite_median(values):
    values = np.asarray(values)
    return (
        float(np.median(values[np.isfinite(values)]))
        if np.isfinite(values).any()
        else None
    )


def run(pilot_snapshot, quiet_snapshot, output_root, docs=None):
    guard = InputGuard()
    snapshot = guard.json(pilot_snapshot)
    pilot_path = guard.protect(
        snapshot["full_pilot_path"], snapshot["full_pilot_sha256"]
    )
    pilot = guard.json(pilot_path)
    for path, expected in pilot["input_sha256"].items():
        guard.protect(path, expected)
    quiet_copy = guard.json(quiet_snapshot)
    quiet = guard.json(quiet_copy["files"]["result"])
    if quiet_copy["events"] != quiet["events"] or quiet["station"] != pilot["station"]:
        raise ValueError("Quiet summary is stale or belongs to a different station")
    epochs = {r["Seed_id"]: r for r in quiet["epoch_keys"]}
    ids = pilot["channel_order"]
    if set(epochs) != set(ids):
        raise ValueError("Quiet fit and pilot channel IDs differ")
    if quiet["settings"].get("tilt_method") != "tiskitpy" or not all(
        seed.rsplit(".", 1)[1].endswith(component)
        for seed, component in zip(ids[1:], ("1", "2", "Z"), strict=True)
    ):
        raise ValueError("Expected the recorded TiSKitPy H1,H2,Z basis")
    pair = guard.json(guard.json(snapshot["gain_band_audit"])["full_audit_path"])
    matched = [
        w
        for w in pair["weeks"]
        if w["station"] == pilot["station"]
        and w["interval_utc"] == pilot["source_week"]
    ]
    if len(matched) != 1:
        raise ValueError("Pilot week is not uniquely bound to the gain-band audit")
    for name in ("pressure", "vertical"):
        expected = matched[0]["response_audit"][name]["epoch"]
        if epochs[expected["Seed_id"]] != expected:
            raise ValueError("Quiet and test P/Z epochs differ")
    for key in ("source_inventory_sha256", "response_snapshot_sha256"):
        for path, expected in quiet[key].items():
            guard.protect(path, expected)
            if pilot["input_sha256"].get(str(Path(path).resolve())) != expected:
                raise ValueError("Quiet and pilot response/inventory hashes differ")
    for path, expected in quiet["inputs"].items():
        guard.protect(ROOT / path, expected)
    for metadata in quiet["source_files"]:
        path = Path(metadata["path"])
        stat = path.stat()
        if (stat.st_size, stat.st_mtime_ns) != (metadata["size"], metadata["mtime_ns"]):
            raise ValueError("Historical quiet/event raw metadata changed")
        guard.protect(path)
    test_start = min(map(UTCDateTime, pilot["starts_utc"]))
    candidates = prior_projections(quiet["events"], test_start)
    for path in (
        __file__,
        ROOT / "../compy_pressure/microseism_fixed_projection.py",
        ROOT / "docs/microseism_table_staging.csv",
    ):
        guard.protect(path)
    rows = []
    for interval, file in zip(pilot["starts_utc"], pilot["files"], strict=True):
        path = guard.protect(file["path"], file["sha256"])
        with np.load(path) as saved:
            f, x = saved["frequency_hz"], saved["window_coefficients"]
            expected = np.einsum("wif,wjf->fij", x, x.conj()) / len(x)
            if not np.allclose(expected, saved["cross_matrix"], rtol=1e-12, atol=0):
                raise ValueError("Saved four-channel cross matrix does not reproduce")
        for candidate in candidates:
            products = projection_products(x, candidate["weights_h1_h2_z"])
            for band in matched[0]["bands"]:
                lo, hi = band["frequency_hz"]
                m = (f >= lo) & (f < hi)
                ratio = products["complex_transfer_ratio"][m]
                selected = products["nominal_coherence_squared"][m] >= 0.5
                rows.append(
                    {
                        "start_utc": interval,
                        "training_start_utc": candidate["training_start_utc"],
                        "training_age_days": (
                            UTCDateTime(interval)
                            - UTCDateTime(candidate["training_end_utc"])
                        )
                        / 86400,
                        "role": candidate["role"],
                        "band_low_hz": lo,
                        "band_high_hz": hi,
                        "windows": len(x),
                        "bins": int(m.sum()),
                        "median_vertical_power_ratio": finite_median(
                            products["vertical_power_ratio"][m]
                        ),
                        "median_transfer_amplitude_ratio": finite_median(abs(ratio)),
                        "median_nominal_coherence_squared": finite_median(
                            products["nominal_coherence_squared"][m]
                        ),
                        "median_projected_coherence_squared": finite_median(
                            products["projected_coherence_squared"][m]
                        ),
                        "coherent_nominal_bins": int(selected.sum()),
                        "coherent_bin_phase_change_rms_deg": float(
                            np.sqrt(np.mean(np.angle(ratio[selected], deg=True) ** 2))
                        )
                        if selected.any()
                        else None,
                        "interpretation": "relative-basis test; not gain correction",
                    }
                )
    guard.verify()
    output = (
        Path(output_root).resolve()
        / f"run_{datetime.now(UTC):%Y%m%dT%H%M%SZ}_{uuid4().hex[:8]}"
    )
    if any(output.is_relative_to(Path(p).parent) for p in guard.hashes):
        raise ValueError("Output would be beside protected inputs")
    output.mkdir(parents=True, exist_ok=False)
    primary = candidates[-1]
    result = {
        "rule_id": RULE_ID,
        "station": pilot["station"],
        "primary_selection": (
            "latest completed prior quiet interval; "
            "no target-data fit or event-QC selection"
        ),
        "primary": primary,
        "candidates": candidates,
        "primary_training_test_disjoint": True,
        "test_intervals_previously_inspected": True,
        "geographic_rotation_performed": False,
        "idealized_rayleigh_injection": ideal_rayleigh_distortion(
            primary["weights_h1_h2_z"]
        ),
        "rows": rows,
        "input_sha256": guard.hashes,
        "protected_inputs_unchanged": True,
        "gain_applied": False,
        "master_or_scores_changed": False,
        "status": "sensitivity_only_no_automatic_tilt_correction_adopted",
        "limitations": [
            "Separate historical training is not a preregistered blind test",
            "Historical catalogue quiet screen does not guarantee no transients",
            "Geographic orientation and native channel handedness mapping not repaired",
            "Long training-test separation; tilt/response can vary with time",
            "Idealized injection does not bound arbitrary mixtures/field distortion",
            "No new gain, waveform, score or frozen-rule selection is produced",
        ],
    }
    atomic_json(output / "review.json", result)
    write_csv(output / "bands.csv", rows)
    if docs:
        compact = {k: v for k, v in result.items() if k != "input_sha256"}
        compact |= {
            "full_review_path": str(output / "review.json"),
            "full_review_sha256": InputGuard.digest(output / "review.json"),
            "protected_file_count": len(guard.hashes),
        }
        atomic_json(Path(docs) / "rr52_microseism_fixed_projection.json", compact)
    guard.verify()
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--quiet-summary", type=Path, required=True)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "results/REVIEWS/MICROSEISM_FIXED_PROJECTION",
    )
    parser.add_argument("--docs", type=Path)
    args = parser.parse_args()
    print(run(args.pilot, args.quiet_summary, args.output_root, args.docs))


if __name__ == "__main__":
    main()
