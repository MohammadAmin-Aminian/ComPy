"""Read-only phase/model audit of reportable stations in an explicit v1 review.

Run: python -m scripts.audit_microseism_gain_bands --review PATH
Writes fresh JSON/CSV/Markdown under results/REVIEWS/MICROSEISM_GAIN_BANDS/.
No waveform, modal solver, calibration, gain application or master update runs.
"""

import argparse
import csv
import statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np
from obspy import UTCDateTime, read_inventory

from compy_pressure.godin_pilot import write_csv
from compy_pressure.microseism_band_audit import (
    AUDIT_ID,
    COHERENT,
    FITS_DEG,
    MINIMUM_BINS,
    UNEXPLAINED_DEG,
    model_band_context,
    phase_band,
)
from compy_pressure.microseism_check import check_paired_run
from compy_pressure.microseism_validation import (
    RULE_ID,
    frequency_station,
    frequency_week,
)
from compy_pressure.ppsd import atomic_json
from scripts.combine_microseism_weeks import InputGuard
from scripts.compare_microseism_methods import channel_poles

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

ROOT = Path(__file__).resolve().parents[1]


def response_audit(paired, data, guard):
    """Reproduce saved complex responses using exact-epoch XML, not new metadata.

    Agreement verifies processing consistency, not sensor absolute calibration.
    """
    inventories = []
    for path, expected in paired["response_snapshots"].items():
        inventories.append(read_inventory(str(guard.protect(path, expected))))
    result = {}
    for name, response_key, output in (
        ("pressure", "pressure_response", "DEF"),
        ("vertical", "displacement_response", "DISP"),
    ):
        epoch = paired[f"{name}_epoch"]
        net, sta, loc, cha = epoch["Seed_id"].split(".")
        channels = []
        for inventory in inventories:
            selected = inventory.select(
                network=net, station=sta, location=loc, channel=cha
            )
            for n in selected:
                for s in n:
                    for c in s:
                        if str(c.start_date) == epoch["channel_start_time"] and str(
                            c.end_date
                        ) == str(epoch["channel_end_time"]):
                            channels.append(c)
        if not channels:
            raise ValueError(f"No exact response epoch for {epoch['Seed_id']}")
        maximum = 0.0
        for c in channels:
            if not (
                c.start_date <= UTCDateTime(paired["start_utc"])
                and (c.end_date is None or c.end_date >= UTCDateTime(paired["end_utc"]))
            ):
                raise ValueError("Paired interval is outside response epoch")
            actual = c.response.get_evalresp_response_for_frequencies(
                data["frequency_hz"], output=output
            )
            saved = data[response_key]
            if not np.allclose(actual, saved, rtol=1e-10, atol=0):
                raise ValueError(f"Saved complex response differs: {epoch['Seed_id']}")
            maximum = max(maximum, float(np.max(abs(actual - saved) / abs(saved))))
        sensitivity = channels[0].response.instrument_sensitivity
        result[name] = {
            "epoch": epoch,
            "maximum_complex_relative_error": maximum,
            "sensitivity": float(sensitivity.value),
            "sensitivity_frequency_hz": float(sensitivity.frequency),
            "input_units": sensitivity.input_units,
            "dip_deg": float(channels[0].dip),
            "stage_count": len(channels[0].response.response_stages),
            "absolute_response_independently_verified": False,
        }
    return result


def audit_week(row, week, guard):
    run = Path(week["run"])
    summary = guard.json(run / "summary.json")
    comparison = guard.json(run / "method_comparison.json")
    guard.protect(run / "summary.json", comparison["run_summary_sha256"])
    records = [
        r
        for r in summary["stations"]
        if f"{r['network']}.{r['station']}" == row["station"]
    ]
    if len(records) != 1 or records[0]["status"] != "checked":
        raise ValueError("Review station is not uniquely checked in its source run")
    record = records[0]
    check = guard.json(record["check_summary"])
    for path, metadata in check["inputs"].items():
        guard.protect(path, metadata["sha256"])
    if (
        check["station"] != row["station"]
        or check["interval_utc"] != week["interval_utc"]
    ):
        raise ValueError("Check station or interval differs from review")
    if frequency_week(check, week["label"]) != week["frequency_check"]:
        raise ValueError("Saved frozen-rule selection no longer reproduces")
    paired_path = guard.protect(record["paired_result"])
    if str(paired_path) not in {str(Path(p).resolve()) for p in check["inputs"]}:
        raise ValueError("Paired comparison is not bound to check inputs")
    paired, data = check_paired_run(paired_path.parent)
    if paired["pressure_epoch"]["Seed_id"].rsplit(".", 2)[0] != row["station"]:
        raise ValueError("Paired spectra belong to another station")
    if [paired["start_utc"], paired["end_utc"]] != week["interval_utc"]:
        raise ValueError("Paired spectra cover a different interval")
    responses = response_audit(paired, data, guard)
    poles = channel_poles(paired)
    bins_path = guard.protect(check["files"]["bins"])
    with bins_path.open() as stream:
        records = list(csv.DictReader(stream))
    bins = {key: np.array([float(r[key]) for r in records]) for key in records[0]}
    f = data["frequency_hz"]
    indices = np.searchsorted(f, bins["frequency_hz"])
    indices = np.clip(indices, 1, len(f) - 1)
    indices -= abs(f[indices - 1] - bins["frequency_hz"]) < abs(
        f[indices] - bins["frequency_hz"]
    )
    if not np.allclose(f[indices], bins["frequency_hz"], rtol=1e-8, atol=0):
        raise ValueError("Cached model bins do not match paired frequency grid")
    results = []
    for band in week["frequency_check"]["eligible_sub_bands"]:
        edges = band["frequency_hz"]
        if edges[1] >= check["water"]["first_resonance_hz"]:
            raise ValueError("Gain-band phase audit requires the first water branch")
        mask = (f[indices] >= edges[0]) & (f[indices] < edges[1])
        context = model_band_context(bins, mask, row.get("rayleigh_baseline_g"))
        for key, expected in (("water_h1", band["g_h1"]), ("water_h2", band["g_h2"])):
            if not np.isclose(context[key], expected, rtol=1e-5, atol=0):
                raise ValueError("Cached model band differs from frozen review")
        phase = phase_band(
            f,
            np.angle(data["transfer_dp"], deg=True),
            data["coherence_squared"],
            edges,
            week["label"],
            poles["pressure"],
            poles["vertical"],
        )
        days = np.array([s[:10] for s in data["retained_window_starts"]])
        daily = []
        for day in sorted(set(days)):
            keep = days == day
            pp = data["pressure_window_psd"][keep].mean(axis=0)
            dd = data["displacement_window_psd"][keep].mean(axis=0)
            dp = data["cross_window_dp"][keep].mean(axis=0)
            ph = phase_band(
                f,
                np.angle(dp / pp, deg=True),
                abs(dp) ** 2 / (pp * dd),
                edges,
                week["label"],
                poles["pressure"],
                poles["vertical"],
            )
            if keep.sum() < 2:
                ph["status"] = "insufficient_averaged_windows"
            daily.append({"day": day, "windows": int(keep.sum()), **ph})
        results.append(
            {**phase, "frozen_band": band, "model_context": context, "daily": daily}
        )
    return {
        "station": row["station"],
        "shift_weeks": week["shift"],
        "interval_utc": week["interval_utc"],
        "prescribed_hypothesis": week["label"],
        "paired_comparison": str(paired_path),
        "check_summary": record["check_summary"],
        "response_audit": responses,
        "positive_real_poles_rad_s": poles,
        "bands": results,
    }


def table_export(review, audited_rows, guard):
    """Epoch-bound, read-only staging table; never merge or assign pass scores.

    These are research execution/repeatability statuses, not an adopted final
    station classification. Missing and nonreportable gains remain blank.
    """
    audited = {r["station"]: r for r in audited_rows}
    cache, result, seen = {}, [], set()
    for row in review["rows"]:
        base = [w for w in row["weeks"] if w["shift"] == 0]
        if len(base) != 1:
            raise ValueError("Each staging row needs exactly one base run")
        run = Path(base[0]["run"])
        if str(run) not in cache:
            cache[str(run)] = guard.json(run / "summary.json")["stations"]
        records = [
            r
            for r in cache[str(run)]
            if f"{r['network']}.{r['station']}" == row["station"]
        ]
        if len(records) != 1:
            raise ValueError("Staging row has an ambiguous base record")
        record = records[0]
        pair = guard.json(record["paired_result"])
        epoch = pair["pressure_epoch"]
        if (
            epoch["Seed_id"] != record["seed_id"]
            or epoch["Seed_id"].rsplit(".", 2)[0] != row["station"]
        ):
            raise ValueError("Staging row pressure epoch differs from source station")
        key = tuple(
            epoch[k]
            for k in (
                "Seed_id",
                "channel_start_time",
                "channel_end_time",
                "inventory_file",
            )
        )
        if key in seen:
            raise ValueError("Duplicate full epoch key in staging table")
        seen.add(key)
        frozen = row["frequency_consistency"]
        phase = audited.get(row["station"])
        status = (
            "conditional_repeatable_unvalidated"
            if frozen["conditional_reportable"]
            else frozen["status"]
        )
        if row["weeks_checked"] == 0:
            status = "unresolved_pairing"
        result.append(
            {
                **epoch,
                "calibration_check_pressure_vertical_microseisms_score": None,
                "calibration_check_pressure_vertical_microseisms_status": status,
                "microseism_rule_id": RULE_ID,
                "microseism_conditional_g": frozen["median_g"]
                if frozen["conditional_reportable"]
                else None,
                "microseism_passing_weeks": frozen["passing_weeks"],
                "microseism_checked_weeks": row["weeks_checked"],
                "microseism_weekly_phase_bands_consistent": phase["consistent_bands"]
                if phase
                else None,
                "microseism_weekly_phase_bands_audited": phase["audited_bands"]
                if phase
                else None,
                "microseism_phase_audit_status": (
                    "weekly_consistent_daily_review_open"
                    if phase and phase["consistent_bands"] == phase["audited_bands"]
                    else "weekly_phase_review_required"
                    if phase
                    else "not_audited"
                ),
                "microseism_tilt_status": "not_assessed",
                "microseism_absolute_vertical_response_status": (
                    "not_independently_verified"
                ),
                "microseism_base_pair_status": record["status"],
                "microseism_base_pair_reason": record.get("paired_reason"),
                "microseism_base_paired_comparison": record["paired_result"],
            }
        )
    return result


def report(review_path, output_root, stations=None, docs=None):
    guard = InputGuard()
    review = guard.json(review_path)
    if (
        review["frequency_rule"]["id"] != RULE_ID
        or not review["protected_inputs_unchanged"]
    ):
        raise ValueError("Expected a protected-input-verified frozen v1 review")
    for path, expected in review["input_sha256"].items():
        guard.protect(path, expected)
    for path in (
        __file__,
        ROOT / "../compy_pressure/microseism_band_audit.py",
        ROOT / "../compy_pressure/microseism_validation.py",
    ):
        guard.protect(path)
    selected = [
        r
        for r in review["rows"]
        if r["frequency_consistency"]["conditional_reportable"]
        and (not stations or r["station"] in stations)
    ]
    if not selected or (stations and set(stations) != {r["station"] for r in selected}):
        raise ValueError(
            "Select reportable stations explicitly, or omit --station for all"
        )
    weeks = []
    station_rows = []
    for row in selected:
        reproduced = frequency_station(row["weeks"])
        saved = {
            k: v
            for k, v in row["frequency_consistency"].items()
            if k != "independent_weeks"
        }
        if reproduced != saved:
            raise ValueError("Frozen station selection no longer reproduces")
        intervals = sorted(
            w["interval_utc"] for w in row["weeks"] if w.get("interval_utc")
        )
        if any(
            UTCDateTime(a[1]) > UTCDateTime(b[0])
            for a, b in zip(intervals, intervals[1:])
        ):
            raise ValueError("Source weeks overlap")
        audited = [
            audit_week(row, w, guard)
            for w in row["weeks"]
            if w.get("frequency_check", {}).get("passes")
        ]
        weeks.extend(audited)
        bands = [b for w in audited for b in w["bands"]]
        daily = [d for b in bands for d in b["daily"]]
        reference = row.get("rayleigh_baseline_g")
        station_rows.append(
            {
                "station": row["station"],
                "passing_weeks": len(audited),
                "audited_bands": len(bands),
                "consistent_bands": sum(
                    b["status"] == "consistent_with_prescribed" for b in bands
                ),
                "worst_weekly_rms_deg": max(
                    b["prescribed_rms_deg"]
                    for b in bands
                    if b["prescribed_rms_deg"] is not None
                ),
                "daily_band_count": len(daily),
                "consistent_daily_bands": sum(
                    d["status"] == "consistent_with_prescribed" for d in daily
                ),
                "daily_status_counts": dict(Counter(d["status"] for d in daily)),
                "frozen_conditional_g": row["frequency_consistency"]["median_g"],
                "fixed_rayleigh_g": reference,
                "difference_pct": 100
                * (row["frequency_consistency"]["median_g"] / reference - 1)
                if reference
                else None,
                "median_dry_vs_water_h1_pct": statistics.median(
                    b["model_context"]["dry_vs_water_h1_pct"] for b in bands
                ),
                "maximum_dry_vs_water_h1_pct": max(
                    b["model_context"]["dry_vs_water_h1_pct"] for b in bands
                ),
                "reference_inside_h1_h2_bands": sum(
                    b["model_context"].get("reference_inside_model_h1_h2", False)
                    for b in bands
                )
                if reference
                else None,
                "tilt_status": "not_assessed_requires_p_h1_h2_z",
                "absolute_vertical_response_status": "not_independently_verified",
                "scientific_score": None,
            }
        )
    staging_rows = table_export(review, station_rows, guard)
    guard.verify()
    output = (
        Path(output_root).resolve()
        / f"run_{datetime.now(UTC):%Y%m%dT%H%M%SZ}_{uuid4().hex[:8]}"
    )
    # Never allow an output root inside any protected file's parent directory.
    if any(output.is_relative_to(Path(p).parent) for p in guard.hashes):
        raise ValueError("Output would be written beside protected inputs")
    output.mkdir(parents=True, exist_ok=False)
    result = {
        "audit_id": AUDIT_ID,
        "source_review": str(Path(review_path).resolve()),
        "frozen_rule_id": RULE_ID,
        "selection_policy": (
            "reportable stations; passing weeks; exact eligible bands"
        ),
        "thresholds": {
            "coherence_squared": COHERENT,
            "minimum_bins": MINIMUM_BINS,
            "fits_deg": FITS_DEG,
            "unexplained_deg": UNEXPLAINED_DEG,
            "role": "legacy descriptive cutoffs, not a new calibration rule",
        },
        "rows": station_rows,
        "weeks": weeks,
        "input_sha256": guard.hashes,
        "protected_inputs_unchanged": True,
        "gain_applied": False,
        "master_or_scores_changed": False,
        "status": (
            "diagnostic_only; frozen selections unchanged; "
            "tilt/model/absolute response remain open"
        ),
    }
    fields = [
        {k: v for k, v in r.items() if k != "daily_status_counts"} for r in station_rows
    ]
    write_csv(output / "stations.csv", fields)
    write_csv(output / "table_staging.csv", staging_rows)
    atomic_json(output / "audit.json", result)
    lines = [
        "# Saved gain-band phase audit",
        "",
        "Descriptive only; frozen selections, gains, masters and scores unchanged.",
        "",
        "| Station | Weeks | Phase-consistent bands | Worst RMS (deg) "
        "| Consistent daily bands | g vs Rayleigh |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in station_rows:
        delta = "—" if r["difference_pct"] is None else f"{r['difference_pct']:+.2f}%"
        lines.append(
            f"| {r['station']} | {r['passing_weeks']} "
            f"| {r['consistent_bands']}/{r['audited_bands']} "
            f"| {r['worst_weekly_rms_deg']:.2f} "
            f"| {r['consistent_daily_bands']}/{r['daily_band_count']} | {delta} |"
        )
    lines += [
        "",
        "Phase consistency cannot exclude coherent tilt "
        "or validate absolute vertical response.",
        "H1–H2 are descriptive noise envelopes, not confidence intervals.",
        "Cached water/dry model sensitivities "
        "are not an Earth-model uncertainty budget.",
        "",
    ]
    (output / "audit.md").write_text("\n".join(lines))
    guard.verify()
    if docs:
        docs = Path(docs)
        if not docs.is_dir():
            raise ValueError("--docs must be an existing directory")
        compact = {
            k: v for k, v in result.items() if k not in {"weeks", "input_sha256"}
        }
        compact |= {
            "full_audit_path": str(output / "audit.json"),
            "full_audit_sha256": InputGuard.digest(output / "audit.json"),
            "source_review_sha256": InputGuard.digest(review_path),
            "protected_file_count": len(guard.hashes),
            "audited_weeks": len(weeks),
            "weekly_band_status_counts": dict(
                Counter(b["status"] for w in weeks for b in w["bands"])
            ),
            "daily_band_status_counts": dict(
                Counter(
                    d["status"] for w in weeks for b in w["bands"] for d in b["daily"]
                )
            ),
            "staging_status_counts": dict(
                Counter(
                    r["calibration_check_pressure_vertical_microseisms_status"]
                    for r in staging_rows
                )
            ),
            "weekly_phase_flags": [
                {
                    "station": w["station"],
                    "shift_weeks": w["shift_weeks"],
                    **{
                        k: b[k]
                        for k in ("frequency_hz", "status", "prescribed_rms_deg")
                    },
                }
                for w in weeks
                for b in w["bands"]
                if b["status"] != "consistent_with_prescribed"
            ],
        }
        atomic_json(docs / "microseism_gain_band_audit.json", compact)
        write_csv(docs / "microseism_table_staging.csv", staging_rows)
        guard.verify()
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument(
        "--station", action="append", help="Reportable NET.STA; default all"
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "results/REVIEWS/MICROSEISM_GAIN_BANDS",
    )
    parser.add_argument(
        "--docs", type=Path, help="Write curated evidence/staging copies"
    )
    args = parser.parse_args()
    print(report(args.review, args.output_root, args.station, args.docs))


if __name__ == "__main__":
    main()
