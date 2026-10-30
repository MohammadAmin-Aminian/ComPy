"""Validate frozen YV trial gains on unused earthquakes; never refit or apply gains."""

import argparse
import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import matplotlib
import numpy as np
from obspy import UTCDateTime
from obspy.taup import TauPyModel

from compy_pressure.archive import discover_waveforms
from compy_pressure.calibration_experiment import pole_factor, real_pressure_poles
from compy_pressure.heldout_validation import (
    frozen_gain_check,
    select_unused_events,
    summarize_errors,
)
from compy_pressure.inventory import load_master_dataset
from compy_pressure.paired_spectra import response_snapshot_path
from compy_pressure.rayleigh_calibration import (
    catalog_events,
    event_geometry,
    event_spectra,
    load_aligned,
    select_calibration_channels,
    water_ratio,
)
from scripts.run_calibration_experiment import Guard

ROOT = Path(__file__).resolve().parents[1]
TRAINING = (
    ROOT / "results/EXPERIMENTS/SOCRATES/all_20261013T195356Z_00cf3343/experiment.json"
)
STATIONS = ["RR01", "RR06", "RR48", "RR52"]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare(args, guard):
    training = guard.json(args.training_experiment)
    if training["status"] != "complete":
        raise ValueError("Training experiment must be complete")
    directory = args.training_experiment.resolve().parent
    reviewed, trials = [], {}
    for entry in training["stations"]:
        if entry["network"] == "YV":
            trial = guard.json(directory / f"YV_{entry['station']}/trial.json")
            reviewed.extend(trial["events"])
            trials[entry["station"]] = (entry, trial)
    development_file = ROOT / "docs/yv_heldout_development_events.json"
    if development_file.is_file():
        reviewed.extend(guard.json(development_file)["events"])
    workflow = guard.json(
        ROOT / "results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json"
    )
    catalogue_path = Path(workflow["options"]["catalog"])
    if not catalogue_path.is_absolute():
        catalogue_path = ROOT / catalogue_path
    catalogue = catalog_events(guard.protect(catalogue_path))
    heldout, excluded = select_unused_events(
        catalogue, reviewed, args.minimum_magnitude
    )
    master_dir = ROOT / "exports"
    for path in master_dir.glob("master_inventory.*"):
        actual = guard.protect(path)
        expected = training["source_sha256"].get(str(actual))
        if expected is None or guard.hashes[actual] != expected:
            raise ValueError(
                f"Master differs from the frozen training snapshot: {path}"
            )
    master = load_master_dataset(master_dir)
    configurations = []
    for station in args.stations:
        entry, trial = trials[station]
        gain = entry["summaries"]["pole_absolute"]["median_pressure_multiplier"]
        if gain is None or gain <= 0:
            raise ValueError(f"No frozen pole-only training gain at {station}")
        fitted = [
            e["variants"]["pole_absolute"]["pressure_multiplier"]
            for e in trial["events"]
            if e.get("variants", {}).get("pole_absolute", {}).get("status")
            == "candidate"
        ]
        if len(fitted) < 3 or not np.isclose(
            np.median(fitted), gain, rtol=0, atol=1e-12
        ):
            raise ValueError(f"Training gain/count provenance mismatch at {station}")
        original = guard.json(trial["baseline_result"])
        if original["settings"]["tilt_method"] != "none":
            raise ValueError(
                "Held-out runner requires the same no-tilt training method"
            )
        for key in ("source_inventory_sha256", "response_snapshot_sha256"):
            for path, expected in original[key].items():
                if guard.hashes[guard.protect(path)] != expected:
                    raise ValueError(f"Response differs from training: {path}")
        rows, channels = select_calibration_channels(
            master, "YV", station, original["epoch_keys"][0], "none"
        )
        for row in rows:
            guard.protect(row.inventory_file)
            guard.protect(response_snapshot_path(row))
        source = (
            original.get("dispersion_model", {}).get("source_file")
            or workflow["options"]["dispersion_file"]
        )
        model_path = Path(source)
        if not model_path.is_absolute():
            model_path = ROOT / model_path
        actual = guard.protect(model_path)
        expected = training["source_sha256"].get(str(actual))
        if expected is None or guard.hashes[actual] != expected:
            raise ValueError(f"Model differs from frozen training: {model_path}")
        model = np.loadtxt(actual, delimiter=",", skiprows=1)
        sensor = channels[0].sensor
        poles = real_pressure_poles(channels[0].response)
        audit = {
            "station": station,
            "pressure_seed_id": rows[0].Seed_id,
            "sensor": {
                key: getattr(sensor, key, None)
                for key in (
                    "type",
                    "description",
                    "manufacturer",
                    "model",
                    "serial_number",
                )
            },
            "pressure_input_units": channels[
                0
            ].response.instrument_sensitivity.input_units,
            "positive_real_pressure_poles_rad_s": poles,
            "alternative_pressure_poles_rad_s": [-p for p in poles],
            "nominal_responses_modified": False,
            "physical_hardware_identity_status": "physical identity unverified",
            "frozen_gain": float(gain),
            "training_candidate_events": len(fitted),
            "training_variant": "pole_absolute",
            "training_settings": original["settings"],
        }
        configurations.append((audit, rows, channels, model))
    return (
        configurations,
        catalogue,
        heldout,
        excluded,
        sorted({e["id"] for e in reviewed}),
    )


def plot_event(path, spectra, check, reflected, gain, station, event):
    import matplotlib.pyplot as plt

    f = spectra["frequency_hz"]
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    axes[0].plot(
        f, 100 * (gain * spectra["psd_ratio_pa"] / spectra["expected_pa"] - 1), "o-"
    )
    axes[0].axhline(0, color="k", lw=0.7)
    axes[0].set_ylabel("Frozen-gain error (%)")
    axes[1].plot(f, spectra["phase_deg"], "o-", label="nominal")
    axes[1].plot(f, reflected, "o-", label="pole alternative")
    axes[1].axhspan(-20, 20, color="green", alpha=0.1)
    axes[1].set_ylabel("P/a phase (degrees)")
    axes[1].legend()
    axes[2].plot(f, spectra["coherence_squared"], "o-")
    axes[2].set_ylabel("Squared coherence")
    axes[2].set_xlabel("Frequency (Hz)")
    axes[2].set_ylim(0, 1.02)
    for axis in axes:
        axis.grid(alpha=0.3)
    fig.suptitle(
        f"{station} {event['id']} Mw {event['magnitude']}\n"
        f"Frozen g={gain:.6f}; {check['status']}"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def validate_station(configuration, catalogue, heldout, archive, taup, output):
    audit, rows, channels, model = configuration
    station, gain, settings = (
        audit["station"],
        audit["frozen_gain"],
        audit["training_settings"],
    )
    records = [
        record for record in archive if record.seed_id in {r.Seed_id for r in rows}
    ]
    output.mkdir()
    write_json(output / "response_audit.json", audit)
    events = []
    latitude, longitude = float(rows[0].Lat), float(rows[0].Lon)
    packet_windows = {
        e["id"]: event_geometry(e, latitude, longitude, model)
        for e in catalogue
        if e["magnitude"] >= 6
    }
    for index, event in enumerate(heldout, 1):
        print(
            f"{station} unused event {index}/{len(heldout)}: "
            f"{event['id']} Mw {event['magnitude']}",
            flush=True,
        )
        item = dict(event, status="input_review", frozen_gain=gain, gain_refitted=False)
        try:
            start, end, distance = packet_windows[event["id"]]
            origin = UTCDateTime(event["origin_utc"])
            if start < origin or end - origin > 43200 or end - start > 7200:
                raise ValueError(
                    "Predicted window fails original origin/duration guard"
                )
            other = [
                eid
                for eid, (left, right, _) in packet_windows.items()
                if eid != event["id"] and left < end and right > start
            ]
            arrivals = taup.get_travel_times(
                event["depth_km"],
                distance / 111.195,
                phase_list=["P", "Pdiff", "PKP", "S", "Sdiff", "SKS"],
            )
            body = [
                {"phase": a.name, "time_utc": str(origin + a.time)}
                for a in arrivals
                if start <= origin + a.time <= end
            ]
            item.update(
                start_utc=str(start),
                end_utc=str(end),
                distance_km=distance,
                other_packet_ids=other,
                body_arrivals_in_window=body,
                window_review="predicted; visual mode isolation remains necessary",
            )
            arrays, actual, rate = load_aligned(
                records, rows, channels, start, float(end - start)
            )
            spectra = event_spectra(
                arrays,
                rows,
                channels,
                rate,
                np.array([1.0]),
                200,
                tuple(settings["band_hz"]),
                settings["pressure_sign"],
            )
            f = spectra["frequency_hz"]
            spectra["expected_pa"] = water_ratio(
                f,
                np.interp(f, model[:, 0], model[:, 1]),
                settings["depth_m"],
                settings["density_kg_m3"],
                settings["sound_speed_m_s"],
            )
            factor = pole_factor(f, audit["positive_real_pressure_poles_rad_s"])
            nominal = frozen_gain_check(spectra, settings, gain, np.ones_like(factor))
            reflected = frozen_gain_check(spectra, settings, gain, factor)
            item.update(
                actual_start_utc=str(actual),
                nominal=nominal,
                reflected=reflected,
                status="isolated_heldout_usable"
                if reflected["status"] == "qc_usable" and not other and not body
                else "arrival_isolation_review"
                if reflected["status"] == "qc_usable"
                else reflected["status"],
            )
            np.savez_compressed(
                output / f"{event['id']}_spectra.npz",
                **spectra,
                reflected_phase_deg=np.angle(spectra["cross_pa"] * factor, deg=True),
            )
            plot_event(
                output / f"{event['id']}_review.png",
                spectra,
                reflected,
                np.angle(spectra["cross_pa"] * factor, deg=True),
                gain,
                station,
                event,
            )
        except ValueError as exc:
            item["reason"] = str(exc)
        events.append(item)
        write_json(output / f"{event['id']}.json", item)
    report = {
        "station": station,
        "frozen_gain": gain,
        "response_audit": audit,
        "primary_isolated_summary": summarize_errors(events),
        "qc_only_including_body_arrival_review": summarize_errors(
            events, require_isolated=False
        ),
        "events": events,
    }
    write_json(output / "validation.json", report)
    return report


def plot_summary(reports, output):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        len(reports), 1, figsize=(11, 3 * len(reports)), squeeze=False
    )
    for axis, report in zip(axes[:, 0], reports, strict=True):
        for isolated, label in [
            (True, "isolated predicted window"),
            (False, "body-arrival review"),
        ]:
            selected = [
                e
                for e in report["events"]
                if e.get("reflected", {}).get("status") == "qc_usable"
                and not e.get("other_packet_ids")
                and (not bool(e.get("body_arrivals_in_window"))) == isolated
            ]
            axis.scatter(
                [UTCDateTime(e["origin_utc"]).datetime for e in selected],
                [e["reflected"]["median_error_percent"] for e in selected],
                label=label,
            )
        axis.axhline(0, color="k", lw=0.7)
        axis.grid(alpha=0.3)
        axis.set_ylabel("Prediction error (%)")
        axis.set_title(
            f"{report['station']}; frozen g={report['frozen_gain']:.6f}; "
            f"isolated usable={report['primary_isolated_summary']['usable_events']}"
        )
        axis.legend(fontsize=8)
    fig.autofmt_xdate()
    fig.tight_layout()
    path = output / "heldout_summary.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-experiment", type=Path, default=TRAINING)
    parser.add_argument("--stations", nargs="+", default=STATIONS)
    parser.add_argument("--minimum-magnitude", type=float, default=6.0)
    parser.add_argument(
        "--data-dir", type=Path, default=Path.home() / "Desktop/Data_Bruit_FM/Data/YV"
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--open-plots", action="store_true")
    args = parser.parse_args()
    guard = Guard()
    configs, catalogue, heldout, excluded, reviewed_ids = prepare(args, guard)
    for audit, _, _, _ in configs:
        print(
            f"{audit['station']}: frozen g={audit['frozen_gain']:.8f}, "
            f"positive poles={audit['positive_real_pressure_poles_rad_s']}",
            flush=True,
        )
    print(
        f"Recorded YV events excluded: {len(reviewed_ids)}; "
        f"unused mainshocks: {len(heldout)}",
        flush=True,
    )
    if args.dry_run:
        guard.verify()
        return
    if not heldout:
        raise ValueError("No independent eligible events; do not reuse training events")
    matplotlib.use("Agg")
    output = (
        args.output_dir
        or ROOT
        / "results/VALIDATION/YV_HELDOUT"
        / f"run_{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}_{uuid4().hex[:8]}"
    ).resolve()
    output.mkdir(parents=True, exist_ok=False)
    for path in [
        Path(__file__),
        ROOT / "../compy_pressure/heldout_validation.py",
        ROOT / "../compy_pressure/calibration_experiment.py",
        ROOT / "../compy_pressure/rayleigh_calibration.py",
        ROOT / "../compy_pressure/paired_spectra.py",
        ROOT / "../compy_pressure/archive.py",
        ROOT / "../compy_pressure/inventory.py",
        ROOT / "scripts/run_calibration_experiment.py",
    ]:
        guard.protect(path)
    write_json(
        output / "selection.json",
        {
            "reviewed_YV_event_ids": reviewed_ids,
            "minimum_magnitude": args.minimum_magnitude,
            "separation_hours": 24,
            "events": heldout,
            "excluded": excluded,
        },
    )
    archive = discover_waveforms(
        args.data_dir, {r.Seed_id for _, rows, _, _ in configs for r in rows}
    )
    raw_stats = {
        record.path: (record.path.stat().st_size, record.path.stat().st_mtime_ns)
        for record in archive
    }
    taup = TauPyModel(model="iasp91")
    reports = [
        validate_station(
            configuration,
            catalogue,
            heldout,
            archive,
            taup,
            output / configuration[0]["station"],
        )
        for configuration in configs
    ]
    event_rows = [
        {
            "station": r["station"],
            "event_id": e["id"],
            "magnitude": e["magnitude"],
            "status": e["status"],
            "frozen_gain": r["frozen_gain"],
            "gain_refitted": False,
            "phase_valid_bins": e.get("reflected", {}).get("phase_valid_bins"),
            "median_error_percent": e.get("reflected", {}).get("median_error_percent"),
            "other_packets": len(e.get("other_packet_ids", [])),
            "body_arrivals_in_window": len(e.get("body_arrivals_in_window", [])),
            "reason": e.get("reason", ""),
        }
        for r in reports
        for e in r["events"]
    ]
    write_csv(output / "heldout_events.csv", event_rows)
    summary = [
        {
            "station": r["station"],
            "frozen_gain": r["frozen_gain"],
            "usable_isolated_events": r["primary_isolated_summary"]["usable_events"],
            "median_error_percent": r["primary_isolated_summary"].get(
                "median_error_percent"
            ),
            "status": r["primary_isolated_summary"]["status"],
        }
        for r in reports
    ]
    write_csv(output / "heldout_summary.csv", summary)
    plot = plot_summary(reports, output)
    guard.verify()
    for path, expected in raw_stats.items():
        stat = path.stat()
        if (stat.st_size, stat.st_mtime_ns) != expected:
            raise RuntimeError(f"Raw file changed during validation: {path}")
    write_json(
        output / "validation.json",
        {
            "status": "complete",
            "gain_refitted": False,
            "gains_applied_to_source_products": False,
            "source_training": str(args.training_experiment.resolve()),
            "source_sha256": {str(p): h for p, h in guard.hashes.items()},
            "raw_file_stats": {str(p): list(stat) for p, stat in raw_stats.items()},
            "inputs_unchanged": True,
            "selected_unused_events": len(heldout),
            "summary": summary,
            "limitations": [
                "Unused events share the training estimator, responses and model",
                "QC/arrival screening does not prove fundamental-mode isolation",
                "Hardware identity/serial calibration remains pending; no XML changed",
                "Event scatter/bootstrap omits absolute calibration systematics",
                "No glitch removal, tilt correction, noise subtraction or score update",
            ],
        },
    )
    print(
        f"\nCOMPLETE. Frozen gains unchanged. Results: {output}\nPlot: {plot}",
        flush=True,
    )
    if args.open_plots:
        subprocess.run(["open", str(plot)], check=True)


if __name__ == "__main__":
    main()
