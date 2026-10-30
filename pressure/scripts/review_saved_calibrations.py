"""Read-only integration checks against four saved network runs and bounded raw windows.

Writes a new review folder only. Does not run a catalogue batch, change responses,
apply gains or write masters. Import from the checkout being reviewed.
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from obspy import UTCDateTime

from compy_pressure.archive import discover_waveforms
from compy_pressure.inventory import load_master_dataset
from compy_pressure.paired_spectra import response_snapshot_path
from compy_pressure.rayleigh_calibration import (
    event_spectra,
    fit_multiplier,
    load_aligned,
    select_calibration_channels,
    water_ratio,
)

RUNS = {
    "YV": "YV_RAYLEIGH_NO_TILT",
    "YL": "YL_ALL_NO_TILT",
    "7D": "7D_ALL_NO_TILT",
    "XO": "XO_ALL_GAUGES_NO_TILT",
}
PROBES = {
    "YV": ["RR03", "RR52", "RR55"],
    "YL": ["A02W", "C07W", "B01W", "C08W"],
    "7D": ["J30A", "J38A", "J45A", "J54A"],
    "XO": ["LA33", "LA34", "WD66", "WD70", "LT02", "WD48"],
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_saved(root, protected):
    networks, records = {}, {}
    for network, folder in RUNS.items():
        run = root / "results/NETWORKS" / folder
        workflow = run / "_workflow/rayleigh.json"
        protected[workflow] = digest(workflow)
        doc = json.loads(workflow.read_text())
        assert doc["execution_status"] == "complete"
        counts, gains, events, error = Counter(), 0, 0, 0.0
        for entry in doc["entries"].values():
            path = Path(entry["update"]["rayleigh_result_file"])
            protected[path] = digest(path)
            result = json.loads(path.read_text())
            counts[result["status"]] += 1
            station = entry["seed_id"].split(".")[1]
            records[network, station] = (run, doc, result)
            fitted = []
            for event in result["events"]:
                if event["status"] != "candidate":
                    continue
                path = Path(event["files"]["spectra"])
                protected[path] = digest(path)
                with np.load(path) as saved:
                    valid = saved["accepted_bins"].astype(bool)
                    fit, _, _ = fit_multiplier(
                        saved["psd_ratio_pa"][valid], saved["expected_pa"][valid]
                    )
                difference = abs(
                    fit["pressure_multiplier"] - event["pressure_multiplier"]
                )
                assert difference < 1e-10, (network, station, difference)
                error = max(error, difference)
                fitted.append(fit["pressure_multiplier"])
                events += 1
            assert len(fitted) == result["candidate_events"]
            if result.get("median_pressure_multiplier") is not None:
                g = float(np.median(fitted))
                assert abs(g - result["median_pressure_multiplier"]) < 1e-10
                assert np.isclose(
                    g,
                    entry["update"]["rayleigh_pressure_multiplier"],
                    rtol=0,
                    atol=1e-10,
                )
                assert np.isclose(
                    1 / g,
                    entry["update"]["rayleigh_sensitivity_multiplier"],
                    rtol=0,
                    atol=1e-10,
                )
                gains += 1
        networks[network] = dict(
            stations=len(doc["entries"]),
            gains=gains,
            accepted_events=events,
            max_saved_gain_error=error,
            status_counts=dict(counts),
        )
        print(network, networks[network], flush=True)
    return networks, records


def check_raw(root, data_root, records, protected):
    probes, source_stats = [], {}
    for network, stations in PROBES.items():
        for station in stations:
            run, workflow, result = records[network, station]
            master_dir = run / "_master"
            if not master_dir.is_dir():
                master_dir = root / "exports"
            for path in master_dir.glob("master_inventory.*"):
                protected[path] = digest(path)
            master = load_master_dataset(master_dir)
            settings = result["settings"]
            rate = settings.get("event_resample_rate_hz")
            rows, channels = select_calibration_channels(
                master, network, station, result["epoch_keys"][0], "none", rate
            )
            for row in rows:
                for path in (Path(row.inventory_file), response_snapshot_path(row)):
                    protected[path] = digest(path)
            candidates = [e for e in result["events"] if e["status"] == "candidate"]
            if candidates:
                chosen = [candidates[0]]
                if station == "C07W":
                    chosen = [
                        min(candidates, key=lambda e: e["pressure_multiplier"]),
                        max(candidates, key=lambda e: e["pressure_multiplier"]),
                    ]
            else:
                chosen = [next(e for e in result["events"] if "start_utc" in e)]
            raw = run / "_inputs/waveforms"
            if not raw.is_dir():
                raw = data_root / network
                if (raw / station).is_dir():
                    raw /= station
            archive = discover_waveforms(raw, {row.Seed_id for row in rows})
            source = (
                result.get("dispersion_model", {}).get("source_file")
                or workflow.get("dispersion_model", {}).get("source_file")
                or workflow["options"]["dispersion_file"]
            )
            model_path = Path(source)
            if not model_path.is_absolute():
                model_path = root / model_path
            model = np.loadtxt(model_path, delimiter=",", skiprows=1)
            protected[model_path] = digest(model_path)
            for event in chosen:
                start, end = (
                    UTCDateTime(event["start_utc"]),
                    UTCDateTime(event["end_utc"]),
                )
                padding = 32 / rate if rate else 0
                for record in archive:
                    if (
                        record.start < end + padding + 1
                        and record.end >= start - padding
                    ):
                        stat = record.path.stat()
                        source_stats[record.path] = [stat.st_size, stat.st_mtime_ns]
                report = dict(
                    network=network,
                    station=station,
                    event=event["id"],
                    saved_status=event["status"],
                )
                try:
                    arrays, _, actual_rate = load_aligned(
                        archive,
                        rows,
                        channels,
                        start,
                        float(end - start),
                        rate,
                        settings.get("native_rate_tolerance_ppm", 0),
                    )
                    spectra = event_spectra(
                        arrays,
                        rows,
                        channels,
                        actual_rate,
                        np.array([1.0]),
                        200,
                        tuple(settings["band_hz"]),
                        settings["pressure_sign"],
                    )
                    f = spectra["frequency_hz"]
                    expected = water_ratio(
                        f,
                        np.interp(f, model[:, 0], model[:, 1]),
                        settings["depth_m"],
                        settings["density_kg_m3"],
                        settings["sound_speed_m_s"],
                    )
                    coherence = spectra["coherence_squared"]
                    coherent = coherence >= settings["coherence_squared_threshold"]
                    valid = coherent & (
                        abs(spectra["phase_deg"]) <= settings["phase_limit_deg"]
                    )
                    report.update(
                        median_coherence=float(np.median(coherence)),
                        phase_valid_bins=int(valid.sum()),
                    )
                    if event["status"] == "candidate":
                        fit, _, _ = fit_multiplier(
                            spectra["psd_ratio_pa"][valid], expected[valid]
                        )
                        delta = abs(
                            fit["pressure_multiplier"] - event["pressure_multiplier"]
                        )
                        assert delta < 1e-8, (network, station, delta)
                        report.update(
                            outcome="candidate_reproduced",
                            gain=fit["pressure_multiplier"],
                            gain_error=delta,
                        )
                    else:
                        report["outcome"] = "rejection_diagnostic"
                        report["fit_reason"] = None
                        try:
                            fit_multiplier(
                                spectra["psd_ratio_pa"][coherent], expected[coherent]
                            )
                        except ValueError as exc:
                            report["fit_reason"] = str(exc)
                except ValueError as exc:
                    if event["status"] == "candidate":
                        raise
                    report.update(outcome="input_rejection", reason=str(exc))
                probes.append(report)
                print(network, station, report["outcome"], flush=True)
    for path, expected in source_stats.items():
        stat = path.stat()
        assert [stat.st_size, stat.st_mtime_ns] == expected, str(path)
    return probes, len(source_stats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    protected = {}
    networks, records = check_saved(args.source_root, protected)
    probes, raw_count = check_raw(args.source_root, args.data_root, records, protected)
    for path, expected in protected.items():
        assert digest(path) == expected, str(path)
    report = dict(
        networks=networks,
        raw_probes=probes,
        protected_hashed_files=len(protected),
        raw_files_stat_checked=raw_count,
        sources_unchanged=True,
        scope="saved-fit replay and bounded raw windows; no batch or applied gain",
    )
    (args.output_dir / "integration_validation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print("Sources unchanged; report:", args.output_dir / "integration_validation.json")


if __name__ == "__main__":
    main()
