"""Selected-station Socrates trials; preserve baseline data, responses and results.

Use every saved spectral event at the selected stations, plus recoverable legacy
boundary events. Keep input/timing/overlap guards and baseline event windows.
No new catalogue selection, PPSD batch, response/master write or gain application.
"""

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import matplotlib
import numpy as np
from obspy import UTCDateTime, read
from scipy.signal import butter, sosfiltfilt

from compy_pressure.archive import discover_waveforms
from compy_pressure.calibration_experiment import (
    evaluate,
    highpass_diagnostic,
    pole_factor,
    real_pressure_poles,
)
from compy_pressure.inventory import load_master_dataset
from compy_pressure.paired_spectra import checked_response, response_snapshot_path
from compy_pressure.rayleigh_calibration import (
    event_spectra,
    load_aligned,
    select_calibration_channels,
    water_ratio,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULTS = {"YV": ["RR03", "RR22", "RR52"], "YL": ["A02W", "B01W", "C07W", "C08W"]}
RUNS = {
    "YV": "YV_RAYLEIGH_NO_TILT",
    "YL": "YL_ALL_NO_TILT",
    "7D": "7D_ALL_NO_TILT",
    "XO": "XO_ALL_GAUGES_NO_TILT",
}
# Pressure multipliers from Aminian (2025), Table 1; Zha (2016) coefficients=1/g.
YV_PUBLISHED = {
    "RR28": 0.78,
    "RR29": 0.70,
    "RR34": 0.58,
    "RR36": 0.64,
    "RR38": 0.60,
    "RR40": 0.74,
    "RR50": 0.68,
    "RR52": 0.67,
}
YL_PUBLISHED = {
    "A02W": 1.50,
    "A06W": 1.53,
    "A08W": 1.43,
    "A10W": 1.82,
    "A12W": 1.67,
    "A14W": 1.97,
    "B01W": 0.19,
    "B04W": 1.41,
    "B09W": 1.32,
    "B11W": 1.90,
    "B14W": 1.62,
    "C01W": 1.85,
    "C03W": 1.55,
    "C05W": 1.77,
    "C07W": 1.99,
    "C08W": 1.78,
    "C09W": 1.84,
    "C10W": 1.51,
    "C11W": 1.83,
    "C16W": 1.79,
    "C17W": 1.43,
    "F02W": 1.47,
    "N02W": 1.73,
    "N03W": 1.80,
    "S01W": 1.81,
}
VARIANTS = {
    "nominal_absolute": (False, "absolute"),
    "pole_absolute": (True, "absolute"),
    "nominal_weighted_log": (False, "weighted_log"),
    "pole_weighted_log": (True, "weighted_log"),
}


class Guard:
    def __init__(self):
        self.hashes = {}

    def protect(self, path):
        path = Path(path).resolve()
        if path not in self.hashes:
            self.hashes[path] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    def json(self, path):
        return json.loads(self.protect(path).read_text())

    def verify(self):
        for path, expected in self.hashes.items():
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise RuntimeError(f"Protected input changed during trial: {path}")


def microseism_rms(data_root, station, start, end, guard):
    """Raw-count event-window diagnostic, not a waveform gain correction."""
    stream = None
    for day in sorted({start.strftime("%Y%m%d"), (end - 1).strftime("%Y%m%d")}):
        paths = sorted(
            (data_root / "YL" / station / "waveforms").glob(
                f"YL.{station}.*.LHZ.{day}.mseed"
            )
        )
        for path in paths:
            part = read(str(guard.protect(path)))
            stream = part if stream is None else stream + part
    if stream is None:
        raise ValueError(f"Missing microseism window at {station}")
    stream.merge(method=0, fill_value=None)
    stream.trim(start, end - 0.5, nearest_sample=True)
    if (
        len(stream) != 1
        or stream[0].stats.sampling_rate != 1
        or np.ma.is_masked(stream[0].data)
        or abs(len(stream[0].data) - round(end - start)) > 1
    ):
        raise ValueError(f"Incomplete 1 Hz microseism window at {station}")
    x = stream[0].data.astype(float)
    if len(x) < 1200 or not np.isfinite(x).all():
        raise ValueError("Insufficient finite microseism samples")
    sos = butter(4, (0.12, 0.3), fs=1.0, btype="band", output="sos")
    x = sosfiltfilt(sos, x - x.mean())[200:-200]
    rms = float(
        np.median(
            [np.sqrt(np.mean(x[k : k + 300] ** 2)) for k in range(0, len(x) - 299, 300)]
        )
    )
    if rms <= 0:
        raise ValueError("Zero microseism power")
    return rms


def c07_state(data_root, event, guard):
    try:
        start, end = UTCDateTime(event["start_utc"]), UTCDateTime(event["end_utc"])
        values = {
            st: microseism_rms(data_root, st, start, end, guard)
            for st in ("C07W", "C08W", "C09W")
        }
        ratio = values["C07W"] / values["C08W"]
        control = values["C08W"] / values["C09W"]
        state = "unclassified"
        if 0.65 < control < 1.05:
            if 0.65 < ratio < 1.2:
                state = "normal"
            elif 0.3 < ratio < 0.6:
                state = "reduced"
        return {
            "state": state,
            "event_window_c07_over_c08": ratio,
            "event_window_c08_over_c09": control,
            "diagnostic_only": True,
        }
    except ValueError as exc:
        return {"state": "unclassified", "reason": str(exc), "diagnostic_only": True}


def get_spectra(event, rows, channels, archive, dispersion, settings, output, guard):
    path = event.get("files", {}).get("spectra")
    if path:
        with np.load(guard.protect(path), allow_pickle=False) as saved:
            return {k: saved[k].copy() for k in saved.files}
    # SENT-002: older boundary failures did not save spectra. Recover only these;
    # never silently retry events rejected by continuity/timing/selection guards.
    if "boundary" not in (event.get("reason") or ""):
        return None
    start, end = UTCDateTime(event["start_utc"]), UTCDateTime(event["end_utc"])
    rate = settings.get("event_resample_rate_hz")
    padding = 32 / rate if rate else 0
    for record in archive:
        if record.start < end + padding + 1 and record.end >= start - padding:
            guard.protect(record.path)
    arrays, _, rate = load_aligned(
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
        rate,
        np.array([1.0]),
        200,
        tuple(settings["band_hz"]),
        settings["pressure_sign"],
    )
    f = spectra["frequency_hz"]
    spectra["expected_pa"] = water_ratio(
        f,
        np.interp(f, dispersion[:, 0], dispersion[:, 1]),
        settings["depth_m"],
        settings["density_kg_m3"],
        settings["sound_speed_m_s"],
    )
    # New artifact only; never replace an original NPZ.
    np.savez_compressed(output / f"{event['id']}_recovered_spectra.npz", **spectra)
    return spectra


def station_trial(
    network, station, entry, workflow, run, data_root, output, guard, master_cache=None
):
    result = guard.json(entry["update"]["rayleigh_result_file"])
    settings = result["settings"]
    output.mkdir()
    if result["status"] == "reference_band_unusable":
        report = {
            "network": network,
            "station": station,
            "saved_station_gain": None,
            "settings": settings,
            "baseline_result": entry["update"]["rayleigh_result_file"],
            "baseline_reason": result.get("reason"),
            "positive_real_pressure_poles_rad_s": [],
            "events": [],
            "summaries": {
                name: {
                    "candidate_events": 0,
                    "median_pressure_multiplier": None,
                    "diagnostic_pooled_median": None,
                    "status": "reference_band_unusable",
                }
                for name in VARIANTS
            },
        }
        (output / "trial.json").write_text(json.dumps(report, indent=2) + "\n")
        return report
    if settings["tilt_method"] != "none":
        raise ValueError(
            "This controlled trial requires the preserved no-tilt baseline"
        )
    master_dir = run / "_master" if (run / "_master").is_dir() else ROOT / "exports"
    for path in master_dir.glob("master_inventory.*"):
        guard.protect(path)
    if master_cache is None:
        master_cache = {}
    if master_dir not in master_cache:
        master_cache[master_dir] = load_master_dataset(master_dir)
    master = master_cache[master_dir]
    rows, channels = select_calibration_channels(
        master,
        network,
        station,
        result["epoch_keys"][0],
        "none",
        settings.get("event_resample_rate_hz"),
    )
    for row in rows:
        guard.protect(row.inventory_file)
        guard.protect(response_snapshot_path(row))
    # Prevent using stale cached spectra with a different input response.
    for key in ("source_inventory_sha256", "response_snapshot_sha256"):
        for path, expected in result[key].items():
            actual = guard.protect(path)
            if guard.hashes[actual] != expected:
                raise ValueError(f"Baseline response provenance mismatch: {path}")
    poles = real_pressure_poles(channels[0].response) if network == "YV" else []
    source = (
        result.get("dispersion_model", {}).get("source_file")
        or workflow["options"]["dispersion_file"]
    )
    model_path = Path(source)
    if not model_path.is_absolute():
        model_path = ROOT / model_path
    dispersion = np.loadtxt(guard.protect(model_path), delimiter=",", skiprows=1)
    boundary = any("boundary" in (e.get("reason") or "") for e in result["events"])
    raw = run / "_inputs/waveforms"
    if not raw.is_dir():
        raw = data_root / network
        if (raw / station).is_dir():
            raw /= station
    archive = (
        discover_waveforms(
            raw,
            {row.Seed_id for row in rows},
        )
        if boundary
        else []
    )
    events = []
    for index, original in enumerate(result["events"], 1):
        print(
            f"{network}.{station} event {index}/{len(result['events'])}: "
            f"{original['id']}",
            flush=True,
        )
        item = {k: original[k] for k in ("id", "origin_utc", "magnitude", "status")}
        item["saved_pressure_multiplier"] = original.get("pressure_multiplier")
        try:
            spectra = get_spectra(
                original, rows, channels, archive, dispersion, settings, output, guard
            )
        except ValueError as exc:
            if original["status"] == "candidate":
                raise
            item["trial_input_reason"] = str(exc)
            events.append(item)
            continue
        if spectra is None:
            item["trial_input_reason"] = (
                "No saved spectra; original input/selection guard retained"
            )
            events.append(item)
            continue
        f = spectra["frequency_hz"]
        # Verify live complex response matches the baseline's exact snapshot.
        checked_response(rows[0], channels[0], f, "DEF")
        checked_response(rows[-1], channels[-1], f, "VEL")
        reflected = pole_factor(f, poles)
        variants = {}
        for name, (reflect, method) in VARIANTS.items():
            factor = reflected if reflect else np.ones_like(reflected)
            variants[name] = evaluate(
                spectra,
                settings,
                bool(original.get("other_predicted_arrival_ids")),
                factor,
                method,
            )
        baseline = variants["nominal_absolute"]
        if original["status"] == "candidate":
            if baseline["status"] != "candidate" or not np.isclose(
                baseline["pressure_multiplier"],
                original["pressure_multiplier"],
                rtol=0,
                atol=1e-10,
            ):
                raise ValueError(
                    f"Nominal replay failed: {network}.{station} {original['id']}"
                )
        item.update(
            variants=variants,
            reflected_phase_deg=np.angle(
                spectra["cross_pa"] * reflected, deg=True
            ).tolist(),
            frequency_hz=f.tolist(),
            nominal_phase_deg=spectra["phase_deg"].tolist(),
            highpass=highpass_diagnostic(
                spectra, reflected, settings["coherence_squared_threshold"]
            ),
        )
        if station == "C07W":
            item["vertical_state"] = c07_state(data_root, original, guard)
        events.append(item)
    summaries = {}
    nominal_count = sum(
        e.get("variants", {}).get("nominal_absolute", {}).get("status") == "candidate"
        for e in events
    )
    if nominal_count != result["candidate_events"]:
        raise ValueError(f"Baseline candidate count changed at {network}.{station}")
    names = list(VARIANTS)
    if station == "C07W":
        names += [name + "_normal_state" for name in VARIANTS]
    minimum = settings["minimum_candidate_events"]
    for name in names:
        normal_only = name.endswith("_normal_state")
        variant = name.removesuffix("_normal_state")
        fitted = [
            e["variants"][variant]["pressure_multiplier"]
            for e in events
            if e.get("variants", {}).get(variant, {}).get("status") == "candidate"
            and (
                not normal_only or e.get("vertical_state", {}).get("state") == "normal"
            )
        ]
        diagnostic = float(np.median(fitted)) if fitted else None
        mixed_c07 = station == "C07W" and not normal_only
        gain = diagnostic if len(fitted) >= minimum and not mixed_c07 else None
        summaries[name] = {
            "candidate_events": len(fitted),
            "median_pressure_multiplier": gain,
            "diagnostic_pooled_median": diagnostic,
            "status": "mixed_vertical_states_diagnostic_only"
            if mixed_c07
            else "provisional_experiment"
            if gain
            else "insufficient_candidates",
        }
    report = {
        "network": network,
        "station": station,
        "baseline_result": entry["update"]["rayleigh_result_file"],
        "saved_station_gain": result.get("median_pressure_multiplier"),
        "settings": settings,
        "positive_real_pressure_poles_rad_s": poles,
        "summaries": summaries,
        "events": events,
    }
    (output / "trial.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


def comparison_rows(reports):
    rows = []
    for report in reports:
        network, station = report["network"], report["station"]
        published = {"YV": YV_PUBLISHED, "YL": YL_PUBLISHED}.get(network, {}).get(
            station
        )
        values = {
            "saved_before": {
                "median_pressure_multiplier": report["saved_station_gain"],
                "status": "historical_saved",
            }
        } | report["summaries"]
        for variant, summary in values.items():
            gain = summary["median_pressure_multiplier"]
            coefficient = (1 / gain if gain else None) if network == "YL" else gain
            rows.append(
                {
                    "network": network,
                    "station": station,
                    "variant": variant,
                    "pressure_multiplier_g": gain,
                    "comparison_coefficient": coefficient,
                    "coefficient_direction": "1/g" if network == "YL" else "g",
                    "published_coefficient": published,
                    "published_difference_percent": 100 * (coefficient / published - 1)
                    if coefficient and published
                    else None,
                    "change_from_saved_g_percent": 100
                    * (gain / report["saved_station_gain"] - 1)
                    if gain and report["saved_station_gain"]
                    else None,
                    "candidate_events": summary.get("candidate_events"),
                    "status": summary["status"],
                }
            )
    return rows


def plot_comparison(rows, output):
    import matplotlib.pyplot as plt

    networks = list(dict.fromkeys(r["network"] for r in rows))
    fig, axes = plt.subplots(
        len(networks), 1, figsize=(16, 5 * len(networks)), squeeze=False
    )
    for axis, network in zip(axes[:, 0], networks, strict=True):
        stations = list(
            dict.fromkeys(r["station"] for r in rows if r["network"] == network)
        )
        if not stations:
            axis.set_visible(False)
            continue
        for i, variant in enumerate(
            ["saved_before", *VARIANTS, "pole_weighted_log_normal_state"]
        ):
            x, y = [], []
            for j, st in enumerate(stations):
                match = next(
                    (
                        r
                        for r in rows
                        if r["network"] == network
                        and r["station"] == st
                        and r["variant"] == variant
                    ),
                    None,
                )
                if match and match["comparison_coefficient"] is not None:
                    x.append(j + (i - 2.5) * 0.055)
                    y.append(match["comparison_coefficient"])
            if x:
                axis.plot(x, y, "o", label=variant, markersize=5)
        for j, st in enumerate(stations):
            match = next(
                r for r in rows if r["network"] == network and r["station"] == st
            )
            if match["published_coefficient"] is not None:
                axis.plot(
                    j,
                    match["published_coefficient"],
                    "rx",
                    markersize=10,
                    label="published"
                    if not any(line.get_label() == "published" for line in axis.lines)
                    else None,
                )
        axis.set_xticks(range(len(stations)), stations)
        axis.tick_params(axis="x", labelrotation=60 if len(stations) > 15 else 0)
        axis.set_title(
            f"{network}: {'1/g' if network == 'YL' else 'g'}; provisional comparisons"
        )
        axis.grid(alpha=0.3)
        axis.legend(fontsize=7)
    fig.tight_layout()
    path = output / "comparison.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def selected_stations(workflow, network, overrides=None, all_stations=False):
    if overrides:
        return overrides
    if all_stations or network not in DEFAULTS:
        return sorted(
            {e["seed_id"].split(".")[1] for e in workflow["entries"].values()}
        )
    return DEFAULTS[network]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", choices=[*RUNS, "both", "all"], default="both")
    parser.add_argument("--stations", nargs="+", help="Single-network override")
    parser.add_argument(
        "--all-stations",
        action="store_true",
        help="Every station in the selected preserved workflow(s)",
    )
    parser.add_argument(
        "--data-root", type=Path, default=Path.home() / "Desktop/Data_Bruit_FM/Data"
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--open-plots", action="store_true")
    args = parser.parse_args()
    if args.stations and (args.network in {"both", "all"} or args.all_stations):
        parser.error("Use --stations with one network")
    networks = (
        list(RUNS)
        if args.network == "all"
        else list(DEFAULTS)
        if args.network == "both"
        else [args.network]
    )
    guard, jobs = Guard(), []
    for network in networks:
        run = ROOT / "results/NETWORKS" / RUNS[network]
        workflow = guard.json(run / "_workflow/rayleigh.json")
        if workflow["execution_status"] != "complete":
            raise ValueError("Wait for the baseline workflow to complete")
        for station in selected_stations(
            workflow, network, args.stations, args.all_stations
        ):
            matches = [
                e
                for e in workflow["entries"].values()
                if e["seed_id"].split(".")[1] == station
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Expected one saved pressure epoch: {network}.{station}"
                )
            jobs.append((network, station, matches[0], workflow, run))
    if args.dry_run:
        for network, station, entry, _, _ in jobs:
            r = guard.json(entry["update"]["rayleigh_result_file"])
            print(f"{network}.{station}: {len(r['events'])} saved records; dry run")
        print(f"Total: {len(jobs)} station epochs across {len(networks)} networks")
        guard.verify()
        return
    matplotlib.use("Agg")
    output = (
        args.output_dir
        or ROOT
        / "results/EXPERIMENTS/SOCRATES"
        / (
            f"{args.network}_{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}_"
            f"{uuid4().hex[:8]}"
        )
    ).resolve()
    output.mkdir(parents=True, exist_ok=False)
    for path in [
        Path(__file__),
        ROOT / "../compy_pressure/calibration_experiment.py",
        ROOT / "../compy_pressure/rayleigh_calibration.py",
        ROOT / "../compy_pressure/paired_spectra.py",
    ]:
        guard.protect(path)
    for path in (ROOT / "exports").glob("master_inventory.*"):
        guard.protect(path)
    master_cache = {}
    reports = [
        station_trial(
            net,
            st,
            entry,
            workflow,
            run,
            args.data_root.resolve(),
            output / f"{net}_{st}",
            guard,
            master_cache,
        )
        for net, st, entry, workflow, run in jobs
    ]
    rows = comparison_rows(reports)
    with (output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    path = plot_comparison(rows, output)
    guard.verify()
    (output / "experiment.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "protected_inputs_unchanged": True,
                "source_sha256": {str(p): h for p, h in guard.hashes.items()},
                "stations": [
                    {
                        k: r[k]
                        for k in (
                            "network",
                            "station",
                            "saved_station_gain",
                            "summaries",
                        )
                    }
                    for r in reports
                ],
                "comparison": rows,
                "scope": "all stored spectra and recoverable boundary events "
                "at selected stations",
                "limitations": [
                    "Experimental, no applied gains or master writes",
                    "Baseline event sampling/windows/sign/model retained",
                    "No INSU glitch removal, noise subtraction or model replacement",
                    "Approximate weights; coherence cap 0.999; no confidence intervals",
                    "C07W event-window labels are diagnostic; hardware is unconfirmed",
                    "High-pass diagnostics never become accepted scalar gains",
                ],
                "published_sources": {
                    "YV": "https://www.ipgp.fr/~stutz/2025-aminian-crawford-stutzmann-al-GJI.pdf",
                    "YL": "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688",
                },
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    print(
        f"\nCOMPLETE. Inputs unchanged. Results: {output}\n"
        f"Comparison: {output / 'comparison.csv'}\nPlot: {path}",
        flush=True,
    )
    if args.open_plots:
        subprocess.run(["open", str(path)], check=True)


if __name__ == "__main__":
    main()
