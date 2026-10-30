"""All-station microseism pressure–vertical check (diagnostic only).

For each pressure epoch with a saved response snapshot in the YV, YL, 7D and
XO masters: compute paired spectra from raw waveforms into this run's own
folder (never beside PPSDs, never the master), then run the microseism check
with the master depth. Recorded Rayleigh gains are compared afterwards; no
gain is fitted, applied or stored. Masters and waveforms are read in place.

Interval rule (fixed before running): the paired-spectra default 7-day interval
centred on the shared data/epoch coverage, unless --start NETWORK=DATE gives a
common week for that network; if that week does not fit, the centred default.
With --base-run RUN --shift-weeks K each station's week is its week in RUN moved
by K weeks, with no fallback (stations outside coverage are left unresolved), so
repeated weeks are independent. No interval is searched for coherence.
"""

import argparse
import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime

from compy_pressure.inventory import load_master_dataset
from compy_pressure.microseism_check import run_check, sha256
from compy_pressure.paired_spectra import has_saved_response, run_paired_spectra
from compy_pressure.ppsd import atomic_json, is_pressure_unit

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

ROOT = Path(__file__).resolve().parents[1]
MAIN = Path.home() / "Projects/bruit-fm"
DATA = Path.home() / "Desktop/Data_Bruit_FM/Data"
# Masters, waveform roots and resampling as used by the completed Rayleigh runs.
NETWORKS = {
    "YV": dict(master=MAIN / "exports", data=DATA),
    "YL": dict(master=MAIN / "exports", data=DATA),
    "7D": dict(
        master=MAIN / "results/NETWORKS/7D_ALL_NO_TILT/_master",
        data=MAIN / "results/NETWORKS/7D_ALL_NO_TILT/_inputs/waveforms",
    ),
    "XO": dict(
        master=MAIN / "results/NETWORKS/XO_ALL_GAUGES_NO_TILT/_master",
        data=MAIN / "results/NETWORKS/XO_ALL_GAUGES_NO_TILT/_inputs/waveforms",
        resample_hz=1.0,
        tolerance_ppm=0.1,
    ),
}
FREQUENCY_MAX = 0.4  # below the 0.5 Hz Nyquist of 1 Hz channels
LAYERS = ROOT / "../compy_pressure/reference_data/rr52_historical_rayleigh_layers.csv"
DENSITY, SOUND_SPEED = 1028.0, 1500.0  # same as the Rayleigh reference settings
GAIN_TABLE = ROOT / "docs/socrates_all_station_comparison.csv"
BASELINE = ROOT / "docs/calibration_network_summary.json"
GAIN_SOURCES = ("baseline", "trial", "published")
RAYLEIGH_SIGN = {
    n["network"]: n["options"]["pressure_sign"]
    for n in json.loads(BASELINE.read_text())["networks"]
}


def apply_sources(sources):
    """Point networks at other masters/waveform folders (e.g. a from-scratch run).

    ``sources`` maps network -> {"master": dir, "data": dir}. Resampling and
    tolerance settings of the completed runs are kept. Returns the mapping
    recorded in settings.json, so shifted weeks reuse the same inputs.
    """
    recorded = {}
    for network, source in (sources or {}).items():
        if network not in NETWORKS or set(source) != {"master", "data"}:
            raise ValueError(f"Invalid microseism source for {network}: {source}")
        paths = {k: Path(v).expanduser().resolve() for k, v in source.items()}
        missing = [str(v) for v in paths.values() if not v.exists()]
        if missing:
            raise FileNotFoundError(f"Microseism sources missing: {missing}")
        NETWORKS[network] = NETWORKS[network] | paths
        recorded[network] = {k: str(v) for k, v in paths.items()}
    return recorded


def recorded_gains(table=GAIN_TABLE):
    """Baseline/trial/published multipliers keyed (network, station), g direction.

    Published coefficients recorded as 1/g (YL, Zha & Webb 2016) are inverted.
    """
    gains = {}
    with Path(table).open() as stream:
        for row in csv.DictReader(stream):
            published = row["published_coefficient"]
            direction = row["comparison_direction"]
            if direction not in {"g", "1/g"}:
                raise ValueError("Unexpected published comparison direction")
            value = float(published) if published else None
            if value is not None and direction == "1/g":
                value = 1 / value
            gains[row["network"], row["station"]] = {
                "baseline": float(row["previous_g"]) if row["previous_g"] else None,
                "trial": float(row["pole_trial_g"]) if row["pole_trial_g"] else None,
                "published": value,
                "published_direction_recorded": direction,
            }
    summary = json.loads(BASELINE.read_text())
    for network in summary["networks"]:
        for item in network["stations"]:
            key = network["network"], item["station"]
            entry = gains.setdefault(key, dict.fromkeys(GAIN_SOURCES))
            if entry["baseline"] != item["pressure_multiplier"]:
                raise ValueError(f"{key}: baseline gain sources disagree")
            entry["rayleigh_pressure_sign"] = network["options"]["pressure_sign"]
    return gains


def pressure_epochs(master, network):
    rows = master[master.Seed_id.str.startswith(f"{network}.")]
    rows = rows[rows.input_units.map(is_pressure_unit)]
    rows = rows[[has_saved_response(r) for _, r in rows.iterrows()]]
    return rows.drop_duplicates("Seed_id").sort_values("Seed_id")


def bound_fractions(bins_csv, band, gain):
    """Structure-free bound test of one recorded gain over long-wave bins."""
    with Path(bins_csv).open() as stream:
        rows = list(csv.DictReader(stream))
    f = np.array([float(r["frequency_hz"]) for r in rows])
    m = (f >= band[0]) & (f <= band[1])
    h1 = np.array([float(r["upper_bound_h1"]) for r in rows])[m]
    h2 = np.array([float(r["upper_bound_h2"]) for r in rows])[m]
    ok, bad = gain <= h1, gain > h2
    return {
        "consistent": float(ok.mean()),
        "violates": float(bad.mean()),
        "inconclusive": float((~ok & ~bad).mean()),
    }


def station_record(result, gains):
    band = result["bands"]["long_wave"]

    # Conditioning = structure-free H2 bound / model-conditional H1 value: the
    # whole noise (H1-H2) plus structure (c_model..inf) spread. Not reference-based.
    def conditioning(sub):
        values = (sub["upper_bound_h2"], sub["mode0_water_h1"])
        return values[0] / values[1] if None not in values else np.inf

    best = min(result["long_wave_sub_bands"], key=conditioning)
    model = band["model_conditional_multiplier"]
    record = {
        "depth_m": result["water"]["depth_m"],
        "first_resonance_hz": result["water"]["first_resonance_hz"],
        "long_wave_hz": band["frequency_hz"],
        "long_wave_coherence_squared": band["median_coherence_squared"],
        "long_wave_phase_deg": band["median_phase_deg"],
        "phase_summary_method": band.get(
            "phase_summary_method", "legacy_linear_median"
        ),
        "long_wave_fraction_sign_minus1": band["fraction_phase_within_90deg_of_zero"],
        "upper_bound_h1": band["structure_free_upper_bound"]["median_upper_bound_h1"],
        "upper_bound_h2": band["structure_free_upper_bound"]["median_upper_bound_h2"],
        "long_wave_mode0_water_h1": model["mode0_water_h1"]["median"],
        "long_wave_mode0_water_h2": model["mode0_water_h2"]["median"],
        "long_wave_mode0_no_water_h1": model["mode0_no_water_h1"]["median"],
        "long_wave_mode0_no_water_h2": model["mode0_no_water_h2"]["median"],
        "best_sub_band_hz": best["frequency_hz"],
        "best_sub_band_coherence_squared": best["median_coherence_squared"],
        "best_sub_band_mode0_water_h1": best["mode0_water_h1"],
        "best_sub_band_mode0_water_h2": best["mode0_water_h2"],
        "best_sub_band_upper_bound_h2": best["upper_bound_h2"],
        "best_sub_band_conditioning_ratio": conditioning(best),
        "literature_band_coherence_squared": result["bands"]["literature_0p1_0p2"][
            "median_coherence_squared"
        ],
        "literature_band_phase_deg": result["bands"]["literature_0p1_0p2"][
            "median_phase_deg"
        ],
        "check_summary": result["files"]["summary"],
    }
    low, high = best["mode0_water_h1"], best["mode0_water_h2"]
    for source in GAIN_SOURCES:
        g = gains.get(source)
        record[f"{source}_g"] = g
        if g is None or low is None or high is None:
            continue
        record[f"{source}_best_h1_pct"] = 100 * (low / g - 1)
        record[f"{source}_best_h2_pct"] = 100 * (high / g - 1)
        record[f"{source}_within_best_bracket"] = bool(low <= g <= high)
        record[f"{source}_bound"] = bound_fractions(
            result["files"]["bins"], band["frequency_hz"], g
        )
    return record


def shifted_start(base_run, network, station, weeks):
    """The station's week in base_run moved by whole weeks; None if it had none."""
    path = Path(base_run) / "stations" / f"{network}.{station}.json"
    start = json.loads(path.read_text()).get("start_utc") if path.exists() else None
    return None if start is None else str(UTCDateTime(start) + weeks * 7 * 86400)


def process_station(row, gains, network, start, master, days, root, fallback=True):
    station = row.Seed_id.split(".")[1]
    config = NETWORKS[network]
    data = config["data"] / network / station
    if not data.is_dir():
        data = config["data"]

    def paired_run(start):
        return run_paired_spectra(
            config["master"],
            data_dir=data,
            network=network,
            station=station,
            pressure_id=row.Seed_id,
            start=start,
            days=days,
            frequency_max=FREQUENCY_MAX,
            output_directory=root / "paired" / network / station / uuid4().hex[:12],
            event_resample_rate=config.get("resample_hz"),
            native_rate_tolerance_ppm=config.get("tolerance_ppm", 0),
            update_master=False,
            reference_rows=master,
        )

    record = {
        "network": network,
        "station": station,
        "seed_id": row.Seed_id,
        "interval_rule": (
            "centred_default"
            if not start
            else ("common_week" if fallback else "shifted_week")
        ),
    }
    paired = paired_run(start)
    if start and fallback and paired["status"] == "insufficient_common_coverage":
        record.update(common_week_reason=paired.get("reason"))
        record["interval_rule"] = "centred_default"
        paired = paired_run(None)
    record.update(
        paired_status=paired["status"],
        paired_reason=paired.get("reason"),
        paired_result=paired["files"].get("result"),
        start_utc=paired.get("start_utc"),
        end_utc=paired.get("end_utc"),
        retained_windows=paired.get("retained_segments", 0),
        requested_windows=paired.get("requested_segments"),
    )
    if paired["status"] != "comparison_only":
        record["status"] = "unresolved_pairing"
        return record
    depth = float(row.Profondeur)
    result = run_check(
        Path(paired["files"]["result"]).parent,
        None,
        LAYERS,
        root / "checks" / network,
        water=(depth, DENSITY, SOUND_SPEED),
    )
    known = gains.get((network, station), {})
    record.update(status="checked", **station_record(result, known))
    record["rayleigh_pressure_sign"] = RAYLEIGH_SIGN[network]
    return record


def plot_summary(records, destination, show):
    done = [r for r in records if r.get("status") == "checked"]
    if not done:
        return
    x = np.arange(len(done))
    figure, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)
    low = np.array([r["best_sub_band_mode0_water_h1"] for r in done], float)
    high = np.array([r["best_sub_band_mode0_water_h2"] for r in done], float)
    axes[0].vlines(
        x,
        low,
        high,
        color="b",
        lw=4,
        alpha=0.5,
        label="microseism H1–H2 (best sub-band)",
    )
    axes[0].plot(
        x,
        [r["upper_bound_h2"] for r in done],
        "r_",
        ms=12,
        label="structure-free bound (H2)",
    )
    markers = {"baseline": "ko", "trial": "C1x", "published": "g^"}
    for source, style in markers.items():
        axes[0].plot(
            x,
            [r.get(f"{source}_g") or np.nan for r in done],
            style,
            label=f"{source} Rayleigh g",
        )
    axes[0].set(
        ylabel="pressure multiplier g",
        ylim=(0, 2.5),
        title="Microseism check vs recorded Rayleigh gains (diagnostic)",
    )
    axes[0].legend(fontsize=8, ncol=3)
    for source, style in (("baseline", "ko"), ("trial", "C1x")):
        for end, alpha in (("h1", 1.0), ("h2", 0.4)):
            axes[1].plot(
                x,
                [r.get(f"{source}_best_{end}_pct", np.nan) for r in done],
                style,
                alpha=alpha,
                label=f"{source} ({end.upper()})",
            )
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].set(ylabel="microseism / Rayleigh − 1 (%)", ylim=(-50, 50))
    axes[1].legend(fontsize=8, ncol=4)
    axes[2].bar(
        x - 0.2,
        [r["long_wave_fraction_sign_minus1"] for r in done],
        0.4,
        label="fraction of long-wave bins with sign −1 phase",
    )
    axes[2].bar(
        x + 0.2,
        [r["best_sub_band_coherence_squared"] for r in done],
        0.4,
        label="best sub-band coherence²",
    )
    twin = axes[2].twinx()
    twin.semilogy(x, [r["best_sub_band_conditioning_ratio"] for r in done], "kd", ms=4)
    twin.set_ylabel("conditioning ratio (◆, log)")
    axes[2].set(ylim=(0, 1.05), xticks=x)
    axes[2].set_xticklabels(
        [f"{r['network']}.{r['station']}" for r in done], rotation=90, fontsize=6
    )
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(True, alpha=0.3)
    figure.tight_layout()
    figure.savefig(destination, dpi=130)
    if show:
        plt.show()
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--network", action="append", choices=sorted(NETWORKS), help="Default all"
    )
    parser.add_argument("--station", action="append", help="Repeatable; default all")
    parser.add_argument(
        "--start",
        action="append",
        default=[],
        metavar="NETWORK=DATE",
        help="Common-week start for one network; default centred intervals",
    )
    parser.add_argument("--days", type=int, help="Default: base-run duration, else 7")
    parser.add_argument(
        "--base-run",
        type=Path,
        help="Completed run whose per-station weeks --shift-weeks moves",
    )
    parser.add_argument(
        "--shift-weeks",
        type=int,
        help="Whole weeks to move each station's base-run week (nonzero)",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Default: beside --base-run, else results/EXPERIMENTS/MICROSEISM_ALL",
    )
    parser.add_argument("--resume", type=Path, help="Existing run folder to continue")
    parser.add_argument(
        "--gain-table",
        type=Path,
        default=GAIN_TABLE,
        help="Rayleigh/published table with the socrates comparison columns",
    )
    parser.add_argument("--open-plots", action="store_true")
    parser.add_argument(
        "--sources",
        type=Path,
        help="JSON {network: {master, data}} overriding the completed-run inputs",
    )
    args = parser.parse_args()
    if (args.base_run is None) != (args.shift_weeks is None) or args.shift_weeks == 0:
        parser.error("--base-run and a nonzero --shift-weeks go together")
    if args.base_run and args.start:
        parser.error("--start cannot be combined with --base-run")
    base = (
        json.loads((args.base_run / "settings.json").read_text())
        if args.base_run
        else None
    )
    sources = (
        json.loads(args.sources.read_text())
        if args.sources
        else (base.get("sources") if base else None)
    )
    sources = apply_sources(sources)
    networks = args.network or (base["networks"] if base else list(NETWORKS))
    if base and set(networks) - set(base["networks"]):
        parser.error("--network must belong to --base-run")
    if args.days is None:
        args.days = base["days"] if base else 7
    if args.days <= 0 or (base and args.days != base["days"]):
        parser.error("Days must be positive and match --base-run when shifting weeks")
    if args.output_root is None:  # shifted weeks live beside their base run
        args.output_root = (
            args.base_run.resolve().parent
            if args.base_run
            else ROOT / "results/EXPERIMENTS/MICROSEISM_ALL"
        )
    starts = dict(item.split("=", 1) for item in args.start)
    if set(starts) - set(networks):
        raise ValueError("--start names a network that is not selected")

    masters = sorted(
        {
            NETWORKS[n]["master"] / name
            for n in networks
            for name in ("master_inventory.csv", "master_inventory.xlsx")
        }
    )
    masters = [p for p in masters if p.exists()]
    before = {str(p): sha256(p) for p in masters}
    settings = {
        "networks": networks,
        "starts": starts,
        "days": args.days,
        "base_run": str(args.base_run.resolve()) if args.base_run else None,
        "shift_weeks": args.shift_weeks,
        "interval_rule": __doc__.split("Interval rule")[1].split("\n\n")[0],
        "network_inputs": {
            n: {k: str(v) for k, v in NETWORKS[n].items()} for n in networks
        },
        "sources": sources or None,
        "frequency_max_hz": FREQUENCY_MAX,
        "master_sha256": before,
        "layers_file": str(LAYERS),
        "density_kg_m3": DENSITY,
        "sound_speed_m_s": SOUND_SPEED,
        "depth_source": "master Profondeur",
        "rayleigh_pressure_sign": {n: RAYLEIGH_SIGN[n] for n in networks},
        "gain_sources": {
            str(args.gain_table): sha256(args.gain_table),
            str(BASELINE): sha256(BASELINE),
        },
        "best_sub_band_rule": "min H2 bound / mode0 H1, 4 fixed log sub-bands",
        "status": "diagnostic_only; no gain, master or score written",
    }
    if args.resume:
        root = args.resume
        saved = json.loads((root / "settings.json").read_text())
        keys = ("networks", "starts", "days", "base_run", "shift_weeks", "sources")
        if any(saved.get(k) != settings[k] for k in keys):
            raise ValueError("Resume must keep the original networks and intervals")
    else:
        stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
        root = args.output_root / f"run_{stamp}_{uuid4().hex[:8]}"
        root.mkdir(parents=True)
        settings["git_commit"] = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
        ).stdout.strip()
        atomic_json(root / "settings.json", settings)
    (root / "stations").mkdir(exist_ok=True)
    gains = recorded_gains(args.gain_table)
    records = []
    for network in networks:
        master = load_master_dataset(NETWORKS[network]["master"])
        epochs = pressure_epochs(master, network)
        if args.station:
            epochs = epochs[epochs.Seed_id.str.split(".").str[1].isin(args.station)]
        for row in epochs.itertuples():
            station = row.Seed_id.split(".")[1]
            path = root / "stations" / f"{network}.{station}.json"
            if path.exists():
                records.append(json.loads(path.read_text()))
                print(f"{network}.{station}: reused {records[-1]['status']}")
                continue
            start = (
                shifted_start(args.base_run, network, station, args.shift_weeks)
                if args.base_run
                else starts.get(network)
            )
            if args.base_run and start is None:
                record = {
                    "network": network,
                    "station": station,
                    "seed_id": row.Seed_id,
                    "status": "no_base_interval",
                }
                atomic_json(path, record)
                records.append(record)
                print(f"{network}.{station}: {record['status']}", flush=True)
                continue
            try:
                record = process_station(
                    row,
                    gains,
                    network,
                    start,
                    master,
                    args.days,
                    root,
                    fallback=not args.base_run,
                )
            except Exception as exc:  # isolate one station; the batch continues
                record = {
                    "network": network,
                    "station": station,
                    "seed_id": row.Seed_id,
                    "status": "processing_error",
                    "error": repr(exc),
                }
            atomic_json(path, record)
            records.append(record)
            print(f"{network}.{station}: {record['status']}", flush=True)

    after = {str(p): sha256(p) for p in masters}
    if after != before:
        raise RuntimeError("A master inventory changed during the microseism batch")
    columns = sorted({k for r in records for k in r if not isinstance(r[k], dict)})
    with (root / "stations.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=columns, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)
    counts = {}
    for r in records:
        key = f"{r['network']}:{r['status']}"
        counts[key] = counts.get(key, 0) + 1
    atomic_json(
        root / "summary.json",
        {
            "status_counts": counts,
            "stations": records,
            "master_unchanged": True,
            "gain_applied": False,
        },
    )
    plot_summary(records, root / "microseism_check.png", args.open_plots)
    print(f"{counts}: {root}")


if __name__ == "__main__":
    main()
