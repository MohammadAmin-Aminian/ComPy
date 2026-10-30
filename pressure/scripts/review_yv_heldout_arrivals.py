"""Plot flagged native P/Z windows without refitting gains or changing event QC."""

import argparse
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from scipy.signal.windows import tukey

from compy_pressure.archive import discover_waveforms
from compy_pressure.calibration_experiment import pole_factor
from compy_pressure.paired_spectra import checked_response
from compy_pressure.rayleigh_calibration import load_aligned
from scripts.run_calibration_experiment import Guard
from scripts.run_yv_heldout_validation import ROOT, TRAINING, prepare, write_json


def plot_window(event, configuration, records, taup, destination):
    audit, rows, channels, model = configuration
    settings = audit["training_settings"]
    origin = UTCDateTime(event["origin_utc"])
    start, end = UTCDateTime(event["start_utc"]), UTCDateTime(event["end_utc"])
    left = origin - 300
    duration = float(np.ceil((end + 600 - left) / 60) * 60)
    arrays, actual, rate = load_aligned(records, rows, channels, left, duration)
    frequencies = np.fft.rfftfreq(arrays.shape[1], 1 / rate)
    lo, hi = settings["band_hz"]
    mask = (frequencies >= lo) & (frequencies <= hi)
    taper = tukey(arrays.shape[1], alpha=0.1)
    factor = pole_factor(frequencies[mask], audit["positive_real_pressure_poles_rad_s"])
    arrivals = [
        {"phase": a.name, "time_utc": str(origin + a.time)}
        for a in taup.get_travel_times(
            event["depth_km"],
            event["distance_km"] / 111.195,
            phase_list=["P", "Pdiff", "PKP", "S", "Sdiff", "SKS"],
        )
        if actual <= origin + a.time <= actual + duration
    ]
    mf = model[:, 0]
    group = model[(mf >= lo) & (mf <= hi), 2]
    group = np.r_[group, np.interp([lo, hi], mf, model[:, 2])]
    group_times = event["distance_km"] * 1000 / group
    block = max(1, round(rate))
    blocks = arrays.shape[1] // block
    elapsed = float(actual - origin) + np.arange(blocks) * block / rate
    fig, axes = plt.subplots(2, 2, figsize=(14, 7), sharex=True)
    for i, (row, channel) in enumerate(zip(rows, channels, strict=True)):
        data = arrays[i] - arrays[i].mean()
        raw = data[: blocks * block].reshape(blocks, block)
        axes[i, 0].fill_between(elapsed, raw.min(axis=1), raw.max(axis=1))
        axes[i, 0].set_ylabel(f"{row.Seed_id}\nDemeaned counts")
        response = checked_response(
            row, channel, frequencies[mask], "DEF" if i == 0 else "VEL"
        )
        if not np.all(np.isfinite(response) & (abs(response) > 0)):
            raise ValueError("Invalid overview response")
        fourier = np.zeros(len(frequencies), dtype=complex)
        fourier[mask] = np.fft.rfft(data * taper)[mask] / response
        if i == 0:
            fourier[mask] *= settings["pressure_sign"] * factor
        else:
            fourier[mask] *= 2j * np.pi * frequencies[mask]
            if channel.dip > 0:
                fourier *= -1
        wave = np.fft.irfft(fourier, arrays.shape[1])
        axes[i, 1].plot(elapsed, wave[: blocks * block : block], lw=0.8)
        axes[i, 1].set_ylabel("Pressure (Pa)" if i == 0 else "Upward a (m/s²)")
    for ax in axes.flat:
        ax.axvspan(
            start - origin,
            end - origin,
            color="purple",
            alpha=0.12,
            label="Original validation window",
        )
        ax.axvspan(
            group_times.min(),
            group_times.max(),
            color="green",
            alpha=0.12,
            label="Reference group-arrival range",
        )
        for arrival in arrivals:
            ax.axvline(
                UTCDateTime(arrival["time_utc"]) - origin,
                ls="--",
                lw=0.8,
                color="red" if arrival["phase"].startswith("P") else "orange",
                label=f"Predicted {arrival['phase']}",
            )
        ax.grid(alpha=0.2)
        ax.set_xlim(elapsed[0], elapsed[-1])
    axes[0, 0].set_title("Raw counts: one-second minimum/maximum")
    axes[0, 1].set_title(f"{lo:g}–{hi:g} Hz; phase-only pole reflection in memory")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    for ax in axes[1]:
        ax.set_xlabel("Seconds after earthquake origin")
    fig.suptitle(
        f"{audit['station']} {event['id']} Mw {event['magnitude']} — arrival review\n"
        "10% Tukey taper for display only; no gain applied; original QC unchanged"
    )
    fig.tight_layout()
    fig.savefig(destination, dpi=130)
    plt.close(fig)
    return {
        "station": audit["station"],
        "event_id": event["id"],
        "original_status": event["status"],
        "plot": str(destination),
        "overview_start_utc": str(actual),
        "overview_seconds": duration,
        "native_rate_hz": rate,
        "predicted_arrivals": arrivals,
        "original_window_start_utc": str(start),
        "original_window_end_utc": str(end),
        "gains_refitted_or_applied": False,
        "original_qc_changed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--stations", nargs="+", default=["RR01", "RR06"])
    parser.add_argument("--training-experiment", type=Path, default=TRAINING)
    parser.add_argument("--minimum-magnitude", type=float, default=6.0)
    parser.add_argument(
        "--data-dir", type=Path, default=Path.home() / "Desktop/Data_Bruit_FM/Data/YV"
    )
    args = parser.parse_args()
    guard = Guard()
    guard.protect(__file__)
    completed = guard.json(args.validation_dir / "validation.json")
    if completed["status"] != "complete":
        raise ValueError("Arrival review requires a completed validation")
    configs, _, _, _, _ = prepare(args, guard)
    selected = []
    for configuration in configs:
        audit = configuration[0]
        report = guard.json(args.validation_dir / audit["station"] / "validation.json")
        if report["frozen_gain"] != audit["frozen_gain"]:
            raise ValueError("Training and validation frozen gains differ")
        selected.extend(
            (event, configuration)
            for event in report["events"]
            if event["status"] == "arrival_isolation_review"
        )
    records = discover_waveforms(
        args.data_dir, {r.Seed_id for _, rows, _, _ in configs for r in rows}
    )
    stats = {
        r.path: (r.path.stat().st_size, r.path.stat().st_mtime_ns) for r in records
    }
    output = (
        args.output_dir
        or ROOT
        / "results/VALIDATION/YV_ARRIVAL_REVIEW"
        / f"waveforms_{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}_{uuid4().hex[:8]}"
    )
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    taup = TauPyModel(model="iasp91")
    reviews = []
    for event, config in selected:
        station = config[0]["station"]
        print(f"{station} {event['id']}: original window, no refit", flush=True)
        reviews.append(
            plot_window(
                event,
                config,
                records,
                taup,
                output / f"{station}_{event['id']}_waveforms.png",
            )
        )
    guard.verify()
    for path, expected in stats.items():
        if (path.stat().st_size, path.stat().st_mtime_ns) != expected:
            raise RuntimeError(f"Raw file changed: {path}")
    write_json(
        output / "review.json",
        {
            "status": "complete",
            "source_validation": str(args.validation_dir.resolve()),
            "reviews": reviews,
            "inputs_unchanged": True,
            "source_sha256": {str(p): h for p, h in guard.hashes.items()},
            "raw_file_stats": {str(p): list(v) for p, v in stats.items()},
            "limitations": [
                "Visual review is exploratory, not an automatic mode classifier.",
                "Predicted phases do not prove contamination; QC labels unchanged.",
                "Display taper/filter changes amplitude envelopes; no new estimates.",
            ],
        },
    )
    print(f"COMPLETE: {output}", flush=True)


if __name__ == "__main__":
    main()
