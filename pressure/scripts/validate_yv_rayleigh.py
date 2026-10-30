"""Read-only YV event/quiet review and independent raw-count spectral audit.

Run from the project root. Writes new validation products only; no master updates.
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from obspy import Stream, UTCDateTime, read, read_inventory
from scipy.signal import welch

from compy_pressure.rayleigh_validation import count_domain_products


def load_counts(paths, seed, start, seconds, rate):
    """Independent MiniSEED extraction; exact continuity, no filling/interpolation."""
    end = start + seconds
    dates = {(start + i * 86400).strftime("%Y.%j") for i in range(3)}
    dates.add(end.strftime("%Y.%j"))
    stream = Stream()
    for path in paths:
        p = Path(path)
        if seed in p.name and any(p.name.endswith(d) for d in dates):
            stream += read(str(p), starttime=start - 2 / rate, endtime=end + 2 / rate)
    stream = stream.select(id=seed)
    stream.sort(keys=["starttime"])
    for left, right in zip(stream, stream[1:]):
        offset = (right.stats.starttime - left.stats.starttime) * rate
        sample = round(offset)
        if abs(offset - sample) / rate > 1e-6:
            raise ValueError(f"Fractional sample-grid offset in {seed}")
        if sample > len(left):
            raise ValueError(f"Gap in independent read of {seed}")
        overlap = min(len(left) - sample, len(right))
        if overlap > 0 and not np.array_equal(
            left.data[sample : sample + overlap], right.data[:overlap]
        ):
            raise ValueError(f"Conflicting overlap in independent read of {seed}")
    stream.merge(method=0, fill_value=None)
    if len(stream) != 1 or stream[0].stats.sampling_rate != rate:
        raise ValueError(f"Nonunique trace/rate: {seed}")
    tr = stream[0]
    sample = int(round((start - tr.stats.starttime) * rate))
    actual = tr.stats.starttime + sample / rate
    if abs(actual - start) > 1e-6:
        raise ValueError(f"Independent sample grid mismatch: {seed}")
    x = tr.data[sample : sample + round(seconds * rate)]
    if (
        len(x) != round(seconds * rate)
        or np.ma.is_masked(x)
        or not np.isfinite(x).all()
    ):
        raise ValueError(f"Noncontinuous samples: {seed}")
    return x.astype(float)


def snapshot_channels(row, summary):
    station_root = Path(row["result_file"]).parents[2].parent
    channels = []
    for key in summary["epoch_keys"]:
        seed = key["Seed_id"]
        code = seed.split(".")[-1]
        inv = read_inventory(str(next(station_root.glob(f"*{code}*/inventory.xml"))))
        found = [
            c
            for n in inv
            for s in n
            for c in s
            if c.code == code
            and str(c.start_date or "") == key["channel_start_time"]
            and str(c.end_date or "") == key["channel_end_time"]
        ]
        if len(found) != 1:
            raise ValueError("Snapshot epoch ambiguity")
        channels.append((seed, found[0]))
    return channels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary", type=Path, default=Path("docs/yv_rayleigh_completed_summary.json")
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/VALIDATION/YV_rayleigh")
    )
    args = parser.parse_args()
    args.output_dir = args.output_dir / f"run_{uuid4().hex[:12]}"
    args.output_dir.mkdir(parents=True, exist_ok=False)
    source = json.loads(args.summary.read_text())
    before = {
        str(p): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in Path("exports").glob("master_inventory.*")
    }
    audit = {
        "scope": (
            "No response corrections/master updates; independent CSD still shares "
            "response metadata and saved rotation"
        ),
        "candidate_event_reviews": [],
        "raw_estimator_checks": [],
        "response_flags": [],
        "quiet_checks": [],
        "errors": [],
    }
    reference = np.genfromtxt(
        "../compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv",
        delimiter=",",
        names=True,
    )
    representatives = {"RR03", "RR07", "RR09", "RR22", "RR40", "RR52", "RR55"}
    for row in source["stations"]:
        summary = json.loads(Path(row["result_file"]).read_text())
        if not summary.get("epoch_keys") or not summary.get("events"):
            continue
        station = summary["station"].split(".")[-1]
        channels = snapshot_channels(row, summary)
        pressure = channels[0][1]
        poles = [
            complex(p)
            for st in pressure.response.response_stages
            for p in getattr(st, "poles", [])
            if complex(p).real > 0
        ]
        fir = [
            {
                "stage": st.stage_sequence_number,
                "dc_coefficient_sum": float(sum(st.numerator)),
            }
            for st in pressure.response.response_stages
            if getattr(st, "numerator", [])
            and not np.isclose(sum(st.numerator), 1, atol=1e-4)
        ]
        audit["response_flags"].append(
            {
                "station": station,
                "pressure_positive_real_poles": [[p.real, p.imag] for p in poles],
                "unnormalized_digital_coefficient_stages": fir,
                "input_units": pressure.response.instrument_sensitivity.input_units,
                "sensor_label": pressure.sensor.type if pressure.sensor else None,
            }
        )
        events = [e for e in summary["events"] if e["status"] == "candidate"]
        paths = sorted({s["path"] for s in summary["source_files"]})
        settings = summary["settings"]
        rate = float(pressure.sample_rate)
        for event in events:
            start, end = UTCDateTime(event["start_utc"]), UTCDateTime(event["end_utc"])
            body = {
                k: bool(start <= UTCDateTime(event[k]) <= end)
                for k in ["predicted_first_p_utc", "predicted_first_s_utc"]
                if event.get(k)
            }
            rec = {
                "station": station,
                "event": event["id"],
                "body_arrival_inside_window": body,
                "quiet_start_utc": event["quiet_start_utc"],
                "tilt_angle_deg": event["tilt_angle_deg"],
                "quiet_vertical_power_ratio": event["quiet_vertical_power_ratio"],
                "packet_review": (
                    "Inspection required; timing flags do not prove contamination"
                ),
            }
            audit["candidate_event_reviews"].append(rec)
            try:
                quiet_start = UTCDateTime(event["quiet_start_utc"])
                # The saved request may precede the first native sample by <one sample.
                # Determine that first sample independently, then require exact grids.
                first_path = next(
                    Path(p)
                    for p in paths
                    if channels[0][0] in Path(p).name
                    and Path(p).name.endswith(quiet_start.strftime("%Y.%j"))
                )
                hdr = read(str(first_path), headonly=True)[0]
                offset = np.ceil((quiet_start - hdr.stats.starttime) * rate - 1e-5)
                actual = hdr.stats.starttime + offset / rate
                arrays = np.stack(
                    [
                        load_counts(paths, seed, actual, 7200, rate)
                        for seed, _ in channels
                    ]
                )
                time_ratios = []
                impulsive = []
                for x in arrays:
                    blocks = x.reshape(12, -1)
                    variances = np.var(blocks, axis=1)
                    time_ratios.append(float(np.max(variances) / np.median(variances)))
                    dx = np.diff(x)
                    mad = np.median(abs(dx - np.median(dx)))
                    impulsive.append(
                        float(np.max(abs(dx - np.median(dx))) / mad) if mad else None
                    )
                audit["quiet_checks"].append(
                    {
                        "station": station,
                        "event": event["id"],
                        "max_10min_variance_over_median_counts": time_ratios,
                        "largest_first_difference_over_mad_counts": impulsive,
                        "interpretation": (
                            "Screening metrics, not acceptance cutoffs; "
                            "instrument responses differ between channels"
                        ),
                    }
                )
                # Count-domain stationarity plot; no new deconvolution or gain applied.
                if station in representatives:
                    f, psd = welch(
                        arrays, rate, nperseg=round(1800 * rate), noverlap=0, axis=1
                    )
                    indices = (f >= 0.001) & (f <= 0.1)
                    fig, axes = plt.subplots(2, 1, figsize=(9, 7))
                    for j, (seed, _) in enumerate(channels):
                        axes[0].semilogy(
                            np.arange(12) * 10 + 5,
                            np.var(arrays[j].reshape(12, -1), axis=1),
                            label=seed.split(".")[-1],
                        )
                        axes[1].loglog(
                            f[indices], psd[j, indices], label=seed.split(".")[-1]
                        )
                    axes[0].set(
                        xlabel="Minutes in quiet interval",
                        ylabel="Variance [counts squared]",
                    )
                    axes[1].set(
                        xlabel="Frequency [Hz]", ylabel="PSD [counts squared / Hz]"
                    )
                    axes[0].legend()
                    axes[1].legend()
                    fig.suptitle(
                        f"{station} {event['id']} quiet screening: "
                        "raw counts, not physical comparison"
                    )
                    fig.tight_layout()
                    fig.savefig(args.output_dir / f"{station}_{event['id']}_quiet.png")
                    plt.close(fig)
            except Exception as exc:
                audit["errors"].append(
                    {
                        "station": station,
                        "event": event["id"],
                        "stage": "quiet",
                        "reason": str(exc),
                    }
                )
        if station not in representatives:
            continue
        available = events or [
            e for e in summary["events"] if e.get("files", {}).get("spectra")
        ]
        if not available:
            continue
        event = (
            max(available, key=lambda e: e.get("relative_rms", 0))
            if events
            else next(e for e in available if e["status"] == "polarity_review")
        )
        try:
            start = UTCDateTime(event["actual_start_utc"])
            duration = event["segment_count"] * 200
            arrays = np.stack(
                [
                    load_counts(paths, seed, start, duration, rate)
                    for seed, _ in channels
                ]
            )
            weights = np.array(event["rotation_coefficients"])
            if channels[3][1].dip > 0:
                weights[2] *= -1

            def response(f):
                return np.stack(
                    [
                        c.response.get_evalresp_response_for_frequencies(
                            f, output="DEF" if j == 0 else "VEL"
                        )
                        for j, (_, c) in enumerate(channels)
                    ]
                )

            result = count_domain_products(
                arrays,
                rate,
                response,
                weights,
                settings["band_hz"],
                settings["pressure_sign"],
            )
            with np.load(event["files"]["spectra"]) as saved:
                differences = {
                    k: float(np.max(abs(result[k] / saved[k] - 1)))
                    for k in ["pressure_psd", "acceleration_psd", "psd_ratio_pa"]
                }
                differences["coherence_abs"] = float(
                    np.max(
                        abs(result["coherence_squared"] - saved["coherence_squared"])
                    )
                )
                differences["phase_deg_abs"] = float(
                    np.max(
                        abs(
                            np.angle(
                                np.exp(
                                    1j
                                    * np.deg2rad(
                                        result["phase_deg"] - saved["phase_deg"]
                                    )
                                ),
                                deg=True,
                            )
                        )
                    )
                )
            record = {
                "station": station,
                "event": event["id"],
                "status": event["status"],
                "differences": differences,
            }
            f = result["frequency_hz"]
            velocity = np.interp(
                f, reference["frequency_hz"], reference["phase_velocity_m_s"]
            )
            k = (
                2
                * np.pi
                * f
                * np.sqrt(settings["sound_speed_m_s"] ** -2 - velocity**-2)
            )
            theory = settings["density_kg_m3"] * np.tan(k * settings["depth_m"]) / k
            mask = (result["coherence_squared"] >= 0.95) & (
                abs(result["phase_deg"]) <= 20
            )
            if event["status"] == "candidate":
                g = float(
                    np.linalg.lstsq(
                        result["psd_ratio_pa"][mask, None], theory[mask], rcond=None
                    )[0][0]
                )
                record.update(
                    independent_gain=g, saved_gain=event["pressure_multiplier"]
                )
            audit["raw_estimator_checks"].append(record)
            print(station, event["id"], differences, flush=True)
        except Exception as exc:
            audit["errors"].append(
                {
                    "station": station,
                    "event": event["id"],
                    "stage": "raw",
                    "reason": str(exc),
                }
            )
        (args.output_dir / "validation.json").write_text(
            json.dumps(audit, indent=2) + "\n"
        )
    after = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in before}
    if before != after:
        raise RuntimeError("Master changed during read-only validation")
    audit.update(
        master_hashes_preserved=True,
        master_sha256=before,
        event_status_counts=dict(
            Counter(e["status"] for e in audit["raw_estimator_checks"])
        ),
    )
    (args.output_dir / "validation.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(
        "Candidates reviewed:",
        len(audit["candidate_event_reviews"]),
        "quiet intervals:",
        len(audit["quiet_checks"]),
        "errors:",
        len(audit["errors"]),
    )


if __name__ == "__main__":
    main()
