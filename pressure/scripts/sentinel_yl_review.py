"""Read-only Sentinel review of YL Rayleigh outliers; never writes scientific products.

Independent of the production spectral path: raw MiniSEED reads, SciPy Welch
cross-spectra, response division of the cross-spectrum and a separately written
water-column formula. Shared with production: StationXML responses, the
dispersion table and the saved event windows. Outputs go to a new directory.
"""

import argparse
import csv
import glob
import json
from pathlib import Path

import numpy as np
from obspy import UTCDateTime, read, read_inventory
from scipy.signal import butter, csd, sosfiltfilt


def station_counts(root, station, channel, start, end):
    stream = None
    for day in sorted({start.strftime("%Y%m%d"), (end - 1).strftime("%Y%m%d")}):
        pattern = root / station / "waveforms" / f"YL.{station}.*.{channel}.{day}.mseed"
        for path in glob.glob(str(pattern)):
            part = read(path)
            stream = part if stream is None else stream + part
    if stream is None:
        raise ValueError(f"{station}.{channel}: no files")
    stream.merge(method=0, fill_value=None)
    stream.trim(start, end - 0.5 / stream[0].stats.sampling_rate, nearest_sample=True)
    if len(stream) != 1 or np.ma.is_masked(stream[0].data):
        raise ValueError(f"{station}.{channel}: gapped window")
    return stream[0].data.astype(float), stream[0].stats.sampling_rate


def water(frequency, phase_velocity, depth, density=1028.0, sound=1500.0):
    """z down, P(0)=0, dP/dz(H)=rho*a_up: P/a_up = rho*tan(kz*H)/kz."""
    kz = 2 * np.pi * frequency * np.sqrt(1 / sound**2 - 1 / phase_velocity**2 + 0j)
    return np.real(density * np.tan(kz * depth) / kz)


def pressure_over_acceleration(root, station, start, end, dispersion, depth, band):
    inventory = read_inventory(str(root / station / "inventory.xml"))
    channels = {c.code: c for c in inventory[0][0]}
    pressure, rate = station_counts(root, station, "LDH", start, end)
    vertical, vertical_rate = station_counts(root, station, "LHZ", start, end)
    if rate != vertical_rate:
        raise ValueError("Unequal sample rates")
    n = min(len(pressure), len(vertical))
    pressure, vertical = pressure[:n], vertical[:n]
    options = dict(fs=rate, window="hann", nperseg=round(200 * rate), noverlap=0)
    f, spp = csd(pressure, pressure, detrend="constant", **options)
    _, szz = csd(vertical, vertical, detrend="constant", **options)
    _, szp = csd(vertical, pressure, detrend="constant", **options)  # conj(Z)*P
    keep = (f >= band[0]) & (f <= band[1])
    f, spp, szz, szp = f[keep], spp[keep].real, szz[keep].real, szp[keep]
    rp = channels["LDH"].response.get_evalresp_response_for_frequencies(f, "DEF")
    rz = channels["LHZ"].response.get_evalresp_response_for_frequencies(f, "VEL")
    upward = -1.0 if channels["LHZ"].dip > 0 else 1.0
    omega = 2 * np.pi * f
    pp = spp / abs(rp) ** 2
    aa = omega**2 * szz / abs(rz) ** 2
    pa = upward * szp / (rp * np.conj(rz)) * (-1j * omega)  # <P conj(A)>
    ratio = np.sqrt(pp / aa)
    theory = water(f, np.interp(f, dispersion[:, 0], dispersion[:, 1]), depth)
    return f, ratio, theory, abs(pa) ** 2 / (pp * aa), np.angle(pa, deg=True)


def daily_microseism_rms(root, station, channel, band=(0.12, 0.3)):
    sos = butter(4, band, btype="band", fs=1.0, output="sos")
    result = {}
    pattern = root / station / "waveforms" / f"YL.{station}.*.{channel}.*.mseed"
    for path in sorted(glob.glob(str(pattern))):
        stream = read(path)
        stream.merge(method=0, fill_value=None)
        hourly = []
        for trace in stream.split():
            if trace.stats.sampling_rate != 1.0 or trace.stats.npts < 7200:
                continue
            data = sosfiltfilt(sos, trace.data.astype(float) - trace.data.mean())
            data = data[600:-600]
            hourly += [
                np.sqrt(np.mean(data[k : k + 3600] ** 2))
                for k in range(0, len(data) - 3599, 3600)
            ]
        if hourly:
            result[path.split(".")[-2]] = float(np.median(hourly))
    return result


def vertical_state(ratio):
    """Classify C07W/C08W LHZ ratio; extremes (e.g. deployment settling) are other."""
    if 0.65 < ratio < 1.2:
        return "normal"
    if 0.3 < ratio < 0.6:
        return "reduced"
    return "other"


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True, help="Data/YL")
    parser.add_argument("--workflow", type=Path, required=True, help="rayleigh.json")
    parser.add_argument("--dispersion-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    dispersion = np.loadtxt(args.dispersion_file, delimiter=",", skiprows=1)
    workflow = json.loads(args.workflow.read_text())
    saved = {
        entry["seed_id"].split(".")[1]: json.loads(Path(entry["files"][0]).read_text())
        for entry in workflow["entries"].values()
    }
    summary = {}

    # 1. C07W vertical amplitude state against two neighbours.
    rms = {
        f"{s}.{c}": daily_microseism_rms(args.data_root, s, c)
        for s in ("C07W", "C08W", "C09W")
        for c in ("LHZ", "LDH")
    }
    days = sorted(set(rms["C07W.LHZ"]) & set(rms["C08W.LHZ"]) & set(rms["C09W.LHZ"]))
    rows = []
    for day in days:
        z78 = rms["C07W.LHZ"][day] / rms["C08W.LHZ"][day]
        rows.append(
            dict(
                day=day,
                z_c07w_over_c08w=z78,
                z_c07w_over_c09w=rms["C07W.LHZ"][day] / rms["C09W.LHZ"][day],
                z_c08w_over_c09w=rms["C08W.LHZ"][day] / rms["C09W.LHZ"][day],
                p_c07w_over_c08w=rms["C07W.LDH"].get(day, np.nan)
                / rms["C08W.LDH"].get(day, np.nan),
                state=vertical_state(z78),
            )
        )
    write_csv(args.output_dir / "c07w_daily_vertical_ratio.csv", rows)
    state = {row["day"]: row["state"] for row in rows}
    events = []
    for event in saved["C07W"]["events"]:
        gain = event.get(
            "pressure_multiplier", event.get("diagnostic_pressure_amplitude_multiplier")
        )
        if gain is not None:
            day = UTCDateTime(event["origin_utc"]).strftime("%Y%m%d")
            events.append(
                dict(
                    event_id=event["id"],
                    origin_utc=event["origin_utc"],
                    status=event["status"],
                    pressure_multiplier=gain,
                    vertical_state=state.get(day, "no_ratio"),
                )
            )
    write_csv(args.output_dir / "c07w_event_gain_by_state.csv", events)
    ratios = np.array([row["z_c07w_over_c08w"] for row in rows])
    summary["c07w"] = {
        "days": len(rows),
        "reduced_days": int(sum(r["state"] == "reduced" for r in rows)),
        "normal_days": int(sum(r["state"] == "normal" for r in rows)),
        "median_reduced_ratio": float(
            np.median(ratios[(ratios > 0.3) & (ratios < 0.6)])
        ),
        "median_normal_ratio": float(
            np.median(ratios[(ratios > 0.65) & (ratios < 1.2)])
        ),
        "median_c08w_over_c09w": float(
            np.median([r["z_c08w_over_c09w"] for r in rows])
        ),
        "candidate_gain_by_state": {
            s: [
                e["pressure_multiplier"]
                for e in events
                if e["vertical_state"] == s and e["status"] == "candidate"
            ]
            for s in ("reduced", "normal", "other", "no_ratio")
        },
    }

    # 2. B01W events that reached the gain-grid boundary, with neighbours.
    b01w = []
    for event in saved["B01W"]["events"]:
        if "boundary" not in (event.get("reason") or ""):
            continue
        for station in ("B01W", "B04W", "A02W"):
            match = next(
                (
                    e
                    for e in saved[station]["events"]
                    if e["id"] == event["id"] and "start_utc" in e
                ),
                None,
            )
            if match is None:
                continue
            f, ratio, theory, coherence, phase = pressure_over_acceleration(
                args.data_root,
                station,
                UTCDateTime(match["start_utc"]),
                UTCDateTime(match["end_utc"]),
                dispersion,
                saved[station]["settings"]["depth_m"],
                (0.02, 0.085),
            )
            b01w += [
                dict(
                    event_id=event["id"],
                    station=station,
                    frequency_hz=float(f[i]),
                    ratio_over_theory=float(ratio[i] / theory[i]),
                    coherence_squared=float(coherence[i]),
                    phase_deg=float(phase[i]),
                )
                for i in range(len(f))
            ]
    write_csv(args.output_dir / "b01w_boundary_events_spectra.csv", b01w)

    # 3. Timing-guard exclusions across the workflow.
    spreads = []
    for result in saved.values():
        for event in result["events"]:
            review = event.get("timing_review", {})
            if "channels" in review:
                spreads.append(
                    max(c["offset_from_earliest_seconds"] for c in review["channels"])
                )
    summary["timing"] = {
        "between_channel_exclusions": len(spreads),
        "spread_seconds_min_median_max": [
            float(np.min(spreads)),
            float(np.median(spreads)),
            float(np.max(spreads)),
        ]
        if spreads
        else None,
        "max_phase_deg_at_band_max": float(
            360 * workflow["options"]["band"][1] * max(spreads)
        )
        if spreads
        else None,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(args.output_dir)


if __name__ == "__main__":
    main()
