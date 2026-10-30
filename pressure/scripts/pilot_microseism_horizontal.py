"""Bounded horizontal-coupling pilot, bound to one week in a saved phase audit.

No physical tilt estimate, cleaned trace, fitted gain or score is produced.
The two diagnostic intervals must be inside one saved week; at most 24 h each.
"""

import argparse
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
from obspy import UTCDateTime

from compy_pressure.archive import SampleGridError, index_waveform_files, read_window
from compy_pressure.godin_pilot import complete_samples, write_csv
from compy_pressure.microseism_check import check_paired_run
from compy_pressure.microseism_horizontal import conditional_products, corrected_fourier
from compy_pressure.paired_spectra import checked_response, response_snapshot_path
from compy_pressure.ppsd import atomic_json, select_epoch
from scripts.combine_microseism_weeks import InputGuard

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

ROOT = Path(__file__).resolve().parents[1]


def median(values):
    values = np.asarray(values)
    return (
        float(np.median(values[np.isfinite(values)]))
        if np.isfinite(values).any()
        else None
    )


def window_rows(decisions):
    """Keep a fixed CSV schema for successful and refused windows alike."""
    return [
        {key: d.get(key, "") for key in ("start", "status", "actual_start", "reason")}
        for d in decisions
    ]


def run(
    audit_path,
    master_csv,
    data_dir,
    station,
    shift,
    starts,
    hours,
    horizontal_ids,
    output_root,
):
    if not 6 <= hours <= 24 or len(starts) not in (1, 2):
        raise ValueError("Pilot requires one/two intervals, each 6..24 hourly windows")
    guard = InputGuard()
    compact = guard.json(audit_path)
    full_path = guard.protect(compact["full_audit_path"], compact["full_audit_sha256"])
    audit = guard.json(full_path)
    for path, expected in audit["input_sha256"].items():
        guard.protect(path, expected)
    weeks = [
        w
        for w in audit["weeks"]
        if w["station"] == station and w["shift_weeks"] == shift
    ]
    if len(weeks) != 1:
        raise ValueError("Select exactly one audited passing week")
    week = weeks[0]
    pair_path = guard.protect(week["paired_comparison"])
    paired, saved = check_paired_run(pair_path.parent)
    master = pd.read_csv(guard.protect(master_csv), keep_default_na=False, dtype=str)
    ids = [
        paired["pressure_epoch"]["Seed_id"],
        *horizontal_ids,
        paired["vertical_seed_id"],
    ]
    if len(set(ids)) != 4 or any(
        ".".join(seed.split(".")[:2]) != station for seed in ids
    ):
        raise ValueError(
            "Four distinct channel IDs must belong to the selected station"
        )
    left, right = map(
        UTCDateTime,
        paired["interval_utc"]
        if "interval_utc" in paired
        else [paired["start_utc"], paired["end_utc"]],
    )
    starts = list(map(UTCDateTime, starts))
    if any(t < left or t + hours * 3600 > right for t in starts):
        raise ValueError("Pilot intervals must be contained in the saved week")
    intervals = sorted((t, t + hours * 3600) for t in starts)
    if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
        raise ValueError("Pilot intervals must not overlap")
    rows, channels = [], []
    for i, seed in enumerate(ids):
        candidates = master.loc[master.Seed_id.eq(seed)]
        matches = []
        for _, row in candidates.iterrows():
            _, channel = select_epoch(row)
            if channel.start_date <= min(starts) and (
                channel.end_date is None
                or channel.end_date >= max(starts) + hours * 3600
            ):
                matches.append((row, channel))
        if len(matches) != 1:
            raise ValueError(f"Need one exact epoch covering the pilot: {seed}")
        row, channel = matches[0]
        guard.protect(row.inventory_file)
        guard.protect(response_snapshot_path(row))
        if i in (0, 3):
            expected = paired["pressure_epoch" if i == 0 else "vertical_epoch"]
            if any(str(row[k]) != str(v) for k, v in expected.items()):
                raise ValueError("Pilot P/Z epochs differ from saved pair")
        elif channel.dip is None or abs(channel.dip) > 0.1:
            raise ValueError(
                "Both supplied H channels must have horizontal metadata dip"
            )
        if i and str(
            channel.response.instrument_sensitivity.input_units
        ).upper() not in {"M", "M/S", "M/S**2", "M/S^2", "M/S/S"}:
            raise ValueError("Seismic channels need known motion units")
        rows.append(row)
        channels.append(channel)
    rates = {c.sample_rate for c in channels}
    # Statistical conditioning does not require geographical azimuth, but
    # ambiguous metadata prohibit interpreting this as a physical rotation.
    # Spectrally collinear observed channels are masked by conditional_products.
    horizontal_azimuths_distinct = (
        all(c.azimuth is not None for c in channels[1:3])
        and abs(np.sin(np.radians(channels[1].azimuth - channels[2].azimuth))) >= 0.1
    )
    if len(rates) != 1 or next(iter(rates)) != paired["sample_rate_hz"]:
        raise ValueError("Native rates must match; no implicit resampling")
    rate = rates.pop()
    npts = round(3600 * rate)
    full_f = np.fft.rfftfreq(npts, 1 / rate)
    selected = (full_f >= 0.02) & (full_f <= 0.1)
    f = full_f[selected]
    saved_indices = np.flatnonzero(
        (saved["frequency_hz"] >= 0.02) & (saved["frequency_hz"] <= 0.1)
    )
    np.testing.assert_allclose(
        f, saved["frequency_hz"][saved_indices], rtol=1e-12, atol=0
    )
    responses = np.array(
        [
            checked_response(r, c, f, "DEF" if i == 0 else "DISP")
            for i, (r, c) in enumerate(zip(rows, channels))
        ]
    )
    paths = set()
    # Restrict discovery to named daily files and neighboring days for midnight
    # crossings; index actual headers rather than trusting filename coverage.
    for t in starts:
        for day in (t - 86400, t, t + 86400):
            for seed in ids:
                paths.update(
                    Path(data_dir).rglob(f"{seed}.D.{day.year}.{day.julday:03d}")
                )
    for p in paths:
        guard.protect(p)
    records = index_waveform_files(paths)
    for item in paired["source_files"]:
        if Path(item["path"]) in paths:
            stat = Path(item["path"]).stat()
            if (stat.st_size, stat.st_mtime_ns) != (item["size"], item["mtime_ns"]):
                raise ValueError("P/Z raw source differs from saved paired metadata")
    guard.protect(__file__)
    guard.protect(ROOT / "../compy_pressure/microseism_horizontal.py")
    decisions, summaries, products = [], [], []
    for start in starts:
        coefficients = []
        retained = []
        for hour in range(hours):
            requested = start + hour * 3600
            arrays, actuals = [], []
            refusal = None
            for seed in ids:
                rr = [r for r in records if r.seed_id == seed]
                try:
                    stream = read_window(
                        rr,
                        requested,
                        requested + 3600 + 2 / rate,
                        strict_sample_grid=True,
                    )
                except SampleGridError as exc:
                    refusal = str(exc)
                    break
                values, actual = complete_samples(stream, requested, rate, npts)
                arrays.append(values)
                actuals.append(actual)
            if refusal:
                decisions.append(
                    {
                        "start": str(requested),
                        "status": "sample_grid_refused",
                        "reason": refusal,
                    }
                )
                continue
            if any(x is None for x in arrays):
                decisions.append(
                    {"start": str(requested), "status": "incomplete_or_flat"}
                )
                continue
            if max(actuals) - min(actuals) > 1e-6:
                raise ValueError("Four native sample grids do not align")
            if max(actuals) + 3600 > right:
                raise ValueError("Aligned samples cross the saved week boundary")
            x = corrected_fourier(
                arrays, responses, rate, selected, paired["vertical_upward_sign"]
            )
            matches = np.flatnonzero(saved["retained_window_starts"] == str(actuals[0]))
            if len(matches) != 1:
                raise ValueError(
                    "Pilot hour is not uniquely present in saved P/Z windows"
                )
            old = int(matches[0])
            for actual, key in (
                (abs(x[0]) ** 2, "pressure_window_psd"),
                (abs(x[3]) ** 2, "displacement_window_psd"),
                (x[3] * x[0].conj(), "cross_window_dp"),
            ):
                if not np.allclose(
                    actual, saved[key][old, saved_indices], rtol=1e-10, atol=0
                ):
                    raise ValueError(
                        "New P/Z hourly spectrum does not reproduce saved input"
                    )
            coefficients.append(x)
            retained.append(str(requested))
            decisions.append(
                {
                    "start": str(requested),
                    "status": "retained",
                    "actual_start": str(actuals[0]),
                }
            )
            print(
                f"{station} {requested}: four channels verified ({hour + 1}/{hours})",
                flush=True,
            )
        result = conditional_products(coefficients)
        products.append((str(start), np.array(coefficients), result))
        for band in week["bands"]:
            lo, hi = band["frequency_hz"]
            mask = (f >= lo) & (f < hi)
            summaries.append(
                {
                    "start_utc": str(start),
                    "hours_requested": hours,
                    "retained_windows": len(retained),
                    "band_hz": [lo, hi],
                    "bins": int(mask.sum()),
                    "conditional_valid_bins": int(
                        result["valid_conditional"][mask].sum()
                    ),
                    "median_coherence_squared_pz": median(
                        result["coherence_squared_pz"][mask]
                    ),
                    "median_p_power_predictable_from_h": median(
                        result["p_power_predictable_from_h"][mask]
                    ),
                    "median_z_power_predictable_from_h": median(
                        result["z_power_predictable_from_h"][mask]
                    ),
                    "median_partial_coherence_squared_pz_given_h": median(
                        result["partial_coherence_squared_pz_given_h"][mask]
                    ),
                    "median_residual_transfer_amplitude_ratio": median(
                        abs(result["residual_z_p_transfer_ratio"][mask])
                    ),
                }
            )
    guard.verify()
    output = (
        Path(output_root).resolve()
        / f"run_{datetime.now(UTC):%Y%m%dT%H%M%SZ}_{uuid4().hex[:8]}"
    )
    if any(output.is_relative_to(Path(p).parent) for p in guard.hashes):
        raise ValueError("Output must not be beside protected inputs")
    output.mkdir(parents=True, exist_ok=False)
    files = []
    for index, (_, x, result) in enumerate(products):
        p = output / f"interval_{index}.npz"
        np.savez_compressed(p, frequency_hz=f, window_coefficients=x, **result)
        files.append({"path": str(p), "sha256": InputGuard.digest(p)})
    write_csv(output / "windows.csv", window_rows(decisions))
    write_csv(output / "bands.csv", summaries)
    result = {
        "status": "diagnostic_only_not_a_tilt_estimate_or_pressure_calibration",
        "station": station,
        "channel_order": ids,
        "channel_metadata": [
            {
                "seed_id": seed,
                "dip_deg": c.dip,
                "azimuth_deg": c.azimuth,
                "sample_rate_hz": c.sample_rate,
            }
            for seed, c in zip(ids, channels)
        ],
        "horizontal_azimuths_distinct_in_metadata": bool(horizontal_azimuths_distinct),
        "source_week": week["interval_utc"],
        "starts_utc": list(map(str, starts)),
        "hours_per_interval": hours,
        "segment_seconds": 3600,
        "taper": "periodic Hann; segment mean removed; one-sided density",
        "convention": "S_ij = mean(X_i * conjugate(X_j)); P,H1,H2,Z",
        "horizontal_condition_limit": 1e6,
        "saved_pz_spectra_reproduced": True,
        "bands": summaries,
        "files": files,
        "input_sha256": guard.hashes,
        "protected_inputs_unchanged": True,
        "gain_applied": False,
        "master_or_scores_changed": False,
        "limitations": [
            "Targeted diagnostic intervals, not an independent validation set",
            "Few hourly averages; neighboring bins correlated; no confidence intervals",
            "H projection also removes real Rayleigh motion; cannot identify tilt",
            "Response, orientation, clock and water/Earth models remain conditional",
            "Missing/duplicate H azimuths prohibit physical rotation inference",
            "No cleaned waveform, gain correction or scientific score is produced",
        ],
    }
    atomic_json(output / "pilot.json", result)
    guard.verify()
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gain-audit", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--station", required=True)
    parser.add_argument("--shift", type=int, default=0)
    parser.add_argument("--start", action="append", required=True)
    parser.add_argument("--hours", type=int, default=6)
    parser.add_argument("--horizontal-ids", nargs=2, required=True)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "results/REVIEWS/MICROSEISM_HORIZONTAL",
    )
    args = parser.parse_args()
    print(
        run(
            args.gain_audit,
            args.master_csv,
            args.data_dir,
            args.station,
            args.shift,
            args.start,
            args.hours,
            args.horizontal_ids,
            args.output_root,
        )
    )


if __name__ == "__main__":
    main()
