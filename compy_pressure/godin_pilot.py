"""Calculate a screened, full-window Godin-style pressure spectrum from raw data."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime, read_inventory
from obspy.signal import PPSD
from scipy.signal import periodogram

from compy_pressure.archive import discover_waveforms, read_window
from compy_pressure.common import complete_samples, write_csv  # noqa: F401 (re-exported)
from compy_pressure.infragravity_test import (
    MODEL_COLOR,
    NO_TANH_COLOR,
    NO_TANH_LABEL,
    REFERENCE_URL,
    godin_pressure_db,
    positive,
)
from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    OUTPUT_DIRECTORY,
    load_master_dataset,
    update_infragravity_test_results,
)
from compy_pressure.ppsd import response_handling, select_epoch
from compy_pressure.pressure_test import (
    response_band,
    shade_excluded_band,
    validate_pressure,
)

DATA_DIRECTORY = Path.home() / "Desktop/Data_Bruit_FM/RR52"
WINDOW_SECONDS = 12.42 * 3600
DAY_SECONDS = 86400
SAVE_MAX_HZ = 0.06


def daily_windows(start, days):
    """Beginning, middle and end of each 24-hour interval, with half-open bounds."""
    if not isinstance(days, int) or days < 1:
        raise ValueError("Days must be a positive integer")
    start = UTCDateTime(start)
    offsets = (0, (DAY_SECONDS - WINDOW_SECONDS) / 2, DAY_SECONDS - WINDOW_SECONDS)
    return [
        (
            start + day * DAY_SECONDS + offset,
            start + day * DAY_SECONDS + offset + WINDOW_SECONDS,
        )
        for day in range(days)
        for offset in offsets
    ]


def count_spectrum(data, sample_rate, maximum_hz=SAVE_MAX_HZ):
    """One-sided full-length periodic Hann periodogram, normalized counts^2/Hz.

    Remove the segment mean before tapering. This explicit implementation
    choice is recorded; the source specifies Hann but not its detrending rule.
    No Welch subsegments, octave averaging or FFT padding are used.
    """
    frequency, power = periodogram(
        data,
        fs=sample_rate,
        window="hann",
        detrend="constant",
        return_onesided=True,
        scaling="density",
        nfft=len(data),
    )
    keep = (frequency > 0) & (frequency <= maximum_hz)
    return frequency[keep], power[keep]


def screen_windows(notch_power, reject_fraction=0.05):
    """Reject values strictly above the 95th percentile; ties remain retained."""
    values = np.asarray(notch_power, dtype=float)
    if not values.size or not np.all(np.isfinite(values) & (values >= 0)):
        raise ValueError("Notch powers must be finite and nonnegative")
    if not 0 <= reject_fraction < 1:
        raise ValueError("Reject fraction must satisfy 0 <= fraction < 1")
    threshold = float(np.percentile(values, 100 * (1 - reject_fraction)))
    return values <= threshold, threshold


def power_db(power):
    return 10 * np.log10(np.maximum(power, np.finfo(float).tiny))


def matched_ppsd_rows(
    ppsd,
    values,
    times,
    start,
    end,
    frequency,
    pilot_db,
    retained,
    depth,
    band,
    model_min,
    model_max,
):
    """Compare the same interval and the exact existing PPSD smoothing masks.

    Hann spectra are averaged in dB within each PPSD frequency bin, then
    converted back to linear power and averaged across retained windows.
    This is a separate matched-estimator diagnostic, not the raw pilot mean.
    """
    selected = np.array([t >= start and t + ppsd.ppsd_length <= end for t in times])
    if not selected.any():
        return []
    rows = []
    for index, center in enumerate(1 / ppsd.period_bin_centers):
        lower = float(1 / ppsd.period_bin_right_edges[index])
        upper = float(1 / ppsd.period_bin_left_edges[index])
        if lower < model_min or upper > model_max:
            continue
        mask = (frequency >= lower) & (frequency <= upper)
        if not mask.any():
            continue
        hann = float(
            power_db(np.mean(10 ** (pilot_db[retained][:, mask].mean(axis=1) / 10)))
        )
        previous = float(power_db(np.mean(10 ** (values[selected, index] / 10))))
        rows.append(
            {
                "frequency_hz": float(center),
                "bin_min_hz": lower,
                "bin_max_hz": upper,
                "hann_matched_db_pa2_hz": hann,
                "same_week_ppsd_mean_db_pa2_hz": previous,
                "hann_minus_ppsd_db": hann - previous,
                "model_mean_db_pa2_hz": float(
                    godin_pressure_db(frequency[mask], depth).mean()
                ),
                "comparison": "compared"
                if lower >= band["minimum_hz"] and upper <= band["maximum_hz"]
                else "excluded_sensor_band",
                "ppsd_windows": int(selected.sum()),
            }
        )
    return sorted(rows, key=lambda r: r["frequency_hz"])


def pilot_plots(
    seed_id,
    frequency,
    corrected,
    retained,
    window_rows,
    threshold,
    depth,
    band,
    model_min,
    model_max,
    matched,
):
    """Pressure/model, raw-count screening and matched PPSD/residual figures."""
    keep = (frequency >= model_min) & (frequency <= model_max)
    f = frequency[keep]
    x = f * 1000
    model = godin_pressure_db(f, depth)
    mean = power_db(corrected[retained].mean(axis=0))[keep]
    unscreened = power_db(corrected.mean(axis=0))[keep]
    percentiles = np.percentile(power_db(corrected[retained][:, keep]), [5, 95], axis=0)
    figures = []
    figure, axes = plt.subplots(
        2, 1, figsize=(11, 8), sharex=True, layout="constrained"
    )
    axes[0].fill_between(
        x,
        *percentiles,
        color="#0072B2",
        alpha=0.18,
        label="Retained window 5th–95th spread",
    )
    axes[0].plot(
        x, unscreened, color="#999999", linewidth=1, label="All-window linear mean"
    )
    axes[0].plot(x, mean, color="#0072B2", linewidth=1.5, label="Retained linear mean")
    axes[0].plot(x, model, color=MODEL_COLOR, linewidth=2, label="Godin model")
    axes[0].plot(
        x,
        godin_pressure_db(f, depth, tanh_factor=False),
        color=NO_TANH_COLOR,
        linewidth=2,
        linestyle=":",
        label=NO_TANH_LABEL,
    )
    axes[0].set_ylabel(r"Pressure PSD [dB re 1 Pa$^2$/Hz]")
    axes[0].legend(fontsize=8)
    axes[1].plot(
        x,
        mean - model,
        color="#0072B2",
        linewidth=1.2,
        label="Retained mean minus model",
    )
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set_ylabel("Mean minus model [dB]")
    axes[1].set_xlabel("Frequency [mHz]")
    for axis in axes:
        axis.set_xscale("log")
        shade_excluded_band(
            axis,
            model_min * 1000,
            model_max * 1000,
            {
                "minimum_hz": band["minimum_hz"] * 1000,
                "maximum_hz": band["maximum_hz"] * 1000,
            },
        )
        axis.grid(which="both", alpha=0.2)
    figure.suptitle(
        f"{seed_id}: full 12.42-hour Hann spectra\n"
        f"{retained.sum()}/{len(retained)} windows retained; depth {depth:g} m "
        "(unverified), gray outside selected sensor band",
        fontsize=11,
    )
    figures.append(figure)

    figure, (timeline, distribution) = plt.subplots(
        1, 2, figsize=(12, 5.5), layout="constrained"
    )
    valid = [r for r in window_rows if r["notch_mean_counts2_hz"] is not None]
    scores = np.array([r["notch_mean_counts2_hz"] for r in valid])
    dates = [UTCDateTime(r["actual_start"]).datetime for r in valid]
    for selection, color, label in [
        (retained, "#0072B2", "Retained"),
        (~retained, "#D55E00", "Rejected"),
    ]:
        timeline.scatter(
            np.array(dates)[selection],
            power_db(scores[selection]),
            color=color,
            label=label,
        )
    timeline.axhline(
        float(power_db(threshold)),
        color="black",
        linestyle="--",
        label="95th percentile",
    )
    timeline.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    timeline.set_xlabel("Window start (UTC)")
    timeline.set_ylabel(r"50–60 mHz mean PSD [dB re 1 count$^2$/Hz]")
    timeline.legend(fontsize=8)
    distribution.hist(
        power_db(scores), bins=min(12, len(scores)), color="#0072B2", alpha=0.7
    )
    distribution.axvline(
        float(power_db(threshold)),
        color="#D55E00",
        linestyle="--",
        label="95th percentile threshold",
    )
    distribution.set_xlabel(r"50–60 mHz mean PSD [dB re 1 count$^2$/Hz]")
    distribution.set_ylabel("Number of windows")
    distribution.legend(fontsize=8)
    for axis in (timeline, distribution):
        axis.grid(alpha=0.2)
    figure.suptitle(
        f"{seed_id}: transient screening before response correction\n"
        "Highest notch powers excluded; threshold ties retained",
        fontsize=11,
    )
    figures.append(figure)

    figure, axes = plt.subplots(
        2, 1, figsize=(11, 7), sharex=True, layout="constrained"
    )
    if matched:
        x = np.array([r["frequency_hz"] for r in matched]) * 1000
        for field, color, label in [
            ("hann_matched_db_pa2_hz", "#0072B2", "Screened Hann, matched PPSD bins"),
            (
                "same_week_ppsd_mean_db_pa2_hz",
                "#009E73",
                "Same-week 3-hour PPSD, linear mean",
            ),
            ("model_mean_db_pa2_hz", MODEL_COLOR, "Model, matched Hann bins"),
        ]:
            axes[0].plot(x, [r[field] for r in matched], color=color, label=label)
        axes[1].plot(x, [r["hann_minus_ppsd_db"] for r in matched], color="#0072B2")
        axes[0].legend(fontsize=8)
    else:
        axes[0].text(
            0.5,
            0.5,
            "No matching PPSD comparison bins",
            ha="center",
            transform=axes[0].transAxes,
        )
    axes[0].set_ylabel(r"Pressure PSD [dB re 1 Pa$^2$/Hz]")
    axes[1].set_ylabel("Hann minus PPSD [dB]")
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Frequency [mHz]")
    for axis in axes:
        axis.set_xscale("log")
        shade_excluded_band(
            axis,
            model_min * 1000,
            model_max * 1000,
            {
                "minimum_hz": band["minimum_hz"] * 1000,
                "maximum_hz": band["maximum_hz"] * 1000,
            },
        )
        axis.grid(which="both", alpha=0.2)
    figure.suptitle(
        f"{seed_id}: same-interval estimator comparison\n"
        "dB averaged within octave bins; linear power averaged across windows",
        fontsize=11,
    )
    figures.append(figure)
    return figures


def run_godin_pilot(
    master_dir=OUTPUT_DIRECTORY,
    *,
    data_dir=DATA_DIRECTORY,
    network="YV",
    station="RR52",
    start=None,
    days=7,
    depth_m=None,
    sensor_min=None,
    sensor_max=None,
    model_min=0.001,
    model_max=0.02,
    update_master=True,
    show=True,
    selected_rows=None,
    records=None,
    return_updates=False,
):
    """Process one pressure epoch with bounded reads and separate saved products."""
    master = load_master_dataset(master_dir) if selected_rows is None else selected_rows
    if not {"ppsd_units", "ppsd_status", "ppsd_npz_file"}.issubset(master):
        raise ValueError("Run ppsd.py first to populate saved pressure results")
    rows = master[
        master.Seed_id.str.startswith(f"{network}.{station}.")
        & master.get("ppsd_units", "").eq("Pa^2/Hz")
        & master.get("ppsd_status", "").isin(["complete", "partial"])
    ].drop_duplicates(EPOCH_COLUMNS + ["inventory_file"])
    if len(rows) != 1:
        raise ValueError("Pilot requires exactly one saved pressure channel epoch")
    row = rows.iloc[0]
    source = Path(row.ppsd_npz_file).expanduser().resolve()
    ppsd = PPSD.load_npz(str(source))
    previous_values = validate_pressure(ppsd, row)
    _, channel = select_epoch(row)
    if response_handling(channel)[1] != "Pa^2/Hz":
        raise ValueError("Pilot requires a pressure-to-counts response")
    snapshot = source.parent / "inventory.xml"
    saved_response = read_inventory(str(snapshot)).get_response(
        ppsd.id, ppsd.times_processed[0]
    )
    band, _, _ = response_band(
        ppsd, snapshot, sensor_min=sensor_min, sensor_max=sensor_max
    )
    depth = positive(depth_m if depth_m is not None else row.Profondeur, "Water depth")
    model_min, model_max = (
        positive(model_min, "Minimum model frequency"),
        positive(model_max, "Maximum model frequency"),
    )
    if not model_min < model_max <= SAVE_MAX_HZ:
        raise ValueError("Model range must satisfy 0 < min < max <= 0.06 Hz")
    if records is None:
        records = discover_waveforms(data_dir, {ppsd.id})
    records = [
        r
        for r in records
        if r.seed_id == ppsd.id
        and (channel.start_date is None or r.end >= channel.start_date)
        and (channel.end_date is None or r.start < channel.end_date)
    ]
    if not records:
        raise ValueError("No waveform data in selected pressure response epoch")
    earliest = max(
        min(r.start for r in records),
        channel.start_date or min(r.start for r in records),
    )
    latest = min(
        max(r.end + 1 / r.sample_rate for r in records),
        channel.end_date or max(r.end + 1 / r.sample_rate for r in records),
    )
    if not isinstance(days, int) or days < 1:
        raise ValueError("Days must be a positive integer")
    if start is None:
        # Center the pilot within coverage, avoiding deployment startup data.
        middle_start = earliest + (latest - earliest - days * DAY_SECONDS) / 2
        start = UTCDateTime(middle_start.date)
        if start < earliest or start + days * DAY_SECONDS > latest:
            start = middle_start
    start = UTCDateTime(start)
    windows = daily_windows(start, days)
    end = start + days * DAY_SECONDS
    if start < earliest or end > latest:
        raise ValueError("Requested pilot interval is outside waveform coverage")
    if channel.end_date is not None and end > channel.end_date:
        raise ValueError("Pilot interval crosses the response epoch")
    sample_rate = float(channel.sample_rate)
    if sample_rate != ppsd.sampling_rate:
        raise ValueError("Saved PPSD sample rate differs from selected inventory epoch")
    npts = round(WINDOW_SECONDS * sample_rate)
    if not np.isclose(npts / sample_rate, WINDOW_SECONDS, rtol=0, atol=1e-8):
        raise ValueError("12.42 hours must contain an integer number of samples")
    if sample_rate / 2 <= SAVE_MAX_HZ:
        raise ValueError("Sampling rate cannot cover the screening notch")
    used_paths = sorted(
        {r.path for r in records if r.start < end + 2 / sample_rate and r.end >= start}
    )
    source_metadata = [
        {"path": str(p), "size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
        for p in used_paths
    ]
    spectra, window_rows = [], []
    frequency = np.fft.rfftfreq(npts, 1 / sample_rate)
    frequency = frequency[(frequency > 0) & (frequency <= SAVE_MAX_HZ)]
    notch = (frequency >= 0.05) & (frequency <= 0.06)
    if not notch.any():
        raise ValueError("FFT cannot resolve screening notch")
    response = channel.response.get_evalresp_response_for_frequencies(
        frequency, output="DEF"
    )
    saved = saved_response.get_evalresp_response_for_frequencies(
        frequency, output="DEF"
    )
    if not np.allclose(response, saved, rtol=1e-10, atol=0):
        raise ValueError("Source response differs from saved PPSD response")
    response_power = np.abs(response) ** 2
    if not np.all(np.isfinite(response_power) & (response_power > 0)):
        raise ValueError("Pressure response is zero or nonfinite")
    for index, (window_start, window_end) in enumerate(windows):
        stream = read_window(records, window_start, window_end + 2 / sample_rate)
        data, actual_start = complete_samples(stream, window_start, sample_rate, npts)
        item = {
            "requested_start": str(window_start),
            "requested_end": str(window_end),
            "actual_start": str(actual_start) if actual_start is not None else "",
            "actual_end_exclusive": str(actual_start + npts / sample_rate)
            if actual_start is not None
            else "",
            "notch_mean_counts2_hz": None,
            "decision": "missing_gapped_or_flat",
        }
        if data is not None:
            if (
                channel.end_date is not None
                and actual_start + npts / sample_rate > channel.end_date
            ):
                raise ValueError("Aligned samples cross response epoch")
            _, power = count_spectrum(data, sample_rate)
            item["notch_mean_counts2_hz"] = float(power[notch].mean())
            item["decision"] = "pending_screening"
            spectra.append(power)
        window_rows.append(item)
        print(
            f"{ppsd.id}: window {index + 1}/{len(windows)} {item['decision']}",
            flush=True,
        )
    if not spectra:
        raise ValueError("No complete continuous windows in selected interval")
    for item in source_metadata:
        current = Path(item["path"]).stat()
        if (current.st_size, current.st_mtime_ns) != (item["size"], item["mtime_ns"]):
            raise ValueError("Source waveform changed during pilot processing")
    counts = np.array(spectra)
    valid = [item for item in window_rows if item["notch_mean_counts2_hz"] is not None]
    retained, threshold = screen_windows(
        [item["notch_mean_counts2_hz"] for item in valid]
    )
    for item, accepted in zip(valid, retained, strict=True):
        item["decision"] = "retained" if accepted else "rejected_high_notch_power"
    corrected = counts / response_power
    mean = corrected[retained].mean(axis=0)
    model = godin_pressure_db(frequency, depth)
    eligible = (frequency >= max(model_min, band["minimum_hz"])) & (
        frequency <= min(model_max, band["maximum_hz"])
    )
    residual = power_db(mean) - model
    matched = matched_ppsd_rows(
        ppsd,
        previous_values,
        ppsd.times_processed,
        start,
        end,
        frequency,
        power_db(corrected),
        retained,
        depth,
        band,
        model_min,
        model_max,
    )
    output = (
        source.parent / "godin_pilot" / f"{start.strftime('%Y%m%dT%H%M%S')}__{days}days"
    )
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output / "spectra.npz",
        frequency_hz=frequency,
        counts_psd=counts,
        pressure_psd=corrected,
        retained=retained,
        response_power_counts2_pa2=response_power,
        window_start_utc=np.array([r["actual_start"] for r in valid]),
        mean_retained_pa2_hz=mean,
        mean_all_pa2_hz=corrected.mean(axis=0),
    )
    write_csv(output / "windows.csv", window_rows)
    write_csv(
        output / "ppsd_comparison.csv",
        matched,
        fieldnames=[
            "frequency_hz",
            "bin_min_hz",
            "bin_max_hz",
            "hann_matched_db_pa2_hz",
            "same_week_ppsd_mean_db_pa2_hz",
            "hann_minus_ppsd_db",
            "model_mean_db_pa2_hz",
            "comparison",
            "ppsd_windows",
        ],
    )
    write_csv(
        output / "frequency_comparison.csv",
        [
            {
                "frequency_hz": float(f),
                "mean_retained_pa2_hz": float(p),
                "mean_retained_db_pa2_hz": float(power_db(p)),
                "model_db_pa2_hz": float(m),
                "mean_minus_model_db": float(r),
                "comparison": "compared" if ok else "excluded_band",
            }
            for f, p, m, r, ok in zip(
                frequency, mean, model, residual, eligible, strict=True
            )
        ],
    )
    figures = pilot_plots(
        ppsd.id,
        frequency,
        corrected,
        retained,
        window_rows,
        threshold,
        depth,
        band,
        model_min,
        model_max,
        matched,
    )
    plot_names = ("pressure_model", "screening", "ppsd_comparison")
    for figure, name in zip(figures, plot_names, strict=True):
        figure.savefig(output / f"{name}.png", dpi=160)
    summary = {
        "seed_id": ppsd.id,
        "channel_start_time": row.channel_start_time,
        "channel_end_time": row.channel_end_time,
        "status": "comparison_only" if eligible.any() else "no_comparable_bins",
        "score": None,
        "reference_url": REFERENCE_URL,
        "start_utc": str(start),
        "end_utc_exclusive": str(end),
        "days": days,
        "window_seconds": WINDOW_SECONDS,
        "nfft": npts,
        "frequency_resolution_hz": sample_rate / npts,
        "hann_equivalent_noise_bandwidth_hz": 1.5 * sample_rate / npts,
        "sample_rate_hz": sample_rate,
        "windows_per_day": 3,
        "detrend": "constant",
        "window": "periodic_hann",
        "fft_padding": False,
        "averaging": "arithmetic mean in linear power after counts-based screening",
        "saved_frequency_max_hz": SAVE_MAX_HZ,
        "units": "Pa^2/Hz",
        "requested_windows": len(windows),
        "complete_windows": len(valid),
        "missing_windows": len(windows) - len(valid),
        "retained_windows": int(retained.sum()),
        "rejected_windows": int((~retained).sum()),
        "screening_band_hz": [0.05, 0.06],
        "screening_metric": "mean linear counts^2/Hz",
        "screening_rule": "strictly above NumPy linear-interpolated "
        "95th percentile; retain ties",
        "screening_threshold_counts2_hz": threshold,
        "actual_rejected_fraction": float((~retained).mean()),
        "depth_m": depth,
        "depth_basis": "user_supplied"
        if depth_m is not None
        else "inventory_elevation_estimate",
        "depth_independently_verified": False,
        "sensor_band": band,
        "model_band_hz": [model_min, model_max],
        "compared_fft_bins": int(eligible.sum()),
        "median_mean_minus_model_db": float(np.median(residual[eligible]))
        if eligible.any()
        else None,
        "same_interval_ppsd_windows": sum(
            t >= start and t + ppsd.ppsd_length <= end for t in ppsd.times_processed
        ),
        "response_snapshot_sha256": hashlib.sha256(snapshot.read_bytes()).hexdigest(),
        "source_ppsd_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_files": source_metadata,
        "files": {
            name: str(output / name)
            for name in [
                "spectra.npz",
                "windows.csv",
                "frequency_comparison.csv",
                "ppsd_comparison.csv",
                *[f"{name}.png" for name in plot_names],
            ]
        },
        "limitations": [
            "Short pilot differs from the original yearlong mean; "
            "no calibration acceptance test is assigned.",
            "Depth, instrument band and regional model applicability "
            "remain unverified.",
            "Overlapping Hann windows and adjacent FFT bins are correlated.",
            "Finite-sample screening fraction can differ from 5%; "
            "ties remain retained.",
            "Detrending and periodic-Hann convention are explicit "
            "implementation choices.",
            "Same-week comparison changes window length, taper, "
            "internal averaging and screening together.",
        ],
    }
    json_path = output / "comparison.json"
    json_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    update = row[EPOCH_COLUMNS + ["inventory_file"]].to_dict() | {
        "infragravity_pilot_status": summary["status"],
        "infragravity_pilot_result_file": str(json_path),
        "infragravity_pilot_npz_file": str(output / "spectra.npz"),
        "infragravity_pilot_start_utc": str(start),
        "infragravity_pilot_end_utc": str(end),
        "infragravity_pilot_window_seconds": WINDOW_SECONDS,
        "infragravity_pilot_retained_windows": int(retained.sum()),
        "infragravity_pilot_rejected_windows": int((~retained).sum()),
        "infragravity_pilot_compared_fft_bins": int(eligible.sum()),
        "infragravity_pilot_median_mean_minus_model_db": summary[
            "median_mean_minus_model_db"
        ],
        **{
            f"infragravity_pilot_{name}_file": str(output / f"{name}.png")
            for name in plot_names
        },
    }
    if update_master:
        update_infragravity_test_results(master_dir, [update])
    print(
        f"Saved {output}; retained {retained.sum()}/{len(retained)} windows; "
        f"median mean-minus-model {summary['median_mean_minus_model_db']} dB",
        flush=True,
    )
    if show:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)
    return (summary, [update]) if return_updates else summary


def show_saved_pilot(
    master_dir=OUTPUT_DIRECTORY, network="YV", station="RR52", *, show=True
):
    """Open saved PNGs without recalculating spectra or changing master files."""
    master = load_master_dataset(master_dir)
    if "infragravity_pilot_result_file" not in master:
        raise ValueError("Run godin_pilot.py first to save pilot results")
    rows = master[
        master.Seed_id.str.startswith(f"{network}.{station}.")
        & master.infragravity_pilot_result_file.ne("")
    ]
    if len(rows) != 1:
        raise ValueError("Expected exactly one saved pilot for this station")
    path = Path(rows.iloc[0].infragravity_pilot_result_file)
    summary = json.loads(path.read_text())
    figures = []
    for name in ("pressure_model.png", "screening.png", "ppsd_comparison.png"):
        image_path = Path(summary["files"][name])
        picture = plt.imread(image_path)
        figure, axis = plt.subplots(figsize=(12, 8))
        axis.imshow(picture)
        axis.axis("off")
        figure.subplots_adjust(left=0, right=1, top=1, bottom=0)
        figures.append(figure)
        print(f"Opened {image_path}", flush=True)
    if show:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-dir", type=Path, default=OUTPUT_DIRECTORY)
    parser.add_argument("--data-dir", type=Path, default=DATA_DIRECTORY)
    parser.add_argument("--network", default="YV")
    parser.add_argument("--station", default="RR52")
    parser.add_argument(
        "--start", help="First interval start in UTC; default a week near mid-recording"
    )
    parser.add_argument("--days", type=int, default=7)
    parser.add_argument("--depth-m", type=float)
    parser.add_argument("--sensor-min-frequency", type=float)
    parser.add_argument("--sensor-max-frequency", type=float)
    parser.add_argument("--model-min-frequency", type=float, default=0.001)
    parser.add_argument("--model-max-frequency", type=float, default=0.02)
    parser.add_argument("--no-update-master", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Open latest saved pilot plots without waveform processing",
    )
    args = parser.parse_args()
    if args.plot_only:
        show_saved_pilot(
            args.master_dir, args.network, args.station, show=not args.no_show
        )
        return
    run_godin_pilot(
        args.master_dir,
        data_dir=args.data_dir,
        network=args.network,
        station=args.station,
        start=args.start,
        days=args.days,
        depth_m=args.depth_m,
        sensor_min=args.sensor_min_frequency,
        sensor_max=args.sensor_max_frequency,
        model_min=args.model_min_frequency,
        model_max=args.model_max_frequency,
        update_master=not args.no_update_master,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()
