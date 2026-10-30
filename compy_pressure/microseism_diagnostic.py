"""Model-free microseism pressure/vertical diagnostic from saved paired spectra."""

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory

from compy_pressure.phase import (
    PHASE_CODE_SHA256,
    PHASE_SUMMARY_METHOD,
    circular_median_degrees,
)

MICROSEISM_BAND_HZ = (0.1, 0.2)
RUAN_COHERENCE_REFERENCE = 0.95
SAIKIA_DAILY_COHERENCE_REFERENCE = 0.8
RR52_FIXED_PRESSURE_MULTIPLIER = 0.6787278140712512
REQUIRED_ARRAYS = {
    "frequency_hz",
    "pressure_psd",
    "displacement_psd",
    "cross_dp",
    "transfer_dp",
    "coherence_squared",
    "valid",
    "pressure_window_psd",
    "displacement_window_psd",
    "cross_window_dp",
    "retained_window_starts",
}


def sha256(path):
    """Return a streaming SHA-256 without modifying the input."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_and_check_spectra(path, expected_segments=None):
    """Load saved paired spectra and reproduce their aggregate products."""
    path = Path(path)
    with np.load(path) as saved:
        missing = REQUIRED_ARRAYS - set(saved.files)
        if missing:
            raise ValueError(f"Saved spectra lack {sorted(missing)}")
        arrays = {name: saved[name].copy() for name in saved.files}

    f = arrays["frequency_hz"]
    valid = arrays["valid"].astype(bool)
    ppw = arrays["pressure_window_psd"]
    ddw = arrays["displacement_window_psd"]
    dpw = arrays["cross_window_dp"]
    starts = arrays["retained_window_starts"]
    if f.ndim != 1 or len(f) < 2 or not np.all(np.diff(f) > 0):
        raise ValueError("Frequency grid must be one-dimensional and increasing")
    if ppw.shape != ddw.shape or ppw.shape != dpw.shape:
        raise ValueError("Saved window auto/cross spectra have different shapes")
    if ppw.ndim != 2 or ppw.shape[1] != len(f) or len(starts) != ppw.shape[0]:
        raise ValueError("Saved window spectra do not match frequency/time axes")
    if expected_segments is not None and len(starts) != expected_segments:
        raise ValueError("Retained segment count differs from comparison metadata")
    if valid.shape != f.shape:
        raise ValueError("Validity mask does not match the frequency grid")

    pp = np.mean(ppw, axis=0)
    dd = np.mean(ddw, axis=0)
    dp = np.mean(dpw, axis=0)
    calculated_valid = (
        np.isfinite(pp) & np.isfinite(dd) & np.isfinite(dp) & (pp > 0) & (dd > 0)
    )
    h = np.divide(dp, pp, out=np.zeros_like(dp), where=calculated_valid)
    coherence = np.divide(
        abs(dp) ** 2,
        pp * dd,
        out=np.zeros_like(pp),
        where=calculated_valid,
    )
    reproduced = {
        "pressure_psd": pp,
        "displacement_psd": dd,
        "cross_dp": dp,
        "transfer_dp": h,
        "coherence_squared": np.clip(coherence, 0, 1),
    }
    for name, calculated in reproduced.items():
        if not np.allclose(arrays[name], calculated, rtol=1e-12, atol=0):
            raise ValueError(f"Saved {name} does not reproduce from window spectra")
    if np.any(valid & ~calculated_valid):
        raise ValueError("Saved validity mask includes nonfinite/nonpositive spectra")
    return arrays


def daily_products(arrays, band=MICROSEISM_BAND_HZ):
    """Re-average the saved hourly spectra by UTC day, without new waveforms."""
    f = arrays["frequency_hz"]
    valid = arrays["valid"].astype(bool)
    band_mask = valid & (f >= band[0]) & (f <= band[1])
    if not band_mask.any():
        raise ValueError("No valid frequencies in the requested microseism band")
    days = np.array([str(value)[:10] for value in arrays["retained_window_starts"]])
    products = []
    for day in sorted(set(days)):
        selected = days == day
        pp = np.mean(arrays["pressure_window_psd"][selected], axis=0)
        dd = np.mean(arrays["displacement_window_psd"][selected], axis=0)
        dp = np.mean(arrays["cross_window_dp"][selected], axis=0)
        good = valid & np.isfinite(pp) & np.isfinite(dd) & (pp > 0) & (dd > 0)
        h = np.divide(dp, pp, out=np.zeros_like(dp), where=good)
        coherence = np.divide(
            abs(dp) ** 2,
            pp * dd,
            out=np.zeros_like(pp),
            where=good,
        )
        coherence = np.clip(coherence, 0, 1)
        average = float(np.mean(coherence[band_mask]))
        products.append(
            {
                "date_utc": day,
                "segments": int(selected.sum()),
                "transfer": h,
                "coherence_squared": coherence,
                "band_mean_coherence_squared": average,
                "band_median_coherence_squared": float(np.median(coherence[band_mask])),
                "band_median_transfer_amplitude_m_pa": float(
                    np.median(abs(h[band_mask]))
                ),
                "passes_saikia_reference": bool(
                    average > SAIKIA_DAILY_COHERENCE_REFERENCE
                ),
            }
        )
    return products


def diagnostic_metrics(
    arrays,
    daily,
    pressure_multiplier=RR52_FIXED_PRESSURE_MULTIPLIER,
    band=MICROSEISM_BAND_HZ,
):
    """Describe nominal and fixed-reference transfer functions; never fit a gain."""
    if not np.isfinite(pressure_multiplier) or pressure_multiplier <= 0:
        raise ValueError("Fixed pressure multiplier must be finite and positive")
    f = arrays["frequency_hz"]
    h = arrays["transfer_dp"]
    coherence = arrays["coherence_squared"]
    valid = arrays["valid"].astype(bool)
    band_mask = valid & (f >= band[0]) & (f <= band[1])
    if not band_mask.any():
        raise ValueError("No valid frequencies in the requested microseism band")
    ruan = band_mask & (coherence > RUAN_COHERENCE_REFERENCE)
    selected_days = [item for item in daily if item["passes_saikia_reference"]]
    selected_stack = (
        np.mean([item["transfer"] for item in selected_days], axis=0)
        if selected_days
        else np.full(h.shape, np.nan + 0j)
    )
    phase = np.angle(h, deg=True)
    fixed_h = h / pressure_multiplier
    return {
        "band_mask": band_mask,
        "ruan_mask": ruan,
        "fixed_transfer": fixed_h,
        "selected_daily_stack": selected_stack,
        "summary": {
            "frequency_band_hz": list(band),
            "valid_bins": int(band_mask.sum()),
            "median_coherence_squared": float(np.median(coherence[band_mask])),
            "median_nominal_transfer_amplitude_m_pa": float(
                np.median(abs(h[band_mask]))
            ),
            "median_nominal_transfer_phase_deg": circular_median_degrees(
                phase[band_mask]
            ),
            "phase_summary_method": PHASE_SUMMARY_METHOD,
            "phase_summary_code_sha256": PHASE_CODE_SHA256,
            "fixed_pressure_multiplier": float(pressure_multiplier),
            "fixed_reference_relation": "H_fixed = H_nominal / g_pressure",
            "median_fixed_reference_transfer_amplitude_m_pa": float(
                np.median(abs(fixed_h[band_mask]))
            ),
            "ruan_high_coherence_reference": {
                "coherence_squared_greater_than": RUAN_COHERENCE_REFERENCE,
                "bins": int(ruan.sum()),
                "fraction_of_band": float(ruan.sum() / band_mask.sum()),
                "frequency_min_hz": float(f[ruan].min()) if ruan.any() else None,
                "frequency_max_hz": float(f[ruan].max()) if ruan.any() else None,
                "median_nominal_transfer_amplitude_m_pa": (
                    float(np.median(abs(h[ruan]))) if ruan.any() else None
                ),
                "median_fixed_reference_transfer_amplitude_m_pa": (
                    float(np.median(abs(fixed_h[ruan]))) if ruan.any() else None
                ),
            },
            "saikia_daily_reference": {
                "band_mean_coherence_squared_greater_than": (
                    SAIKIA_DAILY_COHERENCE_REFERENCE
                ),
                "selected_days": len(selected_days),
                "available_days": len(daily),
                "selected_dates_utc": [item["date_utc"] for item in selected_days],
                "selected_stack_median_nominal_transfer_amplitude_m_pa": (
                    float(np.median(abs(selected_stack[band_mask])))
                    if selected_days
                    else None
                ),
                "selected_stack_median_fixed_reference_amplitude_m_pa": (
                    float(
                        np.median(abs(selected_stack[band_mask])) / pressure_multiplier
                    )
                    if selected_days
                    else None
                ),
            },
        },
    }


def verify_saved_provenance(summary_path, summary):
    """Verify recorded source stats/hashes and hash each saved pilot artifact."""
    problems = []
    for item in summary.get("source_files", []):
        path = Path(item["path"])
        if not path.exists():
            problems.append(f"missing source {path}")
            continue
        stat = path.stat()
        if stat.st_size != item["size"] or stat.st_mtime_ns != item["mtime_ns"]:
            problems.append(f"changed source metadata {path}")
    for group in ("response_snapshots", "inventory_sha256"):
        for name, expected in summary.get(group, {}).items():
            path = Path(name)
            if not path.exists() or sha256(path) != expected:
                problems.append(f"changed or missing protected input {path}")
    if problems:
        raise ValueError("; ".join(problems))

    paths = {"comparison": Path(summary_path)}
    for key in ("spectra", "frequency_table", "windows", "plot"):
        path = Path(summary["files"][key])
        if not path.exists():
            raise ValueError(f"Missing saved pilot artifact {path}")
        paths[key] = path
    return {
        name: {"path": str(path.resolve()), "sha256": sha256(path)}
        for name, path in paths.items()
    }


def response_metadata(summary):
    """Report response-snapshot support without claiming a field-verified band."""
    seeds = [summary["pressure_epoch"]["Seed_id"], summary["vertical_seed_id"]]
    snapshots = [Path(path) for path in summary["response_snapshots"]]
    frequencies = np.array([0.1, 0.15, 0.2])
    result = {}
    for seed in seeds:
        channel = None
        for snapshot in snapshots:
            selected = read_inventory(snapshot).select(
                network=seed.split(".")[0],
                station=seed.split(".")[1],
                location=seed.split(".")[2],
                channel=seed.split(".")[3],
            )
            if len(selected) and len(selected[0]) and len(selected[0][0]):
                channel = selected[0][0][0]
                break
        if channel is None:
            raise ValueError(f"No response snapshot contains {seed}")
        sensitivity = channel.response.instrument_sensitivity
        response = channel.response.get_evalresp_response_for_frequencies(
            frequencies, output="DEF"
        )
        reference = channel.response.get_evalresp_response_for_frequencies(
            [sensitivity.frequency], output="DEF"
        )[0]
        result[seed] = {
            "sensor_description": channel.sensor.description,
            "input_units": sensitivity.input_units,
            "sample_rate_hz": float(channel.sample_rate),
            "sensitivity": float(sensitivity.value),
            "sensitivity_frequency_hz": float(sensitivity.frequency),
            "response_amplitude_relative_to_reference": {
                str(float(frequency)): float(abs(value / reference))
                for frequency, value in zip(frequencies, response, strict=True)
            },
            "stationxml_elevation_m": float(channel.elevation),
        }
    return result


def plot_diagnostic(arrays, daily, metrics, pressure_multiplier, station):
    """Plot the saved transfer, coherence, and fixed comparison reference."""
    f = arrays["frequency_hz"]
    h = arrays["transfer_dp"]
    coherence = arrays["coherence_squared"]
    valid = arrays["valid"].astype(bool)
    ruan = metrics["ruan_mask"]
    selected_stack = metrics["selected_daily_stack"]
    view = valid & (f >= 0.075) & (f <= 0.225)
    figure, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)

    for item in daily:
        axes[0].semilogy(
            f[view], abs(item["transfer"][view]), color="0.78", linewidth=0.7
        )
        axes[1].plot(
            f[view],
            np.angle(item["transfer"][view], deg=True),
            color="0.82",
            linewidth=0.7,
        )
        axes[2].plot(
            f[view], item["coherence_squared"][view], color="0.82", linewidth=0.7
        )
    axes[0].semilogy(f[view], abs(h[view]), color="navy", label="saved weekly H")
    if np.isfinite(selected_stack).any():
        axes[0].semilogy(
            f[view],
            abs(selected_stack[view]),
            color="seagreen",
            label="mean of daily-reference selections",
        )
        axes[1].plot(
            f[view], np.angle(selected_stack[view], deg=True), color="seagreen"
        )
    axes[0].semilogy(
        f[view],
        abs(metrics["fixed_transfer"][view]),
        color="firebrick",
        linestyle="--",
        label=f"fixed gP={pressure_multiplier:.6g}",
    )
    axes[0].scatter(
        f[ruan], abs(h[ruan]), s=12, color="black", zorder=3, label="coherence² > 0.95"
    )
    axes[1].plot(f[view], np.angle(h[view], deg=True), color="navy")
    axes[1].scatter(f[ruan], np.angle(h[ruan], deg=True), s=12, color="black", zorder=3)
    axes[1].axhline(0, color="0.35", linewidth=0.8)
    axes[2].plot(f[view], coherence[view], color="navy")
    axes[2].axhline(
        RUAN_COHERENCE_REFERENCE,
        color="black",
        linestyle="--",
        linewidth=0.9,
        label="Ruan >0.95 reference",
    )
    for axis in axes:
        axis.axvspan(*MICROSEISM_BAND_HZ, color="gold", alpha=0.12)
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("|D/P| (m/Pa)")
    axes[1].set_ylabel("D/P phase (degrees)")
    axes[2].set_ylabel("Squared coherence")
    axes[2].set_xlabel("Frequency (Hz)")
    axes[2].set_ylim(0, 1.03)
    axes[0].legend(fontsize=8, loc="best")
    axes[2].legend(fontsize=8, loc="best")
    figure.suptitle(
        f"{station} saved microseism pressure–vertical diagnostic\n"
        "No fitted gain, modal model, tilt correction, or score"
    )
    figure.tight_layout()
    return figure


def write_csv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_diagnostic(
    summary_path,
    output_directory,
    *,
    pressure_multiplier,
    curated_json=None,
    curated_plot=None,
    comparison_context=None,
):
    """Create a new diagnostic from saved products, without waveform/master writes."""
    summary_path = Path(summary_path).resolve()
    output = Path(output_directory).resolve()
    output.mkdir(parents=True, exist_ok=False)
    summary = json.loads(summary_path.read_text())
    if summary.get("status") != "comparison_only":
        raise ValueError("Input paired-spectra run is not a completed comparison")
    if summary.get("transfer_direction") != "pressure_to_upward_displacement":
        raise ValueError("Unexpected transfer direction")
    if summary.get("transfer_units") != "m/Pa":
        raise ValueError("Unexpected transfer units")
    provenance = verify_saved_provenance(summary_path, summary)
    arrays = load_and_check_spectra(
        summary["files"]["spectra"], summary.get("retained_segments")
    )
    daily = daily_products(arrays)
    metrics = diagnostic_metrics(arrays, daily, pressure_multiplier)
    response = response_metadata(summary)
    station = ".".join(summary["pressure_epoch"]["Seed_id"].split(".")[:2])

    figure_path = output / "microseism_pressure_vertical_diagnostic.png"
    figure = plot_diagnostic(arrays, daily, metrics, pressure_multiplier, station)
    figure.savefig(figure_path, dpi=160)
    plt.close(figure)

    f = arrays["frequency_hz"]
    h = arrays["transfer_dp"]
    frequency_path = output / "frequency_diagnostic.csv"
    write_csv(
        frequency_path,
        [
            {
                "frequency_hz": float(f[i]),
                "valid": bool(arrays["valid"][i]),
                "coherence_squared": float(arrays["coherence_squared"][i]),
                "nominal_transfer_real_m_pa": float(h[i].real),
                "nominal_transfer_imag_m_pa": float(h[i].imag),
                "nominal_transfer_amplitude_m_pa": float(abs(h[i])),
                "nominal_transfer_phase_deg": float(np.angle(h[i], deg=True)),
                "fixed_reference_transfer_amplitude_m_pa": float(
                    abs(metrics["fixed_transfer"][i])
                ),
                "ruan_high_coherence_reference_bin": bool(metrics["ruan_mask"][i]),
            }
            for i in range(len(f))
        ],
        [
            "frequency_hz",
            "valid",
            "coherence_squared",
            "nominal_transfer_real_m_pa",
            "nominal_transfer_imag_m_pa",
            "nominal_transfer_amplitude_m_pa",
            "nominal_transfer_phase_deg",
            "fixed_reference_transfer_amplitude_m_pa",
            "ruan_high_coherence_reference_bin",
        ],
    )
    daily_path = output / "daily_diagnostic.csv"
    daily_rows = [
        {key: value for key, value in item.items() if not isinstance(value, np.ndarray)}
        for item in daily
    ]
    write_csv(daily_path, daily_rows, list(daily_rows[0]))

    result = {
        "status": "diagnostic_only",
        "station": station,
        "source_comparison": str(summary_path),
        "pressure_seed_id": summary["pressure_epoch"]["Seed_id"],
        "vertical_seed_id": summary["vertical_seed_id"],
        "input_interval_utc": [summary["start_utc"], summary["end_utc"]],
        "retained_segments": int(summary["retained_segments"]),
        "estimator": summary["cross_convention"],
        "transfer_direction": summary["transfer_direction"],
        "transfer_units": summary["transfer_units"],
        "pressure_gain_fit_performed": False,
        "rayleigh_gain_refit_performed": False,
        "fixed_reference_applied_to_original_data": False,
        "identifiability": {
            "equation": (
                "H_nominal(f) = (g_pressure/g_vertical) * eta(f; structure, water)"
            ),
            "conclusion": "not_identifiable_from_saved_PZ_spectra_alone",
            "missing_constraints": [
                "independently verified vertical amplitude response",
                "station-specific coupled water/sediment/crust modal model",
                "independently verified water depth",
            ],
            "interpretation": (
                "The fixed Rayleigh pressure multiplier only rescales the displayed "
                "transfer by its reciprocal; it is not estimated by this diagnostic."
            ),
        },
        "method_limits": [
            "saved pressure/vertical spectra only; no waveform recalculation",
            "no horizontal spectra or tilt correction",
            "no local-event/coda exclusion",
            "bounded interval rather than a seasonal or deployment-scale sample",
            "published coherence rules are comparison references, not adopted scores",
            "metadata response shape is not an independent field calibration",
            (
                "no modal prediction, gain estimate, uncertainty interval, "
                "or classification"
            ),
        ],
        "response_snapshot_metadata": response,
        "results": metrics["summary"],
        "daily": daily_rows,
        "provenance": {
            "saved_input_artifacts": provenance,
            "analyzer_sha256": sha256(Path(__file__)),
            "source_metadata_and_recorded_hashes_verified": True,
        },
        "comparison_context": comparison_context,
        "output_files": {
            "frequency_table": {
                "path": str(frequency_path),
                "sha256": sha256(frequency_path),
            },
            "daily_table": {"path": str(daily_path), "sha256": sha256(daily_path)},
            "plot": {"path": str(figure_path), "sha256": sha256(figure_path)},
        },
    }
    result_path = output / "diagnostic.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    if curated_json:
        curated_json = Path(curated_json)
        curated_json.parent.mkdir(parents=True, exist_ok=True)
        if curated_json.exists():
            raise FileExistsError(f"Refusing to replace curated result {curated_json}")
        shutil.copyfile(result_path, curated_json)
    if curated_plot:
        curated_plot = Path(curated_plot)
        curated_plot.parent.mkdir(parents=True, exist_ok=True)
        if curated_plot.exists():
            raise FileExistsError(f"Refusing to replace curated plot {curated_plot}")
        shutil.copyfile(figure_path, curated_plot)
    return result_path


def main():
    root = Path(__file__).resolve().parents[1] / "pressure"
    default_summary = root / (
        "results/NETWORKS/TEST/YV/RR52/"
        "YV.RR52.00.BDH__20121019T090000__58afc2d6f005/"
        "paired_spectra/run_bfbef17f2c13/comparison.json"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, default=default_summary)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--fixed-pressure-multiplier",
        type=float,
    )
    parser.add_argument("--curated-json", type=Path)
    parser.add_argument("--curated-plot", type=Path)
    args = parser.parse_args()
    pressure_multiplier = args.fixed_pressure_multiplier
    if pressure_multiplier is None:
        if args.comparison.resolve() != default_summary.resolve():
            parser.error(
                "--fixed-pressure-multiplier is required for a non-RR52 comparison"
            )
        pressure_multiplier = RR52_FIXED_PRESSURE_MULTIPLIER
    input_summary = json.loads(args.comparison.read_text())
    station = ".".join(
        input_summary["pressure_epoch"]["Seed_id"].split(".")[:2]
    ).replace(".", "_")
    output = args.output_dir or (
        root / f"results/DIAGNOSTICS/{station}_MICROSEISM" / f"run_{uuid4().hex[:12]}"
    )
    result = run_diagnostic(
        args.comparison,
        output,
        pressure_multiplier=pressure_multiplier,
        curated_json=args.curated_json,
        curated_plot=args.curated_plot,
    )
    print(result)


if __name__ == "__main__":
    main()
