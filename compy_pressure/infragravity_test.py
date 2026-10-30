"""Compare saved pressure spectra with the depth-dependent Godin IG-wave model."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal import PPSD

from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    OUTPUT_DIRECTORY,
    TEST_COLUMNS,
    load_master_dataset,
    update_infragravity_test_results,
)
from compy_pressure.pressure_test import (
    response_band,
    shade_excluded_band,
    validate_pressure,
)

GRAVITY = 9.81
PRESSURE_NORMALIZATION = 2000.0  # Equation 6: Pa^2/Hz, not per radian/second.
REFERENCE_URL = "https://doi.org/10.1002/grl.50418"
SCORE_COLUMN = TEST_COLUMNS[2]
STATUS_COLUMN = SCORE_COLUMN.removesuffix("_score") + "_status"
MODEL_COLOR = "#D55E00"
NO_TANH_COLOR = "#7B3294"  # Eq. 6 without the tanh(kH)^(-3/2) factor
NO_TANH_LABEL = r"Godin Eq. 6 without tanh$^{-3/2}$ (2000 cosh$^{-2}$kH)"


def positive(value, name):
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite and positive") from exc
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def dispersion_kh(frequency_hz, depth_m, *, gravity=GRAVITY):
    """Solve x*tanh(x)=(2*pi*f)^2*H/g for x=kH by vectorized bisection."""
    depth_m = positive(depth_m, "Water depth")
    gravity = positive(gravity, "Gravity")
    frequencies = np.asarray(frequency_hz, dtype=float)
    if not frequencies.size or not np.all(np.isfinite(frequencies) & (frequencies > 0)):
        raise ValueError("Model frequencies must be finite and positive")
    a = (2 * np.pi * frequencies) ** 2 * depth_m / gravity
    if not np.isfinite(a).all():
        raise ValueError("Frequency/depth combination exceeds numerical support")
    lower, upper = a.copy(), a + np.sqrt(a) + 1
    for _ in range(80):
        middle = (lower + upper) / 2
        below = middle * np.tanh(middle) < a
        lower = np.where(below, middle, lower)
        upper = np.where(below, upper, middle)
    return (lower + upper) / 2


def godin_pressure_db(frequency_hz, depth_m, *, gravity=GRAVITY, tanh_factor=True):
    """Equation 6 in dB re 1 Pa^2/Hz; stable even for large kH.

    S_p = 2000 * tanh(kH)^(-3/2) * cosh(kH)^(-2).
    Source: Godin et al., DOI 10.1002/grl.50418, Equation 6.
    Angular frequency enters dispersion through omega=2*pi*f only; the
    published normalization already uses per-Hz power density.
    With tanh_factor=False the tanh(kH)^(-3/2) factor is dropped, giving
    2000 * cosh(kH)^(-2), a display variant; kH still uses full dispersion.
    """
    x = dispersion_kh(frequency_hz, depth_m, gravity=gravity)
    log_cosh = np.logaddexp(x, -x) - np.log(2)
    log_power = np.log(PRESSURE_NORMALIZATION) - 2 * log_cosh
    if tanh_factor:
        log_power -= 1.5 * np.log(np.tanh(x))
    return 10 / np.log(10) * log_power


def compare_spectra(
    ppsd,
    values,
    depth_m,
    sensor_band,
    *,
    model_min=0.001,
    model_max=0.02,
    gravity=GRAVITY,
):
    """Apply identical frequency-bin averaging to theory and saved dB spectra."""
    model_min = positive(model_min, "Minimum model frequency")
    model_max = positive(model_max, "Maximum model frequency")
    if model_min >= model_max or model_max > ppsd.sampling_rate / 2:
        raise ValueError("Model band must satisfy minimum < maximum <= Nyquist")
    frequencies = 1 / ppsd.period_bin_centers
    lower_edges, upper_edges = (
        1 / ppsd.period_bin_right_edges,
        1 / ppsd.period_bin_left_edges,
    )
    percentiles = np.percentile(values, [5, 50, 95], axis=0)
    # Match ObsPy's actual period mask, including floating-point boundary rules.
    fft_frequencies = np.fft.rfftfreq(ppsd.nfft, d=1 / ppsd.sampling_rate)[1:]
    fft_periods = 1 / fft_frequencies
    rows = []
    for index, frequency in enumerate(frequencies):
        lower, upper = float(lower_edges[index]), float(upper_edges[index])
        reason = ""
        model = None
        if lower < model_min or upper > model_max:
            reason = "smoothing_bin_outside_model_band"
        else:
            mask = (fft_periods >= ppsd.period_bin_left_edges[index]) & (
                fft_periods <= ppsd.period_bin_right_edges[index]
            )
            bin_frequencies = fft_frequencies[mask]
            if not len(bin_frequencies):
                reason = "no_fft_frequencies_in_bin"
            else:
                model = float(
                    np.mean(
                        godin_pressure_db(bin_frequencies, depth_m, gravity=gravity)
                    )
                )
        if not reason and (
            lower < sensor_band["minimum_hz"] or upper > sensor_band["maximum_hz"]
        ):
            reason = "smoothing_bin_outside_sensor_band"
        p05, median, p95 = map(float, percentiles[:, index])
        rows.append(
            {
                "frequency_hz": float(frequency),
                "smoothing_minimum_hz": lower,
                "smoothing_maximum_hz": upper,
                "model_db_pa2_hz": model,
                "p05_db_pa2_hz": p05,
                "median_db_pa2_hz": median,
                "p95_db_pa2_hz": p95,
                "diagnostic_median_minus_model_db": median - model
                if model is not None
                else None,
                "diagnostic_p05_minus_model_db": p05 - model
                if model is not None
                else None,
                "diagnostic_p95_minus_model_db": p95 - model
                if model is not None
                else None,
                "median_minus_model_db": median - model if not reason else None,
                "comparison": "compared" if not reason else "excluded",
                "exclusion_reason": reason,
            }
        )
    return sorted(rows, key=lambda row: row["frequency_hz"]), percentiles


def residual_plot(ppsd, rows, sensor_band, lower, upper):
    """Separate all diagnostic residuals from eligible, overlapping PSD bins."""
    available = [row for row in rows if row["model_db_pa2_hz"] is not None]
    included = [row for row in available if row["comparison"] == "compared"]
    figure, (overview, selected) = plt.subplots(1, 2, figsize=(12, 5.8))
    x = np.array([row["frequency_hz"] for row in available]) * 1000
    overview.fill_between(
        x,
        [row["diagnostic_p05_minus_model_db"] for row in available],
        [row["diagnostic_p95_minus_model_db"] for row in available],
        color="#888888",
        alpha=0.22,
        label="5th–95th window spread (diagnostic)",
    )
    overview.plot(
        x,
        [row["diagnostic_median_minus_model_db"] for row in available],
        color="#666666",
        linewidth=1.8,
        label="Diagnostic median",
    )
    shade_excluded_band(
        overview,
        lower * 1000,
        upper * 1000,
        {
            "minimum_hz": sensor_band["minimum_hz"] * 1000,
            "maximum_hz": sensor_band["maximum_hz"] * 1000,
        },
    )
    overview.set_xscale("log")
    overview.set_title(
        "Full diagnostic view\nGray region: outside selected sensor band"
    )
    if included:
        x = np.array([row["frequency_hz"] for row in included]) * 1000
        y = np.array([row["median_minus_model_db"] for row in included])
        overview.scatter(
            x,
            y,
            color="#0072B2",
            edgecolor="white",
            zorder=5,
            label="Eligible comparison bins",
        )
        selected.hlines(
            y,
            [r["smoothing_minimum_hz"] * 1000 for r in included],
            [r["smoothing_maximum_hz"] * 1000 for r in included],
            color="#888888",
            linewidth=2.5,
            alpha=0.6,
            label="Horizontal: full smoothing-bin range",
        )
        selected.vlines(
            x,
            [r["diagnostic_p05_minus_model_db"] for r in included],
            [r["diagnostic_p95_minus_model_db"] for r in included],
            color="#0072B2",
            linewidth=2.5,
            label="Vertical: 5th–95th window spread",
        )
        selected.plot(
            x, y, "o", color="#0072B2", markersize=7, label="Eligible median residual"
        )
        minimum = min(r["smoothing_minimum_hz"] for r in included) * 1000
        maximum = max(r["smoothing_maximum_hz"] for r in included) * 1000
        selected.set_xlim(minimum * 0.96, maximum * 1.04)
        selected.set_ylim(
            min(0, min(r["diagnostic_p05_minus_model_db"] for r in included)) - 1,
            max(0, max(r["diagnostic_p95_minus_model_db"] for r in included)) + 1,
        )
    else:
        selected.text(
            0.5,
            0.5,
            "No eligible comparison bins",
            ha="center",
            transform=selected.transAxes,
        )
        selected.set_xlim(lower * 1000, upper * 1000)
    selected.set_title(
        f"{len(included)} eligible bins: detailed view\n"
        "Overlapping smoothing bins are correlated"
    )
    for axis in (overview, selected):
        axis.axhline(
            0,
            color="black",
            linestyle="--",
            linewidth=1,
            label="Observed median = model",
        )
        axis.set_xlabel("Frequency [mHz]")
        axis.set_ylabel("Observed minus PSD-bin-averaged model [dB]")
        axis.grid(which="both", alpha=0.2)
        axis.legend(fontsize=8, loc="upper left" if axis is overview else "lower right")
    figure.suptitle(f"{ppsd.id}: Figure 3 — Godin residuals", fontsize=13)
    figure.subplots_adjust(left=0.08, right=0.98, bottom=0.21, top=0.80, wspace=0.25)
    figure.text(
        0.08,
        0.04,
        "Positive = above prediction; negative = below prediction. "
        "Vertical spread is not a confidence interval.\n"
        "Only blue points contribute to the comparison summary. "
        "The model is a reference, not an acceptance boundary.",
        fontsize=9,
    )
    return figure


def comparison_plots(
    ppsd,
    percentiles,
    rows,
    depth_m,
    sensor_band,
    *,
    model_min=0.001,
    model_max=0.02,
    gravity=GRAVITY,
):
    """Return percentile/model, PPSD/model and residual views in physical units."""
    frequencies = 1 / ppsd.period_bin_centers
    order = np.argsort(frequencies)
    lower = max(model_min, float(min(frequencies)))
    upper = min(model_max, float(max(frequencies)))
    if lower >= upper:
        raise ValueError("Saved spectra have no displayed overlap with the model band")
    dense_f = np.geomspace(lower, upper, 400)
    dense_model = godin_pressure_db(dense_f, depth_m, gravity=gravity)
    dense_no_tanh = godin_pressure_db(
        dense_f, depth_m, gravity=gravity, tanh_factor=False
    )
    available = [row for row in rows if row["model_db_pa2_hz"] is not None]
    included = [row for row in rows if row["comparison"] == "compared"]
    figures = []

    figure, axis = plt.subplots(figsize=(10, 6), layout="constrained")
    axis.fill_between(
        frequencies[order],
        percentiles[0, order],
        percentiles[2, order],
        color="#0072B2",
        alpha=0.2,
        label="Observed 5th–95th percentiles",
    )
    axis.plot(
        frequencies[order],
        percentiles[1, order],
        color="#0072B2",
        label="Observed median",
        linewidth=2,
    )
    axis.plot(
        dense_f,
        dense_model,
        color=MODEL_COLOR,
        alpha=0.55,
        label="Godin equation (unsmoothed)",
    )
    axis.plot(
        dense_f,
        dense_no_tanh,
        color=NO_TANH_COLOR,
        linestyle=":",
        linewidth=2,
        label=NO_TANH_LABEL,
    )
    axis.plot(
        [r["frequency_hz"] for r in available],
        [r["model_db_pa2_hz"] for r in available],
        color=MODEL_COLOR,
        linestyle="--",
        linewidth=2,
        label="Godin averaged over PSD bins",
    )
    if included:
        axis.scatter(
            [r["frequency_hz"] for r in included],
            [r["median_db_pa2_hz"] for r in included],
            color="#0072B2",
            edgecolor="white",
            zorder=5,
            label="Compared median bins",
        )
    shade_excluded_band(axis, lower, upper, sensor_band)
    axis.set_xscale("log")
    axis.set_xlabel("Frequency [Hz]")
    axis.set_ylabel(r"Pressure PSD [dB re 1 Pa$^2$/Hz]")
    axis.set_title(
        f"{ppsd.id}: Figure 1 — Godin comparison at {depth_m:g} m\n"
        "Depth and sensor band require review"
    )
    # Model view must not autoscale to irrelevant seismic/microseism frequencies.
    visible = (frequencies >= lower) & (frequencies <= upper)
    if not visible.any():
        # A narrow requested view can lie between saved period centers. The
        # adjacent displayed spectra still bound the histogram cells; record
        # no eligible bins instead of crashing during axis autoscaling.
        visible[np.argmin(abs(np.log(frequencies / lower)))] = True
        visible[np.argmin(abs(np.log(frequencies / upper)))] = True
    axis.set_ylim(
        min(float(percentiles[0, visible].min()), float(dense_model.min())) - 3,
        max(
            float(percentiles[2, visible].max()),
            float(dense_model.max()),
            float(dense_no_tanh.max()),
        )
        + 3,
    )
    axis.grid(which="both", alpha=0.2)
    axis.legend(fontsize=8)
    figures.append(figure)

    figure = ppsd.plot(
        show=False,
        xaxis_frequency=True,
        show_noise_models=False,
        show_coverage=False,
        period_lim=(lower, upper),
        cmap=plt.get_cmap("viridis"),
    )
    axis = figure.axes[0]
    axis.plot(
        [r["frequency_hz"] for r in available],
        [r["model_db_pa2_hz"] for r in available],
        color=MODEL_COLOR,
        linewidth=2.2,
        linestyle="--",
        label="Godin averaged over PSD bins",
        path_effects=[path_effects.withStroke(linewidth=4, foreground="white")],
    )
    axis.plot(
        dense_f,
        dense_no_tanh,
        color=NO_TANH_COLOR,
        linewidth=2,
        linestyle=":",
        label=NO_TANH_LABEL,
        path_effects=[path_effects.withStroke(linewidth=4, foreground="white")],
    )
    shade_excluded_band(axis, lower, upper, sensor_band)
    axis.set_ylim(figures[0].axes[0].get_ylim())
    axis.set_ylabel(r"Pressure PSD [dB re 1 Pa$^2$/Hz]")
    axis.set_title(
        f"{ppsd.id}: Figure 2 — pressure PPSD and Godin model ({depth_m:g} m)"
    )
    axis.legend(fontsize=8, loc="lower left")
    figure.set_size_inches(11, 7.5)
    axis.set_position([0.10, 0.22, 0.73, 0.66])
    figure.axes[1].set_position([0.86, 0.22, 0.025, 0.66])
    figure.axes[1].set_ylabel("Occurrence [%]")
    figure.text(
        0.1,
        0.055,
        "Orange: Godin model. Color map: occurrence probability.\n"
        "Gray: outside selected sensor band; whole smoothing bins must "
        "fit both bands.\nDiagnostic only; no acceptance tolerance or "
        "calibration correction has been inferred.",
        fontsize=8,
    )
    figures.append(figure)

    figures.append(residual_plot(ppsd, rows, sensor_band, lower, upper))
    return figures


def run_infragravity_test(
    master_dir=OUTPUT_DIRECTORY,
    *,
    network="YV",
    station="RR52",
    depth_m=None,
    model_min=0.001,
    model_max=0.02,
    sensor_min=None,
    sensor_max=None,
    loss_db=3,
    gravity=GRAVITY,
    update_master=True,
    show=True,
    selected_rows=None,
    return_updates=False,
):
    """Save diagnostic predictions and plots; leave calibration scores untouched."""
    master = load_master_dataset(master_dir) if selected_rows is None else selected_rows
    if not {"ppsd_status", "ppsd_npz_file", "ppsd_units"}.issubset(master):
        raise ValueError("Run ppsd.py first to populate saved pressure results")
    rows = master[
        master.Seed_id.str.startswith(f"{network}.{station}.")
        & master.ppsd_status.isin(["complete", "partial"])
        & master.ppsd_units.eq("Pa^2/Hz")
        & master.ppsd_npz_file.ne("")
    ].drop_duplicates(EPOCH_COLUMNS + ["inventory_file"])
    if rows.empty:
        raise ValueError(f"No saved pressure PPSDs for {network}.{station}")
    summaries, updates, figures = [], [], []
    for _, row in rows.iterrows():
        depth = positive(
            depth_m if depth_m is not None else row.get("Profondeur", ""), "Water depth"
        )
        source = Path(row.ppsd_npz_file).expanduser().resolve()
        ppsd = PPSD.load_npz(str(source))
        values = validate_pressure(ppsd, row)
        snapshot = source.parent / "inventory.xml"
        band, _, _ = response_band(
            ppsd,
            snapshot,
            sensor_min=sensor_min,
            sensor_max=sensor_max,
            loss_db=loss_db,
        )
        diagnostics, percentiles = compare_spectra(
            ppsd,
            values,
            depth,
            band,
            model_min=model_min,
            model_max=model_max,
            gravity=gravity,
        )
        included = [r for r in diagnostics if r["comparison"] == "compared"]
        residuals = np.array([r["median_minus_model_db"] for r in included])
        status = "comparison_only" if included else "no_comparable_bins"
        effective_min = max(model_min, band["minimum_hz"])
        effective_max = min(model_max, band["maximum_hz"])
        if effective_min >= effective_max:
            effective_min = effective_max = None
        output = source.parent / "infragravity"
        output.mkdir(exist_ok=True)
        csv_path, json_path = (
            output / "frequency_comparison.csv",
            output / "comparison.json",
        )
        with csv_path.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(diagnostics[0]), lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(diagnostics)
        plot_paths = {
            name: output / f"{name}.png"
            for name in ("model_comparison", "ppsd_with_model", "model_residuals")
        }
        new_figures = comparison_plots(
            ppsd,
            percentiles,
            diagnostics,
            depth,
            band,
            model_min=model_min,
            model_max=model_max,
            gravity=gravity,
        )
        for figure, path in zip(new_figures, plot_paths.values(), strict=True):
            figure.savefig(path, dpi=160)
        figures.extend(new_figures)
        summary = {
            "seed_id": ppsd.id,
            "status": status,
            "score": None,
            "reference": {
                "doi": "10.1002/grl.50418",
                "equation": 6,
                "url": REFERENCE_URL,
                "normalization_pa2_hz": PRESSURE_NORMALIZATION,
            },
            "depth_m": depth,
            "depth_basis": "user_supplied"
            if depth_m is not None
            else "inventory_elevation_estimate",
            "depth_independently_verified": False,
            "gravity_m_s2": gravity,
            "psd_units": "Pa^2/Hz",
            "psd_windows": len(values),
            "model_band_min_hz": model_min,
            "model_band_max_hz": model_max,
            "sensor_band": band,
            "effective_band_min_hz": effective_min,
            "effective_band_max_hz": effective_max,
            "compared_bins": len(included),
            "diagnostic_bins": sum(
                r["model_db_pa2_hz"] is not None for r in diagnostics
            ),
            "residual_definition": (
                "Observed window median minus PSD-bin-averaged model, in dB."
            ),
            "spectral_spread_definition": (
                "Empirical 5th–95th window percentiles; not confidence intervals."
            ),
            "smoothing_definition": (
                "Exact ObsPy period mask over the FFT grid; average model dB values."
            ),
            "median_residual_db": float(np.median(residuals)) if included else None,
            "median_absolute_residual_db": float(np.median(abs(residuals)))
            if included
            else None,
            "source_npz_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "response_snapshot_sha256": hashlib.sha256(
                snapshot.read_bytes()
            ).hexdigest(),
            "files": {
                "frequency_csv": str(csv_path),
                **{k: str(v) for k, v in plot_paths.items()},
            },
            "limitations": [
                "Depth and sensor band need independent validation.",
                "Default 0.001–0.02 Hz band is a conservative workflow choice, "
                "not a universal validity limit.",
                "No acceptance tolerance, pass/fail score or correction factor "
                "is inferred.",
                "Existing 3-hour PPSD windows differ from the original "
                "12.42-hour Hann-window and transient-screening procedure.",
                "Model normalization fitted a different ocean region; "
                "local departures need interpretation.",
                "Overlapping windows and frequency smoothing produce "
                "correlated samples.",
            ],
        }
        json_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        update = row[EPOCH_COLUMNS + ["inventory_file"]].to_dict() | {
            "infragravity_status": status,
            "infragravity_result_file": str(json_path),
            "infragravity_csv_file": str(csv_path),
            "infragravity_depth_m": depth,
            "infragravity_depth_basis": summary["depth_basis"],
            "infragravity_compared_bins": len(included),
            "infragravity_median_residual_db": summary["median_residual_db"],
            "infragravity_median_absolute_residual_db": summary[
                "median_absolute_residual_db"
            ],
            "infragravity_band_min_hz": effective_min,
            "infragravity_band_max_hz": effective_max,
            **{
                f"infragravity_{name}_file": str(path)
                for name, path in plot_paths.items()
            },
        }
        if not row.get(SCORE_COLUMN, ""):
            update[STATUS_COLUMN] = status
        updates.append(update)
        summaries.append(summary)
        print(
            f"{ppsd.id}: Godin model at {depth:g} m; {len(included)} compared bins; "
            f"median residual={summary['median_residual_db']} dB; status={status}",
            flush=True,
        )
        for path in plot_paths.values():
            print(f"Saved {path}", flush=True)
    if update_master:
        update_infragravity_test_results(master_dir, updates)
    if show:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)
    return (summaries, updates) if return_updates else summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-dir", type=Path, default=OUTPUT_DIRECTORY)
    parser.add_argument("--network", default="YV")
    parser.add_argument("--station", default="RR52")
    parser.add_argument("--depth-m", type=float)
    parser.add_argument("--model-min-frequency", type=float, default=0.001)
    parser.add_argument("--model-max-frequency", type=float, default=0.02)
    parser.add_argument("--sensor-min-frequency", type=float)
    parser.add_argument("--sensor-max-frequency", type=float)
    parser.add_argument("--response-loss-db", type=float, default=3)
    parser.add_argument("--gravity", type=float, default=GRAVITY)
    parser.add_argument("--no-update-master", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    run_infragravity_test(
        args.master_dir,
        network=args.network,
        station=args.station,
        depth_m=args.depth_m,
        model_min=args.model_min_frequency,
        model_max=args.model_max_frequency,
        sensor_min=args.sensor_min_frequency,
        sensor_max=args.sensor_max_frequency,
        loss_db=args.response_loss_db,
        gravity=args.gravity,
        update_master=not args.no_update_master,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()
