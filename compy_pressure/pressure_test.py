"""Descriptive pressure comparisons with Webb (1998, ≤0.2 Hz) and Brown (≥0.2 Hz)."""

import argparse
import csv
import hashlib
import json
from importlib.resources import files
from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from obspy import read_inventory
from obspy.signal import PPSD

from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    OUTPUT_DIRECTORY,
    TEST_COLUMNS,
    load_master_dataset,
    update_pressure_test_results,
)

SCORE_COLUMN = TEST_COLUMNS[0]
STATUS_COLUMN = SCORE_COLUMN.removesuffix("_score") + "_status"
# Red below and green above the 0.2 Hz split between the two sources.
REFERENCE_COLORS = {"webb_instructions": "#C62828", "brown_figure8": "#2E7D32"}
REFERENCE_NAMES = {
    "webb_instructions": "Webb 1998",
    "brown_figure8": "Brown et al. 2014",
}
STATUS_MARKERS = {
    "within": ("o", "#222222"),
    "below": ("v", "#6A3D9A"),
    "above": ("^", "#E66101"),
}


def reference_limits(axis, diagnostics):
    """Display all supplied reference points separately from numerical inclusion."""
    outline = [path_effects.withStroke(linewidth=4, foreground="white")]
    for reference, color in REFERENCE_COLORS.items():
        rows = [
            row
            for row in diagnostics
            if row["reference"] == reference
            and row["reference_low_db_pa2_hz"] is not None
        ]
        if not rows:
            continue
        for bound, marker, style in [("low", "v", "--"), ("high", "^", "-")]:
            x = [
                row["reference_frequency_hz"]
                if reference == "webb_instructions"
                else row["frequency_hz"]
                for row in rows
            ]
            y = [row[f"reference_{bound}_db_pa2_hz"] for row in rows]
            axis.plot(
                x,
                y,
                color=color,
                linestyle="none" if reference == "webb_instructions" else style,
                marker=marker if reference == "webb_instructions" else None,
                markersize=7,
                linewidth=2.2,
                zorder=6,
                label=(
                    f"{REFERENCE_NAMES[reference]} "
                    f"{'lower' if bound == 'low' else 'upper'} limit"
                ),
                path_effects=outline,
            )


def shade_excluded_band(axis, lower, upper, band):
    """Shade frequencies outside the band while keeping the saved-data axis."""
    if band["minimum_hz"] > lower:
        axis.axvspan(lower, min(upper, band["minimum_hz"]), color="gray", alpha=0.15)
    if band["maximum_hz"] < upper:
        axis.axvspan(max(lower, band["maximum_hz"]), upper, color="gray", alpha=0.15)
    axis.set_xlim(lower, upper)


def ppsd_limits_plot(ppsd, diagnostics, band):
    """Overlay two color-coded references on the observed probability histogram."""
    frequencies = 1 / ppsd.period_bin_centers
    lower, upper = float(min(frequencies)), float(max(frequencies))
    figure = ppsd.plot(
        show=False,
        xaxis_frequency=True,
        show_noise_models=False,
        show_coverage=False,
        period_lim=(lower, upper),
        cmap=plt.get_cmap("viridis"),
    )
    axis = figure.axes[0]
    reference_limits(axis, diagnostics)
    shade_excluded_band(axis, lower, upper, band)
    amplitudes = np.asarray(ppsd.psd_values)
    limits = [
        row[f"reference_{bound}_db_pa2_hz"]
        for row in diagnostics
        for bound in ("low", "high")
        if row[f"reference_{bound}_db_pa2_hz"] is not None
    ]
    axis.set_ylim(
        min(float(amplitudes.min()), min(limits)) - 5,
        max(float(amplitudes.max()), max(limits)) + 5,
    )
    axis.set_ylabel(r"Pressure PSD [dB re 1 Pa$^2$/Hz]")
    axis.set_xlabel("Frequency [Hz]")
    axis.set_title(f"{ppsd.id}: PPSD with Webb and Brown limits")
    axis.legend(loc="upper right", fontsize=8, framealpha=0.96)
    # ObsPy creates these axes with add_axes, so subplots_adjust cannot reserve
    # caption space. Set both positions explicitly to keep labels unclipped.
    figure.set_size_inches(11, 7.5)
    axis.set_position([0.10, 0.22, 0.73, 0.66])
    figure.axes[1].set_position([0.86, 0.22, 0.025, 0.66])
    figure.axes[1].set_ylabel("Occurrence [%]")
    figure.text(
        0.1,
        0.055,
        "Color map: occurrence probability. Red: Webb 1998 Fig. 2 (≤0.2 Hz). "
        "Green: Brown et al. 2014 (≥0.2 Hz).\n"
        "Webb: discrete figure-read bounds. Brown: digitized, PSD-bin averaged.\n"
        "Gray: outside selected sensor band; only fully supported bins are "
        "compared. Screening only.",
        fontsize=8,
    )
    return figure


def test_details_plot(ppsd, diagnostics):
    """Show distance to bounds and within-range window percentages separately."""
    figure, axes = plt.subplots(
        2, 2, figsize=(12, 8), layout="constrained", sharex="col", sharey="row"
    )
    markers = {"within": "o", "below": "v", "above": "^"}
    for column, (reference, color) in enumerate(REFERENCE_COLORS.items()):
        rows = [
            row
            for row in diagnostics
            if row["reference"] == reference and row["comparison"] != "excluded"
        ]
        margin_axis, occupancy_axis = axes[:, column]
        margin_axis.axhline(0, color="black", linewidth=1)
        if reference == "brown_figure8":
            margin_axis.axhspan(-2, 2, color="gray", alpha=0.12)
        for row in rows:
            # Positive clearance means inside both bounds; negative clearance
            # is the dB exceedance below the lower or above the upper limit.
            margin = min(
                row["median_offset_from_low_db"], -row["median_offset_from_high_db"]
            )
            marker = markers[row["comparison"]]
            x = row["frequency_hz"]
            margin_axis.vlines(x, 0, margin, color=color, linewidth=1.6)
            margin_axis.scatter(x, margin, color=color, marker=marker, s=55, zorder=3)
            occupancy_axis.scatter(
                x,
                100 * row["window_fraction_within_bounds"],
                color=color,
                marker=marker,
                s=55,
                zorder=3,
            )
            if row["near_digitization_boundary"]:
                margin_axis.scatter(
                    x,
                    margin,
                    s=120,
                    facecolor="none",
                    edgecolor="black",
                    linewidth=1.2,
                    zorder=4,
                )
        count = sum(row["comparison"] == "within" for row in rows)
        margin_axis.set_title(
            f"{REFERENCE_NAMES[reference]}: {count}/{len(rows)} "
            "compared medians within limits",
            color=color,
        )
        if not rows:
            margin_axis.text(
                0.5,
                0.5,
                "No supported reference bins",
                transform=margin_axis.transAxes,
                ha="center",
            )
        for panel in (margin_axis, occupancy_axis):
            panel.set_xscale("log")
            panel.grid(which="both", alpha=0.2)
        occupancy_axis.set_ylim(-5, 105)
        occupancy_axis.set_xlabel("Compared frequency [Hz]")
    axes[0, 0].set_ylabel("Median clearance to nearest bound [dB]\n+ inside; − outside")
    axes[1, 0].set_ylabel("PSD windows within both limits [%]")
    legend = [
        Line2D(
            [],
            [],
            color="#555555",
            marker=marker,
            linestyle="none",
            label=name.capitalize(),
        )
        for name, marker in markers.items()
    ]
    legend.append(
        Line2D(
            [],
            [],
            color="black",
            marker="o",
            markerfacecolor="none",
            linestyle="none",
            markersize=9,
            label="Digitization review",
        )
    )
    axes[0, 1].legend(handles=legend, fontsize=8, loc="best")
    figure.suptitle(
        f"{ppsd.id}: pressure-level comparison details\n"
        "Excluded bins omitted; percentages are descriptive, "
        "not calibration acceptance",
        fontsize=11,
    )
    return figure


def load_references():
    directory = files("compy_pressure").joinpath("reference_data")
    metadata = json.loads(directory.joinpath("pressure_references.json").read_text())
    with directory.joinpath(metadata["brown"]["data_file"]).open() as stream:
        brown = list(csv.DictReader(stream))
    return metadata, np.array(
        [
            [
                float(row[key])
                for key in ("frequency_hz", "low_db_pa2_hz", "high_db_pa2_hz")
            ]
            for row in brown
        ]
    )


def validate_pressure(ppsd, row):
    if row["ppsd_units"] != "Pa^2/Hz" or ppsd.special_handling != "hydrophone":
        raise ValueError("This comparison requires pressure PSDs in Pa^2/Hz")
    if ppsd.id != row["Seed_id"]:
        raise ValueError("Saved PPSD identifier does not match its master row")
    values = np.asarray(ppsd.psd_values, dtype=float)
    if values.shape != (len(ppsd.times_processed), len(ppsd.period_bin_centers)):
        raise ValueError("PSD values and frequency/time coordinates disagree")
    if not values.size or not np.isfinite(values).all():
        raise ValueError("Pressure spectra are empty or contain non-finite values")
    if len(set(map(str, ppsd.times_processed))) != len(values):
        raise ValueError("Duplicate PSD start times would bias this comparison")
    return values


def response_band(ppsd, snapshot, *, sensor_min=None, sensor_max=None, loss_db=3):
    """Estimate a nominal response band; user bounds override this estimate.

    The inferred band is the contiguous range around the sensitivity reference
    frequency with amplitude attenuation no greater than loss_db. It is not an
    independently verified sensor operating band or calibration result.
    """
    if not np.isfinite(loss_db) or loss_db <= 0:
        raise ValueError("Response-loss threshold must be finite and positive")
    if (sensor_min is None) != (sensor_max is None):
        raise ValueError(
            "Provide both --sensor-min-frequency and --sensor-max-frequency"
        )
    nyquist = ppsd.sampling_rate / 2
    if sensor_min is not None and not (
        np.isfinite(sensor_min)
        and np.isfinite(sensor_max)
        and 0 < sensor_min < sensor_max <= nyquist
    ):
        raise ValueError("Sensor bounds must satisfy 0 < minimum < maximum <= Nyquist")
    inventory = read_inventory(str(snapshot))
    response = inventory.get_response(ppsd.id, ppsd.times_processed[0])
    reference = float(response.instrument_sensitivity.frequency)
    if not 0 < reference < nyquist:
        raise ValueError("Invalid response sensitivity reference frequency")
    frequencies = np.unique(
        np.r_[
            np.geomspace(1 / max(ppsd.period_bin_right_edges), nyquist * 0.999, 2000),
            reference,
        ]
    )
    amplitudes = np.abs(
        response.get_evalresp_response_for_frequencies(frequencies, output="DEF")
    )
    ref_index = int(np.flatnonzero(frequencies == reference)[0])
    if not np.isfinite(amplitudes).all() or amplitudes[ref_index] <= 0:
        raise ValueError("Invalid pressure-response amplitudes")
    relative_db = 20 * np.log10(
        np.maximum(amplitudes, np.finfo(float).tiny) / amplitudes[ref_index]
    )
    keep = relative_db >= -loss_db
    left = right = ref_index
    while left > 0 and keep[left - 1]:
        left -= 1
    while right < len(keep) - 1 and keep[right + 1]:
        right += 1
    minimum, maximum = float(frequencies[left]), float(frequencies[right])
    if minimum >= maximum:
        raise ValueError("No usable interval around the response reference frequency")
    band = {
        "minimum_hz": sensor_min if sensor_min is not None else minimum,
        "maximum_hz": sensor_max if sensor_max is not None else maximum,
        "basis": "user_supplied" if sensor_min is not None else "nominal_response",
        "nominal_minimum_hz": minimum,
        "nominal_maximum_hz": maximum,
        "reference_frequency_hz": reference,
        "maximum_attenuation_db": loss_db,
        "independently_verified": False,
    }
    return band, frequencies, relative_db


def summarize(values, frequencies, lower_edges, upper_edges, band, metadata, brown):
    """Compare only supported bins; preserve skipped points and missing coverage."""
    frequencies = np.asarray(frequencies)
    percentiles = np.percentile(values, [5, 50, 95], axis=0)
    diagnostics = []

    def add(model, index, low, high, target, reason="", tolerance=0):
        p05, median, p95 = map(float, percentiles[:, index])
        usable = (
            lower_edges[index] >= band["minimum_hz"]
            and upper_edges[index] <= band["maximum_hz"]
        )
        if not usable:
            reason = "smoothing_bin_outside_sensor_band"
        if model == "webb_instructions" and not metadata["webb"].get(
            "use_in_comparison", False
        ):
            reason = "unverified_project_instruction_bounds"
        within = low <= median <= high
        diagnostics.append(
            {
                "reference": model,
                "reference_frequency_hz": float(target),
                "frequency_hz": float(frequencies[index]),
                "smoothing_minimum_hz": float(lower_edges[index]),
                "smoothing_maximum_hz": float(upper_edges[index]),
                "reference_low_db_pa2_hz": float(low),
                "reference_high_db_pa2_hz": float(high),
                "p05_db_pa2_hz": p05,
                "median_db_pa2_hz": median,
                "p95_db_pa2_hz": p95,
                "window_fraction_within_bounds": float(
                    np.mean((values[:, index] >= low) & (values[:, index] <= high))
                )
                if not reason
                else None,
                "median_offset_from_low_db": median - low,
                "median_offset_from_high_db": median - high,
                "comparison": "excluded"
                if reason
                else ("within" if within else "below" if median < low else "above"),
                "exclusion_reason": reason,
                "near_digitization_boundary": bool(
                    min(abs(median - low), abs(median - high)) <= tolerance
                )
                if tolerance and not reason
                else False,
            }
        )

    for point in metadata["webb"]["points"]:
        target = point["frequency_hz"]
        index = int(np.argmin(abs(np.log(frequencies / target))))
        # A discrete handout point must match within half a period-center step.
        steps = np.diff(np.sort(np.log(frequencies)))
        reason = (
            "reference_point_outside_saved_bins"
            if (abs(np.log(frequencies[index] / target)) > np.median(steps) / 2 + 1e-8)
            else ""
        )
        add(
            "webb_instructions",
            index,
            point["low_db_pa2_hz"],
            point["high_db_pa2_hz"],
            target,
            reason,
        )

    # Brown is the reference at and above the split; Webb covers lower bins.
    # A bin already compared with a Webb point (the 0.2 Hz bin) is not
    # compared again, so every bin is judged against one source only.
    webb_bins = {
        int(np.argmin(abs(np.log(frequencies / row["reference_frequency_hz"]))))
        for row in diagnostics
        if row["reference"] == "webb_instructions" and row["comparison"] != "excluded"
    }
    brown_minimum = max(brown[0, 0], metadata["brown"].get("minimum_comparison_hz", 0))
    for index, frequency in enumerate(frequencies):
        if frequency < brown_minimum or index in webb_bins:
            continue
        if lower_edges[index] < brown[0, 0] or upper_edges[index] > brown[-1, 0]:
            # Show why the band has no comparison instead of extrapolating.
            add(
                "brown_figure8",
                index,
                0,
                0,
                frequency,
                "smoothing_bin_outside_digitized_reference",
            )
            diagnostics[-1]["reference_low_db_pa2_hz"] = None
            diagnostics[-1]["reference_high_db_pa2_hz"] = None
            diagnostics[-1]["median_offset_from_low_db"] = None
            diagnostics[-1]["median_offset_from_high_db"] = None
            continue
        # ObsPy smooths dB values over FFT frequencies. Smooth the digitized
        # curves over those same linearly spaced FFT bins before comparison.
        fft_step = band["fft_frequency_step_hz"]
        fft_frequencies = (
            np.arange(
                np.ceil(lower_edges[index] / fft_step),
                np.floor(upper_edges[index] / fft_step) + 1,
            )
            * fft_step
        )
        low, high = [
            float(
                np.mean(
                    np.interp(
                        np.log(fft_frequencies), np.log(brown[:, 0]), brown[:, column]
                    )
                )
            )
            for column in (1, 2)
        ]
        add(
            "brown_figure8",
            index,
            low,
            high,
            frequency,
            tolerance=metadata["brown"]["digitization_tolerance_db"],
        )

    included = [row for row in diagnostics if row["comparison"] != "excluded"]
    fraction = (
        float(np.mean([row["comparison"] == "within" for row in included]))
        if included
        else None
    )
    diagnostics.sort(
        key=lambda row: (row["reference"] != "webb_instructions", row["frequency_hz"])
    )
    return diagnostics, percentiles, fraction


def comparison_plot(
    ppsd, percentiles, diagnostics, band, response_f, response_db, metadata, brown
):
    figure, (axis, response_axis) = plt.subplots(
        2,
        1,
        figsize=(11, 8),
        sharex=True,
        layout="constrained",
        gridspec_kw={"height_ratios": [3, 1]},
    )
    frequencies = 1 / ppsd.period_bin_centers
    order = np.argsort(frequencies)
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
    webb = metadata["webb"]["points"] if metadata["webb"].get("use_in_plot") else []
    if webb:
        # Shading between the points guides the eye; only the points are compared.
        axis.fill_between(
            [p["frequency_hz"] for p in webb],
            [p["low_db_pa2_hz"] for p in webb],
            [p["high_db_pa2_hz"] for p in webb],
            color=REFERENCE_COLORS["webb_instructions"],
            alpha=0.12,
            linewidth=0,
        )
        axis.vlines(
            [p["frequency_hz"] for p in webb],
            [p["low_db_pa2_hz"] for p in webb],
            [p["high_db_pa2_hz"] for p in webb],
            color=REFERENCE_COLORS["webb_instructions"],
            linewidth=5,
            alpha=0.75,
            label="Webb 1998 Fig. 2 limits (≤0.2 Hz)",
        )
    used = brown[:, 0] >= metadata["brown"].get("minimum_comparison_hz", 0)
    axis.fill_between(
        brown[used, 0],
        brown[used, 1],
        brown[used, 2],
        color=REFERENCE_COLORS["brown_figure8"],
        alpha=0.25,
        label="Brown et al. 2014 limits (≥0.2 Hz), digitized",
    )
    included = [row for row in diagnostics if row["comparison"] != "excluded"]
    for name, (marker, color) in STATUS_MARKERS.items():
        selected = [row for row in included if row["comparison"] == name]
        if selected:
            axis.scatter(
                [r["frequency_hz"] for r in selected],
                [r["median_db_pa2_hz"] for r in selected],
                color=color,
                marker=marker,
                edgecolor="white",
                zorder=5,
                label=f"Median {name}",
            )
    axis.set_ylabel(r"Pressure PSD [dB re 1 Pa$^2$/Hz]")
    axis.set_title(
        f"{ppsd.id}: expected pressure levels ({len(ppsd.times_processed):,} windows)"
    )
    axis.legend(fontsize=8, loc="best")
    response_axis.plot(response_f, response_db, color="#555555")
    response_axis.axhline(
        -band["maximum_attenuation_db"], linestyle="--", color="black"
    )
    response_axis.set_ylabel("Response relative\nto reference [dB]")
    response_axis.set_ylim(-40, 5)
    response_axis.set_xlabel("Frequency [Hz]")
    for panel in (axis, response_axis):
        panel.set_xscale("log")
        panel.set_xlim(min(frequencies), max(frequencies))
        panel.axvspan(min(frequencies), band["minimum_hz"], color="gray", alpha=0.12)
        panel.axvspan(band["maximum_hz"], max(frequencies), color="gray", alpha=0.12)
        panel.grid(which="both", alpha=0.2)
    figure.suptitle(
        "Pressure-level screening; sensor band and reference bounds need review",
        fontsize=10,
    )
    return figure


def run_pressure_test(
    master_dir=OUTPUT_DIRECTORY,
    *,
    network="YV",
    station="RR52",
    sensor_min=None,
    sensor_max=None,
    loss_db=3,
    write_score=False,
    update_master=True,
    show=True,
    selected_rows=None,
    return_updates=False,
):
    """Reuse NPZ data; write diagnostics and patch only matching pressure rows."""
    if write_score and (sensor_min is None or sensor_max is None):
        raise ValueError("--write-score requires explicit sensor minimum and maximum")
    master = load_master_dataset(master_dir) if selected_rows is None else selected_rows
    required = {"ppsd_status", "ppsd_npz_file", "ppsd_units"}
    if not required.issubset(master):
        raise ValueError("Run ppsd.py first to populate saved pressure results")
    rows = master[
        master.Seed_id.str.startswith(f"{network}.{station}.")
        & master.ppsd_status.isin(["complete", "partial"])
        & master.ppsd_units.eq("Pa^2/Hz")
        & master.ppsd_npz_file.ne("")
    ].drop_duplicates(EPOCH_COLUMNS + ["inventory_file"])
    if rows.empty:
        raise ValueError(f"No saved pressure PPSDs for {network}.{station}")
    metadata, brown = load_references()
    summaries, updates, figures = [], [], []
    for _, row in rows.iterrows():
        source = Path(row.ppsd_npz_file).expanduser().resolve()
        ppsd = PPSD.load_npz(str(source))
        values = validate_pressure(ppsd, row)
        snapshot = source.parent / "inventory.xml"
        band, response_f, response_db = response_band(
            ppsd,
            snapshot,
            sensor_min=sensor_min,
            sensor_max=sensor_max,
            loss_db=loss_db,
        )
        band["fft_frequency_step_hz"] = ppsd.sampling_rate / ppsd.nfft
        diagnostics, percentiles, fraction = summarize(
            values,
            1 / ppsd.period_bin_centers,
            1 / ppsd.period_bin_right_edges,
            1 / ppsd.period_bin_left_edges,
            band,
            metadata,
            brown,
        )
        destination = source.parent / "pressure_expected"
        destination.mkdir(exist_ok=True)
        histogram_path = destination / "ppsd_with_limits.png"
        details_path = destination / "test_details.png"
        csv_path, json_path, png_path = [
            destination / name
            for name in (
                "frequency_comparison.csv",
                "comparison.json",
                "comparison.png",
            )
        ]
        with csv_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(diagnostics[0]))
            writer.writeheader()
            writer.writerows(diagnostics)
        included = [d for d in diagnostics if d["comparison"] != "excluded"]
        score = fraction if write_score else None
        summary = {
            "seed_id": ppsd.id,
            "psd_units": "Pa^2/Hz",
            "psd_windows": len(values),
            "sensor_band": band,
            "source_npz_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "response_snapshot_sha256": hashlib.sha256(
                snapshot.read_bytes()
            ).hexdigest(),
            "references": metadata,
            "tested_frequency_bins": len(included),
            "median_in_range_fraction": fraction,
            "score": score,
            "status": "screened" if score is not None else "comparison_only",
            "score_definition": (
                "Equal-weight fraction of included median bins inside bounds; "
                "no statistical significance or calibration acceptance implied."
            ),
            "reference_coverage": {
                name: sum(d["reference"] == name for d in included)
                for name in ("webb_instructions", "brown_figure8")
            },
            "near_digitization_boundaries": sum(
                d["near_digitization_boundary"] for d in included
            ),
            "limitations": [
                "Nominal response band is an operational estimate; user-supplied "
                "bounds are not independently validated by this code.",
                "Unverified handout points are excluded from current comparisons. "
                "Brown Figure 8 is digitized, not official numerical model data "
                "or a universal seafloor acceptance envelope.",
                "Overlapping PSD windows are correlated; "
                "occupancy fractions are descriptive.",
                "Three-hour smoothed spectra retain the existing processing "
                "settings, not Brown's original windowing.",
                "No complete reference coverage is claimed "
                "between discrete Webb points.",
            ],
            "files": {
                "frequency_csv": str(csv_path),
                "plot": str(png_path),
                "ppsd_with_limits_plot": str(histogram_path),
                "test_details_plot": str(details_path),
            },
        }
        json_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        figure = comparison_plot(
            ppsd,
            percentiles,
            diagnostics,
            band,
            response_f,
            response_db,
            metadata,
            brown,
        )
        figure.savefig(png_path, dpi=160)
        figures.append(figure)
        for factory, path, arguments in [
            (ppsd_limits_plot, histogram_path, (ppsd, diagnostics, band)),
            (test_details_plot, details_path, (ppsd, diagnostics)),
        ]:
            figure = factory(*arguments)
            figure.savefig(path, dpi=160)
            figures.append(figure)
        update = row[EPOCH_COLUMNS + ["inventory_file"]].to_dict() | {
            "pressure_expected_status": summary["status"],
            "pressure_expected_result_file": str(json_path),
            "pressure_expected_csv_file": str(csv_path),
            "pressure_expected_plot_file": str(png_path),
            "pressure_expected_ppsd_plot_file": str(histogram_path),
            "pressure_expected_details_plot_file": str(details_path),
            "pressure_expected_median_in_range_fraction": fraction,
            "pressure_expected_tested_bins": len(included),
            "pressure_expected_band_min_hz": band["minimum_hz"],
            "pressure_expected_band_max_hz": band["maximum_hz"],
            "pressure_expected_band_basis": band["basis"],
        }
        if write_score:
            update[SCORE_COLUMN] = score
            update[STATUS_COLUMN] = summary["status"]
        elif not row.get(SCORE_COLUMN, ""):
            update[STATUS_COLUMN] = "comparison_only"
        updates.append(update)
        summaries.append(summary)
        print(
            f"{ppsd.id}: {len(included)} bins compared; median in-range fraction "
            f"{fraction}; status={summary['status']}\n"
            f"Saved {png_path}\nSaved {histogram_path}\nSaved {details_path}",
            flush=True,
        )
    if update_master:
        update_pressure_test_results(master_dir, updates)
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
    parser.add_argument("--sensor-min-frequency", type=float)
    parser.add_argument("--sensor-max-frequency", type=float)
    parser.add_argument("--response-loss-db", type=float, default=3)
    parser.add_argument(
        "--write-score",
        action="store_true",
        help="Store the descriptive fraction as score; requires explicit sensor bounds",
    )
    parser.add_argument("--no-update-master", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    run_pressure_test(
        args.master_dir,
        network=args.network,
        station=args.station,
        sensor_min=args.sensor_min_frequency,
        sensor_max=args.sensor_max_frequency,
        loss_db=args.response_loss_db,
        write_score=args.write_score,
        update_master=not args.no_update_master,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()
