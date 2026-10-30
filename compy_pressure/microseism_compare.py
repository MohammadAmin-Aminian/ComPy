"""Compare one explicit completed baseline microseism campaign.

This module summarizes diagnostics that already exist.  It does not read new
waveforms, estimate pressure gains, apply the frozen Rayleigh gains, update a
master, or write a scientific score.
"""

import argparse
import csv
import json
import math
import shutil
from collections import Counter
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np

from compy_pressure.microseism_campaign import DEFAULT_MANIFEST, ROOT, load_manifest
from compy_pressure.microseism_diagnostic import (
    load_and_check_spectra,
    sha256,
    verify_saved_provenance,
)

NETWORK_ORDER = ("YV", "YL", "7D", "XO")
NETWORK_COLORS = {
    "YV": "#3b6fb6",
    "YL": "#d17a22",
    "7D": "#4c956c",
    "XO": "#9c5aa5",
}
CSV_FIELDS = [
    "network",
    "station",
    "pressure_seed_id",
    "vertical_seed_id",
    "start_date_utc",
    "end_date_utc",
    "retained_segments",
    "available_days",
    "selected_days",
    "selected_day_fraction",
    "median_coherence_squared",
    "ruan_high_coherence_bins",
    "valid_band_bins",
    "ruan_high_coherence_fraction",
    "median_nominal_transfer_amplitude_m_pa",
    "median_fixed_reference_transfer_amplitude_m_pa",
    "median_nominal_transfer_phase_deg",
    "phase_summary_method",
    "phase_family",
    "phase_distance_to_family_center_deg",
    "fixed_pressure_multiplier_g",
    "comparison_direction",
    "first_method_comparison_coefficient",
    "published_reference",
    "published_coefficient",
    "publication_comparison_status",
    "first_method_publication_difference_percent",
    "source_comparison",
    "diagnostic",
]


def resolve_project_path(path):
    """Resolve an absolute path or a path recorded relative to the project."""
    path = Path(path).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def parse_trailing_json(path):
    """Parse the final JSON object from one campaign station log."""
    text = Path(path).read_text()
    decoder = json.JSONDecoder()
    candidates = [0]
    candidates.extend(index + 1 for index in range(len(text)) if text[index] == "\n")
    parsed = []
    for start in candidates:
        if start >= len(text) or text[start] != "{":
            continue
        try:
            value, end = decoder.raw_decode(text, start)
        except json.JSONDecodeError:
            continue
        if not text[end:].strip() and isinstance(value, dict):
            parsed.append(value)
    if len(parsed) != 1:
        raise ValueError(
            f"Expected one trailing JSON object in {path}, found {len(parsed)}"
        )
    return parsed[0]


def wrapped_phase_degrees(value):
    """Wrap degrees into [-180, 180)."""
    return float((float(value) + 180) % 360 - 180)


def phase_family(value):
    """Nearest-centre descriptor only; not a coherence or hardware verdict."""
    if value is None or not np.isfinite(value):
        return "indeterminate", None
    wrapped = wrapped_phase_degrees(value)
    zero_distance = abs(wrapped)
    reversed_distance = abs(180 - abs(wrapped))
    if zero_distance <= reversed_distance:
        return "near_0_deg", zero_distance
    return "near_180_deg", reversed_distance


def circular_mean_degrees(values):
    """Return the circular mean in [-180, 180)."""
    radians = np.deg2rad(np.asarray(values, dtype=float))
    return wrapped_phase_degrees(
        np.rad2deg(np.arctan2(np.mean(np.sin(radians)), np.mean(np.cos(radians))))
    )


def geometric_spread_factor(values):
    """Return exp(population standard deviation of log(values))."""
    values = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("Transfer amplitudes must be finite and positive")
    return float(np.exp(np.std(np.log(values))))


def _check_hash(path, expected, protected):
    path = resolve_project_path(path)
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"Protected artifact hash changed: {path}")
    previous = protected.setdefault(str(path), actual)
    if previous != actual:
        raise ValueError(f"Conflicting protected hashes for {path}")
    return path


def _station_row(diagnostic, diagnostic_path):
    results = diagnostic["results"]
    ruan = results["ruan_high_coherence_reference"]
    daily = results["saikia_daily_reference"]
    context = diagnostic["comparison_context"]
    first = context["first_method"]
    publication = context["publication"]
    raw_phase = results["median_nominal_transfer_phase_deg"]
    phase = None if raw_phase is None else float(raw_phase)
    family, distance = phase_family(phase)
    selected = int(daily["selected_days"])
    available = int(daily["available_days"])
    return {
        "network": diagnostic["station"].split(".")[0],
        "station": diagnostic["station"],
        "pressure_seed_id": diagnostic["pressure_seed_id"],
        "vertical_seed_id": diagnostic["vertical_seed_id"],
        "start_date_utc": diagnostic["input_interval_utc"][0][:10],
        "end_date_utc": diagnostic["input_interval_utc"][1][:10],
        "retained_segments": int(diagnostic["retained_segments"]),
        "available_days": available,
        "selected_days": selected,
        "selected_day_fraction": selected / available,
        "median_coherence_squared": float(results["median_coherence_squared"]),
        "ruan_high_coherence_bins": int(ruan["bins"]),
        "valid_band_bins": int(results["valid_bins"]),
        "ruan_high_coherence_fraction": float(ruan["fraction_of_band"]),
        "median_nominal_transfer_amplitude_m_pa": float(
            results["median_nominal_transfer_amplitude_m_pa"]
        ),
        "median_fixed_reference_transfer_amplitude_m_pa": float(
            results["median_fixed_reference_transfer_amplitude_m_pa"]
        ),
        "median_nominal_transfer_phase_deg": phase,
        "phase_summary_method": results.get(
            "phase_summary_method", "legacy_linear_median"
        ),
        "phase_family": family,
        "phase_distance_to_family_center_deg": distance,
        "fixed_pressure_multiplier_g": float(first["pressure_multiplier_g"]),
        "comparison_direction": first["comparison_direction"],
        "first_method_comparison_coefficient": float(first["comparison_coefficient"]),
        "published_reference": publication["reference"],
        "published_coefficient": publication["coefficient"],
        "publication_comparison_status": publication["status"],
        "first_method_publication_difference_percent": publication[
            "first_method_difference_percent"
        ],
        "source_comparison": diagnostic["source_comparison"],
        "diagnostic": str(Path(diagnostic_path).resolve()),
    }


def load_campaign(log_root, manifest_path=DEFAULT_MANIFEST, verify_spectra=True):
    """Load one explicit campaign and verify its station/artifact bindings."""
    log_root = Path(log_root).resolve()
    summary_path = log_root / "summary.tsv"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Campaign summary does not exist: {summary_path}")
    manifest_path = Path(manifest_path).resolve()
    manifest = load_manifest(manifest_path)
    expected = {f"{row['network']}.{row['station']}": row for row in manifest}
    with summary_path.open(newline="") as handle:
        entries = list(csv.DictReader(handle, delimiter="\t"))
    names = [entry["station"] for entry in entries]
    if len(names) != len(set(names)):
        raise ValueError("Campaign summary contains duplicate stations")
    missing = sorted(set(expected) - set(names))
    extra = sorted(set(names) - set(expected))
    if missing or extra:
        raise ValueError(
            f"Campaign/manifest station mismatch; missing={missing}, extra={extra}"
        )
    incomplete = [
        entry["station"] for entry in entries if entry["status"] != "complete"
    ]
    if incomplete:
        raise ValueError(f"Campaign has incomplete stations: {incomplete}")

    protected = {
        str(manifest_path): sha256(manifest_path),
        str(summary_path): sha256(summary_path),
    }
    by_station = {entry["station"]: entry for entry in entries}
    rows = []
    for name, manifest_row in expected.items():
        log_path = resolve_project_path(by_station[name]["log"])
        protected[str(log_path)] = sha256(log_path)
        launch = parse_trailing_json(log_path)
        if launch.get("station") != name:
            raise ValueError(f"Log station binding differs for {name}")
        for key in ("pressure_seed_id", "vertical_seed_id"):
            if launch.get(key) != manifest_row[key]:
                raise ValueError(f"Log {key} binding differs for {name}")
        expected_gain = float(manifest_row["fixed_pressure_multiplier"])
        if not math.isclose(launch["fixed_pressure_multiplier"], expected_gain):
            raise ValueError(f"Log fixed gain differs for {name}")

        diagnostic_path = resolve_project_path(launch["diagnostic"])
        protected[str(diagnostic_path)] = sha256(diagnostic_path)
        diagnostic = json.loads(diagnostic_path.read_text())
        if (
            diagnostic.get("station") != name
            or diagnostic.get("status") != "diagnostic_only"
        ):
            raise ValueError(f"Diagnostic binding/status differs for {name}")
        if diagnostic["source_comparison"] != str(
            resolve_project_path(launch["paired_comparison"])
        ):
            raise ValueError(f"Source comparison binding differs for {name}")
        if any(
            diagnostic.get(key) is not False
            for key in (
                "pressure_gain_fit_performed",
                "rayleigh_gain_refit_performed",
                "fixed_reference_applied_to_original_data",
            )
        ):
            raise ValueError(f"Diagnostic mutation/fitting flag differs for {name}")
        if diagnostic["identifiability"]["conclusion"] != (
            "not_identifiable_from_saved_PZ_spectra_alone"
        ):
            raise ValueError(f"Unexpected identifiability claim for {name}")
        context = diagnostic["comparison_context"]
        if context.get("microseism_gain_estimated") is not False:
            raise ValueError(f"Unexpected microseism gain claim for {name}")
        _check_hash(
            context["manifest"]["path"], context["manifest"]["sha256"], protected
        )
        first = context["first_method"]
        if not math.isclose(first["pressure_multiplier_g"], expected_gain):
            raise ValueError(f"Diagnostic fixed gain differs for {name}")
        if first["comparison_direction"] != manifest_row["comparison_direction"]:
            raise ValueError(f"Diagnostic comparison direction differs for {name}")

        for record in diagnostic["provenance"]["saved_input_artifacts"].values():
            _check_hash(record["path"], record["sha256"], protected)
        for record in diagnostic["output_files"].values():
            _check_hash(record["path"], record["sha256"], protected)

        comparison_path = resolve_project_path(diagnostic["source_comparison"])
        comparison = json.loads(comparison_path.read_text())
        verify_saved_provenance(comparison_path, comparison)
        if verify_spectra:
            load_and_check_spectra(
                comparison["files"]["spectra"], diagnostic["retained_segments"]
            )
        rows.append(_station_row(diagnostic, diagnostic_path))
    return rows, protected, summary_path


def summarize_group(rows):
    """Calculate transparent descriptive statistics for one station group."""
    nominal = [row["median_nominal_transfer_amplitude_m_pa"] for row in rows]
    fixed = [row["median_fixed_reference_transfer_amplitude_m_pa"] for row in rows]
    nominal_spread = geometric_spread_factor(nominal)
    fixed_spread = geometric_spread_factor(fixed)
    families = Counter(row["phase_family"] for row in rows)
    selected = sum(row["selected_days"] for row in rows)
    available = sum(row["available_days"] for row in rows)
    return {
        "stations": len(rows),
        "retained_segments": {
            "minimum": min(row["retained_segments"] for row in rows),
            "median": float(np.median([row["retained_segments"] for row in rows])),
            "maximum": max(row["retained_segments"] for row in rows),
        },
        "median_coherence_squared": {
            "minimum": min(row["median_coherence_squared"] for row in rows),
            "median": float(
                np.median([row["median_coherence_squared"] for row in rows])
            ),
            "maximum": max(row["median_coherence_squared"] for row in rows),
        },
        "stations_with_weekly_median_coherence_at_least_0_8": sum(
            row["median_coherence_squared"] >= 0.8 for row in rows
        ),
        "stations_with_any_ruan_gt_0_95_bins": sum(
            row["ruan_high_coherence_bins"] > 0 for row in rows
        ),
        "median_ruan_gt_0_95_fraction": float(
            np.median([row["ruan_high_coherence_fraction"] for row in rows])
        ),
        "saikia_daily_reference": {
            "selected_days": selected,
            "available_days": available,
            "selected_fraction": selected / available,
            "stations_with_none": sum(row["selected_days"] == 0 for row in rows),
            "stations_with_any": sum(row["selected_days"] > 0 for row in rows),
            "stations_with_all": sum(
                row["selected_days"] == row["available_days"] for row in rows
            ),
        },
        "phase": {
            "circular_mean_deg": circular_mean_degrees(
                [row["median_nominal_transfer_phase_deg"] for row in rows]
            ),
            "family_counts": dict(families),
        },
        "amplitude_dispersion": {
            "nominal_geometric_spread_factor": nominal_spread,
            "fixed_reference_geometric_spread_factor": fixed_spread,
            "fixed_vs_nominal_spread_change_percent": 100
            * (fixed_spread / nominal_spread - 1),
            "interpretation": (
                "descriptive cross-station spread only; water depth, structure, "
                "vertical response, and sampled wavefields are not controlled"
            ),
        },
    }


def publication_statistics(rows):
    """Summarize only publication coefficients already frozen in the manifest."""
    matched = [
        row
        for row in rows
        if row["published_coefficient"] is not None
        and row["first_method_publication_difference_percent"] is not None
    ]

    def stats(selected):
        differences = [
            row["first_method_publication_difference_percent"] for row in selected
        ]
        return {
            "stations": len(selected),
            "difference_range_percent": [min(differences), max(differences)],
            "median_difference_percent": float(np.median(differences)),
        }

    yv = [row for row in matched if row["network"] == "YV"]
    yl = [row for row in matched if row["network"] == "YL"]
    yl_peers = [row for row in yl if row["station"] != "YL.C07W"]
    return {
        "YV_Aminian2025_same_g_direction": stats(yv),
        "YL_ZhaWebb2016_reciprocal_direction_all": stats(yl),
        "YL_ZhaWebb2016_reciprocal_direction_excluding_C07W": stats(yl_peers),
        "unmatched_status_counts": dict(
            Counter(
                row["publication_comparison_status"]
                for row in rows
                if row not in matched
            )
        ),
        "scope": (
            "These differences compare the frozen first-method Rayleigh result "
            "with a published coefficient; they are not microseism gain estimates."
        ),
    }


def summarize(rows):
    """Build global, per-network, and paper-reference summaries."""
    networks = {
        network: summarize_group([row for row in rows if row["network"] == network])
        for network in NETWORK_ORDER
    }
    ranked = sorted(rows, key=lambda row: row["median_coherence_squared"])
    return {
        "all_stations": summarize_group(rows),
        "networks": networks,
        "highest_median_coherence": [
            {
                "station": row["station"],
                "median_coherence_squared": row["median_coherence_squared"],
            }
            for row in reversed(ranked[-5:])
        ],
        "lowest_median_coherence": [
            {
                "station": row["station"],
                "median_coherence_squared": row["median_coherence_squared"],
            }
            for row in ranked[:5]
        ],
        "publication": publication_statistics(rows),
    }


def plot_comparison(rows):
    """Plot the station-level descriptive comparison."""
    x = np.arange(len(rows))
    colors = [NETWORK_COLORS[row["network"]] for row in rows]
    figure, axes = plt.subplots(4, 1, figsize=(16, 13), sharex=True)
    axes[0].scatter(
        x, [row["median_coherence_squared"] for row in rows], c=colors, s=32
    )
    axes[0].axhline(0.8, color="0.25", linestyle="--", linewidth=0.9)
    axes[0].set_ylabel("Weekly median\ncoherence²")
    axes[0].set_ylim(0.5, 1.01)

    axes[1].scatter(
        x, [row["ruan_high_coherence_fraction"] for row in rows], c=colors, s=32
    )
    axes[1].set_ylabel("Fraction of band\nwith coherence² > 0.95")
    axes[1].set_ylim(-0.03, 1.03)

    axes[2].scatter(
        x, [row["median_nominal_transfer_phase_deg"] for row in rows], c=colors, s=32
    )
    for reference in (-180, 0, 180):
        axes[2].axhline(reference, color="0.5", linewidth=0.7)
    axes[2].set_ylabel("Median D/P\nphase (degrees)")
    axes[2].set_ylim(-190, 190)

    nominal = [row["median_nominal_transfer_amplitude_m_pa"] for row in rows]
    fixed = [row["median_fixed_reference_transfer_amplitude_m_pa"] for row in rows]
    for index, (before, after) in enumerate(zip(nominal, fixed, strict=True)):
        axes[3].plot([index, index], [before, after], color="0.72", linewidth=0.8)
    axes[3].scatter(
        x,
        nominal,
        facecolors="none",
        edgecolors=colors,
        s=36,
        label="nominal transfer",
    )
    axes[3].scatter(x, fixed, c=colors, marker="x", s=34, label="fixed-g display")
    axes[3].set_yscale("log")
    axes[3].set_ylabel("Median |D/P|\n(m/Pa)")
    axes[3].legend(loc="upper right", fontsize=8)

    boundaries = []
    for index in range(1, len(rows)):
        if rows[index]["network"] != rows[index - 1]["network"]:
            boundaries.append(index - 0.5)
    for axis in axes:
        for boundary in boundaries:
            axis.axvline(boundary, color="0.82", linewidth=0.8)
        axis.grid(axis="y", alpha=0.2)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([row["station"] for row in rows], rotation=90, fontsize=7)
    figure.suptitle(
        "Baseline microseism P–Z diagnostics (0.1–0.2 Hz)\n"
        "Frozen Rayleigh gains shown only as H_fixed = H_nominal / g"
    )
    figure.tight_layout()
    return figure


def markdown_report(rows, statistics, source_summary):
    """Render the concise scientific comparison report."""
    overall = statistics["all_stations"]
    daily = overall["saikia_daily_reference"]
    phase_counts = overall["phase"]["family_counts"]
    publication = statistics["publication"]
    yv_paper = publication["YV_Aminian2025_same_g_direction"]
    yl_paper = publication["YL_ZhaWebb2016_reciprocal_direction_excluding_C07W"]
    top = ", ".join(
        f"{item['station']} ({item['median_coherence_squared']:.3f})"
        for item in statistics["highest_median_coherence"]
    )
    bottom = ", ".join(
        f"{item['station']} ({item['median_coherence_squared']:.3f})"
        for item in statistics["lowest_median_coherence"]
    )
    rr52 = next(row for row in rows if row["station"] == "YV.RR52")
    rr52_nominal = rr52["median_nominal_transfer_amplitude_m_pa"]
    rr52_fixed = rr52["median_fixed_reference_transfer_amplitude_m_pa"]
    weekly_count = overall["stations_with_weekly_median_coherence_at_least_0_8"]
    lines = [
        "# Baseline microseism pressure–vertical comparison",
        "",
        f"Explicit campaign: `{source_summary}`.",
        "",
        "## Result",
        "",
        (
            f"All {len(rows)} frozen-baseline stations completed. The saved P/Z "
            "diagnostics provide coherence, phase and transfer-amplitude checks, "
            "but **no independent microseism pressure gain**: "
            "`H_nominal=(g_pressure/g_vertical)*eta(f; water, structure)`. "
            "The Rayleigh multiplier is fixed and only changes the displayed "
            "transfer through `H_fixed=H_nominal/g`."
        ),
        "",
        (
            f"Across the campaign, {weekly_count}/{len(rows)} stations have weekly "
            "median coherence² at least 0.8, and "
            f"{overall['stations_with_any_ruan_gt_0_95_bins']}/{len(rows)} contain "
            "at least one 0.1–0.2 Hz bin above the Ruan et al. 0.95 reference. "
            f"The later Saikia daily reference selects {daily['selected_days']}/"
            f"{daily['available_days']} station-days across "
            f"{daily['stations_with_any']}/{len(rows)} stations. These are descriptive "
            "comparisons, not BRUIT-FM pass scores."
        ),
        "",
        f"Highest weekly medians: {top}. Lowest: {bottom}.",
        "",
        (
            f"RR52, the pilot reference, is close to the YV network median: "
            f"coherence² {rr52['median_coherence_squared']:.3f}, "
            f"{rr52['ruan_high_coherence_bins']}/{rr52['valid_band_bins']} "
            f"Ruan-reference bins, {rr52['selected_days']}/"
            f"{rr52['available_days']} selected days and phase "
            f"{rr52['median_nominal_transfer_phase_deg']:+.1f}°. Its nominal/fixed "
            f"display amplitudes remain {rr52_nominal:.3e}/{rr52_fixed:.3e} m/Pa."
        ),
        "",
        "## Network comparison",
        "",
        (
            "| Network | Stations | Segments min/median/max | Median coherence² "
            "(network median) | Weekly median ≥0.8 | Any bins >0.95 | Selected "
            "days | Phase center | Nominal → fixed amplitude spread |"
        ),
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for network in NETWORK_ORDER:
        item = statistics["networks"][network]
        segments = item["retained_segments"]
        item_daily = item["saikia_daily_reference"]
        dispersion = item["amplitude_dispersion"]
        lines.append(
            f"| {network} | {item['stations']} | {segments['minimum']}/"
            f"{segments['median']:.0f}/{segments['maximum']} | "
            f"{item['median_coherence_squared']['median']:.3f} | "
            f"{item['stations_with_weekly_median_coherence_at_least_0_8']} | "
            f"{item['stations_with_any_ruan_gt_0_95_bins']} | "
            f"{item_daily['selected_days']}/{item_daily['available_days']} | "
            f"{item['phase']['circular_mean_deg']:+.1f}° | "
            f"×{dispersion['nominal_geometric_spread_factor']:.3f} → "
            f"×{dispersion['fixed_reference_geometric_spread_factor']:.3f} |"
        )
    lines.extend(
        [
            "",
            (
                f"The phase separates into {phase_counts.get('near_0_deg', 0)} "
                f"near-0° and {phase_counts.get('near_180_deg', 0)} near-180° "
                "stations. The split is by network: positive scalar gains cannot "
                "cause it. Treat it as an unresolved pressure/vertical sign or "
                "response-convention difference until metadata and processing "
                "conventions are audited; do not pool signed admittance curves yet."
            ),
            "",
            (
                "The nominal-to-fixed spread factors are descriptive only. A lower "
                "cross-station spread after displaying the frozen gains is not an "
                "independent validation because water depth, seafloor structure, "
                "vertical response and wavefield differ among stations."
            ),
            "",
            "## First method and papers",
            "",
            (
                f"The already-frozen first-method values retain their prior paper "
                f"comparisons: {yv_paper['stations']} YV values differ from "
                "[Aminian et al. (2025)]"
                "(https://www.ipgp.fr/~stutz/"
                "2025-aminian-crawford-stutzmann-al-GJI.pdf) by "
                f"{yv_paper['difference_range_percent'][0]:+.2f}% "
                f"to {yv_paper['difference_range_percent'][1]:+.2f}% (median "
                f"{yv_paper['median_difference_percent']:+.2f}%) in the same `g` "
                f"direction. The {yl_paper['stations']} YL peers excluding C07W "
                "differ from [Zha & Webb (2016)]"
                "(https://doi.org/10.1002/2015JB012688) by "
                f"{yl_paper['difference_range_percent'][0]:+.2f}% to "
                f"{yl_paper['difference_range_percent'][1]:+.2f}% (median "
                f"{yl_paper['median_difference_percent']:+.2f}%) in reciprocal "
                "`1/g`; C07W remains the known mixed-state outlier at +86.76%. "
                "These numbers come from method 1, not from a microseism fit."
            ),
            "",
            (
                "[Ruan et al. (2014)](https://doi.org/10.1002/2014JB011162) "
                "is represented only by its >0.95 coherence "
                "reference: reproducing its 7D scale factors requires the coupled "
                "modal structure/water model and joint inversion. "
                "[Deng et al. (2022)](https://doi.org/10.1029/2021GL096952) "
                "uses frequency-dependent body-wave bands for XO, so no scalar "
                "microseism/paper difference is reported. The "
                "[Saikia et al. (2020)](https://doi.org/10.1007/s11001-020-09402-0) "
                "daily >0.8 rule is shown only as a selection reference."
            ),
            "",
            "Primary-source interpretation and links are retained in "
            "[the RR52 method note](rr52_microseism_diagnostic.md); the frozen "
            "publication-direction rules are in "
            "[the campaign note](microseism_baseline_run.md).",
            "",
            "## Limits and preservation",
            "",
            "- Campaign intervals were selected independently and are not "
            "synchronized; "
            "network/source-season differences remain.",
            "- Coverage ranges from 24 to 168 hourly segments, so short records should "
            "not be ranked as equivalent week-long estimates.",
            "- The diagnostics contain no horizontal/tilt correction, local-event/coda "
            "screen, modal prediction, gain uncertainty or classification.",
            "- Every explicit log, diagnostic, saved artifact and stored aggregate was "
            "verified. No waveform, response, master, prior result or score was "
            "changed; "
            "no fixed gain was applied to original data.",
            "",
            "Machine-readable station values: "
            "[CSV](microseism_baseline_comparison.csv) and "
            "[JSON](microseism_baseline_comparison.json). "
            "[Comparison figure](figures/microseism_baseline_comparison.png).",
            "",
        ]
    )
    return "\n".join(lines)


def write_csv(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def run_comparison(
    log_root,
    output_directory,
    *,
    manifest_path=DEFAULT_MANIFEST,
    curated_prefix=None,
    curated_plot=None,
):
    """Verify, summarize and plot one explicit campaign."""
    output = Path(output_directory).resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows, protected, summary_path = load_campaign(log_root, manifest_path)
    statistics = summarize(rows)

    csv_path = output / "comparison.csv"
    plot_path = output / "comparison.png"
    report_path = output / "comparison.md"
    json_path = output / "comparison.json"
    write_csv(csv_path, rows)
    figure = plot_comparison(rows)
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    report_path.write_text(markdown_report(rows, statistics, summary_path))
    audit = {
        "status": "diagnostic_comparison_only",
        "source_campaign_summary": str(summary_path),
        "station_count": len(rows),
        "microseism_gain_estimated": False,
        "fixed_gains_refit": False,
        "fixed_gains_applied_to_original_data": False,
        "master_or_score_updated": False,
        "identifiability": (
            "H_nominal(f)=(g_pressure/g_vertical)*eta(f; structure, water); "
            "pressure gain is not identifiable from saved P/Z spectra alone"
        ),
        "threshold_scope": (
            "Ruan >0.95 bin and Saikia daily mean >0.8 are literature comparison "
            "references, not adopted BRUIT-FM scores"
        ),
        "statistics": statistics,
        "stations": rows,
        "protected_sha256": protected,
        "protected_inputs_preserved": True,
        "output_files": {
            "csv": {"path": str(csv_path), "sha256": sha256(csv_path)},
            "plot": {"path": str(plot_path), "sha256": sha256(plot_path)},
            "report": {"path": str(report_path), "sha256": sha256(report_path)},
        },
    }
    json_path.write_text(json.dumps(audit, indent=2) + "\n")

    changed = [path for path, digest in protected.items() if sha256(path) != digest]
    if changed:
        raise RuntimeError(f"Protected inputs changed during comparison: {changed}")
    if curated_prefix:
        curated_prefix = Path(curated_prefix)
        curated_prefix.parent.mkdir(parents=True, exist_ok=True)
        targets = {
            json_path: curated_prefix.with_suffix(".json"),
            csv_path: curated_prefix.with_suffix(".csv"),
            report_path: curated_prefix.with_suffix(".md"),
        }
        for source, target in targets.items():
            if target.exists():
                raise FileExistsError(f"Refusing to replace curated output {target}")
            shutil.copyfile(source, target)
    if curated_plot:
        curated_plot = Path(curated_plot)
        curated_plot.parent.mkdir(parents=True, exist_ok=True)
        if curated_plot.exists():
            raise FileExistsError(f"Refusing to replace curated output {curated_plot}")
        shutil.copyfile(plot_path, curated_plot)
    return json_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--log-root",
        type=Path,
        required=True,
        help="Exact campaign log directory containing summary.tsv",
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--curated-prefix", type=Path)
    parser.add_argument("--curated-plot", type=Path)
    args = parser.parse_args()
    output = args.output_dir or (
        ROOT
        / "results/VALIDATION/MICROSEISM_BASELINE_COMPARISON"
        / f"run_{uuid4().hex[:12]}"
    )
    result = run_comparison(
        args.log_root,
        output,
        manifest_path=args.manifest,
        curated_prefix=args.curated_prefix,
        curated_plot=args.curated_plot,
    )
    audit = json.loads(result.read_text())
    print(json.dumps(audit["statistics"], indent=2))
    print(f"Output: {result.parent}")


if __name__ == "__main__":
    main()
