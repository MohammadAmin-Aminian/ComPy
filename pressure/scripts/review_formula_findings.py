"""Re-check the 2026-10-19 formula-review findings on saved products (read-only).

Writes CSVs to --output-dir only; never changes saved results, masters or scores.

1. Linear vs circular band phase in saved microseism products: the 46-station
   baseline diagnostics and every YV station in each completed all-methods week.
2. Godin-predicted depth bias inside saved nearby-station comparisons.
3. Godin comparison bins reaching below 2 mHz (saved PPSD comparison and pilot).

    PYTHONPATH=.. python scripts/review_formula_findings.py
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from compy_pressure.infragravity_test import godin_pressure_db

ROOT = Path(__file__).resolve().parents[1]
GRAVITY = 9.81
ALL_METHODS = "EXPERIMENTS/THREE_METHODS/run_20261016T144248Z_2515/microseism"
# Completed intervals of the interrupted all-methods run (local_p_integration_closeout).
INTERVALS = {
    "base": "run_20261016T144432Z_536542f9",
    "-12w": "run_20261016T145753Z_90e7e1bb",
    "-6w": "run_20261016T151040Z_c36ae13d",
    "-2w": "run_20261016T152342Z_601f267d",
}


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def wrap(degrees):
    return (np.asarray(degrees) + 180.0) % 360.0 - 180.0


def circular_median(phase_deg):
    """Median of deviations from the circular mean; also the resultant length."""
    z = np.mean(np.exp(1j * np.radians(phase_deg)))
    centre = np.degrees(np.angle(z))
    return float(wrap(centre + np.median(wrap(phase_deg - centre)))), float(abs(z))


def family(value):
    return "near_0" if abs(wrap(value)) <= 90 else "near_180"


def phase_row(label, station, band, phase, coherence, saved):
    linear = float(np.median(phase))
    circular, resultant = circular_median(phase)
    return {
        "source": label,
        "station": station,
        "band": band,
        "bins": len(phase),
        "median_coherence_squared": round(float(np.median(coherence)), 3),
        "saved_linear_median_deg": None if saved is None else round(saved, 2),
        "recomputed_linear_median_deg": round(linear, 2),
        "circular_median_deg": round(circular, 2),
        "resultant_length": round(resultant, 3),
        "linear_minus_circular_deg": round(float(wrap(linear - circular)), 2),
        "family_linear": family(linear),
        "family_circular": family(circular),
        "family_changed": family(linear) != family(circular),
        "fraction_bins_positive": round(float(np.mean(phase > 0)), 3),
    }


def load_spectra(path):
    with np.load(path) as data:
        valid = data["valid"].astype(bool)
        return (
            data["frequency_hz"],
            np.angle(data["transfer_dp"], deg=True),
            data["coherence_squared"],
            valid,
        )


def baseline_phases(baseline_table):
    """Codex microseism diagnostics behind docs/microseism_baseline_comparison.csv."""
    rows = []
    with Path(baseline_table).open() as stream:
        for item in csv.DictReader(stream):
            diagnostic = json.loads(Path(item["diagnostic"]).read_text())
            spectra = diagnostic["provenance"]["saved_input_artifacts"]["spectra"]
            f, phase, coherence, valid = load_spectra(spectra["path"])
            mask = valid & (f >= 0.1) & (f <= 0.2)
            rows.append(
                phase_row(
                    "baseline_46",
                    item["station"],
                    "0.1-0.2Hz",
                    phase[mask],
                    coherence[mask],
                    float(item["median_nominal_transfer_phase_deg"]),
                )
            )
    return rows


def all_methods_phases(microseism, network):
    rows = []
    for label, run in INTERVALS.items():
        for path in sorted((microseism / run / "stations").glob(f"{network}.*.json")):
            item = json.loads(path.read_text())
            if item.get("status") != "checked" or not item.get("paired_result"):
                continue
            f, phase, coherence, valid = load_spectra(
                Path(item["paired_result"]).with_name("spectra.npz")
            )
            low, high = item["long_wave_hz"]
            for band, mask, saved in (
                (
                    "0.1-0.2Hz",
                    (f >= 0.1) & (f <= 0.2),
                    item.get("literature_band_phase_deg"),
                ),
                (
                    "long_wave",
                    (f >= low) & (f <= high),
                    item.get("long_wave_phase_deg"),
                ),
            ):
                mask &= valid
                if mask.any():
                    rows.append(
                        phase_row(
                            f"all_methods_{label}",
                            item["station"],
                            band,
                            phase[mask],
                            coherence[mask],
                            saved,
                        )
                    )
    return rows


def pressure_folders(results, network):
    return sorted(Path(results, "NETWORKS/TEST", network).glob("*/*DH__*"))


def bin_model(frequency, depth, fft_step):
    """ObsPy-style dB average of Godin over FFT bins of a one-octave period bin."""
    low, high = frequency / np.sqrt(2), frequency * np.sqrt(2)
    grid = np.arange(np.ceil(low / fft_step), np.floor(high / fft_step) + 1)
    grid = grid[grid > 0] * fft_step
    return float(np.mean(godin_pressure_db(grid, depth)))


def nearby_depth_bias(results, network):
    """Godin depth offset expected in compared IG-band bins (f <= sqrt(g/2piH))."""
    depths = {}
    for folder in pressure_folders(results, network):
        summary = folder / "infragravity/comparison.json"
        if summary.is_file():
            depths[folder.name.split("__")[0]] = json.loads(summary.read_text())[
                "depth_m"
            ]
    rows = []
    for folder in pressure_folders(results, network):
        path = folder / "nearby/comparison.json"
        if not path.is_file():
            continue
        summary = json.loads(path.read_text())
        target = depths.get(summary["seed_id"])
        references = [
            r
            for r in summary.get("paired_references", [])
            if r["used_in_reference_group"]
        ]
        row = {"station": folder.parent.name, "status": summary["status"]}
        if not references or target is None:
            rows.append(row | {"reference_stations": 0})
            continue
        reference_depths = [depths.get(r["epoch"]["Seed_id"]) for r in references]
        step = summary["response_band"]["fft_minimum_hz"]
        cutoff = float(np.sqrt(GRAVITY / (2 * np.pi * target)))
        with (folder / "nearby/frequency_comparison.csv").open() as stream:
            compared = [r for r in csv.DictReader(stream) if r["difference_median_db"]]
        in_band = [r for r in compared if float(r["frequency_hz"]) <= cutoff]
        predicted = [
            bin_model(float(r["frequency_hz"]), target, step)
            - float(
                np.median(
                    [
                        bin_model(float(r["frequency_hz"]), d, step)
                        for d in reference_depths
                        if d
                    ]
                )
            )
            for r in in_band
        ]
        rows.append(
            row
            | {
                "reference_stations": len(references),
                "target_depth_m": target,
                "reference_depths_m": ";".join(str(d) for d in reference_depths),
                "lowest_compared_hz": round(
                    min(float(r["frequency_hz"]) for r in compared), 5
                ),
                "ig_cutoff_hz": round(cutoff, 5),
                "compared_bins": len(compared),
                "compared_bins_in_ig_band": len(in_band),
                "godin_depth_bias_db_min": round(min(predicted), 2)
                if predicted
                else None,
                "godin_depth_bias_db_max": round(max(predicted), 2)
                if predicted
                else None,
            }
        )
    return rows


def godin_low_bins(results, network):
    rows = []
    for folder in pressure_folders(results, network):
        tables = [("saved_ppsd", folder / "infragravity/frequency_comparison.csv")]
        tables += [
            ("raw_pilot", path)
            for path in sorted(folder.glob("godin_pilot/*/frequency_comparison.csv"))
        ]
        for product, table in tables:
            if not table.is_file():
                continue
            with table.open() as stream:
                compared = [
                    r for r in csv.DictReader(stream) if r["comparison"] == "compared"
                ]
            if product == "saved_ppsd":
                low = [float(r["smoothing_minimum_hz"]) < 0.002 for r in compared]
            else:
                low = [float(r["frequency_hz"]) < 0.002 for r in compared]
            rows.append(
                {
                    "station": folder.parent.name,
                    "product": product,
                    "compared_bins": len(compared),
                    "bins_reaching_below_2mHz": int(sum(low)),
                }
            )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=ROOT / "results")
    parser.add_argument(
        "--baseline-table", default=ROOT / "docs/microseism_baseline_comparison.csv"
    )
    parser.add_argument("--network", default="YV")
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "results/REVIEWS/FORMULA_REVIEW"
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    baseline = baseline_phases(args.baseline_table)
    weekly = all_methods_phases(args.results_dir / ALL_METHODS, args.network)
    nearby = nearby_depth_bias(args.results_dir, args.network)
    godin = godin_low_bins(args.results_dir, args.network)
    write_csv(args.output_dir / "phase_baseline_46_stations.csv", baseline)
    write_csv(args.output_dir / f"{args.network}_phase_all_methods_weeks.csv", weekly)
    write_csv(args.output_dir / f"{args.network}_nearby_depth_bias.csv", nearby)
    write_csv(args.output_dir / f"{args.network}_godin_below_2mhz.csv", godin)
    for name, rows in (("baseline 46", baseline), (f"{args.network} weeks", weekly)):
        error = np.abs([r["linear_minus_circular_deg"] for r in rows])
        coherent = [r for r in rows if r["median_coherence_squared"] >= 0.5]
        off = sum(abs(r["linear_minus_circular_deg"]) > 5 for r in coherent)
        changed = sum(r["family_changed"] for r in rows)
        print(
            f"{name}: {len(rows)} rows; max |linear-circular| {error.max():.1f} deg; "
            f"family changes {changed}; coherent rows >5 deg {off}/{len(coherent)}"
        )
    ig_bins = sum(r.get("compared_bins_in_ig_band", 0) for r in nearby)
    low_bins = sum(r["bins_reaching_below_2mHz"] for r in godin)
    print(
        f"{args.network} nearby: {ig_bins} compared IG-band bins; "
        f"Godin bins below 2 mHz: {low_bins}"
    )
    print(f"Saved to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
