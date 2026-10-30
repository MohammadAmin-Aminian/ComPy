"""Build descriptive station-balanced pressure noise percentiles from approved PPSDs.

The --plan mode requires no approvals and writes a selection audit only. Actual
model calculation requires an adopted policy, fingerprint-bound decisions and
explicit usable frequency/time intervals. Original spectra and masters are
never changed. No calibration multiplier is applied.
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime
from obspy.signal import PPSD

from compy_pressure.plot_network_ppsd import (  # noqa: F401 (station_distribution re-exported)
    accumulate_histogram,
    spectrum,
    station_distribution,
)
from compy_pressure.quality_review import (
    KEYS,
    ROOT,
    Inputs,
    fresh_output,
    validate_review,
    write_csv,
)


def quantile_curves(probability, db_edges, coverage, percentiles):
    if (
        len(percentiles) != 2
        or not np.isfinite(percentiles).all()
        or not 0 < percentiles[0] < percentiles[1] < 100
    ):
        raise ValueError("Require two increasing percentiles strictly inside (0, 100)")
    cdf = probability.cumsum(axis=0)
    levels = (db_edges[:-1] + db_edges[1:]) / 2
    curves = []
    for percentile in percentiles:
        index = np.argmax(cdf >= percentile / 100, axis=0)
        values = levels[index].astype(float)
        values[coverage == 0] = np.nan
        curves.append(values)
    return np.asarray(curves)


def selected_samples(ppsd, row):
    """Require full time windows and full smoothing bins within reviewed bounds."""
    native_row = {"Seed_id": row["Seed_id"], "ppsd_units": row["ppsd"]["units"]}
    f, values, pressure = spectrum(ppsd, native_row)
    if not pressure:
        raise ValueError("This noise model requires pressure in Pa^2/Hz")
    decision = row["review_decision"]
    start, end = map(UTCDateTime, decision["approved_time_range_utc"])
    times = np.asarray([float(t) for t in ppsd.times_processed])
    retained = (times >= float(start)) & (times + ppsd.ppsd_length <= float(end))
    low, high = decision["approved_frequency_band_hz"]
    order = np.argsort(1 / np.asarray(ppsd.period_bin_centers))
    left = (1 / np.asarray(ppsd.period_bin_right_edges))[order]
    right = (1 / np.asarray(ppsd.period_bin_left_edges))[order]
    keep = (left >= low) & (right <= high)
    indices = np.flatnonzero(keep)
    if not retained.any() or not len(indices):
        raise ValueError(f"No fully supported approved windows/bins: {row['Seed_id']}")
    if np.any(np.diff(indices) != 1):
        raise ValueError("Approved frequency interval is not contiguous")
    edges = (1 / np.asarray(ppsd.period_xedges))[::-1]
    edges = edges[indices[0] : indices[-1] + 2]
    return (
        f[keep],
        values[retained][:, keep],
        edges,
        (
            float(times[retained].min()),
            float(times[retained].max() + ppsd.ppsd_length),
        ),
    )


def grouping(row, group_by):
    net, _, _, channel = row["Seed_id"].split(".")
    return {
        "all": "all",
        "network": net,
        "channel": channel,
        "network_channel": net + "_" + channel,
    }[group_by]


def build_noise_model(
    review,
    output,
    *,
    root=ROOT,
    plan=False,
    group_by="network_channel",
    percentiles=(5.0, 95.0),
):
    if group_by not in {"all", "network", "channel", "network_channel"}:
        raise ValueError("Unknown noise grouping")
    # Validate percentile choices even in a planning run.
    quantile_curves(np.array([[1.0]]), np.array([0.0, 1.0]), np.array([1]), percentiles)
    inputs = Inputs(root)
    for code in (Path(__file__), ROOT / "../compy_pressure/plot_network_ppsd.py"):
        inputs.protect(code)
    report = inputs.document(review)
    rows = validate_review(report)
    for path, expected in report["protected_sha256"].items():
        inputs.protect(path, expected)
    selected = [row for row in rows if row["accepted_for_noise_model"]]
    selection = [
        {
            **{k: row[k] for k in KEYS},
            "row_id": row["row_id"],
            "accepted_for_noise_model": row["accepted_for_noise_model"],
            "reason": "explicit_scientific_approval"
            if row["accepted_for_noise_model"]
            else "not_approved_for_noise_model",
        }
        for row in rows
    ]
    summary = {
        "mode": "selection_plan" if plan else "descriptive_noise_model",
        "review": str(inputs.protect(review)),
        "policy_version": report["policy"]["version"],
        "policy_status": report["policy"]["status"],
        "channels_reviewed": len(rows),
        "channels_approved": len(selected),
        "percentiles": list(percentiles),
        "weighting": "equal stations per frequency; epochs pooled within station",
        "group_by": group_by,
        "gain_applied": False,
        "master_updated": False,
        "validated_global_noise_model": False,
        "models": [],
        "protected_inputs_unchanged": True,
    }
    if not plan and not selected:
        raise ValueError(
            "No approved channels. Use --plan; never fabricate scientific approval."
        )
    out = fresh_output(output, root)
    # Prepare all model data before creating outputs, so a bad approval fails closed.
    groups, ranges = defaultdict(dict), defaultdict(list)
    if not plan:
        for row in selected:
            p = row["ppsd"]
            inputs.protect(p["path"], p["sha256"])
            ppsd = PPSD.load_npz(p["path"])
            f, values, edges, interval = selected_samples(ppsd, row)
            group, station = grouping(row, group_by), row["station"]
            for previous in ranges[group, station]:
                if max(previous[0], interval[0]) < min(previous[1], interval[1]):
                    raise ValueError(
                        "Overlapping approved station records would double count data"
                    )
            ranges[group, station].append(interval)
            state = groups[group].setdefault(station, {})
            accumulate_histogram(state, f, values, edges)
    inputs.verify()
    out.mkdir(parents=True)
    write_csv(
        out / "selection.csv",
        selection,
        list(KEYS) + ["row_id", "accepted_for_noise_model", "reason"],
    )
    for group, stations in sorted(groups.items()):
        f, edges, db, probability, coverage = station_distribution(
            list(stations.values())
        )
        low, high = quantile_curves(probability, db, coverage, percentiles)
        curves = [
            {
                "frequency_hz": float(hz),
                "low_db_pa2_hz": float(lo) if n else "",
                "high_db_pa2_hz": float(hi) if n else "",
                "stations": int(n),
            }
            for hz, lo, hi, n in zip(f, low, high, coverage, strict=True)
        ]
        write_csv(
            out / f"{group}_noise_curves.csv",
            curves,
            ["frequency_hz", "low_db_pa2_hz", "high_db_pa2_hz", "stations"],
        )
        np.savez_compressed(
            out / f"{group}_noise_distribution.npz",
            frequency_edges_hz=edges,
            db_edges=db,
            station_balanced_probability=probability,
            stations_per_frequency=coverage,
        )
        fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
        ax.semilogx(f, low, label=f"{percentiles[0]:g}th percentile")
        ax.semilogx(f, high, label=f"{percentiles[1]:g}th percentile")
        ax.set(
            xlabel="Frequency [Hz]",
            ylabel=r"Pressure PSD [dB re 1 Pa$^2$/Hz]",
            title=f"{group}: descriptive approved-station noise percentiles\n"
            "Equal station weights; nominal responses; 1-dB histogram resolution",
        )
        ax.grid(which="both", alpha=0.2)
        ax.legend()
        fig.savefig(out / f"{group}_noise_curves.png", dpi=160)
        plt.close(fig)
        summary["models"].append(
            {
                "group": group,
                "stations": len(stations),
                "retained_windows": sum(s["windows"] for s in stations.values()),
                "frequency_bins_with_coverage": int((coverage > 0).sum()),
            }
        )
    inputs.verify()
    summary["protected_sha256"] = inputs.hashes
    (out / "noise_model_report.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "protected_sha256"}, indent=2
        )
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--plan", action="store_true", help="Selection audit only, no model fitting"
    )
    parser.add_argument(
        "--group-by",
        choices=["all", "network", "channel", "network_channel"],
        default="network_channel",
    )
    parser.add_argument("--percentiles", type=float, nargs=2, default=(5.0, 95.0))
    args = parser.parse_args()
    build_noise_model(
        args.review,
        args.output_dir,
        plan=args.plan,
        group_by=args.group_by,
        percentiles=args.percentiles,
    )


if __name__ == "__main__":
    main()
