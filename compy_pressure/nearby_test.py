"""Compare saved pressure spectra with nearby stations at similar water depth."""

import argparse
import csv
import hashlib
import json
import warnings
from dataclasses import asdict, dataclass
from io import BytesIO
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy.signal import PPSD

from compy_pressure.common import PROJECT_DIRECTORY
from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    OUTPUT_DIRECTORY,
    RESULT_COLUMNS,
    TEST_COLUMNS,
    load_master_dataset,
    update_nearby_test_results,
)
from compy_pressure.pressure_test import response_band, validate_pressure

KEYS = EPOCH_COLUMNS + ["inventory_file"]
STATUS_COLUMN = TEST_COLUMNS[1].removesuffix("_score") + "_status"
FINISHED = {
    "comparison_only",
    "exploratory_only",
    "no_suitable_neighbors",
    "insufficient_common_days",
    "insufficient_reference_stations",
    "no_comparable_frequencies",
    "insufficient_joint_coverage",
}
REFERENCES = [
    {
        "doi": "10.1029/2019EA000783",
        "role": "Compact-array calibration context; not a distance cutoff",
    },
    {
        "doi": "10.1002/grl.50418",
        "role": "Regional infragravity and depth-dependence context; "
        "not a distance cutoff",
    },
]


@dataclass(frozen=True)
class Criteria:
    preferred_distance_km: float = 50
    maximum_distance_km: float = 100
    preferred_depth_fraction: float = 0.05
    maximum_depth_fraction: float = 0.10
    minimum_common_days: int = 30
    minimum_daily_coverage_hours: float = 12
    minimum_reference_stations: int = 2
    response_loss_db: float = 3

    def __post_init__(self):
        scalars = [
            self.preferred_distance_km,
            self.maximum_distance_km,
            self.preferred_depth_fraction,
            self.maximum_depth_fraction,
            self.minimum_daily_coverage_hours,
            self.response_loss_db,
        ]
        if not np.isfinite(scalars).all() or min(scalars) <= 0:
            raise ValueError("Criteria must be finite and positive")
        if self.preferred_distance_km > self.maximum_distance_km:
            raise ValueError("Preferred distance must not exceed maximum distance")
        if not self.preferred_depth_fraction <= self.maximum_depth_fraction < 1:
            raise ValueError("Depth fractions must satisfy preferred <= maximum < 1")
        if self.minimum_daily_coverage_hours > 24:
            raise ValueError("Daily coverage cannot exceed 24 hours")
        if (
            not isinstance(self.minimum_common_days, int)
            or self.minimum_common_days < 1
        ):
            raise ValueError("Common-day minimum must be a positive integer")
        if (
            not isinstance(self.minimum_reference_stations, int)
            or self.minimum_reference_stations < 2
        ):
            raise ValueError("Require at least two distinct reference stations")


def row_key(row):
    return tuple(str(row[name]) for name in KEYS)


def station_key(row):
    return ".".join(row["Seed_id"].split(".")[:2])


def pressure_rows(master):
    needed = {"ppsd_status", "ppsd_units", "ppsd_npz_file"}
    if not needed.issubset(master):
        raise ValueError("Run the PPSD stage before comparing nearby stations")
    return (
        master[
            master.ppsd_status.isin(["complete", "partial"])
            & master.ppsd_units.eq("Pa^2/Hz")
            & master.ppsd_npz_file.ne("")
        ]
        .drop_duplicates(KEYS)
        .copy()
    )


def metadata(row):
    latitude, longitude, depth = [
        float(row[name]) for name in ("Lat", "Lon", "Profondeur")
    ]
    if not np.isfinite([latitude, longitude, depth]).all():
        raise ValueError("Coordinates and water-depth estimate must be finite")
    if not -90 <= latitude <= 90 or not -180 <= longitude <= 180 or depth <= 0:
        raise ValueError("Invalid coordinates or nonpositive water-depth estimate")
    return latitude, longitude, depth


def declared_sensor(row):
    """Compare inventory labels; do not infer actual hardware from channel codes."""
    label = str(row.get("sensor_type", "")).strip().casefold()
    model = str(row.get("sensor_model", "")).strip().casefold()
    return (label, model) if label else None


def known_failure(row):
    for name in RESULT_COLUMNS:
        try:
            if str(row.get(name, "")).strip() and float(row[name]) == 0:
                return True
        except (ValueError, TypeError):
            continue
    return False


def candidate(target, reference, criteria):
    """Make a geographic/metadata decision before loading any reference spectra."""
    audit = {
        "reference_seed_id": reference["Seed_id"],
        "reference_epoch_start": reference["channel_start_time"],
        "reference_epoch_end": reference["channel_end_time"],
        "reference_inventory_file": reference["inventory_file"],
        "distance_km": None,
        "azimuth_degrees": None,
        "depth_difference_m": None,
        "depth_difference_fraction": None,
        "sensor_match": "unknown",
        "tier": "excluded",
        "reason": "",
        "common_days": 0,
        "compared_frequency_bins": 0,
        "median_difference_db": None,
        "used_in_reference_group": False,
    }
    if station_key(target) == station_key(reference):
        audit["reason"] = "same_station"
        return audit
    try:
        lat, lon, depth = metadata(target)
        rlat, rlon, rdepth = metadata(reference)
    except (ValueError, TypeError):
        audit["reason"] = "missing_or_invalid_metadata"
        return audit
    distance, azimuth, _ = gps2dist_azimuth(lat, lon, rlat, rlon)
    fraction = abs(depth - rdepth) / ((depth + rdepth) / 2)
    sensor = declared_sensor(target)
    same_sensor = sensor is not None and sensor == declared_sensor(reference)
    audit.update(
        distance_km=distance / 1000,
        azimuth_degrees=azimuth,
        depth_difference_m=abs(depth - rdepth),
        depth_difference_fraction=fraction,
        sensor_match="matching_declared_labels"
        if same_sensor
        else "unknown"
        if sensor is None or declared_sensor(reference) is None
        else "different_declared_labels",
    )
    if distance / 1000 > criteria.maximum_distance_km:
        audit["reason"] = "outside_maximum_distance"
    elif fraction > criteria.maximum_depth_fraction:
        audit["reason"] = "outside_maximum_depth_difference"
    elif known_failure(reference):
        audit["reason"] = "reference_has_recorded_failed_check"
    else:
        audit["tier"] = (
            "preferred"
            if (
                same_sensor
                and distance / 1000 <= criteria.preferred_distance_km
                and fraction <= criteria.preferred_depth_fraction
            )
            else "exploratory"
        )
    return audit


def union_seconds(intervals):
    """Measure coverage without counting overlapped PPSD windows twice."""
    total = 0.0
    end = None
    for left, right in sorted(intervals):
        if end is None or left > end:
            total += right - left
        elif right > end:
            total += right - end
        end = right if end is None else max(end, right)
    return total


def daily_spectra(times, values, duration, minimum_hours):
    """Daily median dB spectra from contained windows with adequate coverage."""
    groups = {}
    for index, start in enumerate(times):
        day = UTCDateTime(start.date)
        # Midnight-spanning windows are excluded; no partial-window PSD is invented.
        if start + duration <= day + 86400:
            groups.setdefault(str(start.date), []).append((index, float(start - day)))
    days, spectra, counts, coverage = [], [], [], []
    for day, windows in sorted(groups.items()):
        seconds = union_seconds([(s, s + duration) for _, s in windows])
        if seconds + 1e-6 < minimum_hours * 3600:
            continue
        days.append(day)
        spectra.append(np.median(values[[i for i, _ in windows]], axis=0))
        counts.append(len(windows))
        coverage.append(seconds / 3600)
    return (
        np.array(days),
        np.array(spectra).reshape(len(days), values.shape[1]),
        counts,
        coverage,
    )


@dataclass
class DailyPressure:
    row: dict
    frequencies: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    days: np.ndarray
    values: np.ndarray
    band: dict
    processing: tuple
    hashes: dict
    window_counts: list
    coverage_hours: list


def load_daily(row, criteria):
    source = Path(row["ppsd_npz_file"])
    blob = source.read_bytes()
    ppsd = PPSD.load_npz(BytesIO(blob))
    values = validate_pressure(ppsd, row)
    if not np.isclose(ppsd.sampling_rate, float(row["sample_rate_hz"])):
        raise ValueError("PPSD sample rate differs from the master")
    band, _, _ = response_band(
        ppsd, source.parent / "inventory.xml", loss_db=criteria.response_loss_db
    )
    band["fft_minimum_hz"] = ppsd.sampling_rate / ppsd.nfft
    band["fft_maximum_hz"] = ppsd.sampling_rate / 2
    days, daily, counts, coverage = daily_spectra(
        ppsd.times_processed,
        values,
        ppsd.ppsd_length,
        criteria.minimum_daily_coverage_hours,
    )
    order = np.argsort(1 / ppsd.period_bin_centers)
    return DailyPressure(
        dict(row),
        (1 / ppsd.period_bin_centers)[order],
        (1 / ppsd.period_bin_right_edges)[order],
        (1 / ppsd.period_bin_left_edges)[order],
        days,
        daily[:, order],
        band,
        (
            ppsd.ppsd_length,
            ppsd.overlap,
            round(
                float(
                    np.median(
                        np.log2(
                            ppsd.period_bin_right_edges / ppsd.period_bin_left_edges
                        )
                    )
                ),
                12,
            ),
        ),
        {
            "ppsd_sha256": hashlib.sha256(blob).hexdigest(),
            "response_sha256": hashlib.sha256(
                (source.parent / "inventory.xml").read_bytes()
            ).hexdigest(),
        },
        counts,
        coverage,
    )


def align_pair(target, reference, criteria):
    """Interpolate already smoothed peer spectra in log-frequency; never extrapolate."""
    if target.processing != reference.processing:
        return None, "incompatible_processing_settings"
    days, ti, ri = np.intersect1d(target.days, reference.days, return_indices=True)
    if len(days) < criteria.minimum_common_days:
        return {"days": days}, "insufficient_common_days"
    minimum = max(
        target.band["minimum_hz"],
        reference.band["minimum_hz"],
        target.band.get("fft_minimum_hz", 0),
        reference.band.get("fft_minimum_hz", 0),
    )
    maximum = min(
        target.band["maximum_hz"],
        reference.band["maximum_hz"],
        target.band.get("fft_maximum_hz", np.inf),
        reference.band.get("fft_maximum_hz", np.inf),
    )
    tm = (target.lower >= minimum) & (target.upper <= maximum)
    rm = (reference.lower >= minimum) & (reference.upper <= maximum)
    if rm.sum() < 2:
        return {"days": days}, "no_comparable_frequencies"
    rf = reference.frequencies[rm]
    tm &= (target.frequencies >= rf.min()) & (target.frequencies <= rf.max())
    indices = np.flatnonzero(tm)
    if not len(indices):
        return {"days": days}, "no_comparable_frequencies"
    x = np.log(target.frequencies[indices])
    aligned = np.array(
        [np.interp(x, np.log(rf), row[rm]) for row in reference.values[ri]]
    )
    tv = target.values[np.ix_(ti, indices)]
    return {
        "days": days,
        "target_indices": ti,
        "frequency_indices": indices,
        "target": tv,
        "reference": aligned,
        "difference": tv - aligned,
    }, ""


def reference_group(target, pairs, criteria):
    """One epoch/location per station, equal weights and joint daily coverage."""
    unique = {}
    for pair in sorted(
        pairs,
        key=lambda p: (
            p["audit"]["tier"] != "preferred",
            -len(p["aligned"]["days"]),
            p["audit"]["distance_km"],
            p["audit"]["reference_seed_id"],
            p["audit"]["reference_epoch_start"],
        ),
    ):
        unique.setdefault(station_key(pair["data"].row), pair)
    selected = list(unique.values())
    shape = target.values.shape
    stack = np.full((len(selected), *shape), np.nan)
    for index, pair in enumerate(selected):
        aligned = pair["aligned"]
        stack[index][
            np.ix_(aligned["target_indices"], aligned["frequency_indices"])
        ] = aligned["reference"]
    counts = np.sum(np.isfinite(stack), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        peer = np.nanmedian(stack, axis=0) if selected else np.full(shape, np.nan)
    peer[counts < criteria.minimum_reference_stations] = np.nan
    residual = target.values - peer
    keep = np.sum(np.isfinite(residual), axis=0) >= criteria.minimum_common_days
    return {
        "pairs": selected,
        "peer": peer,
        "residual": residual,
        "counts": counts,
        "frequency_mask": keep,
    }


def finite_percentiles(values):
    finite = np.asarray(values)[np.isfinite(values)]
    return (
        tuple(map(float, np.percentile(finite, [5, 50, 95])))
        if len(finite)
        else (None, None, None)
    )


def finite_rms(values):
    """Root mean square of finite dB differences: offset and scatter together."""
    finite = np.asarray(values)[np.isfinite(values)]
    return float(np.sqrt(np.mean(finite**2))) if len(finite) else None


def comparison_plots(target, audits, pairs, group, status, criteria):
    figure, axes = plt.subplots(
        2, 1, figsize=(11, 8), sharex=True, layout="constrained"
    )
    if pairs:
        for index, pair in enumerate(pairs):
            aligned = pair["aligned"]
            f = target.frequencies[aligned["frequency_indices"]]
            axes[0].plot(
                f,
                np.median(aligned["reference"], axis=0),
                alpha=0.45,
                color="#888888",
                label="Pair-matched reference medians" if index == 0 else None,
            )
    keep = group["frequency_mask"]
    if keep.any():
        f = target.frequencies[keep]
        good = np.isfinite(group["residual"][:, keep])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            tv = np.nanmedian(np.where(good, target.values[:, keep], np.nan), axis=0)
            rv = np.nanmedian(group["peer"][:, keep], axis=0)
        percentiles = np.array(
            [finite_percentiles(v) for v in group["residual"][:, keep].T]
        ).T
        axes[0].plot(
            f, tv, color="#0072B2", linewidth=2, label="Target, group-matched days"
        )
        axes[0].plot(
            f, rv, color="#E69F00", linewidth=2, label="Equal-station peer median"
        )
        axes[1].fill_between(
            f,
            percentiles[0],
            percentiles[2],
            alpha=0.2,
            color="#0072B2",
            label="5th–95th daily spread",
        )
        axes[1].plot(
            f, percentiles[1], color="#0072B2", label="Median daily difference"
        )
        rms = [finite_rms(v) for v in group["residual"][:, keep].T]
        axes[1].plot(
            f,
            [np.nan if r is None else r for r in rms],
            color="black",
            linestyle="--",
            linewidth=1.5,
            label="RMS daily difference",
        )
    elif pairs:
        pair = pairs[0]["aligned"]
        f = target.frequencies[pair["frequency_indices"]]
        axes[0].plot(
            f,
            np.median(pair["target"], axis=0),
            color="#0072B2",
            label="Target, first pair's common days",
        )
        axes[1].plot(
            f,
            np.median(pair["difference"], axis=0),
            color="#888888",
            label="Single-pair context only",
        )
        axes[1].plot(
            f,
            np.sqrt(np.mean(np.asarray(pair["difference"]) ** 2, axis=0)),
            color="black",
            linestyle="--",
            linewidth=1.5,
            label="Single-pair RMS difference",
        )
    else:
        for axis in axes:
            axis.text(
                0.5,
                0.5,
                "No eligible paired spectra",
                ha="center",
                transform=axis.transAxes,
            )
    axes[0].set_ylabel(r"Pressure PSD [dB re 1 Pa$^2$/Hz]")
    axes[1].set_ylabel("Target minus reference [dB]")
    axes[1].set_xlabel("Frequency [Hz]")
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
    for axis in axes:
        axis.set_xscale("log")
        axis.grid(which="both", alpha=0.2)
        if axis.get_legend_handles_labels()[0]:
            axis.legend(fontsize=8)
    figure.suptitle(
        f"{target.row['Seed_id']}: nearby pressure comparison\n{status}; "
        "daily spread is not a confidence interval",
        fontsize=11,
    )
    figures = [figure]

    figure, axis = plt.subplots(figsize=(11, 5.8), layout="constrained")
    heat = None
    mask = keep.copy()
    days = target.days
    title = "Joint-reference daily differences"
    if mask.any():
        heat = group["residual"][:, mask]
    elif pairs:
        aligned = pairs[0]["aligned"]
        mask[:] = False
        mask[aligned["frequency_indices"]] = True
        days, heat = aligned["days"], aligned["difference"]
        title = "Single-pair daily differences: context only"
    if heat is not None and heat.size:
        x = np.array([mdates.date2num(UTCDateTime(d).datetime) for d in days])
        all_days = np.arange(x[0], x[-1] + 1)
        z = np.full((len(all_days), heat.shape[1]), np.nan)
        z[(x - x[0]).astype(int)] = heat
        f = target.frequencies[mask]
        if len(f) == 1:
            edges = f[0] * np.array([1 / np.sqrt(2), np.sqrt(2)])
        else:
            mid = np.sqrt(f[:-1] * f[1:])
            edges = np.r_[f[0] ** 2 / mid[0], mid, f[-1] ** 2 / mid[-1]]
        limit = max(1, float(np.nanpercentile(abs(heat), 98)))
        mesh = axis.pcolormesh(
            np.r_[all_days, all_days[-1] + 1],
            edges,
            np.ma.masked_invalid(z.T),
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            shading="flat",
        )
        figure.colorbar(mesh, ax=axis, label="Target minus reference [dB]")
        axis.set_yscale("log")
        axis.xaxis_date()
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    else:
        axis.text(
            0.5,
            0.5,
            "No eligible daily differences",
            ha="center",
            transform=axis.transAxes,
        )
    axis.set_xlabel("UTC day; missing or unsupported cells remain blank")
    axis.set_ylabel("Frequency [Hz]")
    axis.set_title(f"{target.row['Seed_id']}: {title}\n{status}", fontsize=11)
    figures.append(figure)

    figure, (map_axis, depth_axis) = plt.subplots(
        1, 2, figsize=(12, 6), layout="constrained"
    )
    colors = {"preferred": "#009E73", "exploratory": "#E69F00", "excluded": "#999999"}
    for tier, color in colors.items():
        rows = [
            a
            for a in audits
            if a["tier"] == tier
            and a["distance_km"] is not None
            and a["distance_km"] <= criteria.maximum_distance_km * 1.15
        ]
        if rows:
            distances = np.array([a["distance_km"] for a in rows])
            angle = np.deg2rad([a["azimuth_degrees"] for a in rows])
            map_axis.scatter(
                distances * np.sin(angle),
                distances * np.cos(angle),
                color=color,
                label=tier.capitalize(),
                edgecolor="white",
                s=55,
            )
            depth_axis.scatter(
                distances,
                [a["depth_difference_fraction"] * 100 for a in rows],
                color=color,
                s=55,
                label=tier.capitalize(),
            )
    map_axis.scatter(0, 0, marker="*", color="#0072B2", s=140, label="Target")
    for radius, style in [
        (criteria.preferred_distance_km, "--"),
        (criteria.maximum_distance_km, ":"),
    ]:
        map_axis.add_patch(
            plt.Circle(
                (0, 0), radius, fill=False, color="black", linestyle=style, alpha=0.5
            )
        )
    extent = criteria.maximum_distance_km * 1.15
    map_axis.set(
        xlim=(-extent, extent),
        ylim=(-extent, extent),
        xlabel="East from target [km]",
        ylabel="North from target [km]",
        title="Geodesic distance and azimuth",
    )
    map_axis.set_aspect("equal")
    map_axis.legend(fontsize=8)
    depth_axis.axvline(
        criteria.preferred_distance_km, color="black", linestyle="--", alpha=0.5
    )
    depth_axis.axhline(
        criteria.preferred_depth_fraction * 100,
        color="black",
        linestyle="--",
        alpha=0.5,
    )
    depth_axis.axvline(
        criteria.maximum_distance_km, color="black", linestyle=":", alpha=0.5
    )
    depth_axis.axhline(
        criteria.maximum_depth_fraction * 100, color="black", linestyle=":", alpha=0.5
    )
    depth_axis.set(
        xlim=(0, extent),
        ylim=(0, criteria.maximum_depth_fraction * 120),
        xlabel="Distance [km]",
        ylabel="Depth difference / pair mean depth [%]",
        title="Preferred and exploratory selection",
    )
    for axis in (map_axis, depth_axis):
        axis.grid(alpha=0.2)
    figure.suptitle(
        f"{target.row['Seed_id']}: neighbour selection\n"
        "Display is limited to the local range; "
        "all exclusions are saved in candidates.csv",
        fontsize=11,
    )
    figures.append(figure)
    return figures


def write_csv(path, rows, fieldnames=None):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(rows[0]) if rows else fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def ensure_master_idle(results_directory):
    """An unfinished batch may still write the same CSV/Excel pair."""
    path = Path(results_directory) / "_workflow/ppsd.json"
    if path.exists():
        report = json.loads(path.read_text())
        keys = report.get("selected_epoch_keys", [])
        if keys and sum(report.get("status_counts", {}).values()) != len(keys):
            raise ValueError(
                "PPSD batch is unfinished: finish/resume it before updating "
                "nearby diagnostics in the same master. "
                "--no-update-master permits a read-only pilot."
            )


def pool_fingerprint(rows):
    """Include every available reference, so new peers invalidate unresolved caches."""
    digest = hashlib.sha256()
    for _, row in rows.sort_values(KEYS).iterrows():
        names = (
            KEYS
            + [
                "Lat",
                "Lon",
                "Profondeur",
                "sensor_type",
                "sensor_model",
                "ppsd_npz_file",
                "ppsd_status",
                "ppsd_units",
                "sample_rate_hz",
            ]
            + RESULT_COLUMNS
        )
        digest.update(json.dumps([row.get(n, "") for n in names]).encode())
        for path in [
            Path(row["ppsd_npz_file"]),
            Path(row["ppsd_npz_file"]).parent / "inventory.xml",
        ]:
            digest.update(
                hashlib.sha256(path.read_bytes()).digest()
                if path.is_file()
                else b"missing"
            )
    return digest.hexdigest()


def run_nearby_test(
    master_dir=OUTPUT_DIRECTORY,
    *,
    network=None,
    station=None,
    selected_rows=None,
    reference_rows=None,
    criteria=None,
    dataset_cache=None,
    update_master=True,
    return_updates=False,
    show=True,
    results_directory=None,
):
    """Save pair and ensemble diagnostics without assigning calibration scores."""
    criteria = criteria or Criteria()
    master = (
        load_master_dataset(master_dir)
        if selected_rows is None or reference_rows is None
        else None
    )
    targets = pressure_rows(master) if selected_rows is None else selected_rows.copy()
    references = (
        pressure_rows(master)
        if reference_rows is None
        else pressure_rows(reference_rows)
    )
    if network:
        targets = targets[targets.Seed_id.str.split(".").str[0].eq(network)]
    if station:
        targets = targets[targets.Seed_id.str.split(".").str[1].eq(station)]
    targets = targets.drop_duplicates(KEYS)
    if targets.empty:
        raise ValueError("No saved pressure PPSDs for the selected target scope")
    if update_master:
        if (
            results_directory is None
            and Path(master_dir).resolve() == Path(OUTPUT_DIRECTORY).resolve()
        ):
            results_directory = PROJECT_DIRECTORY / "results/NETWORKS/TEST"
        if results_directory is not None:
            ensure_master_idle(results_directory)
    cache = dataset_cache if dataset_cache is not None else {}

    def load(row):
        key = (
            row_key(row),
            criteria.minimum_daily_coverage_hours,
            criteria.response_loss_db,
        )
        if key not in cache:
            try:
                cache[key] = load_daily(row, criteria)
            except Exception as exc:
                cache[key] = exc
        if isinstance(cache[key], Exception):
            raise ValueError(str(cache[key]))
        return cache[key]

    summaries, updates, figures = [], [], []
    for _, target_row in targets.iterrows():
        metadata(target_row)
        target = load(target_row)
        audits, pairs = [], []
        for _, reference_row in references.iterrows():
            audit = candidate(target_row, reference_row, criteria)
            audits.append(audit)
            if audit["reason"]:
                continue
            try:
                reference = load(reference_row)
                aligned, reason = align_pair(target, reference, criteria)
                audit["common_days"] = (
                    len(aligned["days"]) if aligned is not None else 0
                )
                if reason:
                    audit["reason"] = reason
                    audit["tier"] = "excluded"
                    continue
                audit["compared_frequency_bins"] = len(aligned["frequency_indices"])
                audit["median_difference_db"] = float(
                    np.median(np.median(aligned["difference"], axis=0))
                )
                pairs.append({"audit": audit, "data": reference, "aligned": aligned})
            except Exception as exc:
                audit["reason"] = f"unreadable_reference: {exc}"
                audit["tier"] = "excluded"
        preferred = [p for p in pairs if p["audit"]["tier"] == "preferred"]
        group = reference_group(target, preferred, criteria)
        status = "comparison_only"
        if not group["frequency_mask"].any():
            group = reference_group(target, pairs, criteria)
            status = "exploratory_only"
        keep = group["frequency_mask"]
        if not keep.any():
            if not pairs:
                reasons = {a["reason"] for a in audits}
                status = (
                    "insufficient_common_days"
                    if "insufficient_common_days" in reasons
                    else (
                        "no_comparable_frequencies"
                        if "no_comparable_frequencies" in reasons
                        else "no_suitable_neighbors"
                    )
                )
            elif len(group["pairs"]) < criteria.minimum_reference_stations:
                status = "insufficient_reference_stations"
            else:
                status = "insufficient_joint_coverage"
        else:
            for pair in group["pairs"]:
                pair["audit"]["used_in_reference_group"] = True
        frequency_rows = []
        for index in np.flatnonzero(keep):
            valid = np.isfinite(group["residual"][:, index])
            p05, median, p95 = finite_percentiles(group["residual"][:, index])
            frequency_rows.append(
                {
                    "frequency_hz": float(target.frequencies[index]),
                    "paired_days": int(valid.sum()),
                    "target_median_db_pa2_hz": float(
                        np.median(target.values[valid, index])
                    ),
                    "peer_median_db_pa2_hz": float(
                        np.median(group["peer"][valid, index])
                    ),
                    "difference_p05_db": p05,
                    "difference_median_db": median,
                    "difference_p95_db": p95,
                    "difference_rms_db": finite_rms(group["residual"][:, index]),
                }
            )
        daily_rows = []
        for index, day in enumerate(target.days):
            values = group["residual"][index, keep]
            p05, median, p95 = finite_percentiles(values)
            count = int(np.isfinite(values).sum())
            daily_rows.append(
                {
                    "utc_day": str(day),
                    "frequency_bins": count,
                    "median_difference_db": median,
                    "rms_difference_db": finite_rms(values),
                    "frequency_p05_db": p05,
                    "frequency_p95_db": p95,
                    "consistent_frequency_support": bool(
                        keep.any() and count == int(keep.sum())
                    ),
                    "target_window_count": target.window_counts[index],
                    "target_coverage_hours": target.coverage_hours[index],
                }
            )
        pair_rows = []
        for pair in pairs:
            aligned = pair["aligned"]
            for index, frequency_index in enumerate(aligned["frequency_indices"]):
                p05, median, p95 = finite_percentiles(aligned["difference"][:, index])
                pair_rows.append(
                    {
                        "reference_seed_id": pair["audit"]["reference_seed_id"],
                        "reference_epoch_start": pair["audit"]["reference_epoch_start"],
                        "reference_inventory_file": pair["audit"][
                            "reference_inventory_file"
                        ],
                        "tier": pair["audit"]["tier"],
                        "frequency_hz": float(target.frequencies[frequency_index]),
                        "paired_days": len(aligned["days"]),
                        "difference_p05_db": p05,
                        "difference_median_db": median,
                        "difference_p95_db": p95,
                        "difference_rms_db": finite_rms(
                            aligned["difference"][:, index]
                        ),
                    }
                )
        destination = Path(target_row.ppsd_npz_file).parent / "nearby"
        destination.mkdir(exist_ok=True)
        csv_paths = {
            name: destination / f"{name}.csv"
            for name in (
                "candidates",
                "frequency_comparison",
                "daily_comparison",
                "pair_comparisons",
            )
        }
        write_csv(
            csv_paths["candidates"],
            audits,
            list(candidate(target_row, target_row, criteria)),
        )
        write_csv(
            csv_paths["frequency_comparison"],
            frequency_rows,
            [
                "frequency_hz",
                "paired_days",
                "target_median_db_pa2_hz",
                "peer_median_db_pa2_hz",
                "difference_p05_db",
                "difference_median_db",
                "difference_p95_db",
                "difference_rms_db",
            ],
        )
        write_csv(
            csv_paths["daily_comparison"],
            daily_rows,
            [
                "utc_day",
                "frequency_bins",
                "median_difference_db",
                "rms_difference_db",
                "frequency_p05_db",
                "frequency_p95_db",
                "consistent_frequency_support",
                "target_window_count",
                "target_coverage_hours",
            ],
        )
        write_csv(
            csv_paths["pair_comparisons"],
            pair_rows,
            [
                "reference_seed_id",
                "reference_epoch_start",
                "reference_inventory_file",
                "tier",
                "frequency_hz",
                "paired_days",
                "difference_p05_db",
                "difference_median_db",
                "difference_p95_db",
                "difference_rms_db",
            ],
        )
        names = ("spectra_comparison", "daily_residuals", "neighbor_selection")
        plot_paths = {name: destination / f"{name}.png" for name in names}
        plots = comparison_plots(target, audits, pairs, group, status, criteria)
        for name, figure in zip(names, plots, strict=True):
            figure.savefig(plot_paths[name], dpi=160)
            if not show:
                plt.close(figure)
        figures.extend(plots)
        json_path = destination / "comparison.json"
        summary = {
            "seed_id": target_row.Seed_id,
            "status": status,
            "score": None,
            "criteria": asdict(criteria),
            "criteria_basis": "Proposed screening settings; "
            "not published acceptance thresholds",
            "references": REFERENCES,
            "depth_basis": "inventory_elevation_estimate",
            "depth_independently_verified": False,
            "response_band": target.band,
            "qualifying_target_days": len(target.days),
            "paired_reference_epochs": len(pairs),
            "reference_stations": len(group["pairs"]) if keep.any() else 0,
            "preferred_reference_stations": len(
                {station_key(p["data"].row) for p in preferred}
            ),
            "compared_frequency_bins": int(keep.sum()),
            "joint_common_days": int(
                np.isfinite(group["residual"][:, keep]).any(axis=1).sum()
            )
            if keep.any()
            else 0,
            "median_difference_db": float(
                np.median([r["difference_median_db"] for r in frequency_rows])
            )
            if frequency_rows
            else None,
            "target_hashes": target.hashes,
            "processing": {
                "ppsd_length_seconds": target.processing[0],
                "overlap_fraction": target.processing[1],
                "smoothing_width_octaves": target.processing[2],
            },
            "paired_references": [
                {
                    "epoch": {n: p["data"].row[n] for n in KEYS},
                    "hashes": p["data"].hashes,
                    "response_band": p["data"].band,
                    "used_in_reference_group": p["audit"]["used_in_reference_group"],
                }
                for p in pairs
            ],
            "files": {
                **{n: str(p) for n, p in csv_paths.items()},
                **{n: str(p) for n, p in plot_paths.items()},
            },
            "definitions": {
                "depth_difference": "abs(H1-H2) / mean(H1,H2)",
                "spectra": "Daily median of saved dB spectra; "
                "no response correction reapplied",
                "frequency_alignment": "Linear interpolation of smoothed dB spectra "
                "in log frequency, with full-bin band masks and no extrapolation",
                "reference_group": "One epoch/location per distinct network.station, "
                "equal weights; preferred group first, explicitly exploratory group "
                "if needed",
                "residual": "Target daily median minus the same-day peer median; "
                "daily differences summarized per frequency",
                "spread": "Empirical daily spread, not a confidence interval",
            },
            "limitations": [
                "Nearby stations are not independently calibrated references.",
                "Matching declared sensor labels does not verify hardware type.",
                "Depth, coastal setting and seafloor environment require review.",
                "Overlapping windows, smoothed bins and nearby stations "
                "are correlated.",
                "A level difference can be environmental; "
                "no correction or pass/fail is inferred.",
            ],
        }
        json_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        update = {n: target_row[n] for n in KEYS} | {
            "nearby_status": status,
            "nearby_result_file": str(json_path),
            "nearby_reference_stations": summary["reference_stations"],
            "nearby_compared_frequency_bins": int(keep.sum()),
            "nearby_joint_common_days": summary["joint_common_days"],
            "nearby_median_difference_db": summary["median_difference_db"],
            **{f"nearby_{n}_file": str(p) for n, p in (csv_paths | plot_paths).items()},
        }
        if not target_row.get(TEST_COLUMNS[1], ""):
            update[STATUS_COLUMN] = status
        summaries.append(summary)
        updates.append(update)
        print(
            f"{target_row.Seed_id}: nearby={status}; "
            f"{summary['reference_stations']} references; "
            f"{keep.sum()} frequency bins; saved {destination}",
            flush=True,
        )
    if update_master:
        update_nearby_test_results(master_dir, updates)
    if show:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)
    return (summaries, updates) if return_updates else summaries


def add_criteria_arguments(parser):
    for name, default in [
        ("preferred-distance-km", 50),
        ("maximum-distance-km", 100),
        ("preferred-depth-fraction", 0.05),
        ("maximum-depth-fraction", 0.10),
        ("minimum-daily-coverage-hours", 12),
        ("response-loss-db", 3),
    ]:
        parser.add_argument(f"--{name}", type=float, default=default)
    parser.add_argument("--minimum-common-days", type=int, default=30)
    parser.add_argument("--minimum-reference-stations", type=int, default=2)


def criteria_from_args(args):
    return Criteria(
        **{name: getattr(args, name) for name in Criteria.__dataclass_fields__}
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-dir", type=Path, default=OUTPUT_DIRECTORY)
    parser.add_argument("--network")
    parser.add_argument("--station")
    parser.add_argument("--no-update-master", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    add_criteria_arguments(parser)
    args = parser.parse_args()
    run_nearby_test(
        args.master_dir,
        network=args.network,
        station=args.station,
        criteria=criteria_from_args(args),
        update_master=not args.no_update_master,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()
