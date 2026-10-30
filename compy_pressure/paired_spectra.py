"""Synchronized pressure/upward-displacement spectra, coherence and admittance."""

import argparse
import hashlib
import json
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime, read_inventory
from scipy.signal.windows import hann

from compy_pressure.archive import discover_waveforms, read_window
from compy_pressure.common import complete_samples, write_csv
from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    OUTPUT_DIRECTORY,
    load_master_dataset,
    update_paired_results,
)
from compy_pressure.ppsd import (
    atomic_json,
    is_pressure_unit,
    response_handling,
    select_epoch,
)

KEYS = EPOCH_COLUMNS + ["inventory_file"]
FINISHED = {
    "comparison_only",
    "insufficient_common_coverage",
    "no_complete_pairs",
    "no_valid_frequencies",
    "no_vertical_pair",
    "ambiguous_pair",
}


class PairingError(ValueError):
    """Unresolved metadata pairing, distinct from a processing exception."""

    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


def select_pair(pressure_row, master, vertical_id=None, *, allow_unequal_rates=False):
    """Use response units and orientation, and refuse ambiguous overlapping epochs."""
    _, pressure = select_epoch(pressure_row)
    if response_handling(pressure)[1] != "Pa^2/Hz":
        raise ValueError("Selected channel is not pressure")
    net, sta, loc, _ = pressure_row["Seed_id"].split(".")
    candidates = master[master.Seed_id.str.startswith(f"{net}.{sta}.")]
    if vertical_id:
        candidates = candidates[candidates.Seed_id.eq(vertical_id)]
    else:
        candidates = candidates[candidates.Seed_id.str.split(".").str[2].eq(loc)]
    matches = []
    for _, row in candidates.drop_duplicates(KEYS).iterrows():
        units = str(row.get("input_units", "")).upper().replace(" ", "")
        if units not in {"M", "M/S", "M/S**2", "M/S^2", "M/S/S"}:
            continue
        _, channel = select_epoch(row)
        if response_handling(channel)[1] != "m^2/s^4/Hz":
            continue
        if channel.dip is None or not np.isclose(abs(channel.dip), 90, atol=0.1):
            continue
        starts = [c.start_date for c in (pressure, channel) if c.start_date]
        ends = [c.end_date for c in (pressure, channel) if c.end_date]
        left, right = max(starts, default=None), min(ends, default=None)
        if left is None or right is None or left < right:
            matches.append((row, channel))
    if len(matches) != 1:
        raise PairingError(
            "no_vertical_pair" if not matches else "ambiguous_pair",
            f"Expected one overlapping vertical epoch, found {len(matches)}; "
            "select a vertical SEED ID or resolve epoch ambiguity explicitly",
        )
    row, vertical = matches[0]
    if pressure.sample_rate != vertical.sample_rate and not allow_unequal_rates:
        raise ValueError("Unequal sample rates: no implicit resampling is permitted")
    return pressure, row, vertical


def spectral_terms(p, z, sample_rate, pressure_response, vertical_response, sign=1):
    """One-sided Hann densities after complex response correction exactly once."""
    if len(p) != len(z) or len(p) < 4:
        raise ValueError("Paired arrays must have the same length, at least four")
    if not np.isfinite(p).all() or not np.isfinite(z).all():
        raise ValueError("Nonfinite paired samples")
    taper = hann(len(p), sym=False)
    pf = np.fft.rfft((p - np.mean(p)) * taper)
    zf = np.fft.rfft((z - np.mean(z)) * taper)
    if (
        np.shape(pressure_response) != pf.shape
        or np.shape(vertical_response) != zf.shape
    ):
        raise ValueError("Responses must match the full rFFT grid")
    rp, rz = np.asarray(pressure_response), np.asarray(vertical_response)
    valid = np.isfinite(rp) & np.isfinite(rz) & (abs(rp) > 0) & (abs(rz) > 0)
    valid[0] = False  # DISP response is singular at DC; no DC output is interpreted.
    scale = np.full(pf.size, 2 / (sample_rate * np.sum(taper**2)))
    scale[0] /= 2
    if len(p) % 2 == 0:
        scale[-1] /= 2
    pc = np.divide(pf, rp, out=np.zeros_like(pf), where=valid)
    zc = sign * np.divide(zf, rz, out=np.zeros_like(zf), where=valid)
    return (
        abs(pc) ** 2 * scale,
        abs(zc) ** 2 * scale,
        zc * pc.conj() * scale,
        valid,
    )


def averaged_products(pp, zz, zp):
    """Average complex cross-power first; no single-window coherence result."""
    if len(pp) < 2 or np.shape(pp) != np.shape(zz) or np.shape(pp) != np.shape(zp):
        raise ValueError("At least two matching segments are required")
    pp, zz, zp = np.mean(pp, axis=0), np.mean(zz, axis=0), np.mean(zp, axis=0)
    valid = np.isfinite(pp) & np.isfinite(zz) & np.isfinite(zp) & (pp > 0) & (zz > 0)
    h = np.divide(zp, pp, out=np.zeros_like(zp), where=valid)
    coherence = np.divide(abs(zp) ** 2, pp * zz, out=np.zeros_like(pp), where=valid)
    return pp, zz, zp, h, np.clip(coherence, 0, 1), valid


def source_metadata(records, start, end):
    return [
        {"path": str(p), "size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
        for p in sorted({r.path for r in records if r.start < end and r.end >= start})
    ]


def response_snapshot_path(row):
    explicit = row.get("rayleigh_response_snapshot_file", "")
    return (
        Path(explicit)
        if explicit
        else Path(row["ppsd_npz_file"]).parent / "inventory.xml"
    )


def has_saved_response(row):
    """Return whether an exact saved response snapshot is available."""
    try:
        return response_snapshot_path(row).is_file()
    except (KeyError, TypeError, ValueError):
        return False


def checked_response(row, channel, frequency, output):
    """Compare a live response to its exact saved snapshot on this grid."""
    response = channel.response.get_evalresp_response_for_frequencies(
        frequency, output=output
    )
    snapshot = response_snapshot_path(row)
    inv = read_inventory(str(snapshot))
    _, _, loc, code = row["Seed_id"].split(".")
    saved = [
        c
        for n in inv
        for s in n
        for c in s
        if c.code == code
        and c.location_code == loc
        and str(c.start_date or "") == row["channel_start_time"]
        and str(c.end_date or "") == row["channel_end_time"]
    ]
    if len(saved) != 1:
        raise ValueError("Saved response snapshot does not match the exact epoch")
    if saved[0].dip != channel.dip:
        raise ValueError("Live orientation differs from saved PPSD snapshot")
    old = saved[0].response.get_evalresp_response_for_frequencies(
        frequency, output=output
    )
    if not np.allclose(response, old, rtol=1e-10, atol=0, equal_nan=True):
        raise ValueError("Live response differs from saved PPSD snapshot")
    return response


def plot_products(f, pp, zz, h, coherence, valid, title):
    fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    axes[0].loglog(f[valid], pp[valid], label="Pressure (Pa²/Hz)")
    twin = axes[0].twinx()
    twin.loglog(f[valid], zz[valid], color="tab:orange", label="Upward D (m²/Hz)")
    axes[0].set_ylabel("Pa²/Hz")
    twin.set_ylabel("m²/Hz", color="tab:orange")
    axes[1].semilogx(f[valid], coherence[valid])
    axes[1].set_ylabel("Squared coherence")
    axes[1].set_ylim(0, 1.02)
    axes[2].loglog(f[valid], abs(h[valid]))
    axes[2].set_ylabel("|D/P| (m/Pa)")
    axes[3].semilogx(f[valid], np.angle(h[valid], deg=True))
    axes[3].set_ylabel("Phase (degrees)")
    axes[3].set_ylim(-180, 180)
    axes[3].set_xlabel("Frequency (Hz)")
    for axis in axes:
        axis.axvspan(0.1, 0.2, alpha=0.12, color="green")
        axis.grid(which="both", alpha=0.2)
    fig.suptitle(
        title + "\nUncorrected for tilt; shaded microseism band; diagnostic only"
    )
    fig.tight_layout()
    return fig


def run_paired_spectra(
    master_dir=OUTPUT_DIRECTORY,
    *,
    data_dir=None,
    network="YV",
    station="RR52",
    start=None,
    days=7,
    segment_seconds=3600,
    frequency_min=0.01,
    frequency_max=0.5,
    vertical_id=None,
    pressure_id=None,
    output_directory=None,
    event_resample_rate=None,
    native_rate_tolerance_ppm=0,
    update_master=False,
    show=False,
    selected_rows=None,
    reference_rows=None,
    records=None,
    return_updates=False,
):
    """Bounded raw reads into separate run folders; scores are never assigned."""
    master = (
        load_master_dataset(master_dir) if reference_rows is None else reference_rows
    )
    scope = master if selected_rows is None else selected_rows
    rows = scope[scope.Seed_id.str.startswith(f"{network}.{station}.")]
    rows = rows[rows.input_units.map(is_pressure_unit)]
    if pressure_id:
        rows = rows[rows.Seed_id.eq(pressure_id)]
    rows = rows.drop_duplicates(KEYS)
    if len(rows) != 1:
        raise ValueError("Select exactly one pressure epoch")
    row = rows.iloc[0]
    if not has_saved_response(row):
        raise ValueError("Pressure epoch needs an exact saved response snapshot")
    if update_master and not row.get("ppsd_npz_file", ""):
        raise ValueError("Master updates require a completed PPSD-backed epoch")
    if (
        not isinstance(days, int)
        or days < 1
        or not np.isfinite(segment_seconds)
        or segment_seconds <= 0
    ):
        raise ValueError("Days and segment duration must be positive")
    if not 0 < frequency_min < frequency_max:
        raise ValueError("Require 0 < minimum frequency < maximum frequency")
    output = (
        Path(output_directory)
        if output_directory is not None
        else response_snapshot_path(row).parent
        / "paired_spectra"
        / f"run_{uuid4().hex[:12]}"
    )
    output.mkdir(parents=True, exist_ok=False)
    summary = {
        "pressure_epoch": {n: str(row[n]) for n in KEYS},
        "status": "no_complete_pairs",
        "files": {},
        "segment_seconds": segment_seconds,
        "days": days,
        "cross_convention": "S_dp = mean(D * conjugate(P))",
        "transfer_direction": "pressure_to_upward_displacement",
        "transfer_units": "m/Pa",
        "pressure_psd_units": "Pa^2/Hz",
        "vertical_psd_units": "m^2/Hz",
        "code_sha256": {
            name: hashlib.sha256(
                Path(__file__).with_name(name).read_bytes()
            ).hexdigest()
            for name in (
                "paired_spectra.py",
                "archive.py",
                "common.py",
                "ppsd.py",
                "rayleigh_calibration.py",
            )
        },
        "numpy_version": np.__version__,
        "taper": "periodic Hann",
        "detrend": "segment mean",
        "overlap": 0,
        "frequency_min_hz": frequency_min,
        "frequency_max_hz": frequency_max,
        "event_resampling": None,
        "limitations": [
            "No tilt correction or event selection",
            "Sensor bands and depth unverified",
            "Segments may be correlated; no confidence interval or calibration score",
            "Timing/polarity follow metadata; no independent clock verification",
        ],
    }
    updates = []

    def finish(status, reason=None):
        summary["status"] = status
        if reason:
            summary["reason"] = reason
        path = output / "comparison.json"
        summary["files"]["result"] = str(path)
        atomic_json(path, summary)
        patch = {n: row[n] for n in KEYS} | {
            "paired_status": status,
            "paired_result_file": str(path),
            "paired_segments": summary.get("retained_segments", 0),
            "paired_vertical_seed_id": summary.get("vertical_seed_id", ""),
            "paired_plot_file": summary["files"].get("plot", ""),
            "paired_npz_file": summary["files"].get("spectra", ""),
            "paired_start_utc": summary.get("start_utc", ""),
            "paired_end_utc": summary.get("end_utc", ""),
            "paired_microseism_median_coherence": summary.get(
                "microseism_median_coherence"
            ),
            "paired_error": None,
        }
        updates.append(patch)
        if update_master:
            assert_workflows_finished(Path(row.ppsd_npz_file).parents[3])
            update_paired_results(master_dir, updates)
        return (summary, updates) if return_updates else summary

    try:
        pressure, vrow, vertical = select_pair(
            row,
            master,
            vertical_id,
            allow_unequal_rates=event_resample_rate is not None,
        )
    except PairingError as exc:
        return finish(exc.status, str(exc))
    if not has_saved_response(vrow):
        raise ValueError("Vertical epoch needs an exact saved response snapshot")
    summary["vertical_epoch"] = {n: str(vrow[n]) for n in KEYS}
    summary["vertical_seed_id"] = vrow.Seed_id
    native_rates = [float(pressure.sample_rate), float(vertical.sample_rate)]
    rate = float(event_resample_rate or pressure.sample_rate)
    if event_resample_rate is not None:
        if not np.isfinite(event_resample_rate) or not 0 < event_resample_rate <= min(
            native_rates
        ):
            raise ValueError("Event resample rate must be positive and not upsample")
        if (
            not np.isfinite(native_rate_tolerance_ppm)
            or not 0 <= native_rate_tolerance_ppm <= 0.1
        ):
            raise ValueError("Native-rate tolerance must be between 0 and 0.1 ppm")
        summary["event_resampling"] = {
            "method": "polyphase Kaiser beta=12; UTC Lanczos a=12",
            "target_rate_hz": rate,
            "native_rates_hz": native_rates,
            "padding_seconds": 32 / rate,
            "native_rate_tolerance_ppm": native_rate_tolerance_ppm,
            "response_correction": "native complex response after resampling counts",
            "clock_correction_applied": False,
            "gaps_filled": False,
        }
    if frequency_max >= rate / 2:
        raise ValueError("Maximum frequency must be below Nyquist")
    npts = round(segment_seconds * rate)
    if npts < 4 or not np.isclose(npts / rate, segment_seconds, rtol=0, atol=1e-8):
        raise ValueError("Segment duration must contain an integer number of samples")
    if records is None:
        if data_dir is None:
            data_dir = Path.home() / "Desktop/Data_Bruit_FM/Data"
        records = discover_waveforms(data_dir, {row.Seed_id, vrow.Seed_id})
    pr = [r for r in records if r.seed_id == row.Seed_id]
    vr = [r for r in records if r.seed_id == vrow.Seed_id]
    if not pr or not vr:
        return finish("insufficient_common_coverage", "Raw paired channels unavailable")
    left = max(
        min(r.start for r in pr),
        min(r.start for r in vr),
        *[c.start_date for c in (pressure, vertical) if c.start_date],
    )
    right = min(
        max(r.end + 1 / r.sample_rate for r in pr),
        max(r.end + 1 / r.sample_rate for r in vr),
        *[c.end_date for c in (pressure, vertical) if c.end_date],
    )
    duration = days * 86400
    if right - left < duration:
        return finish(
            "insufficient_common_coverage",
            "Requested days exceed shared epoch coverage",
        )
    if start is None:
        middle = left + (right - left - duration) / 2
        start = UTCDateTime(middle.date)
        if start < left or start + duration > right:
            start = middle
    start = UTCDateTime(start)
    end = start + duration
    if start < left or end > right:
        return finish(
            "insufficient_common_coverage",
            "Requested interval outside shared response/data coverage",
        )
    summary.update(
        start_utc=str(start),
        end_utc=str(end),
        sample_rate_hz=rate,
        vertical_dip_deg=vertical.dip,
        vertical_upward_sign=-1 if vertical.dip > 0 else 1,
    )
    f = np.fft.rfftfreq(npts, 1 / rate)
    # Evaluate only positive frequencies; DISP response at DC is undefined.
    rp, rz = np.zeros(f.size, complex), np.zeros(f.size, complex)
    band = (f >= frequency_min) & (f <= frequency_max)
    rp[band] = checked_response(row, pressure, f[band], "DEF")
    rz[band] = checked_response(vrow, vertical, f[band], "DISP")
    if not band.any():
        raise ValueError("No FFT bins in selected frequency range")
    padding = 32 / rate if event_resample_rate is not None else 0
    metadata = source_metadata(
        pr + vr,
        start - padding,
        end + padding + 2 / min(native_rates),
    )
    summary["source_files"] = metadata
    summary["response_snapshots"] = {
        str(p): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in {response_snapshot_path(row), response_snapshot_path(vrow)}
    }
    summary["input_units"] = {
        "pressure": pressure.response.instrument_sensitivity.input_units,
        "vertical": vertical.response.instrument_sensitivity.input_units,
    }
    summary["inventory_sha256"] = {
        str(p): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in {Path(row.inventory_file), Path(vrow.inventory_file)}
    }
    powers_p, powers_z, cross, decisions = [], [], [], []
    count = int(np.floor(duration / segment_seconds))
    for i in range(count):
        ws = start + i * segment_seconds
        decision = "missing_gapped_nonfinite_or_flat"
        arrays, actuals = [], []
        if event_resample_rate is not None:
            # Imported lazily because rayleigh_calibration uses the pairing helpers.
            from compy_pressure.rayleigh_calibration import load_aligned

            try:
                loaded, actual, _ = load_aligned(
                    records,
                    [row, vrow],
                    [pressure, vertical],
                    ws,
                    segment_seconds,
                    event_resample_rate,
                    native_rate_tolerance_ppm,
                )
                arrays = list(loaded)
                actuals = [actual, actual]
            except ValueError as exc:
                decision = f"resampling_refused: {exc}"
                arrays = [None, None]
                actuals = [None, None]
        else:
            for rr in (pr, vr):
                stream = read_window(rr, ws, ws + segment_seconds + 2 / rate)
                array, actual = complete_samples(stream, ws, rate, npts)
                arrays.append(array)
                actuals.append(actual)
        if all(a is not None for a in arrays):
            if abs(actuals[0] - actuals[1]) > 1e-6:
                decision = "unaligned_sample_grids"
            elif (
                actuals[0] + segment_seconds > min(end, right) + 1e-6
                or actuals[1] + segment_seconds > min(end, right) + 1e-6
            ):
                decision = "aligned_window_outside_interval_or_epoch"
            else:
                pp, zz, zp, mask = spectral_terms(
                    *arrays, rate, rp, rz, summary["vertical_upward_sign"]
                )
                powers_p.append(pp[band])
                powers_z.append(zz[band])
                cross.append(zp[band])
                decision = "retained"
        decisions.append(
            {
                "requested_start": str(ws),
                "decision": decision,
                "pressure_actual_start": str(actuals[0] or ""),
                "vertical_actual_start": str(actuals[1] or ""),
            }
        )
        if (i + 1) % 24 == 0:
            print(f"{row.Seed_id}: paired segments {i + 1}/{count}", flush=True)
    for item in metadata:
        stat = Path(item["path"]).stat()
        if (stat.st_size, stat.st_mtime_ns) != (item["size"], item["mtime_ns"]):
            raise ValueError("Source waveform changed during paired processing")
    for path, digest in (
        summary["inventory_sha256"] | summary["response_snapshots"]
    ).items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError("Source inventory changed during paired processing")
    audit = output / "windows.csv"
    write_csv(
        audit,
        decisions,
        fieldnames=[
            "requested_start",
            "decision",
            "pressure_actual_start",
            "vertical_actual_start",
        ],
    )
    summary["files"]["windows"] = str(audit)
    summary.update(requested_segments=count, retained_segments=len(powers_p))
    if len(powers_p) < 2:
        return finish(
            "no_complete_pairs", "Fewer than two synchronized, nonflat segments"
        )
    pp, zz, zp, h, coherence, valid = averaged_products(powers_p, powers_z, cross)
    valid &= mask[band]
    f = f[band]
    if not valid.any():
        return finish(
            "no_valid_frequencies", "No finite positive paired spectral support"
        )
    npz = output / "spectra.npz"
    np.savez_compressed(
        npz,
        frequency_hz=f,
        pressure_psd=pp,
        displacement_psd=zz,
        cross_dp=zp,
        transfer_dp=h,
        coherence_squared=coherence,
        valid=valid,
        pressure_window_psd=powers_p,
        displacement_window_psd=powers_z,
        cross_window_dp=cross,
        pressure_response=rp[band],
        displacement_response=rz[band],
        retained_window_starts=[
            d["pressure_actual_start"] for d in decisions if d["decision"] == "retained"
        ],
    )
    table = output / "spectra.csv"
    write_csv(
        table,
        [
            {
                "frequency_hz": float(f[i]),
                "pressure_pa2_hz": float(pp[i]),
                "displacement_m2_hz": float(zz[i]),
                "cross_dp_real": float(zp[i].real),
                "cross_dp_imag": float(zp[i].imag),
                "coherence_squared": float(coherence[i]),
                "transfer_real_m_pa": float(h[i].real),
                "transfer_imag_m_pa": float(h[i].imag),
                "transfer_amplitude_m_pa": float(abs(h[i])),
                "transfer_phase_deg": float(np.angle(h[i], deg=True)),
                "valid": bool(valid[i]),
            }
            for i in range(len(f))
        ],
    )
    summary["files"]["frequency_table"] = str(table)
    micro = valid & (f >= 0.1) & (f <= 0.2)
    summary["microseism_bins"] = int(micro.sum())
    summary["microseism_median_coherence"] = (
        float(np.median(coherence[micro])) if micro.any() else None
    )
    summary["valid_bins"] = int(valid.sum())
    figure = plot_products(
        f,
        pp,
        zz,
        h,
        coherence,
        valid,
        f"{row.Seed_id} / {vrow.Seed_id}: {len(powers_p)} paired segments",
    )
    plot = output / "spectra_coherence_transfer.png"
    figure.savefig(plot, dpi=140)
    summary["files"].update(spectra=str(npz), plot=str(plot))
    if show:
        plt.show()
    plt.close(figure)
    return finish("comparison_only")


def assert_workflows_finished(results_directory, *, ignore_stage=None):
    """Refuse master writes when any saved network batch is unfinished (not a lock)."""
    for path in (Path(results_directory) / "_workflow").glob("*.json"):
        if path.stem == ignore_stage:
            continue
        report = json.loads(path.read_text())
        keys = report.get("selected_epoch_keys", [])
        if keys and sum(report.get("status_counts", {}).values()) != len(keys):
            raise ValueError(
                f"Unfinished {path.stem} workflow: run without master updates"
            )


def view_saved(summary_path, *, show=True):
    """Display a saved PNG without loading raw waveforms or updating the master."""
    summary = json.loads(Path(summary_path).read_text())
    plot = summary.get("files", {}).get("plot")
    if not plot or not Path(plot).is_file():
        raise ValueError("This paired result has no saved plot")
    figure, axis = plt.subplots(figsize=(10, 11))
    axis.imshow(plt.imread(plot))
    axis.axis("off")
    figure.tight_layout()
    if show:
        plt.show()
    plt.close(figure)
    return plot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-dir", type=Path, default=OUTPUT_DIRECTORY)
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--network", default="YV")
    parser.add_argument("--station", default="RR52")
    parser.add_argument("--pressure-id")
    parser.add_argument("--vertical-id")
    parser.add_argument("--output-directory", type=Path)
    parser.add_argument("--event-resample-rate", type=float)
    parser.add_argument("--native-rate-tolerance-ppm", type=float, default=0)
    parser.add_argument("--start")
    parser.add_argument("--days", type=int, default=7)
    parser.add_argument("--segment-seconds", type=float, default=3600)
    parser.add_argument("--frequency-min", type=float, default=0.01)
    parser.add_argument("--frequency-max", type=float, default=0.5)
    parser.add_argument("--update-master", action="store_true")
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--plot-only", type=Path, metavar="COMPARISON_JSON")
    args = vars(parser.parse_args())
    plot_only = args.pop("plot_only")
    if plot_only is not None:
        print(view_saved(plot_only))
        return
    summary = run_paired_spectra(**args)
    print(f"{summary['status']}: {summary['files']['result']}")


if __name__ == "__main__":
    main()
