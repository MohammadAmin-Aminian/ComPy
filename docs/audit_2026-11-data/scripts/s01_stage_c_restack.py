"""Stage C of the whole-archive run, recomputed from cached per-window spectra.

Pass ``cache`` (ObsPy): for each station, read the cleaned data exactly as
``pipeline/c_compliance.py`` does (``load_windows``) and save the per-window
Welch spectra of ``compy.compliance.spectra`` up to ``FMAX`` (at least the
station's f_c), with window start times. Nothing in the archive is written.

Pass ``restack`` (no new waveform reads): run stage C's own ``run_gain`` on the
cache, through ``compy.compliance.stacked_compliance_from_spectra``, for four
variants:

* ``gji_inclusive``: GJI gates, the published gate (must equal 06_compliance);
* ``wide_inclusive``: the wide gates with the 15 dB outlier cut of
  ``scripts/stage_c_wide.py`` (must equal 06_compliance_wide, used by stage E);
* ``gji_exclusive``, ``wide_exclusive``: the same with ``gate="exclusive"``.

Outputs (``<OUT>/stage_c/<variant>/<NET>/results/<STA>/``) have the stage C
format, so ``compy.results.load_station`` reads them, plus
``<STA>_window_levels.npy`` (accepted windows: time, log10 median |S_pz|/S_pp
in the band, pressure dB, acceleration dB) for ``compy.anomaly.scan_levels``,
and ``<STA>_window_table.npz`` (every window: time, band coherence, levels,
band-median Z/P phase) for the acceptance and tide analyses.

Usage (conda env bruitfm, from the repository root)::

    python docs/audit_2026-11-data/scripts/s01_stage_c_restack.py --workers 10
    python docs/audit_2026-11-data/scripts/s01_stage_c_restack.py --stations YV.RR52

Stations: every station of ``12_inversion_all_stations/tables/station_plan.csv``
(at least 100 windows in the wide pass). Stations with a stable epoch
(``epoch_end``) are also restacked with stage C's ``--end`` rule.
"""

import argparse
import csv
import json
import multiprocessing as mp
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = Path(os.environ.get("COMPY_AUDIT_OUT", ARCHIVE / "13_audit"))
CACHE = OUT / "spectra_cache"
STAGE_C = OUT / "stage_c"
OVERRIDE = ARCHIVE / "00_inventory/configs/compliance_override.py"
FMAX = 0.05  # Hz; raised to 1.05 f_c where needed
# scripts/stage_c_wide.py
WIDE_P, WIDE_A, OUTLIER_DB = (10.0, 90.0), (-170.0, -100.0), 15.0
VARIANTS = ("gji_inclusive", "wide_inclusive", "gji_exclusive", "wide_exclusive")
REFERENCE = {"gji": ARCHIVE / "06_compliance", "wide": ARCHIVE / "06_compliance_wide"}


def plan():
    """(network, station, epoch_end or None) of the stations to process."""
    path = ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv"
    with open(path) as fh:
        return [
            (r["network"], r["station"], r["epoch_end"] or None)
            for r in csv.DictReader(fh)
        ]


def _stage_c_module(net, sta):
    """Import pipeline/c_compliance.py configured for one archive station (read only)."""
    os.environ.update(P_NET=net, COMPY_STATION=sta, P_ROOT=str(ARCHIVE / "06_compliance"),
                      COMPY_RUN_DIR=str(ARCHIVE / "06_compliance" / net),
                      COMPY_PIPELINE_CONFIG=str(OVERRIDE), MPLBACKEND="Agg")  # fmt: skip
    sys.path[:0] = [str(REPO / "pipeline"), str(REPO)]
    import c_compliance

    return c_compliance


# --------------------------------------------------------------------------
# Pass 1: per-window spectra
# --------------------------------------------------------------------------
def cache_station(net, sta):
    path = CACHE / net / f"{sta}.npz"
    if path.exists():
        return "cached"
    c = _stage_c_module(net, sta)
    from obspy import read_inventory

    from compy.compliance import spectra

    inv = read_inventory(str(c.INVENTORY)).select(station=sta)
    depth = -inv[0][0].elevation
    fc = float(np.sqrt(c.GRAVITY / (2 * np.pi * depth)))
    windows, dropped = c.load_windows(sta)
    if not windows:
        raise RuntimeError("no windows")
    fmax = max(FMAX, 1.05 * fc)
    starts, zz, pp, zp, coh, nsamp = [], [], [], [], [], []
    keep = None
    for w in windows:
        f, a, b, x, k = spectra(w, c.NSEG)
        if keep is None:
            keep = f <= fmax
            f_keep = f[keep]
        starts.append(float(w[0].stats.starttime))
        zz.append(a[keep])
        pp.append(b[keep])
        zp.append(x[keep])
        coh.append(k[keep])
        nsamp.append(len(w[0]))
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, f=f_keep, starts=np.array(starts), zz=np.array(zz), pp=np.array(pp),
             zp=np.array(zp), coh=np.array(coh), n_samples=np.array(nsamp), depth=depth,
             fc=fc, dropped=dropped, gains=json.dumps(c.gains_for(sta)),
             nseg=c.NSEG, window_s=c.WINDOW_S)  # fmt: skip
    tmp.rename(path)
    return f"{len(starts)} windows"


# --------------------------------------------------------------------------
# Pass 2: restack
# --------------------------------------------------------------------------
def _month(t):
    d = datetime.fromtimestamp(t, tz=timezone.utc)
    return f"{d.year}-{d.month:02d}"


def restack_station(net, sta, epoch_end=None):
    c = _stage_c_module(net, sta)
    from compy.compliance import SelectionCriteria, stacked_compliance_from_spectra

    d = np.load(CACHE / net / f"{sta}.npz")
    f, n_samples = d["f"], int(d["n_samples"][0])
    if np.any(d["n_samples"] != n_samples):
        raise RuntimeError("windows of different length")
    depth, fc = float(d["depth"]), float(d["fc"])
    rows_all = list(zip(d["starts"], d["zz"], d["pp"], d["zp"], d["coh"]))
    gains = json.loads(str(d["gains"]))
    epochs = [(None, rows_all)]
    if epoch_end:
        from obspy import UTCDateTime

        end = float(UTCDateTime(epoch_end))
        width = n_samples / 1.25  # stage C keeps windows with endtime <= end
        epochs.append(
            (epoch_end, [r for r in rows_all if r[0] + (n_samples - 1) / 1.25 <= end])
        )
        del width
    report = {}
    for variant in VARIANTS:
        gates, gate = variant.split("_")
        for end, rows in epochs:
            g = dict(gains)
            if end:
                cal = json.loads(
                    (c.WORK / sta / "calibration" / "compy_gain.json").read_text()
                )
                from obspy import UTCDateTime

                before = [x for t, x in zip(cal["events_used"], cal["event_gains"])
                          if UTCDateTime(t) < UTCDateTime(end)]  # fmt: skip
                if before:
                    g["compy"] = float(np.median(before))
            report[(variant, end)] = _restack_variant(c, SelectionCriteria, stacked_compliance_from_spectra,
                                                      net, sta, f, n_samples, depth, fc, rows, g, gates, gate,
                                                      end, len(d["starts"]), int(d["dropped"]))  # fmt: skip
    return report


def _restack_variant(c, SelectionCriteria, from_spectra, net, sta, f, n_samples, depth, fc,
                     rows, gains, gates, gate, end, n_all, dropped):  # fmt: skip
    tightened = []

    def stack(windows, depth_, gain, band, criteria, nseg, *, output_band=None, gate="inclusive",
              gravity_correction=True):  # fmt: skip
        if gates == "wide" and not any(criteria is x for x in tightened):
            # scripts/stage_c_wide.py: median levels of the windows the wide criteria pass
            # (inclusive semantics), then upper limits at most OUTLIER_DB above them.
            sel = (f >= band[0]) & (f <= band[1])
            levels = []
            for _, zz, pp, _zp, coh in windows:
                pg = pp * gain**2
                if criteria.check(f, zz, pg, coh, sel) is None:
                    levels.append((10 * np.log10(np.mean(pg[sel])),
                                   10 * np.log10(np.mean(zz[sel] * (2 * np.pi * f[sel]) ** 4))))  # fmt: skip
            if not levels:
                raise ValueError("no windows passed compliance quality thresholds")
            p, a = np.median(levels, axis=0)
            criteria = SelectionCriteria(criteria.min_median_coherence,
                                         (WIDE_P[0], min(WIDE_P[1], p + OUTLIER_DB)),
                                         (WIDE_A[0], min(WIDE_A[1], a + OUTLIER_DB)), criteria.name)  # fmt: skip
            tightened.append(criteria)
        return from_spectra(f, windows, depth_, gain, band, criteria, nseg, n_samples=n_samples,
                            output_band=output_band, gate=gate,
                            gravity_correction=gravity_correction)  # fmt: skip

    c.stacked_compliance = stack
    c.GATE = gate
    if gates == "wide":
        c.SelectionCriteria = lambda thr, _p, _a, name: SelectionCriteria(
            thr, WIDE_P, WIDE_A, name.replace("gji2025", "wide")
        )
    else:
        c.SelectionCriteria = SelectionCriteria
    out = STAGE_C / f"{gates}_{gate}" / net / "results" / sta
    if end:
        out = out / "epoch"
    out.mkdir(parents=True, exist_ok=True)
    summary = {"station": sta, "depth_m": depth, "f_c_hz": fc, "windows_total": len(rows),
               "windows_dropped_events": dropped, "window_s": c.WINDOW_S, "nseg": c.NSEG,
               "gains": gains, "gate": gate, "gates": gates, "epoch_end": end, "results": {},
               "source": "13_audit restack of cached stage C spectra"}  # fmt: skip
    for label, gain in gains.items():
        result, threshold, attempts = c.run_gain(rows, depth, gain, fc)
        entry = {"gain": gain, "attempts": {str(k): v for k, v in attempts.items()}}
        if result is None:
            summary["results"][label] = entry
            continue
        in_band = (result.freqs >= result.band[0]) & (result.freqs <= result.band[1])
        accepted = {float(t) for t in result.window_starts}
        months = {}
        for r in rows:
            if float(r[0]) in accepted:
                months.setdefault(_month(float(r[0])), []).append(r)
        monthly = {}
        for key, ws in sorted(months.items()):
            if len(ws) < 20:
                continue
            m = stack(ws, depth, gain, result.band, result.criteria, c.NSEG,
                      output_band=result.band, gate=gate)  # fmt: skip
            monthly[key] = {"n": m.n_windows, "compliance": m.compliance.tolist(),
                            "uncertainty": m.uncertainty.tolist()}  # fmt: skip
        entry.update({
            "threshold": threshold, "band_hz": result.band, "n_windows": result.n_windows,
            "n_averages": result.n_averages, "rejected": result.rejected,
            "median_coherence_in_band": float(np.median(result.coherence[in_band])),
            "median_phase_deg_in_band": float(np.median(np.abs(result.phase_deg[in_band]))),
            "median_relative_uncertainty": float(np.median(
                result.uncertainty[in_band] / result.compliance[in_band])),
            "monthly": monthly,
            "criteria": {"pressure_db": list(result.criteria.pressure_psd_db),
                         "acceleration_db": list(result.criteria.acceleration_psd_db)},
        })  # fmt: skip
        summary["results"][label] = entry
        np.savez(out / f"{sta}_compliance_{label}.npz", freqs=result.freqs,
                 compliance=result.compliance, uncertainty=result.uncertainty,
                 coherence=result.coherence, phase_deg=result.phase_deg,
                 window_mean=result.window_mean, window_std=result.window_std,
                 band=np.array(result.band), all_freqs=result.all_freqs,
                 mean_coherence_all=result.mean_coherence_all, S_zz=result.stack["S_zz"],
                 S_pp=result.stack["S_pp"], S_pz=result.stack["S_pz"], gain=gain, depth=depth,
                 window_starts=np.array([float(t) for t in result.window_starts]),
                 bin_counts=result.bin_counts)  # fmt: skip
        if gate == "inclusive":
            sel = (f >= result.band[0]) & (f <= result.band[1])
            levels, table = [], []
            for t, zz, pp, zp, coh in rows:
                pg, zg = pp * gain**2, zp * gain
                tf = np.abs(zg[sel]) / pg[sel]
                lev = (float(t), float(np.log10(np.median(tf))), float(10 * np.log10(np.mean(pg[sel]))),
                       float(10 * np.log10(np.mean(zz[sel] * (2 * np.pi * f[sel]) ** 4))))  # fmt: skip
                table.append(lev + (float(np.median(coh[sel])),
                                    float(np.degrees(np.angle(np.sum(zg[sel]))))))  # fmt: skip
                if float(t) in accepted:
                    levels.append(lev)
            np.save(out / f"{sta}_window_levels.npy", np.array(levels))
            tab = np.array(table)
            np.savez(out / f"{sta}_window_table_{label}.npz", time=tab[:, 0], log10_tf=tab[:, 1],
                     pressure_db=tab[:, 2], acceleration_db=tab[:, 3], coherence=tab[:, 4],
                     phase_deg=tab[:, 5], accepted=np.isin(tab[:, 0], list(accepted)),
                     band=np.array(result.band))  # fmt: skip
    (out / f"{sta}_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    return compare(net, sta, gates, gate, end, out)


def compare(net, sta, gates, gate, end, out):
    """Maximum relative difference from the archived stage C result (inclusive gate only)."""
    if gate != "inclusive":
        return None
    ref = REFERENCE[gates] / net / "results" / sta
    if end:
        ref = ref / "epoch"
    diffs = {}
    for path in sorted(out.glob(f"{sta}_compliance_*.npz")):
        r = ref / path.name
        if not r.exists():
            diffs[path.stem] = "no reference"
            continue
        a, b = np.load(path), np.load(r)
        if len(a["freqs"]) != len(b["freqs"]) or not np.array_equal(
            a["window_starts"], b["window_starts"]
        ):
            diffs[path.stem] = (
                f"windows {len(a['window_starts'])} vs {len(b['window_starts'])}"
            )
            continue
        keys = (
            "compliance",
            "uncertainty",
            "coherence",
            "phase_deg",
            "window_mean",
            "window_std",
        )
        diffs[path.stem] = max(float(np.max(np.abs(a[k] - b[k]) / np.maximum(np.abs(b[k]), 1e-300)))
                               for k in keys)  # fmt: skip
    return diffs


# --------------------------------------------------------------------------
def worker(args):
    net, sta, end, stage = args
    t0 = time.time()
    log = OUT / "logs" / f"{net}.{sta}.{stage}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    try:
        if stage == "cache":
            msg = cache_station(net, sta)
        else:
            rep = restack_station(net, sta, end)
            msg = json.dumps(
                {f"{v}|{e or ''}": r for (v, e), r in rep.items()}, default=str
            )
        log.write_text(msg + "\n")
        return net, sta, stage, "ok", msg, time.time() - t0
    except Exception:
        log.write_text(traceback.format_exc())
        return (
            net,
            sta,
            stage,
            "error",
            traceback.format_exc().splitlines()[-1],
            time.time() - t0,
        )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--stage", choices=("cache", "restack", "all"), default="all")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument(
        "--stations", nargs="*", help="NET.STA codes (default: the stage E plan)"
    )
    args = ap.parse_args()
    todo = plan()
    if args.stations:
        todo = [x for x in todo if f"{x[0]}.{x[1]}" in args.stations]
    stages = ("cache", "restack") if args.stage == "all" else (args.stage,)
    ctx = mp.get_context("spawn")
    status = OUT / "logs" / "status.jsonl"
    status.parent.mkdir(parents=True, exist_ok=True)
    for stage in stages:
        tasks = [(n, s, e, stage) for n, s, e in todo]
        with (
            ctx.Pool(args.workers, maxtasksperchild=1) as pool,
            open(status, "a") as fh,
        ):
            for i, res in enumerate(pool.imap_unordered(worker, tasks), 1):
                fh.write(
                    json.dumps(
                        dict(
                            zip(
                                ("net", "sta", "stage", "status", "msg", "seconds"), res
                            )
                        )
                    )
                    + "\n"
                )
                fh.flush()
                print(f"[{i}/{len(tasks)}] {res[0]}.{res[1]} {stage} {res[3]} {res[5]:.0f}s {res[4][:160]}",
                      flush=True)  # fmt: skip


if __name__ == "__main__":
    main()
