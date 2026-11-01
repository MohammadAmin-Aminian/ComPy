"""Stage B: clean one station's 1.25 Hz data for compliance (ComPy 2.2).

Order (Aminian et al. 2025, Fig. 2, adapted for unattended processing):
  1. Common coverage of the Z, 1, 2 and pressure channels; segments shorter than MIN_SEGMENT_H
     are skipped.
  2. Earthquake spans from the USGS catalogue, M >= 5.85, 1.5 days per
     magnitude unit (Deen et al. 2017; thesis 3.1.3).
  3. Periodic-transient (glitch) detection per DEGLITCH_BLOCK_D-day block on
     each channel; a stacked template is subtracted only when the fold
     strength reaches DEGLITCH_MIN_STRENGTH.
  4. Hourly tilt angles (TiSKitPy CleanRotator, hours with earthquakes
     skipped), median smoothed over ROT_SMOOTH_H hours and applied per sample.
  5. Local events by recursive STA/LTA on the rotated vertical.
  6. Removal of horizontal-coherent noise from the vertical in equal windows of at most
     CLEAN_WINDOW_H hours, transfer functions estimated outside event spans.
Pressure is never modified except by step 3 (if a transient is detected).

Outputs in work/<STA>/clean/: segment MiniSEED files, spans.json,
tilt_hourly.csv, deglitch.json, summary.json.
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from obspy import Stream, UTCDateTime, read

from config import CATALOG, CHANNELS, COMPY, DEGLITCH, P, WORK, Z

sys.path.insert(0, str(COMPY))

import compy  # noqa: E402
import compy_cleaning as cleaning  # noqa: E402
import tiskitpy as tiskit  # noqa: E402

MIN_SEGMENT_H = 13
DEGLITCH_BLOCK_D = 15
DEGLITCH_MIN_STRENGTH = 2.0
ROT_SMOOTH_H = 24
ROT_LIMIT_DEG = 20.0
CLEAN_WINDOW_H = 24
EQ_MINMAG = 5.85
EQ_DAYS_PER_MAG = 1.5


def load_segments(station):
    files = sorted((WORK / station / "phys").glob(f"{station}.*.mseed"))
    st = Stream()
    for f in files:
        st += read(str(f))
    st.merge(method=0)
    intervals = None
    for channel in CHANNELS:
        parts = [
            (tr.stats.starttime, tr.stats.endtime)
            for tr in st.select(channel=channel).split()
        ]
        if intervals is None:
            intervals = parts
            continue
        intervals = [
            (max(a0, b0), min(a1, b1))
            for a0, a1 in intervals
            for b0, b1 in parts
            if min(a1, b1) - max(a0, b0) > 0
        ]
    segments = []
    for start, end in intervals or []:
        if end - start < MIN_SEGMENT_H * 3600:
            continue
        seg = st.slice(start, end).split()
        seg.traces = [tr for tr in seg if tr.stats.npts > 0]
        if sorted(tr.stats.channel for tr in seg) == sorted(CHANNELS):
            for tr in seg:
                tr.data = tr.data.astype(np.float64)
            segments.append(seg)
    return segments


def overlaps(t0, t1, spans):
    return any(s < t1 and e > t0 for s, e in spans)


def deglitch(seg, eq_spans):
    report = []
    t_start, t_end = seg[0].stats.starttime, seg[0].stats.endtime
    block = DEGLITCH_BLOCK_D * 86400
    for tr in seg:
        t0 = t_start
        while t0 < t_end:
            t1 = min(t0 + block, t_end)
            part = tr.slice(t0, t1)
            if part.stats.npts * part.stats.delta < 3 * 86400:
                break
            detect = part.copy()
            detect.detrend("linear")
            if tr.stats.channel != P:
                detect.differentiate()
                detect.differentiate()
            detect.filter("bandpass", freqmin=0.002, freqmax=0.3, zerophase=True)
            found = cleaning.periodic_template(
                detect,
                period_range=(3596.0, 3604.0),
                period_step=0.02,
                exclude=eq_spans,
            )
            entry = {
                "channel": tr.stats.channel,
                "start": str(t0),
                "end": str(t1),
                "period_s": found.period,
                "strength": found.strength,
                "applied": False,
            }
            if found.strength >= DEGLITCH_MIN_STRENGTH:
                disp = part.copy()
                disp.detrend("linear")
                disp.filter("highpass", freq=0.001, zerophase=True)
                template = cleaning.periodic_template(
                    disp,
                    period_range=(found.period - 0.05, found.period + 0.05),
                    period_step=0.01,
                    exclude=eq_spans,
                )
                fs = tr.stats.sampling_rate
                i0 = int(round((t0 - t_start) * fs))
                times = t0.timestamp + np.arange(part.stats.npts) / fs
                tr.data[i0 : i0 + part.stats.npts] -= template.evaluate(times)
                entry.update(applied=True, template_rms=float(np.std(template.template)),
                             data_rms=float(np.std(disp.data)))  # fmt: skip
            report.append(entry)
            t0 = t1
    return report


def hourly_tilt(seg, eq_spans):
    rows = []
    for window in compy.split_stream(seg.select(component="[Z12]"), 3600):
        t0 = window[0].stats.starttime
        if overlaps(t0, t0 + 3600, eq_spans):
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rot = tiskit.CleanRotator(window, remove_eqs=False, verbose=False)
        except Exception:
            continue
        if np.isfinite(rot.angle) and rot.angle <= ROT_LIMIT_DEG:
            rows.append((t0 + 1800, float(rot.angle), float(rot.azimuth % 360),
                         float(rot.variance_reduction)))  # fmt: skip
    return rows


def as_timespans(spans):
    if not spans:
        return None
    return tiskit.TimeSpans(
        start_times=[s for s, _ in spans], end_times=[e for _, e in spans]
    )


def main():
    import logging

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("station")
    args = parser.parse_args()
    logging.disable(logging.INFO)
    station = args.station
    out_dir = WORK / station / "clean"
    out_dir.mkdir(parents=True, exist_ok=True)
    segments = load_segments(station)
    if not segments:
        raise SystemExit(f"{station}: no usable segments")
    first, last = segments[0][0].stats.starttime, segments[-1][0].stats.endtime
    eqs = tiskit.TimeSpans.from_eqs(
        (first, last), minmag=EQ_MINMAG, days_per_magnitude=EQ_DAYS_PER_MAG,
        eq_file=str(CATALOG), save_eq_file=False,
    )  # fmt: skip
    eq_spans = list(zip(eqs.start_times, eqs.end_times))
    all_local, tilt_rows, glitch_rows, seg_info = [], [], [], []
    for k, seg in enumerate(segments):
        start, end = seg[0].stats.starttime, seg[0].stats.endtime
        print(f"{station} segment {k}: {start} - {end}", flush=True)
        clock = UTCDateTime()
        if DEGLITCH:
            glitch_rows += deglitch(seg, eq_spans)
        print(f"  deglitch {'done' if DEGLITCH else 'skipped (config.DEGLITCH)'} "
              f"{UTCDateTime() - clock:.0f} s", flush=True)  # fmt: skip
        rows = hourly_tilt(seg, eq_spans)
        print(
            f"  tilt done {UTCDateTime() - clock:.0f} s ({len(rows)} hours)", flush=True
        )
        tilt_rows += rows
        if len(rows) >= 3:
            seg = cleaning.rotate_smoothly(
                seg, [r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows],
                ROT_SMOOTH_H,
            )  # fmt: skip
        local = cleaning.sta_lta_spans(seg.select(channel=Z)[0])
        all_local += local
        avoid = as_timespans(sorted(eq_spans + local))
        print(f"  rotation + STA/LTA done {UTCDateTime() - clock:.0f} s", flush=True)
        cleaned = compy.remove_coherent_noise(
            seg, CLEAN_WINDOW_H, avoid_spans=avoid, cover_all=True
        )
        print(
            f"  coherent-noise removal done {UTCDateTime() - clock:.0f} s", flush=True
        )
        for tr in cleaned:
            tr.data = tr.data.astype(np.float32)
        path = out_dir / f"{station}_seg{k:02d}.mseed"
        cleaned.write(str(path), format="MSEED", encoding="FLOAT32")
        seg_info.append({"file": path.name, "start": str(start), "end": str(end),
                         "tilt_estimates": len(rows), "local_events": len(local)})  # fmt: skip
    (out_dir / "spans.json").write_text(json.dumps({
        "earthquakes": [[str(s), str(e)] for s, e in eq_spans],
        "local": [[str(s), str(e)] for s, e in all_local],
        "rule": f"USGS M>={EQ_MINMAG}, {EQ_DAYS_PER_MAG} d per magnitude unit; STA/LTA 60/1800 s",
    }, indent=1))  # fmt: skip
    with open(out_dir / "tilt_hourly.csv", "w") as fh:
        fh.write("centre_utc,angle_deg,azimuth_deg,variance_reduction\n")
        for t, a, az, vr in tilt_rows:
            fh.write(f"{t},{a:.5f},{az:.3f},{vr:.5f}\n")
    (out_dir / "deglitch.json").write_text(json.dumps(glitch_rows, indent=1))
    summary = {"station": station, "segments": seg_info,
               "deglitch_applied": sum(r["applied"] for r in glitch_rows),
               "deglitch_blocks": len(glitch_rows), "deglitch_enabled": DEGLITCH,
               "rotation_smoothing_h": ROT_SMOOTH_H, "clean_window_h": CLEAN_WINDOW_H,
               "created": str(UTCDateTime())}  # fmt: skip
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
