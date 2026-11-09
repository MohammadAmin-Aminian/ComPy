"""Residual Z/P phase of every inverted station against its instrument metadata.

For each station of the stage E plan, the stacked phase of the wide-gate stage C
result (the one stage E inverted, stable epoch where used) is compared with 0/180
degrees (``compy.anomaly.phase_models``), and joined with:

* the station's instruments (``00_inventory/stations.csv``: seismometer family,
  pressure-sensor type, network);
* the analog (poles-and-zeros) part of the metadata responses that stage A
  deconvolved: corner frequencies and the phase they impose in the compliance
  band, for the vertical and for the pressure channel.

If a metadata corner is wrong by a factor k, the residual phase is a fixed
function of the metadata corner and k. A residual that follows one instrument
type, and is absent where that instrument has no corner in the band, points to
that instrument. Output: ``<OUT>/phase/phase_metadata.csv`` and ``.json``.
"""

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from compy import anomaly  # noqa: E402
from compy.results import load_station  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = Path(os.environ.get("COMPY_AUDIT_OUT", ARCHIVE / "13_audit")) / "phase"
RESULTS = Path(os.environ.get("PHASE_RESULTS", ARCHIVE / "06_compliance_wide"))


def analog_response(channel, f):
    """Product of the channel's poles-and-zeros stages at ``f`` (no gain, no FIR delays)."""
    h = np.ones_like(f, complex)
    corners = []
    for stage in channel.response.response_stages:
        if stage.__class__.__name__ != "PolesZerosResponseStage":
            continue
        rad = "RADIANS" in stage.pz_transfer_function_type.upper()
        s = 1j * 2 * np.pi * f if rad else 1j * f
        z = np.array(stage.zeros, complex)
        p = np.array(stage.poles, complex)
        h *= np.prod(s[:, None] - z[None, :], axis=1) / np.prod(
            s[:, None] - p[None, :], axis=1
        )
        scale = 1 / (2 * np.pi) if rad else 1.0
        corners += [abs(x) * scale for x in p if 0 < abs(x) * scale < 0.2]
    return h, sorted(corners)


def metadata(inventories, net, sta, z_code, p_code, t_mid):
    from obspy import UTCDateTime, read_inventory

    if net not in inventories:
        inventories[net] = read_inventory(
            str(ARCHIVE / "00_inventory/stationxml" / f"{net}.xml")
        )
    # Station epochs of the merged XML do not always cover the channel epochs (XO),
    # so select by station and pick the channel epoch that covers the record.
    inv = inventories[net].select(station=sta)
    t = UTCDateTime(t_mid)
    out = {}
    for key, code in (("Z", z_code), ("P", p_code)):
        chans = [c for n in inv for s in n for c in s if c.code == code]
        inside = [
            c
            for c in chans
            if c.start_date <= t and (c.end_date is None or c.end_date >= t)
        ]
        out[key] = (inside or chans or [None])[0]
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with open(ARCHIVE / "00_inventory/stations.csv") as fh:
        stations = {(r["network"], r["station"]): r for r in csv.DictReader(fh)}
    with open(ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv") as fh:
        plan = list(csv.DictReader(fh))
    with open(ARCHIVE / "12_inversion_all_stations/tables/inversion_results.csv") as fh:
        inv_res = {r["sid"]: r for r in csv.DictReader(fh)}
    inventories, rows, curves = {}, [], {}
    for p in plan:
        net, sta, sid = p["network"], p["station"], p["sid"]
        folder = RESULTS / net / "results" / sta
        if p["epoch_end"]:
            folder = folder / "epoch"
        try:
            s = load_station(folder, sta, p["wide_gain_used"] or None)
        except FileNotFoundError:
            continue
        sel = s.in_band
        if sel.sum() < 3:
            continue
        f = s.freqs[sel]
        m = anomaly.phase_models(f, s.phase_deg[sel])
        res = m["residual_deg"]
        info = stations.get((net, sta), {})
        meta = metadata(inventories, net, sta, info.get("Z", ""), info.get("P", ""),
                        float(np.median(s.window_starts)))  # fmt: skip
        row = {
            "sid": sid, "network": net, "station": sta, "depth_m": s.depth,
            "sensor_family": info.get("sensor_family"), "seismometer": info.get("seismometer"),
            "pressure_type": info.get("pressure_type"), "z_channel": info.get("Z"),
            "p_channel": info.get("P"), "f_lo": f[0], "f_hi": f[-1], "n_bins": len(f),
            "polarity": m["polarity"], "res_lo": res[0], "res_hi": res[-1],
            "res_mean": float(np.mean(res)), "res_slope_deg": float(res[0] - res[-1]),
            "coh_median": float(np.median(s.coherence[sel])),
            "rms_constant": m["constant"]["rms_deg"], "rms_timing": m["timing"]["rms_deg"],
            "dt_s": m["timing"]["params"]["dt_s"],
            "rms_dpg": m["dpg_corner"]["rms_deg"], "rms_seis": m["seismometer_corner"]["rms_deg"],
            "dpg_bias_lo": m["dpg_corner"]["bias_low_pct"][0], "dpg_bias_hi": m["dpg_corner"]["bias_high_pct"][0],
            "seis_bias_lo": m["seismometer_corner"]["bias_low_pct"][0],
            "seis_bias_hi": m["seismometer_corner"]["bias_high_pct"][0],
            "seis_period_true": m["seismometer_corner"]["params"]["period_true_s"],
            "quality": inv_res.get(sid, {}).get("quality"), "tier": p["tier"],
        }  # fmt: skip
        for key in ("Z", "P"):
            ch = meta[key]
            if ch is None:
                row.update({f"meta_{key}_corners_mhz": "", f"meta_{key}_phase_lo": np.nan,
                            f"meta_{key}_phase_hi": np.nan})  # fmt: skip
                continue
            h, corners = analog_response(ch, f)
            ph = np.degrees(np.unwrap(np.angle(h)))
            row[f"meta_{key}_corners_mhz"] = ";".join(f"{1e3 * c:.3g}" for c in corners)
            row[f"meta_{key}_lowest_mhz"] = 1e3 * corners[0] if corners else np.nan
            # phase the deconvolution removed, relative to the band top
            row[f"meta_{key}_phase_lo"] = float(ph[0] - ph[-1])
            row[f"meta_{key}_phase_hi"] = 0.0
            row[f"meta_{key}_sensor"] = (
                (ch.sensor.description or ch.sensor.model or "") if ch.sensor else ""
            )
        rows.append(row)
        curves[sid] = {"f": f.tolist(), "residual": res.tolist(), "phase": s.phase_deg[sel].tolist(),
                       "coherence": s.coherence[sel].tolist()}  # fmt: skip
        print(f"{sid:12s} {row['sensor_family'][:18]:18s} {row['pressure_type'][:10]:10s} "
              f"res {res[0]:+6.1f} -> {res[-1]:+6.1f}  metaP {row.get('meta_P_phase_lo', np.nan):+6.1f} "
              f"metaZ {row.get('meta_Z_phase_lo', np.nan):+6.1f}", flush=True)  # fmt: skip
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "phase_metadata.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    (OUT / "phase_curves.json").write_text(json.dumps(curves))


if __name__ == "__main__":
    main()
