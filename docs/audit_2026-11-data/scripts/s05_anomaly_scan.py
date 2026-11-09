"""Anomaly scan of every inverted station, with cross-checks that need the data.

Inputs (all produced by the earlier scripts or by the whole-archive run):

* ``13_audit/stage_c/wide_inclusive``: the stage C result stage E inverted
  (bit-identical to 06_compliance_wide), now with per-window levels;
* ``13_audit/stage_c/wide_exclusive``: the same with the exclusive gate;
* ``12_inversion_all_stations``: stage E records (variant ``final``);
* ``13_audit/phase``, ``13_audit/dpg_tide``, ``13_audit/level_steps``;
* stage D event gains (``06_compliance/<NET>/work/<STA>/calibration/compy_gain.json``).

Runs ``compy.anomaly.scan_station`` and ``scan_network`` (per network), then adds:

* ``gate_bias``: exclusive minus inclusive band-mean compliance [%];
* ``which_sensor``: for every |Z/P| changepoint, the relative Z and P level steps
  (``s04``) within 10 days, and the stage D event gains before/after;
* ``dpg_response``: phase-fitted leak-corner factor, the tidal ratio R_M2 and the
  compliance bias at the band bottom implied by the leak-only model after a gain
  fitted at 20-60 mHz.

Outputs in ``13_audit/anomaly``: ``anomalies.csv``/``.md`` (ranked),
``stations.csv`` (one row per station), ``network_info.json``.
"""

import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).parent))
from compy import anomaly  # noqa: E402
from compy.results import load_station  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
AUDIT = Path(os.environ.get("COMPY_AUDIT_OUT", ARCHIVE / "13_audit"))
OUT = AUDIT / "anomaly"


def ts(iso):
    return (
        datetime.fromisoformat(iso.replace("Z", ""))
        .replace(tzinfo=timezone.utc)
        .timestamp()
    )


def read_csv(path):
    return list(csv.DictReader(open(path))) if Path(path).exists() else []


def leak_bias(row, f_band):
    """Compliance bias [%] at the band bottom implied by the leak-only model (k_lowest)."""
    import s03_dpg_tide as tide
    from obspy import read_inventory

    net, sta = row["sid"].split(".")
    inv = read_inventory(
        str(ARCHIVE / "00_inventory/stationxml" / f"{net}.xml")
    ).select(station=sta)
    chans = [c for n in inv for s in n for c in s if c.code == row["p_channel"]]
    if not chans:
        return None
    stages = tide.paz(chans[0])
    fcal = np.geomspace(0.02, 0.06, 9)
    out = []
    for k in (row["k_lowest"], row["k_lowest_lo"], row["k_lowest_hi"]):
        k = float(k)
        # compliance ~ 1/|P_est| = |H_meta / H_true|, normalised in the calibration band
        r = lambda f: np.abs(
            tide.analog(stages, f) / tide.analog(stages, f, k, "lowest")
        )  # noqa: E731
        out.append(100 * (r(f_band[:1])[0] / np.mean(r(fcal)) - 1))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plan = read_csv(ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv")
    inv_res = {
        r["sid"]: r
        for r in read_csv(
            ARCHIVE / "12_inversion_all_stations/tables/inversion_results.csv"
        )
    }
    phase = {r["sid"]: r for r in read_csv(AUDIT / "phase/phase_metadata.csv")}
    tide_rows = {r["sid"]: r for r in read_csv(AUDIT / "dpg_tide/dpg_tide.csv")}
    level_steps = read_csv(AUDIT / "level_steps/level_steps.csv")
    stations_by_net, infos_by_net, all_anoms, rows = {}, {}, [], []
    for p in plan:
        net, sta, sid = p["network"], p["station"], p["sid"]
        sub = Path(net) / "results" / sta / ("epoch" if p["epoch_end"] else "")
        inc, exc = (
            AUDIT / "stage_c/wide_inclusive" / sub,
            AUDIT / "stage_c/wide_exclusive" / sub,
        )
        label = p["wide_gain_used"] or None
        if not (inc / f"{sta}_summary.json").exists():
            continue
        s = load_station(
            inc,
            sta,
            label,
            inversion_dirs=[
                ARCHIVE
                / "12_inversion_all_stations"
                / net
                / "results"
                / sta
                / "inversion"
            ],
        )
        found, info = anomaly.scan_station(s)
        for a in found:
            a.station = sid
        stations_by_net.setdefault(net, {})[sid] = s
        infos_by_net.setdefault(net, {})[sid] = info
        row = {"sid": sid, "network": net, "tier": p["tier"], "quality": inv_res.get(sid, {}).get("quality"),
               "depth_m": s.depth, "label": s.label, "gain": s.gain, "n_windows": len(s.window_starts),
               "band_lo_mhz": 1e3 * s.band[0], "band_hi_mhz": 1e3 * s.band[1],
               "epoch_end": p["epoch_end"]}  # fmt: skip
        # exclusive gate
        try:
            e = load_station(exc, sta, label)
            sel_i, sel_e = s.in_band, e.in_band
            common = np.intersect1d(s.freqs[sel_i], e.freqs[sel_e])
            # bins where the exclusive gate keeps at least half the inclusive windows
            counts = np.load(exc / f"{sta}_compliance_{e.label}.npz")["bin_counts"]
            n_bin = np.interp(common, e.all_freqs, counts)
            common = common[n_bin >= 0.5 * len(s.window_starts)]
            if len(common):
                ci = np.interp(common, s.freqs, s.compliance)
                ce = np.interp(common, e.freqs, e.compliance)
                d = 100 * (ce / ci - 1)
                med = float(np.median(d))
                row.update(gate_bias_median_pct=med, gate_bias_lo_pct=float(d[0]),
                           gate_bias_hi_pct=float(d[-1]), n_windows_exclusive=len(e.window_starts),
                           acceptance=s.acceptance, gate_bins=len(common))  # fmt: skip
                if abs(med) >= 1.0:
                    found.append(anomaly.Anomaly(sid, "selection", "gate_bias",
                        f"the published (inclusive) gate is {-med:+.2f} % above the exclusive one "
                        f"(median of {len(common)} bins; {-d[0]:+.2f} % at {1e3 * common[0]:.1f} mHz; "
                        f"acceptance {100 * (s.acceptance or 0):.0f} %)",
                        min(10.0, abs(med) / 0.3), float(-med),
                        evidence={"per_bin_pct": np.round(d, 2).tolist()}))  # fmt: skip
        except FileNotFoundError:
            pass
        # level steps and stage D gains around each |Z/P| changepoint
        cps = (
            info.get("levels", {}).get("changepoints", {}).get("transfer function", [])
        )
        cal_path = (
            ARCHIVE
            / "06_compliance"
            / net
            / "work"
            / sta
            / "calibration"
            / "compy_gain.json"
        )
        cal = json.loads(cal_path.read_text()) if cal_path.exists() else {}
        ev_t = (
            np.array([ts(t) for t in cal.get("events_used", [])])
            if cal.get("events_used")
            else np.array([])
        )
        ev_g = np.array(cal.get("event_gains", []), float)
        which = []
        for c in cps:
            near = [r for r in level_steps if r["network"] == net and r["station"] == sta
                    and abs(ts(r["time"]) - c.time) < 10 * 86400]  # fmt: skip
            z = [float(r["shift_db"]) for r in near if r["channel"] == "Z"]
            pz = [float(r["shift_db"]) for r in near if r["channel"] == "P"]
            g_before, g_after = ev_g[ev_t < c.time], ev_g[ev_t >= c.time]
            gain_ratio = (
                float(np.median(g_after) / np.median(g_before))
                if len(g_before) and len(g_after)
                else None
            )
            sensor = (
                "Z"
                if z and not pz
                else "P"
                if pz and not z
                else "both"
                if z and pz
                else "neither"
            )
            which.append({"time": anomaly._iso(c.time), "tf_shift_pct": anomaly._pct_from_log10(c.shift),
                          "z_level_db": z, "p_level_db": pz, "sensor": sensor,
                          "stage_d_gain_ratio": gain_ratio, "events_before": int(len(g_before)),
                          "events_after": int(len(g_after))})  # fmt: skip
        row["steps"] = json.dumps(which)
        row["n_steps"] = len(which)
        row["epoch_shift_pct"] = info.get("levels", {}).get("stack_shift_pct")
        # DPG response
        ph, td = phase.get(sid, {}), tide_rows.get(sid, {})
        row.update(
            {
                k: ph.get(k)
                for k in (
                    "pressure_type",
                    "sensor_family",
                    "res_lo",
                    "res_hi",
                    "rms_timing",
                    "rms_dpg",
                    "rms_seis",
                )
            }
        )
        row.update(
            {
                k: td.get(k)
                for k in (
                    "k_lowest",
                    "k_lowest_lo",
                    "k_lowest_hi",
                    "R_M2",
                    "R_pred_lowest",
                    "snr",
                )
            }
        )
        if td.get("k_lowest") and ph.get("p_channel"):
            try:
                b = leak_bias(dict(td, p_channel=ph["p_channel"]), s.freqs[s.in_band])
                if b:
                    row.update(
                        leak_bias_pct=b[0],
                        leak_bias_range=f"{min(b[1:]):+.1f}..{max(b[1:]):+.1f}",
                    )
                    if abs(b[0]) >= 3 and td.get("tide_status") == "ok":
                        rm2 = float(td["R_M2"])
                        agree = (rm2 > 1) == (float(td["k_lowest"]) < 1)
                        found.append(anomaly.Anomaly(sid, "instrument", "dpg_response",
                            f"pressure response mismatch: compliance phase gives leak corner x{float(td['k_lowest']):.2f}, "
                            f"tide (M2) amplitude after deconvolution is {20 * np.log10(rm2):+.1f} dB "
                            f"({'same' if agree else 'opposite'} sign); implied compliance bias at {s.freqs[s.in_band][0] * 1e3:.1f} mHz "
                            f"{b[0]:+.1f} % after the 20-60 mHz gain",
                            min(10.0, abs(float(ph.get('res_lo') or 0)) / 1.0), float(b[0]),
                            evidence={"k_lowest": float(td["k_lowest"]), "R_M2": rm2, "agree": agree}))  # fmt: skip
            except Exception as exc:  # noqa: BLE001
                row["leak_bias_error"] = str(exc)
        all_anoms += found
        rows.append(row)
        print(f"{sid:12s} {len(found):3d} anomalies, steps {len(which)}", flush=True)
    network_info = {}
    for net, st in stations_by_net.items():
        if len(st) < 3:
            continue
        found, info = anomaly.scan_network(st, infos_by_net[net], name=f"network {net}")
        all_anoms += found
        network_info[net] = {
            k: v for k, v in info.items() if k in ("phase_sign", "phase_vs_sediment")
        }
        cm = info.get("common_mode")
        if cm:
            network_info[net]["common_mode"] = {"p_value": cm["p_value"], "mean_corr": cm["network_mean_corr"],
                                                "keys": cm["keys"], "median_dev": np.nanmedian(cm["matrix"], axis=0).tolist()}  # fmt: skip
        network_info[net]["synchronous_steps"] = [[(anomaly._iso(t), c, anomaly._pct_from_log10(x)) for t, c, x in g]
                                                  for g in info.get("synchronous_steps", []) if len({c for _, c, _ in g}) >= 2]  # fmt: skip
    ranked = anomaly.rank(all_anoms)
    # The score saturates at 40 (sigma capped at 8, effect at 20 %), so many records tie;
    # break ties by tier (tier 1 = interpreted stations first), then by |effect|.
    tiers = {p["sid"]: p["tier"] for p in plan}
    ranked.sort(key=lambda a: (-round(a.score, 6), tiers.get(a.station, "1") != "1",
                               -(abs(a.effect_pct) if np.isfinite(a.effect_pct) else 0)))  # fmt: skip
    anomaly.to_csv(ranked, OUT / "anomalies.csv")
    (OUT / "anomalies.md").write_text(anomaly.to_markdown(ranked))
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "stations.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    (OUT / "network_info.json").write_text(
        json.dumps(network_info, indent=1, default=str)
    )
    print(f"{len(ranked)} ranked anomalies, {len(rows)} stations")


if __name__ == "__main__":
    main()
