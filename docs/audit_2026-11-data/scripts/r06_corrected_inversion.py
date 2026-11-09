"""Stage E with the DPG response corrected, on the eight YV paper stations.

The compliance stage E inverted (06_compliance_wide, stable epoch where used) is
corrected for the leak-corner factor k of each station (k from the M2 tide, r05,
where the tide is clear; otherwise from the phase, s03) with
``compy.response_check.correct_compliance``: compliance, uncertainty, per-window
statistics and monthly stacks are divided by (1 + bias(f)). Stage E variant
``final`` of ComPy 2.8.0 (the code of 12_inversion_all_stations, fixed seeds) is
then run on a new folder, so any change in Vs comes from the correction alone.

Output: ``13_audit/research/inversion_corrected/YV/results/<STA>/`` and
``inversion_corrected/comparison.csv`` (Vs at 0.5-5 km, sediment, chi2, gain
factor; corrected vs 12_inversion_all_stations).
"""

import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from compy import response_check as rc  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = ARCHIVE / "13_audit/research/inversion_corrected"
CODE = Path(
    os.environ.get("COMPY_STAGE_E_CODE", Path.home() / "Desktop/BRUIT-FM/ComPy")
)  # v2.8.0
STATIONS = ("RR28", "RR29", "RR34", "RR36", "RR38", "RR40", "RR50", "RR52")
DEPTHS_M = (500, 1000, 2000, 3000, 5000)


def factor_k(sid):
    for row in csv.DictReader(open(ARCHIVE / "13_audit/research/response_bracket.csv")):
        if row["sid"] == sid and row.get("k_tide"):
            return float(row["k_tide"]), "tide"
    for row in csv.DictReader(open(ARCHIVE / "13_audit/dpg_tide/dpg_tide.csv")):
        if row["sid"] == sid and row.get("k_lowest"):
            return float(row["k_lowest"]), "phase"
    return None, None


def prepare(sta, plan):
    from obspy import read_inventory

    sid = f"YV.{sta}"
    src = ARCHIVE / "06_compliance_wide/YV/results" / sta
    if plan[sid]["epoch_end"]:
        src = src / "epoch"
    dst = OUT / "YV/results" / sta
    dst.mkdir(parents=True, exist_ok=True)
    k, source = factor_k(sid)
    phase = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/phase/phase_metadata.csv"))
    }
    inv = read_inventory(str(ARCHIVE / "00_inventory/stationxml/YV.xml")).select(
        station=sta
    )
    ch = [c for n in inv for s in n for c in s if c.code == phase[sid]["p_channel"]][0]
    stages = rc.paz_stages(ch)
    summary = json.loads((src / f"{sta}_summary.json").read_text())
    for path in src.glob(f"{sta}_compliance_*.npz"):
        d = dict(np.load(path))
        factor = d["compliance"] / rc.correct_compliance(
            d["compliance"], stages, d["freqs"], k
        )
        for key in ("compliance", "uncertainty", "window_mean", "window_std"):
            d[key] = d[key] / factor
        np.savez(dst / path.name, **d)
        label = path.stem.split("_compliance_")[1]
        band = (d["freqs"] >= d["band"][0]) & (d["freqs"] <= d["band"][1])
        fb = d["freqs"][band]
        for m in summary["results"][label].get("monthly", {}).values():
            if len(m["compliance"]) == len(fb):
                fm = np.asarray(m["compliance"]) / rc.correct_compliance(
                    m["compliance"], stages, fb, k
                )
                m["compliance"] = (np.asarray(m["compliance"]) / fm).tolist()
                m["uncertainty"] = (np.asarray(m["uncertainty"]) / fm).tolist()
    summary["response_correction"] = {"k": k, "k_source": source, "model": "lowest pole scaled",
                                      "script": "docs/audit_2026-11-data/scripts/r06_corrected_inversion.py"}  # fmt: skip
    (dst / f"{sta}_summary.json").write_text(json.dumps(summary, indent=1))
    return k, source


def invert(sta):
    env = dict(os.environ, P_NET="YV", P_ROOT=str(OUT), COMPY_STATION=sta, COMPY_RUN_DIR=str(OUT / "YV"),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MPLBACKEND="Agg",
               COMPY_PIPELINE_CONFIG=str(ARCHIVE / "00_inventory/configs/compliance_override.py"))  # fmt: skip
    log = OUT / f"YV.{sta}.log"
    with open(log, "w") as fh:
        return subprocess.run([sys.executable, "e_invert.py", sta, "final", "--processes", "5"],
                              cwd=CODE / "pipeline", env=env, stdout=fh, stderr=subprocess.STDOUT).returncode  # fmt: skip


def vs_median(folder, sta):
    """Median Vs at DEPTHS_M from the posterior draws, as scripts/all_stations_collect.py does."""
    sys.path.insert(0, str(CODE))
    import compy.inversion.dram as ci

    j = json.loads((folder / f"{sta}_final.json").read_text())
    d = np.load(folder / f"{sta}_final.npz")
    m = j["model"]
    model = ci.CrustalModel(moho_depth=m["moho_depth"], vp_vs_range=tuple(m["vp_vs_range"]),
                            sediment_max=m["sediment_max"], require_moho_jump=True)  # fmt: skip
    names = list(d["names"])
    theta = d["theta"].reshape(-1, len(names))
    rng = np.random.default_rng(0)
    pick = theta[rng.choice(len(theta), min(8000, len(theta)), replace=False)][
        :, : len(model.names)
    ]
    draws = np.array([model.vs_at(t, np.array(DEPTHS_M, float)) for t in pick])
    return np.median(draws, axis=0)


def main():
    plan = {
        r["sid"]: r
        for r in csv.DictReader(
            open(ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv")
        )
    }
    ks = {sta: prepare(sta, plan) for sta in STATIONS}
    with ThreadPoolExecutor(2) as pool:
        codes = dict(zip(STATIONS, pool.map(invert, STATIONS)))
    ref = {
        r["sid"]: r
        for r in csv.DictReader(
            open(ARCHIVE / "12_inversion_all_stations/tables/inversion_results.csv")
        )
    }
    rows = []
    for sta in STATIONS:
        sid = f"YV.{sta}"
        folder = OUT / "YV/results" / sta / "inversion"
        rec_path = folder / f"{sta}_final.json"
        row = {"sid": sid, "k": ks[sta][0], "k_source": ks[sta][1], "exit": codes[sta]}
        if rec_path.exists():
            rec = json.loads(rec_path.read_text())
            par = rec.get("parameters", {})
            row.update(chi2=rec.get("chi2", {}).get("median"), n_data=rec.get("n_data"),
                       gain_factor=par.get("gain", {}).get("median"),
                       sediment_m=par.get("sediment_thickness", {}).get("median"))  # fmt: skip
            vs = vs_median(folder, sta)
            for z, v in zip(DEPTHS_M, vs):
                row[f"vs_{z}"] = float(v)
                row[f"vs_{z}_ref"] = (
                    float(ref[sid][f"vs_{z}"]) if ref[sid].get(f"vs_{z}") else np.nan
                )
                row[f"dvs_{z}_pct"] = 100 * (row[f"vs_{z}"] / row[f"vs_{z}_ref"] - 1)
            row["sediment_ref"] = float(ref[sid]["sediment_m"])
            row["chi2_ref"] = float(ref[sid]["chi2"])
            row["gain_factor_ref"] = float(ref[sid]["gain_factor"])
        rows.append(row)
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "comparison.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    print(json.dumps(rows, indent=0, default=str)[:3000])


if __name__ == "__main__":
    main()
