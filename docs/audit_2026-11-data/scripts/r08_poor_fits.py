"""Why do 7D.J35A and J44A not fit? Band top vs response correction.

The inversion report flags J35A (chi2 246/13) and J44A (179/11). Their residuals
are smooth and largest at the band top (0.82-0.91 f_c, where pressure is strongly
attenuated) and, at J35A, at the band bottom (where the DPG response bias acts).
Stage E (``final``, ComPy 2.8.0, fixed seeds) is rerun in new folders on:

* ``top075``: the same compliance with the band top capped at 0.75 f_c;
* ``top075_corrected``: the same plus the DPG response correction (r06).

Output: ``13_audit/research/poor_fits/<variant>/7D/results/<STA>`` and
``poor_fits/comparison.csv`` (chi2/N per variant vs 12_inversion_all_stations).
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
sys.path.insert(0, str(HERE))
from compy import response_check as rc  # noqa: E402
from r06_corrected_inversion import CODE, factor_k  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = ARCHIVE / "13_audit/research/poor_fits"
STATIONS = ("J35A", "J44A")
TOP = 0.75


def prepare(sta, variant):
    from obspy import read_inventory

    sid = f"7D.{sta}"
    src = ARCHIVE / "06_compliance_wide/7D/results" / sta
    dst = OUT / variant / "7D/results" / sta
    dst.mkdir(parents=True, exist_ok=True)
    summary = json.loads((src / f"{sta}_summary.json").read_text())
    fc = float(summary["f_c_hz"])
    stages = None
    if variant.endswith("corrected"):
        phase = {
            r["sid"]: r
            for r in csv.DictReader(open(ARCHIVE / "13_audit/phase/phase_metadata.csv"))
        }
        inv = read_inventory(str(ARCHIVE / "00_inventory/stationxml/7D.xml")).select(
            station=sta
        )
        ch = [
            c for n in inv for s in n for c in s if c.code == phase[sid]["p_channel"]
        ][0]
        stages = rc.paz_stages(ch)
        k, _ = factor_k(sid)
    for path in src.glob(f"{sta}_compliance_*.npz"):
        d = dict(np.load(path))
        old_band = d["band"].copy()
        top = d["freqs"][d["freqs"] <= TOP * fc].max()
        d["band"] = np.array([old_band[0], min(old_band[1], top)])
        if stages is not None:
            factor = d["compliance"] / rc.correct_compliance(
                d["compliance"], stages, d["freqs"], k
            )
            for key in ("compliance", "uncertainty", "window_mean", "window_std"):
                d[key] = d[key] / factor
        np.savez(dst / path.name, **d)
        label = path.stem.split("_compliance_")[1]
        entry = summary["results"][label]
        keep = (d["freqs"] >= old_band[0]) & (d["freqs"] <= old_band[1])
        new = d["freqs"][keep] <= d["band"][1]
        for m in entry.get("monthly", {}).values():
            if len(m["compliance"]) == keep.sum():
                comp, unc = np.asarray(m["compliance"]), np.asarray(m["uncertainty"])
                if stages is not None:
                    fm = comp / rc.correct_compliance(comp, stages, d["freqs"][keep], k)
                    comp, unc = comp / fm, unc / fm
                m["compliance"], m["uncertainty"] = (
                    comp[new].tolist(),
                    unc[new].tolist(),
                )
        entry["band_hz"] = d["band"].tolist()
    summary["research_variant"] = {"variant": variant, "band_top_fraction_of_fc": TOP}
    (dst / f"{sta}_summary.json").write_text(json.dumps(summary, indent=1))


def invert(args):
    sta, variant = args
    root = OUT / variant
    env = dict(os.environ, P_NET="7D", P_ROOT=str(root), COMPY_STATION=sta, COMPY_RUN_DIR=str(root / "7D"),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MPLBACKEND="Agg",
               COMPY_PIPELINE_CONFIG=str(ARCHIVE / "00_inventory/configs/compliance_override.py"))  # fmt: skip
    with open(root / f"7D.{sta}.log", "w") as fh:
        return subprocess.run([sys.executable, "e_invert.py", sta, "final", "--processes", "5"],
                              cwd=CODE / "pipeline", env=env, stdout=fh, stderr=subprocess.STDOUT).returncode  # fmt: skip


def main():
    tasks = [(s, v) for v in ("top075", "top075_corrected") for s in STATIONS]
    for sta, variant in tasks:
        prepare(sta, variant)
    with ThreadPoolExecutor(2) as pool:
        codes = list(pool.map(invert, tasks))
    ref = {
        r["sid"]: r
        for r in csv.DictReader(
            open(ARCHIVE / "12_inversion_all_stations/tables/inversion_results.csv")
        )
    }
    rows = []
    for (sta, variant), code in zip(tasks, codes):
        path = OUT / variant / "7D/results" / sta / "inversion" / f"{sta}_final.json"
        row = {"sid": f"7D.{sta}", "variant": variant, "exit": code,
               "chi2_ref": float(ref[f"7D.{sta}"]["chi2"]), "n_ref": int(ref[f"7D.{sta}"]["n_data"])}  # fmt: skip
        if path.exists():
            rec = json.loads(path.read_text())
            row.update(
                chi2=rec["chi2"]["median"], n=rec["n_data"], max_rhat=rec["max_rhat"]
            )
        rows.append(row)
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "comparison.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
