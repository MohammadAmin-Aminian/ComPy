"""Stage H: fill report/template.tex from results/tables/report_numbers.json and compile it.

Writes report/XE_young_orca_compliance_report.tex and .pdf (xelatex, twice),
with the figures copied next to it. Placeholders are <<NAME>>; any left
unfilled stop the build.
"""

import json
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
from obspy import read_inventory

from config import FIGURES, INVENTORY, RESULTS, ROOT, STATIONS, WORK

REPORT = ROOT / "report"
FIGS = ("fig01_map.png", "fig03_psd_stages.png", "fig04_tilt.png", "fig05_calibration.png", "fig06_coherogram.png",
        "fig07_windows.png", "fig08_compliance.png", "fig09_heatmaps.png", "fig10_profiles.png")  # fmt: skip
REF = "CC11"


def f1(x, d=1):
    return f"{x:.{d}f}"


def rng(values, d=2, scale=1.0):
    v = [x * scale for x in values]
    return f"{min(v):.{d}f}--{max(v):.{d}f}"


def qtxt(q, d=2, scale=1e-3):
    return f"{q[1] * scale:.{d}f} ({q[0] * scale:.{d}f}--{q[2] * scale:.{d}f})"


def cut_fraction():
    cut = total = 0.0
    for s in STATIONS:
        for p in (WORK / s / "qc").glob("*.json"):
            j = json.loads(p.read_text())
            cut += j["cut_hours"]
            total += 24 * len(j["days"])
    return cut / total if total else 0.0


def bump_db():
    """Median vertical acceleration PSD of cleaned REF data in 8-12 mHz [dB]."""
    from scipy.signal import welch
    from obspy import read

    vals = []
    for path in sorted((WORK / REF / "clean").glob(f"{REF}_seg*.mseed"))[:8]:
        tr = read(str(path)).select(channel="BHZ")[0]
        x = tr.data.astype(float)
        for i in range(0, len(x) - 4500, 4500 * 6):
            f, p = welch(
                x[i : i + 4500] - x[i : i + 4500].mean(), fs=1.25, nperseg=1024
            )
            sel = (f > 0.008) & (f < 0.012)
            vals.append(10 * np.log10(np.median(p[sel] * (2 * np.pi * f[sel]) ** 4)))
    return float(np.median(vals))


def table1(recs):
    rows = []
    for r in recs:
        if r.get("n_windows", 0) < 500 or "band" not in r:
            continue
        if r.get("gain_compy") and r.get("gain_compy_std"):
            g = f"{r['gain_compy']:.3f} $\\pm$ {r['gain_compy_std']:.3f} ({r['gain_events']})"
        elif r.get("gain_compy"):
            g = f"{r['gain_compy']:.3f} ({r['gain_events']})"
        else:
            g = "--"
        gb = f"{r['gain_bruitfm']:.3f}" if r.get("gain_bruitfm") else "--"
        rows.append(f"{r['station']} & {r['lat']:.4f} & {r['lon']:.4f} & {r['depth']:.0f} & {g} & {gb} & "
                    f"{r['band'][0] * 1e3:.1f}--{r['band'][1] * 1e3:.1f} & {r['threshold']} & {r['n_windows']} & "
                    f"{r['coherence']:.3f} \\\\")  # fmt: skip
    return (
        "\\begin{table*}\n\\centering\\small\n\\caption{Stations with a compliance measurement. Gain: ComPy Rayleigh-wave "
        "calibration (event standard deviation, number of events); BRUIT-FM: independent microseism calibration where "
        "available. Frequency: compliance band. Threshold: coherence gate. Coherence: band median of the stack.}\n"
        "\\label{tab:stations}\n\\begin{tabular}{lrrrlrrrrr}\n\\toprule\n"
        "Station & Lat. [$^\\circ$] & Lon. [$^\\circ$] & Depth [m] & Gain (ComPy) & Gain (BRUIT-FM) & Band [mHz] & "
        "Thr. & Windows & Coh. \\\\\n\\midrule\n"
        + "\n".join(rows)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}"
    )


def table2(recs, depths):
    idx = [depths.index(z) for z in (500.0, 1000.0, 2000.0, 4000.0)]
    rows = []
    for r in recs:
        if not r["status"].startswith("inverted"):
            continue
        vs = np.array(r["vs"])
        cells = " & ".join(
            f"{vs[1, i] / 1e3:.2f} ({vs[0, i] / 1e3:.2f}--{vs[2, i] / 1e3:.2f})"
            for i in idx
        )
        vm = (
            f"{np.array(r['vs_main'])[1, idx[-1]] / 1e3:.2f}"
            if r.get("vs_main")
            else "--"
        )
        mark = "$^*$" if "uncalibrated" in r["status"] else ""
        rows.append(f"{r['station']}{mark} & {r['chi2']:.1f}/{r['n_data']} & {r['rhat']:.3f} & "
                    f"{r['sediment'][1]:.0f} & {r['vp_vs'][1]:.2f} & {cells} & {vm} \\\\")  # fmt: skip
    return (
        "\\begin{table*}\n\\centering\\small\n\\caption{Inversion results: fit ($\\chi^2/N$), largest $\\hat R$, "
        "sediment thickness and crustal Vp/Vs (posterior medians), Vs [km\\,s$^{-1}$] as median (5--95 per cent) at "
        "depths below the seafloor, and Vs at 4\\,km with Vp/Vs = 1.85. $^*$Gain not calibrated (sampled over 0.6--1.4).}\n"
        "\\label{tab:vs}\n\\resizebox{\\textwidth}{!}{\\begin{tabular}{lrrrrllllr}\n\\toprule\n"
        "Station & $\\chi^2/N$ & $\\hat R$ & Sed. [m] & Vp/Vs & Vs 0.5\\,km & Vs 1\\,km & Vs 2\\,km & Vs 4\\,km & "
        "Vs 4\\,km (1.85) \\\\\n\\midrule\n"
        + "\n".join(rows)
        + "\n\\bottomrule\n\\end{tabular}}\n\\end{table*}"
    )


def main():
    nums = json.loads((RESULTS / "tables" / "report_numbers.json").read_text())
    depths = nums["depths"]
    recs = nums["stations"]
    di = {z: depths.index(z) for z in depths}
    inv = read_inventory(str(INVENTORY))
    all_depths = [-s.elevation for s in inv[0]]
    ok = [r for r in recs if r["status"].startswith("inverted")]
    cal = [r for r in ok if r["status"] == "inverted"]
    measured = [r for r in recs if r.get("n_windows", 0) >= 500 and "band" in r]
    with_gain = [r for r in measured if r.get("gain_compy")]
    spread = [
        100 * r["gain_compy_std"] / r["gain_compy"]
        for r in with_gain
        if r.get("gain_compy_std")
    ]
    pairs = [
        (r["station"], r["gain_compy"], r["gain_bruitfm"])
        for r in with_gain
        if r.get("gain_bruitfm")
    ]
    ratios = [a / b for _, a, b in pairs]
    agree = [x for x in ratios if 0.8 <= x <= 1.25]
    outliers = [
        f"{n} (ComPy {a:.2f}, BRUIT-FM {b:.2f})"
        for n, a, b in pairs
        if not 0.8 <= a / b <= 1.25
    ]
    pooled = nums["pooled"]
    pv, pm, ps, c1 = (
        np.array(pooled[k]) for k in ("vs", "vs_main", "serp", "crust1_vs")
    )
    uncal = [r["station"] for r in ok if "uncalibrated" in r["status"]]
    nocal = [
        r["station"]
        for r in recs
        if r.get("calibration_error") and r.get("n_windows", 0) >= 500
    ]
    cc07 = next(
        (r for r in recs if r["station"] == "CC07" and r.get("monthly_dev")), None
    )
    others = [
        max(abs(x) for x in r["monthly_dev"])
        for r in cal
        if r.get("monthly_dev") and r["station"] != "CC07"
    ]
    vs2000 = [np.array(r["vs"])[1, di[2000.0]] / 1e3 for r in cal]
    xe2000 = pv[1, di[2000.0]] / 1e3
    slower = [
        name
        for name, v in (("RR36", 3.09), ("RR40", 3.56), ("RR52", 3.59))
        if v < min(vs2000)
    ]
    ridge = (f"Only {', '.join(slower)} {'is' if len(slower) == 1 else 'are'} slower than every Young ORCA station "
             f"({min(vs2000):.2f}--{max(vs2000):.2f}\\,km\\,s$^{{-1}}$ at 2\\,km); the others fall inside the normal-crust range."
             if slower else
             f"None of them is slower than the Young ORCA range ({min(vs2000):.2f}--{max(vs2000):.2f}\\,km\\,s$^{{-1}}$ at 2\\,km), "
             "so at that depth the ridge stations are not distinguishable from normal crust by Vs alone.")  # fmt: skip
    fill = {
        "DATE": time.strftime("%d %B %Y"),
        "COMPY_VERSION": "2.4.0",
        "NSTA_TOTAL": str(len(inv[0])),
        "DEPTH_MIN": f"{min(all_depths):.0f}",
        "DEPTH_MAX": f"{max(all_depths):.0f}",
        "GAIN_RANGE": rng([r["gain_compy"] for r in with_gain]),
        "GAIN_SPREAD": f"{np.median(spread):.1f}",
        "GAIN_SPREAD_RANGE": rng(spread, 1),
        "BAND_LO": f"{np.median([r['band'][0] for r in cal]) * 1e3:.1f}",
        "BAND_HI": f"{np.median([r['band'][1] for r in cal]) * 1e3:.1f}",
        "N_OK": str(len(ok)),
        "N_MEASURED": str(len(measured)),
        "WIN_MEDIAN": f"{np.median([r['n_windows'] for r in ok]):.0f}",
        "CHI2N": f"{np.median([r['chi2'] / r['n_data'] for r in ok]):.2f}",
        "SED_RANGE": f"{min(r['sediment'][1] for r in cal):.0f}--{max(r['sediment'][1] for r in cal):.0f}",
        "C1SED_RANGE": f"{min(r['crust1_sediment'] for r in cal):.0f}--{max(r['crust1_sediment'] for r in cal):.0f}",
        "VS300": f1(pv[1, di[300.0]] / 1e3, 2),
        "VS1000": f1(pv[1, di[1000.0]] / 1e3, 2),
        "VS2000": f1(xe2000, 2),
        "VS2000_MED": f1(xe2000, 2),
        "VS4000": f1(pv[1, di[4000.0]] / 1e3, 2),
        "VS4000_MAIN": f1(pm[1, di[4000.0]] / 1e3, 2),
        "C1_4000": f1(c1[di[4000.0]] / 1e3, 2),
        "SERP500": qtxt(ps[:, di[500.0]], 0, 1),
        "SERP1000": qtxt(ps[:, di[1000.0]], 0, 1),
        "SERP2000": qtxt(ps[:, di[2000.0]], 0, 1),
        "SERP1000_MED": f"{ps[1, di[1000.0]]:.0f}",
        "N_PROCESSED": str(len(recs)),
        "N_NOSEG": str(sum("no usable segments" in r["status"] for r in recs)),
        "NOSEG_LIST": ", ".join(
            r["station"] for r in recs if "no usable segments" in r["status"]
        ),
        "GATES_LIST": ", ".join(
            r["station"] for r in recs if r["status"] == "no compliance (gates)"
        )
        or "none",
        "N_REJECTED": str(sum(r["status"].startswith("rejected") for r in recs)),
        "REJ_LIST": ", ".join(
            f"{r['station']} ({r['n_windows']})"
            for r in recs
            if r["status"].startswith("rejected")
        ),
        "CUT_FRACTION": f"{100 * cut_fraction():.0f}\\,per\\,cent",
        "REF_STA": REF,
        "BUMP_DB": f"{bump_db():.0f}",
        "N_CAL": str(len(with_gain)),
        "CAL_EVENTS": f"{min(r['gain_events'] for r in with_gain)}--{max(r['gain_events'] for r in with_gain)}",
        "GAIN_RATIO": (
            f"{np.median(agree):.3f} ({min(agree):.3f}--{max(agree):.3f}) at {len(agree)} of {len(ratios)} "
            "stations"
            + (
                f"; they disagree at {', '.join(outliers)}, where one of the two methods failed"
                if outliers
                else ""
            )
        ),  # fmt: skip
        "UNCAL_TEXT": (
            f"No event passed the coherence test at {' and '.join(nocal)}; "
            f"{'these stations were' if len(uncal) > 1 else 'this station was'} kept with the nominal response "
            "and a gain sampled over 0.6--1.4, and "
            f"{'are' if len(uncal) > 1 else 'is'} excluded from the pooled results."
            if nocal
            else ""
        ),  # fmt: skip
        "N_CAL_INV": str(len(cal)),
        "COH_RANGE": rng([r["coherence"] for r in cal], 3),
        "REL_ERR": f"{100 * np.median([r['rel_error'] for r in cal]):.1f}",
        "ETA_LO_RANGE": rng([r["eta"]["eta"][0] * 1e9 for r in cal], 3),
        "MONTHLY_TYPICAL": f"{100 * np.median(others):.0f}--{100 * np.percentile(others, 90):.0f}"
        if others
        else "--",
        "CC07_DEV": f"{100 * cc07['monthly_dev'][0]:+.0f} to {100 * cc07['monthly_dev'][1]:+.0f}"
        if cc07
        else "--",
        "RHAT_MAX": f"{max(r['rhat'] for r in ok):.3f}",
        "VS500_RANGE": rng([np.array(r["vs"])[1, di[500.0]] for r in cal], 2, 1e-3),
        "VS1000_RANGE": rng([np.array(r["vs"])[1, di[1000.0]] for r in cal], 2, 1e-3),
        "VS2000_RANGE": rng(vs2000, 2),
        "WR1000": f"{np.median([r['width_ratio'][di[1000.0]] for r in cal if r.get('width_ratio')]):.2f}",
        "WR2000": f"{np.median([r['width_ratio'][di[2000.0]] for r in cal if r.get('width_ratio')]):.2f}",
        "SENS50": f"{np.median([r['sens50'] for r in cal]) / 1e3:.1f}",
        "SENS90": f"{np.median([r['sens90'] for r in cal]) / 1e3:.1f}",
        "VPVS_RANGE": rng([r["vp_vs"][1] for r in cal]),
        "RIDGE_COMPARISON": ridge,
        "TABLE1": table1(recs),
        "TABLE2": table2(recs, depths),
    }
    text = (REPORT / "template.tex").read_text()
    for key, value in fill.items():
        text = text.replace(f"<<{key}>>", value)
    left = sorted(set(re.findall(r"<<([A-Z0-9_]+)>>", text)))
    if left:
        raise SystemExit(f"unfilled placeholders: {left}")
    for name in FIGS:
        shutil.copy(FIGURES / name, REPORT / name)
    tex = REPORT / "XE_young_orca_compliance_report.tex"
    tex.write_text(text)
    for _ in range(2):
        subprocess.run(
            ["xelatex", "-interaction=nonstopmode", tex.name],
            cwd=REPORT,
            capture_output=True,
        )
    pdf = tex.with_suffix(".pdf")
    print(
        f"{pdf} ({pdf.stat().st_size / 1e6:.1f} MB)"
        if pdf.exists()
        else "PDF NOT BUILT"
    )
    (REPORT / "fill.json").write_text(
        json.dumps(
            {k: v for k, v in fill.items() if not k.startswith("TABLE")}, indent=1
        )
    )


if __name__ == "__main__":
    main()
