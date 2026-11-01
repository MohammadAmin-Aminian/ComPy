"""Stage I: web version of the report (summary + the PDF's pages as images).

Usage: python i_html.py OUT_DIR
Writes OUT_DIR/index.html and OUT_DIR/pages/page-N.png from
report/XE_young_orca_compliance_report.pdf and report/fill.json.
"""

import html
import json
import subprocess
import sys
from pathlib import Path

from config import RESULTS, ROOT

PDF = ROOT / "report" / "XE_young_orca_compliance_report.pdf"

CSS = """
:root { --bg:#f4f6f7; --paper:#ffffff; --fg:#1b2329; --muted:#5a6772; --rule:#d6dde2; --accent:#0d5c76; --accent-soft:#e1eef2;
  --serif:"Source Serif 4", Georgia, serif; --sans:"IBM Plex Sans", system-ui, sans-serif; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#111619; --paper:#182026; --fg:#e2e8ec;
  --muted:#9aa8b2; --rule:#2b3741; --accent:#6dc4dc; --accent-soft:#1c3440; color-scheme:dark; } }
:root[data-theme="dark"] { --bg:#111619; --paper:#182026; --fg:#e2e8ec; --muted:#9aa8b2; --rule:#2b3741; --accent:#6dc4dc;
  --accent-soft:#1c3440; color-scheme:dark; }
body { background:var(--bg); color:var(--fg); font-family:var(--serif); font-size:17px; line-height:1.6; }
.wrap { max-width:980px; margin:0 auto; padding-block:36px 72px; padding-left:max(16px, env(safe-area-inset-left, 0px));
  padding-right:max(16px, env(safe-area-inset-right, 0px)); }
.kicker { font-family:var(--sans); font-size:12px; letter-spacing:.12em; text-transform:uppercase; color:var(--accent); font-weight:600; }
h1 { font-size:clamp(26px,3.6vw,36px); line-height:1.2; margin:8px 0; text-wrap:balance; font-weight:600; }
h2 { font-size:22px; margin:36px 0 10px; font-weight:600; }
.meta { font-family:var(--sans); font-size:14px; color:var(--muted); }
p, li { max-width:70ch; }
.facts { display:grid; grid-template-columns:repeat(auto-fit, minmax(170px, 1fr)); gap:10px; margin:18px 0; font-family:var(--sans); }
.fact { background:var(--paper); border:1px solid var(--rule); padding:10px 12px; min-width:0; }
.fact b { display:block; font-size:20px; color:var(--accent); font-variant-numeric:tabular-nums; }
.fact span { font-size:13px; color:var(--muted); }
.abstract { background:var(--paper); border:1px solid var(--rule); border-left:4px solid var(--accent); padding:16px 20px; }
.pages { display:flex; flex-direction:column; gap:18px; margin-top:12px; }
.pages img { width:100%; height:auto; background:#fff; border:1px solid var(--rule); box-shadow:0 1px 3px rgba(0,0,0,.08); }
.note { font-family:var(--sans); font-size:14px; color:var(--muted); }
code { font-family:"IBM Plex Mono", ui-monospace, monospace; font-size:.85em; background:var(--accent-soft); padding:1px 5px; }
"""


def main():
    out = Path(sys.argv[1])
    (out / "pages").mkdir(parents=True, exist_ok=True)
    for old in (out / "pages").glob("page-*.png"):
        old.unlink()
    subprocess.run(
        ["pdftoppm", "-png", "-r", "130", str(PDF), str(out / "pages" / "page")],
        check=True,
    )
    pages = sorted(
        (out / "pages").glob("page-*.png"), key=lambda p: int(p.stem.split("-")[-1])
    )
    fill = json.loads((ROOT / "report" / "fill.json").read_text())
    nums = json.loads((RESULTS / "tables" / "report_numbers.json").read_text())
    facts = [
        (fill["N_OK"], f"stations inverted (of {fill['NSTA_TOTAL']})"),
        (
            f"{fill['BAND_LO']}–{fill['BAND_HI']} mHz",
            "compliance band, coherence > 0.9",
        ),
        (fill["GAIN_RANGE"].replace("--", "–"), "DPG gain (Rayleigh, M ≥ 7)"),
        (f"{fill['VS1000']} km/s", "Vs at 1 km below seafloor (pooled)"),
        (
            f"{fill['VS4000']}–{fill['VS4000_MAIN']} km/s",
            "Vs at 4 km (Vp/Vs free – 1.85)",
        ),
        (f"{fill['SERP1000_MED']} %", "apparent serpentinization at 1 km"),
    ]
    tex = (ROOT / "report" / "XE_young_orca_compliance_report.tex").read_text()
    abstract = tex.split("\\textbf{SUMMARY}\\\\")[1].split("\\\\[4pt]")[0]
    for a, b in (("\\,", " "), ("$\\sim$", "~"), ("$\\pm$", "±"), ("$M\\geq7$", "M ≥ 7"), ("$\\chi^2/N$", "χ²/N"),
                 ("km\\,s$^{-1}$", "km/s"), ("km s$^{-1}$", "km/s"), ("--", "–"), ("\\ ", " "), ("per cent", "%"),
                 ("\n", " ")):  # fmt: skip
        abstract = abstract.replace(a, b)
    rows = "".join(
        f'<div class="fact"><b>{html.escape(str(v))}</b><span>{html.escape(t)}</span></div>'
        for v, t in facts
    )
    imgs = "".join(
        f'<img src="pages/{p.name}" alt="Report page {k + 1}" loading="lazy">'
        for k, p in enumerate(pages)
    )
    page = f"""<title>Young ORCA Compliance</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">ComPy 2.4 · second network · draft report</div>
  <h1>Shallow crustal shear velocity of ~40 Ma Pacific lithosphere from seafloor compliance at the Young ORCA array</h1>
  <div class="meta">For M.-A. Aminian · {html.escape(fill["DATE"])} · network XE (Pacific Array Young ORCA, 2018–2019) ·
   {len(nums["inverted"])} stations inverted</div>
  <div class="facts">{rows}</div>
  <h2>Summary</h2>
  <div class="abstract"><p>{html.escape(abstract.strip())}</p></div>
  <p class="note">Below is the full report in the layout of Aminian et al. (2025, GJI). The PDF and LaTeX source are in
  <code>BRUIT-FM/compliance_xe_2026-10-29/report/</code>, the results in <code>results/</code> and the figures in
  <code>figures/</code>.</p>
  <h2>Report</h2>
  <div class="pages">{imgs}</div>
</div>
"""
    (out / "index.html").write_text(page)
    print(out / "index.html", len(pages), "pages")


if __name__ == "__main__":
    main()
