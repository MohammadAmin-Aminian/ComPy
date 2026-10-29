"""Stage E: HTML report (report/index.html, figures copied to report/img)."""

import json
import shutil
from html import escape

from a_invert import EPOCH, ROOT
from c_figures import STATIONS, VARIANT_LABELS

REPORT = ROOT / "report"


def load_tables():
    tables = ROOT / "results" / "tables"
    main = {r["station"]: r for r in json.loads((tables / "inversion_summary.json").read_text())}
    variants = json.loads((tables / "variants.json").read_text())
    return main, variants


def main_table(main):
    head = ("<tr><th>Station</th><th>Record</th><th>Band [mHz]</th><th>χ² / N</th><th>R-hat</th>"
            "<th>ESS</th><th>Sediment [m]</th><th>Vs 0.5 km</th><th>Vs 1 km</th><th>Vs 2 km</th>"
            "<th>50 % / 90 % sensitivity [km]</th></tr>")  # fmt: skip
    rows = []
    for s in STATIONS:
        r = main[s]
        bad = r["chi2_median"] > 2 * r["n_data"]
        rows.append(
            f"<tr{' class=\"flag\"' if bad else ''}><td>{s}</td><td>{'stable epoch' if r['epoch'] else 'whole'}</td>"
            f"<td>{r['band_mhz']}</td><td>{r['chi2_median']:.1f} / {r['n_data']}</td>"
            f"<td>{r['max_rhat']:.3f}</td><td>{r['min_ess']}</td><td>{r['sediment_m']}</td>"
            f"<td>{r['vs_0.5km']}</td><td>{r['vs_1km']}</td><td>{r['vs_2km']}</td>"
            f"<td>{r['sens50_km']:.1f} / {r['sens90_km']:.1f}</td></tr>"
        )  # fmt: skip
    return f"<table>{head}{''.join(rows)}</table>"


def variant_table(variants):
    head = ("<tr><th>Station</th><th>Variant</th><th>χ² / N</th><th>R-hat</th><th>Sediment [m]</th>"
            "<th>Vs 0.5 km</th><th>Vs 1 km</th><th>Vs 2 km</th><th>Vs 3 km</th></tr>")  # fmt: skip
    rows = []
    for r in variants:
        label = VARIANT_LABELS[r["variant"]]
        if r["vp_vs"]:
            label += f" ({r['vp_vs']:.2f}, 90 %: {r['vp_vs_90'][0]:.2f}–{r['vp_vs_90'][1]:.2f})"
        cls = [] if r["converged"] else ["warn"]
        if r["variant"] == "main":
            cls.append("lead")
        vs = r["vs"]
        rows.append(
            f"<tr class=\"{' '.join(cls)}\"><td>{r['station']}</td><td>{escape(label)}</td>"
            f"<td>{r['chi2']:.1f} / {r['n']}</td><td>{r['max_rhat']:.2f}</td><td>{r['sediment']:.0f}</td>"
            f"<td>{vs['0.5']:.0f}</td><td>{vs['1.0']:.0f}</td><td>{vs['2.0']:.0f}</td><td>{vs['3.0']:.0f}</td></tr>"
        )  # fmt: skip
    return f"<table class=\"compact\">{head}{''.join(rows)}</table>"


def figure(name, caption):
    return (f'<figure><div class="plate"><img src="img/{name}" alt="{escape(caption)}" loading="lazy"></div>'
            f"<figcaption>{caption}</figcaption></figure>")  # fmt: skip


def main():
    (REPORT / "img").mkdir(parents=True, exist_ok=True)
    for png in (ROOT / "figures").glob("*.png"):
        shutil.copy(png, REPORT / "img" / png.name)
    main_rows, variants = load_tables()
    vpvs = {r["station"]: r for r in variants if r["variant"] == "vpvs"}
    fixed = {r["station"]: r for r in variants if r["variant"] == "fixedvp"}
    body = TEMPLATE.format(
        main_table=main_table(main_rows),
        variant_table=variant_table(variants),
        rr38_vpvs=(f"{vpvs['RR38']['vp_vs']:.2f}, 90 %: {vpvs['RR38']['vp_vs_90'][0]:.2f}–"
                   f"{vpvs['RR38']['vp_vs_90'][1]:.2f}") if "RR38" in vpvs else "n/a",
        rr38_vpvs_chi2=f"{vpvs['RR38']['chi2']:.1f}" if "RR38" in vpvs else "n/a",
        vpvs_list=", ".join(f"{s} {vpvs[s]['vp_vs']:.2f}" for s in STATIONS if s in vpvs and s != "RR38"),
        fixed_bad=", ".join(f"{s} (χ² {fixed[s]['chi2']:.0f}/{fixed[s]['n']}, R-hat {fixed[s]['max_rhat']:.2f})"
                            for s in STATIONS if s in fixed and not fixed[s]["converged"]),  # fmt: skip
        figs_main=figure("vs_profiles.png", "Vs posteriors of the main run with CRUST1.0 and the profiles "
                         "described in Aminian et al. (2025, §5.2, approximate). Grey: below the depth holding "
                         "90 % of the lowest-frequency sensitivity.")
        + figure("compliance_fit.png", "Measured compliance (±1σ including the 5 % gain error) and the "
                 "posterior predictive median and 90 % band."),
        figs_sens=figure("sensitivity.png", "Posterior median Vs for each modelling choice.")
        + figure("sediment.png", "Thickness of the slow uppermost layer (Hamilton sediment profile)."),
        figs_station="".join(figure(f"{s}_inversion.png", f"{s}{' (stable epoch)' if s in EPOCH else ''}: "
                                    "(a) data fit, (b) Vs posterior density with prior, CRUST1.0, GJI 2025 "
                                    "and Occam, (c) cumulative Vs sensitivity, (d) sediment thickness.")
                             for s in STATIONS),  # fmt: skip
    )
    (REPORT / "index.html").write_text(body)
    print(REPORT / "index.html")


TEMPLATE = """<title>YV Vs Inversion</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400&display=swap">
<style>
/* Layout: one reading column for text; tables and figures widen to the page. */
:root {{
  --bg: #f6f7f8; --surface: #ffffff; --ink: #14181c; --muted: #555d66; --rule: #d6dbe0;
  --accent: #1f5f99; --warn-bg: #fbeee4; --warn-ink: #8a3b12; --lead-bg: #eaf1f8; --plate: #ffffff;
  --serif: "Source Serif 4", Georgia, serif; --sans: "IBM Plex Sans", system-ui, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #12161a; --surface: #1a2026; --ink: #e8ecef; --muted: #a3adb6; --rule: #2e373f;
  --accent: #7fb3e6; --warn-bg: #3a2418; --warn-ink: #f2b48c; --lead-bg: #1c2a38; --plate: #f4f5f6;
  color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #12161a; --surface: #1a2026; --ink: #e8ecef; --muted: #a3adb6; --rule: #2e373f;
  --accent: #7fb3e6; --warn-bg: #3a2418; --warn-ink: #f2b48c; --lead-bg: #1c2a38; --plate: #f4f5f6;
  color-scheme: dark; }}
body {{ background: var(--bg); color: var(--ink); font: 16px/1.6 var(--sans);
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px)); }}
main {{ max-width: 1120px; margin: 0 auto; padding-block: 40px 64px; display: grid; gap: 28px; }}
.prose {{ max-width: 68ch; }}
h1, h2, h3 {{ font-family: var(--serif); font-weight: 600; text-wrap: balance; line-height: 1.2; margin: 0; }}
h1 {{ font-size: 2.1rem; }} h2 {{ font-size: 1.45rem; margin-top: 12px; }} h3 {{ font-size: 1.1rem; }}
p, ul, ol {{ margin: 0; }} ul, ol {{ padding-left: 1.2em; display: grid; gap: 6px; }}
.meta {{ color: var(--muted); font-size: .92rem; }}
.eyebrow {{ text-transform: uppercase; letter-spacing: .08em; font-size: .78rem; color: var(--accent); font-weight: 600; }}
section {{ display: grid; gap: 14px; }}
code, .mono {{ font-family: var(--mono); font-size: .88em; }}
.box {{ background: var(--surface); border: 1px solid var(--rule); border-radius: 6px; padding: 18px 20px; display: grid; gap: 10px; }}
.wrap {{ overflow-x: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: .9rem; font-variant-numeric: tabular-nums; background: var(--surface); }}
table.compact {{ font-size: .82rem; }}
th, td {{ text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--rule); vertical-align: top; }}
th {{ font-weight: 600; color: var(--muted); white-space: nowrap; }}
tr.flag td, tr.warn td {{ background: var(--warn-bg); color: var(--warn-ink); }}
tr.lead td {{ background: var(--lead-bg); }}
figure {{ margin: 0; display: grid; gap: 8px; }}
.plate {{ background: var(--plate); border-radius: 4px; padding: 8px; overflow-x: auto; }}
.plate img {{ display: block; width: 100%; min-width: 640px; height: auto; }}
figcaption {{ color: var(--muted); font-size: .9rem; max-width: 90ch; }}
.grid2 {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px; }}
.grid2 > * {{ min-width: 0; }}
a {{ color: var(--accent); }}
a:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
</style>
<main>
<header class="prose" style="display:grid;gap:10px">
  <p class="eyebrow">ComPy 2.3 · RHUM-RUM network YV · 29 October 2026</p>
  <h1>Shear velocity beneath the YV broadband stations</h1>
  <p class="meta">Bayesian inversion of the 2.2 compliance results with the new <code>compy_inversion</code> module.
  Branch <code>claude/inversion-2026-10-29</code>. All 121 tests pass.</p>
</header>

<section class="prose">
  <h2>Summary</h2>
  <ul>
    <li>The earlier inversion code had errors that shaped the published Vs models:
      <ul>
        <li>the v1 roughness term did nothing;</li>
        <li>the likelihood used the residual norm instead of its square;</li>
        <li>Vp and density followed Vs, contrary to the text;</li>
        <li>proposals and thickness moves broke detailed balance;</li>
        <li>the errors were single-window standard deviations, tens of times too large for a stacked mean;</li>
        <li>the compliance lacked the wave-attraction term.</li>
      </ul>
      The v2 default <code>alpha=0.25</code> is a 4 m/s Vs-jump prior.</li>
    <li>The new module follows Zha &amp; Webb (2016):
      <ul>
        <li>a Hamilton sediment of free thickness over crustal layers that double in thickness to the CRUST1.0 Moho;</li>
        <li>a covariance with the monthly scatter, a 5 % gain error shared by all frequencies and a 1 % floor;</li>
        <li>a DRAM sampler with a sampled smoothing weight, plus parallel tempering when needed.</li>
      </ul></li>
    <li>All eight main runs converge (R-hat ≤ 1.03) and seven fit to χ² below N. RR38 does not fit with
      Vp/Vs 1.85.</li>
    <li>Compliance hardly constrains Vp/Vs. With one free crustal Vp/Vs in [1.7, 2.4], the 90 % intervals
      span most of the prior at seven stations (medians {vpvs_list}). RR38 alone pushes it to the upper bound
      ({rr38_vpvs}) and still fits worse than the others (χ² {rr38_vpvs_chi2} for 6). Very high Vp/Vs, as in
      serpentinized peridotite, or a low bias in its low-frequency compliance could explain this.</li>
    <li>Vp fixed at CRUST1.0 (the paper's assumption) makes the posterior multimodal. It fits poorly or fails
      to converge at: {fixed_bad}.</li>
    <li>The data sense the top kilometres. Half the lowest-frequency sensitivity lies above 0.2–2.1 km.
      Deeper layers mostly reflect the prior and the layer geometry.</li>
  </ul>
</section>

<section>
  <h2 class="prose">Main results</h2>
  <p class="prose">Main run: Vp/Vs 1.85, gain from ComPy calibration, stable epochs where the pressure/vertical
  sensitivity changed. Vs in m/s, as median (5–95 %). Highlighted row: the fit is poor (χ² above 2N).</p>
  <div class="wrap">{main_table}</div>
  {figs_main}
</section>

<section>
  <h2 class="prose">What changes Vs</h2>
  <div class="prose box">
    <p>Each variant changes one choice relative to the main run:</p>
    <ul>
      <li>Vp from Brocher (2005);</li>
      <li>a free sediment Vs scale;</li>
      <li>the gains.csv pressure gain;</li>
      <li>the whole record instead of the stable epoch;</li>
      <li>a ±1 km Moho;</li>
      <li>Vp and density fixed at CRUST1.0 (tempered);</li>
      <li>one free crustal Vp/Vs in [1.7, 2.4] (tempered).</li>
    </ul>
    <p>The gain and record choices barely matter. The Vp assumption and the sediment model matter most.</p>
  </div>
  {figs_sens}
  <div class="wrap">{variant_table}</div>
  <p class="meta prose">Highlighted rows did not converge (R-hat ≥ 1.05 or ESS ≤ 100); their numbers describe
  the mode most chains reached.</p>
</section>

<section class="prose">
  <h2>Audit of the earlier inversion</h2>
  <ol>
    <li><b>v1 roughness was a no-op.</b> <code>sqrt(sum(np.gradient(vs, 2)))</code> uses 2 as the grid spacing,
      and the sum telescopes to √((V<sub>bottom</sub>−V<sub>top</sub>)/2). Internal structure never changed it.</li>
    <li><b>v1 likelihood</b> was <code>exp(-0.5‖r‖)</code>, using the norm and not its square.</li>
    <li><b>v1 coupled Vp and ρ to Vs</b> (Poisson 0.25, Gardner), while the paper says they were fixed.</li>
    <li><b>v1 proposals broke detailed balance.</b> Out-of-bounds Vs moves were redrawn or reflected. Thickness
      moves added the same step to two layers, so the total thickness was not conserved.</li>
    <li><b>v1 bounds</b> were 0.15–1.2 × the start (0.9–1.1 for the deepest layer), not those stated in the paper.</li>
    <li><b>v1 errors</b> were the standard deviation of single-window estimates (about 30–60 % of η). That is
      roughly √n too large for a mean of several thousand windows.</li>
    <li><b>Pre-2.0 compliance</b> multiplied the wave-attraction term by ω²/(ω²+c), effectively dropping it. It
      also used a free-air constant of <code>3.07*10e-6</code> = 3.07e-5. Both lowered η and made deep Vs too high.</li>
    <li><b>The v2 default alpha = 0.25</b> on a 1-m grid in m/s is a Gaussian prior on Vs jumps with σ = 4 m/s.</li>
    <li><b>Paper errata.</b> Eq. 13 has the acceptance inequality reversed, the misfit (1/N)√Σ is non-standard,
      and the 70–90 % acceptance target misquotes Roberts et al. (1997), whose optimum is 23 %.</li>
  </ol>
  <p class="meta">Details: <code>docs/INVERSION.md</code> in the branch. The legacy functions are kept unchanged.</p>
</section>

<section class="prose">
  <h2>Method</h2>
  <div class="box">
    <h3>Data</h3>
    <p>Gravity-corrected compliance from 2.2, restricted to the coherent band. The covariance has four parts:</p>
    <ul>
      <li>the Bendat &amp; Piersol error (about 0.2 %);</li>
      <li>the monthly-stack covariance divided by the number of months, shrunk halfway to its diagonal. Months
        vary by 2–3 % (5–13 % at RR38), with correlation 0.8–0.9 across frequency except at RR34 (0.1);</li>
      <li>a 5 % gain error along ∂η/∂ln g, which excludes the gain-independent wave attraction;</li>
      <li>a 1 % floor.</li>
    </ul>
    <h3>Model</h3>
    <p>A Hamilton (1976, 1979) silt-clay Vs(z) and Vp(z) sediment of 0–1000 m thickness, in 8 sublayers. Below
    it, crustal layers of 100, 200, 400, 800 and 1600 m and one more down to the CRUST1.0 Moho, then a mantle
    half-space. Vs is uniform in 500–4700 m/s; Vp = 1.85 Vs and ρ from Christensen &amp; Shaw (1970). A
    second-difference roughness prior on crustal Vs has its scale sampled between 20 and 3000 m/s.</p>
    <h3>Sampler</h3>
    <p>Eight DRAM chains of 80 000 iterations, 30 000 of them burn-in, starting from prior draws. The adaptation
    targets 23 % acceptance and is frozen after burn-in. A prior-only run gives the reference widths. The
    fixed-Vp and free-Vp/Vs runs use 8-temperature parallel tempering.</p>
    <h3>Checks</h3>
    <ul>
      <li>The forward model equals <code>raydep</code> to round-off and the static half-space limit.</li>
      <li>DRAM, the hierarchical weight and tempering reproduce analytic targets.</li>
      <li>A synthetic RR52-like model is recovered.</li>
      <li>Draws are reproducible bit for bit.</li>
    </ul>
  </div>
</section>

<section>
  <h2 class="prose">Station figures</h2>
  {figs_station}
</section>

<section class="prose">
  <h2>Decisions for you</h2>
  <ol>
    <li><b>Vp assumption.</b> Vp/Vs 1.85 (Zha &amp; Webb), Brocher, a free Vp/Vs, or fixed CRUST1.0 Vp (paper).
      The data do not choose a Vp/Vs, so the free-Vp/Vs run is the most honest: its Vs intervals include
      that uncertainty. Fixed CRUST1.0 Vp fits RR50 and RR52 poorly.</li>
    <li><b>RR38.</b> Interpret a high Vp/Vs (serpentinization), or treat its low-frequency compliance as biased
      and exclude the station.</li>
    <li><b>Sediment Vs profile.</b> Hamilton is used here; Zha &amp; Webb used Ruan et al. (2014), which can be
      passed in if you have the relation.</li>
    <li><b>Carried over from the compliance run:</b> the RR32 → RR34 substitution, the stable epochs and the
      gain choice. None of these changed the inversion much.</li>
  </ol>
</section>

<footer class="meta prose">
  Code: <code>compy_inversion.py</code>, <code>tests/test_compy_inversion.py</code>,
  <code>_Example/4_Vs_inversion.py</code>, <code>docs/INVERSION.md</code>. Run folder:
  <code>~/Desktop/BRUIT-FM/inversion_yv_2026-10-29</code> (scripts a–e, results, figures).
</footer>
</main>
"""

if __name__ == "__main__":
    main()
