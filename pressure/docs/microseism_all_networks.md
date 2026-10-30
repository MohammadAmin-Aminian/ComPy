# Microseism check at all stations, and the four-column comparison

**Current main integration (15 October):** the legacy 35-station best-band gain
set below is superseded. Use [the frozen-rule review](microseism_frequency_table.md)
and [completed independent weeks](microseism_independent_weeks.md): seven
conditional repeatable estimates, five outside the YL design set. See also the
[integration decision](reviews/microseism_integration.md). Claude's core water
calculation, Codex's model-free diagnostics and the handoff's frequency/week rule
are retained together. All gains remain conditional research results.

Branch `claude/microseism-all`. This branch merges `origin/main` (Codex's
microseism campaign, diagnostic and `paired_spectra` resampling and output
directory) under the independent Claude method. Method:
[microseism_pressure_vertical_check.md](microseism_pressure_vertical_check.md).

## Command

```bash
cd ~/Projects/bruit-fm-claude
bash scripts/run_microseism_check_all.sh        # about 25–35 min, 89 stations
```

Useful options:

* `--network YL` (repeatable) or `--station A02W` limits the run;
* `--start YV=2013-06-10` sets a common week for one network;
* `--resume <run folder>` continues an interrupted run.

The wrapper prints the run folder. The four-column table is
`method_comparison.md`, with the figure in `method_comparison.png` and full
values in `method_comparison.csv`.

## What runs

**Inputs.** Every pressure epoch with a saved response snapshot is used, read in
place from the main checkout:

| Network | Master | Waveforms | Notes |
| --- | --- | --- | --- |
| YV | `exports/` | `~/Desktop/Data_Bruit_FM/Data` | 43 stations |
| YL | `exports/` | `~/Desktop/Data_Bruit_FM/Data` | 25 stations |
| 7D | `7D_ALL_NO_TILT/_master` | its `_inputs/waveforms` | 15 stations |
| XO | `XO_ALL_GAUGES_NO_TILT/_master` | its `_inputs/waveforms` | 6 stations; audited 1 Hz working copies |

Channels pair automatically: each station has one pressure and one vertical
epoch.

**Interval.** The paired-spectra default week centred on each station's shared
coverage. This is the same rule as Codex's campaign, so for its 46 stations the
paired spectra are identical. A smoke test confirmed this byte for byte at
YL.A02W, 7D.J30A, XO.LA33 and XO.WD66.

**Columns.** `compare_microseism_methods.py` builds four columns per station:

1. **Claude microseism check.** Phase screen with positive-real poles read from
   each channel's saved response (pressure and vertical), the implied pressure
   sign, and the g bracket (H1–H2) where informative: conditioning ≤ 1.5 and
   phase explained.
2. **Rayleigh teleseismic calibration.** Saved baseline (`previous_g`) and
   pole-trial g from `docs/socrates_all_station_comparison.csv`, plus the
   pressure sign each Rayleigh run assumed (YV −1, YL/7D/XO +1). Nothing is
   rerun. Use `--gain-table` to supply a regenerated table with the same
   columns.
3. **Codex microseism diagnostic.** `bruitfm.microseism_diagnostic.run_diagnostic`
   runs on the same paired spectra: median coherence and phase over
   0.1–0.2 Hz, the Ruan > 0.95 fraction, Saikia selected days, and nominal
   amplitude. Codex's method reports no gain. The display multiplier is the
   baseline g, otherwise the trial g, otherwise 1. Where the user's 46-station
   Codex campaign exists, its phase is attached as a cross-check.
4. **Published values in the g direction:**
   * Aminian et al. (2025) Table 1 for YV;
   * Zha & Webb (2016) Table 1 for YL, inverted from 1/g;
   * Deng et al. (2022) for XO WD48, WD66 and WD70. These are per band within
     ~0.02–0.07 Hz, and WD66's values are read from their Fig. 4g. They are
     compared only when the check's sub-band overlaps exactly one of their
     bands (otherwise they are marked † in the table);
   * 7D: none, because Ruan et al. (2014) fit scale factors but did not tabulate
     them.

**Writes.** Nothing outside the run folder. Masters are hash-checked before and
after.

## Results: all 89 stations

Four runs, one per network, each `diagnostic_only`, with masters unchanged and no
gain applied:

| Network | Run folder (`results/EXPERIMENTS/MICROSEISM_ALL/`) | Checked | Unpaired |
| --- | --- | --- | --- |
| YV | `run_20261014T183508Z_1f096dca` | 41 | RR12, RR42 |
| YL | `run_20261014T183501Z_0a1ff170` | 25 | — |
| 7D | `run_20261014T183456Z_61c4966a` | 12 | J28A, J53A, J63A |
| XO | `run_20261014T183418Z_736f654a` | 6 | — |

Combined outputs:

* [microseism_method_comparison_table.md](microseism_method_comparison_table.md);
* [microseism_method_comparison.csv](microseism_method_comparison.csv) and
  [.json](microseism_method_comparison.json);
* [figures/microseism_method_comparison.png](figures/microseism_method_comparison.png).

Merge command:
`compare_microseism_methods.py --combine --docs docs <run> <run> ...`

### Same spectra, two codes

Where Codex's campaign used the same week (46 stations), its saved phase equals
column 3 exactly, so the paired spectra are identical. Claude's coherent-bin
mean over 0.05–0.30 Hz and Codex's median over 0.1–0.2 Hz agree within a few
degrees. The exceptions are stations with almost no coherence (RR31, RR53,
RR57, YL.A14W, XO.LT02) and the four 7D stations below.

### Polarity

69 stations fit a phase hypothesis within 20°: 65 within 10°, and 4
approximately, within 10–20° (RR55, B09W, C08W, C09W). Of these, 68 imply the
same pressure sign the Rayleigh run assumed: YV ≈ 0° → −1; YL, 7D and
XO ≈ 180° → +1.

* The YV-versus-rest split therefore matches the per-network sign conventions
  already in use. Both independent methods agree on relative P/Z polarity.
  Neither establishes which channel carries the inversion.
* The exception is **YV.RR47** (171°, sign +1 against the assumed −1).

Phase not explained (rms > 20°; amplitude not interpreted):

* **YV RR17, RR18, RR22, RR56.** RR17 and RR18 have large recorded poles; RR22
  and RR56 sit at ±110°.
* **YL.B01W, at 150°.** Zha & Webb's coefficient (0.19 → g 5.26) is anomalous
  too.
* **7D J35A, J36A, J43A, J44A.** These are near 0° at 0.1–0.2 Hz, which is the
  opposite polarity to the rest of 7D, and drift strongly negative below about
  0.1 Hz. None has a Rayleigh gain, which is consistent with a response problem
  at these four.
* **7D.J52A (133°) and 7D.J38A (rms 21°, marginal).**
* **XO.LT02**, which is incoherent.

### Gain (35 informative stations)

* **YL agrees with Rayleigh.** At 12 of 17 informative stations the microseism
  midpoint is within ±2% of the Rayleigh baseline.
  * The baseline is 3.5–7.6% above Zha & Webb (2016) at every station except
    C07W. The microseism shows the same offset everywhere except C03W.
  * At conditioning ≈ 1.01 the microseism value is almost free of the Earth
    model. It still shares several inputs with the Rayleigh baseline: the
    nominal pressure and vertical responses, the master depths and
    ρw = 1028 kg/m³.
  * The common offset can therefore come from one of those shared inputs or
    from a convention in Zha & Webb that is not reproduced here. These data
    cannot tell the two apart (open question; see the
    [self-review](reviews/claude_microseism_self_review.md)).
* **YL.C07W.** *(retracted 15 Oct)* The base week (0.544–0.595) seemed to
  reject the baseline 0.269. Three of the four repeated weeks, however, give
  0.27–0.31, which agrees with it. See "Repeated weeks: first results".
* **YL C03W, C11W, A12W, S01W.** These differ from baseline by −7%, +8%, +4% and
  +3%.
* **YV.** RR40 and RR52 agree with baseline and published within 1%. RR28 and
  RR29 brackets contain both references. RR38 again sits 3.7% below baseline
  and 2% above published. RR50 is again high (+12%).
  * One hypothesis for RR50 is a crust slower than the RR52 model (Aminian et
    al. 2025, Fig. 9g). It is not shown, because the Rayleigh reference uses
    the same model in a different band.
  * The YV published values (Aminian et al. 2025) come from the same
    teleseismic Rayleigh P/a method family as column 2. They are not an
    independent check of it.
  * This run used each station's default centred week, not the common
    2013-06-10 week of [yv_microseism_check.md](yv_microseism_check.md).
  * As a result RR34 (a March week, conditioning 2.97) is not informative here,
    and the RR29 and RR36 brackets are wider.
  * H1–H2 width therefore depends on the week's long-wave coherence.
  * Reproducibility across independent weeks is **not yet tested**:
    * the RR38, RR40 and RR50 weeks overlap the 2013-06-10 week by 4–6 days;
    * RR52 used the identical week;
    * only RR28 and RR29 (1 day of overlap) and RR34 (March) are near
      independent, and their brackets changed.
* **7D.** J30A, J45A and J54A sit 4–8% below baseline. J29A (≈ 1.00), J55A
  (≈ 0.79) and J67A (≈ 1.00) are new values; no Rayleigh or published gain
  exists for them.
* **XO.**
  * LA33 and LA34 are 5–6% below baseline.
  * WD66 and WD70 are 2–3% above baseline.
  * Deng et al. (2022) comparisons are made only when frequency matched.
    * WD66's sub-band 0.034–0.043 Hz only partially overlaps their 0.02–0.04 Hz
      band (~0.92). The earlier +1.2…+2.7% calculation is historical, not a full
      frequency match. New comparisons require whole-sub-band containment.
    * WD70: their 0.91–1.71 is frequency dependent with no band mapping
      given, so there is **no test**.
    * WD48: the sub-band is 0.115–0.171 Hz because the station is shallow.
      This lies outside Deng's ~0.02–0.07 Hz, so its 1.065–1.122 against their
      0.95 is **not a like-for-like test**.
* **Earth-model dependence.** Every gain bracket is a mode-0 value in the
  historical RR52 layer model. The structure-free maximum (vertical incidence,
  H2) bounds how far the Earth model can move g upward:
  * +1–3% at 12 YL stations;
  * +4–6% at RR52, YL S01W, 7D J45A, J54A, J55A, J67A and XO WD66, WD70;
  * +8–13% at RR38, RR40, RR50, YL C07W, C11W and XO WD48;
  * +19–36% at RR28, RR29, RR36, YL A10W, C03W, 7D J29A, J30A and XO LA33,
    LA34. There the value leans on the RR52 model.
* **Uncertainty.** Brackets are noise brackets from one week, with no
  confidence intervals. The repeated-week runs below give the first
  empirical spread.

## Repeated weeks (empirical uncertainty)

A single week gives only a noise bracket. To measure week-to-week spread, each
station's week in a completed base run is moved by whole weeks. There is no
fallback: stations whose coverage does not reach a shifted week stay
unresolved. Each shift is a new run folder, compared as usual. The thresholds
were fixed before running: stable ≤ 5%, moderate ≤ 10%, otherwise variable
(range of weekly midpoints ÷ median).

```bash
cd ~/Projects/bruit-fm-claude
bash scripts/run_microseism_weeks.sh <base run folder>
```

The default shifts are −12, −6, +6 and +12 weeks; pass others after the folder.
Shifted runs are written beside their base run, where the combiner looks for
them. Handoff for the remaining work: [microseism_claude_handoff.md](microseism_claude_handoff.md).
After all four networks:

```bash
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  scripts/combine_microseism_weeks.py <YV base> <YL base> <7D base> <XO base> \
  --docs docs
```

This writes `microseism_weeks.csv/.json`, `microseism_weeks_table.md` and
`figures/microseism_weeks.png`.

### Repeated weeks: first results (YL only; 15 October)

YL ran all four shifts (−12, −6, +6, +12 weeks). The other networks did not
complete:

* **XO and YV:** no shifted runs were made.
* **7D:** every shifted week fell outside the on-disk data, which covers only
  33–41 days (late November–December 2011). 7D needs shifts of ±1–2 weeks.
* **YL:** the four shifts ran twice with identical results, so the combiner
  now keeps one run per shift.

**Week-to-week spread** (YL stations with ≥ 2 informative weeks):

| Class | Stations |
| --- | --- |
| stable (≤ 5%) | A10W, B09W, C05W, C09W, F02W, N03W |
| moderate (≤ 10%) | A02W, B04W, B14W, C11W, C16W |
| variable (> 10%) | A06W, A12W, B11W, C01W, C03W, C07W, C17W, S01W |

**Method flaw found (not an instrument effect).**

* At YL.C01W, the base week gives 0.58 and the other four weeks 1.12–1.15.
* Codex's model-free 0.1–0.2 Hz P/Z amplitude is approximately flat across all
  five weeks (3.0–3.5 × 10⁻⁷ m/Pa). This supports the sub-band-selection explanation;
  the ratio alone does not establish that both instruments were unchanged.
* In every week the 0.14–0.20 Hz sub-band gives about twice the lower bands.
  The best-sub-band rule picked the low bands when they were coherent and
  the top band when they were not.
* Above about 0.1 Hz the mode-0 RR52 model fails, as already found at RR52
  (sediments and mode content; Ruan et al. 2014). There the true admittance
  can fall *below* the mode-0 value, outside the assumed bracket.
* The conditioning test cannot detect this.

**Scale.** Across the 100 informative station-weeks (four base runs and four
YL shifts):

* 48 used a sub-band reaching above 0.1 Hz;
* 63 either rest on a single sub-band or have sub-bands (conditioning ≤ 1.5)
  that disagree by more than 10%.

**Consistent set.** In the base week only 16 of the 35 informative stations
pass that frequency-consistency check:

* YV: RR40, RR52;
* YL: A02W, B04W, B14W, C03W, C05W, C11W, C16W, C17W, N03W;
* 7D: J45A, J54A, J55A, J67A;
* XO: WD66.

**YL.C07W.** All sub-bands agree within each week, but the weeks differ by
×2: 0.27–0.31 in three weeks and 0.54–0.60 in the base week. The Rayleigh
baseline (0.269) and pole trial (0.526) split the same way. This points to a
real time-varying long-period response or wavefield at C07W; it is a lead, not
a conclusion.

**Consequence.** The 35-station "informative" gain set above is
**superseded**. Gains should be quoted only where coherent sub-bands
≤ 0.1 Hz agree within 10% and the repeated weeks are stable. This rule was
added after seeing these data (post hoc), so it must be applied blind to the
remaining networks' weeks before any value is trusted.

## Rayleigh method (column 2), optional rerun

The saved 89-station Rayleigh results are already complete and reviewed; column
2 reads them. To regenerate the trial independently (historical command; writes
a new root in the main checkout):

```bash
cd ~/Projects/bruit-fm
PYTHONPATH=src conda run --no-capture-output -n bruitfm python \
  scripts/run_calibration_experiment.py --network all --all-stations --open-plots
```
