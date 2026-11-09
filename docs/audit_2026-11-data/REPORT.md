# Whole-archive audit: the open tests of the ComPy 3 and pressure audits, run on the data (November 2026)

Branch `compy3-audit-data`, stacked on `pressure-audit` (PR #22), which is stacked on
`compy3-audit` (PR #21). Those two audits ran without ObsPy or waveforms and left their
discriminating tests "for the IPGP machine". This branch runs them on the whole-archive
run (`~/Desktop/ComPy`, ComPy 2.7–2.8, November 2026):

- **287 inverted stations** (at least 100 windows), in 17 networks;
- cleaned 1.25 Hz data, raw counts, hourly tilt, 3-h PPSDs with cross-spectra;
- stage C, D and E outputs.

**No default changes:** every archived stage C result is reproduced bit-for-bit (§7).
Everything that would change a number is a question under [Decisions](#9-decisions).

**Evidence labels**, as in #21 and #22:

- **Verified:** computed here, with the script named.
- **Hypothesis:** consistent with the evidence, with the test that would settle it.

## 0. What was run

| Script (`scripts/`) | What it does | Inputs |
|---|---|---|
| `s01_stage_c_restack.py` | Per-window stage C spectra of every station, cached once. Then stage C's own `run_gain` restacks them through `stacked_compliance_from_spectra` for GJI and wide gates, inclusive and exclusive. It also saves per-window levels, which stage C does not | `02_cleaned`, stage D gains |
| `s02_phase_metadata.py` | Residual Z/P phase of the stack that stage E inverted, joined with instrument type and the poles and zeros that stage A deconvolved | stage C, StationXML |
| `s03_dpg_tide.py` | Fits a corner factor k to the phase, measures M2 and K1 in hourly counts, deconvolves the metadata and compares with the EOT20 tide | counts, StationXML, EOT20 |
| `s04_level_steps.py` | Z and P levels of each station minus its network median (microseism and compliance bands), with changepoints | `04_ppsd` |
| `s05_anomaly_scan.py` | `compy.anomaly` on every station and network, plus cross-checks: gate bias, which sensor stepped, stage D event gains, the pressure response | all of the above, stage E |
| `s06_tilt_tides.py` | Tidal lines in the hourly tilt angle and the tilt-correction variance reduction | `03_tilt` |
| `s07_figures.py` | Figures and `numbers.json` | — |

The waveform figures of #21 were made on real data with
`examples/compy3_waveform_diagnostics.py` (`figures/waveform/`). Every number below is in
[`numbers.json`](numbers.json) or in the tables of [`data/`](data/README.md).

**Not run:**

- the Rayleigh calibration extended to 10–20 mHz (#21, §5.A test 2);
- the cleaning-window boundary test (#21 M7);
- #22's YR hydrophone and 7K microseism tests;
- anything at the native sampling rate.

## 1. Executive summary

1. **The residual Z/P phase is the pressure gauge's response (verified).** It was #21's
   main open question.
   - **Which gauge type shows it:**
     - absolute gauges (APG, no high-pass): −0.6° (n = 33);
     - differential gauges (DPG): +6.0° at the band bottom (n = 168, 88 % positive);
     - hydrophones: −5.7° (n = 57).
   - **With the seismometer held fixed:** at the deep, good 7D stations with the same
     Trillium Compact, APGs give −0.2° and DPGs +5.4° at 8 mHz (Mann-Whitney p = 0.02).
   - **Ruled out:** the seismometer-corner hypothesis.
2. **The tide confirms it independently (verified, Doran et al. 2019 test).**
   - **The tide is visible.** DPG counts record the M2 tide clearly: amplitude 10–240 times
     the noise at 145 DPGs.
   - **The metadata over-attenuate it.** After deconvolution, M2 is a median +8.6 dB
     above EOT20 (90 % of 145 DPGs above), so the metadata attenuate too much at tidal
     periods.
   - **The phase predicts the tide.** The phase-fitted corner predicts the measured tidal
     ratio (Spearman 0.48, p = 10⁻⁹).
   - **It holds station by station within a network:** −0.67 (7B), −0.66 (XE), −0.87 (XO).
     This is unit-to-unit variation of the DPG time constant.
   - **Update ([RESEARCH.md](RESEARCH.md)):** the tide ratio is flat across O1, K1, M2
     and S2, so the tide measures the leak corner directly. That corner predicts the
     compliance-band phase (Spearman 0.82, n = 119). The +3.4 dB offset of the first
     version came from the weaker phase-only corner estimate.
3. **Implied compliance bias.** After the 20–60 mHz gain:
   - DPG stations: median −3.0 % at the band bottom from the tide (p10 −6.4 %,
     p90 −0.6 %), or −2.0 % from the phase;
   - YV: −3 to −7 %, and RR40 +7 to +9 %.

   This is much narrower than #21's envelope of −1 to −22 %. Stage E's sampled gain does
   not absorb it (ρ = 0.00).
4. **Polarity follows the DPG family, and the compliance phase measures it (verified).**
   - **The rule:** SIO-metadata DPGs sit at 0°, LDEO/WHOI-metadata DPGs at 180°, and both
     families are mixed inside 7D. This explains #22's P7, the 20/12 local-P split in 7D.
   - **It predicts the stage D sign:** at 107 of 108 resolved stations.
   - **It matches the local-P polarity:** at 61 of 64 channels.
   - **BRUIT-FM's `rayleigh_sign` is not a polarity:** it disagrees with the rule at 69 of
     154 stations.
5. **The coherence-gate bias is real but small (verified).** Inclusive (published) minus
   exclusive:
   - overall: median +0.24 %, p90 +2.2 %; at least 1 % at 58 of 280 stations;
   - YV +1.2 %, against the 1.4–4.1 % that #21 modelled; XE +0.2 %;
   - the bias rises as acceptance falls (Spearman −0.77).
6. **Epoch mixing goes beyond the 108 stations stage D could scan (verified).**
   - **How common it is:** 165 of 287 stations have at least one \|Z/P\| changepoint. The
     stack inverted by stage E shifts by at least 1 % at 112 stations and at least 5 % at
     57, 21 of them tier-1 stations rated "good".
   - **The steps are real sensitivity changes:** they agree with the stage D event gains
     (ρ = 0.71, n = 73).
   - **YV, September 2013:** #21's synchronous step is reproduced with the independent
     whole-archive processing.
   - **RR50 and RR36 (correction):** their "pressure" jumps are steps of the **vertical**
     (+4 to +5 dB in both bands, P flat).
   - **8J:** the hydrophones drift together, by 10–50 %, between August 2021 and
     February 2022.
7. **Tidal currents modulate the tilt.** The tilt-correction variance reduction has an M2
   line at 103 of 271 stations: 78 % of 7B, 75 % of 8J and 60 % of YV stations. At 91 % of
   7B stations the tilt angle also has a 24-h (S1) line; that cause is open.
8. **Tools:**
   - `compy.response_check` (tested): corner fit, tide fit, bias, sign from phase;
   - seven scripts;
   - a ranked table of 1,948 anomalies.

   **Bugs found:**
   - **Pressure atlas, APG tide check:** a single NaN dropped every hour, so 7 APGs returned
     "no tidal signal".
   - **`stage_c_wide.py`:** cannot run the exclusive gate.
   - **XO StationXML:** station epochs that do not cover its channels.
   - **Anomaly score:** it saturates.

## 2. Audit findings (severity-ranked)

Severity and status follow #21:

- **High:** can move compliance by at least 1 % or invalidate a station.
- **Status:** Fixed, Option (off by default) or Open.

| # | Sev. | Where | Finding | Evidence | Effect | Status |
|---|---|---|---|---|---|---|
| A1 | High | Stage A responses of every DPG (StationXML), `pipeline/a_preprocess.py` (`remove_response`) | **The metadata DPG response is wrong at low frequency**, station by station: the true leak corner is lower than in the metadata at most stations (k < 1), higher at some (RR40, ZZ) | **Verified** by two independent observables (§5.A) | Compliance −3.0 % median at the band bottom from the tide (p10–p90: −6.4 to −0.6 %; −2.0 % from the phase), YV −3 to −7 %, RR40 +7 to +9 % | Open: D1, D2 |
| A2 | High | `12_inversion_all_stations` inputs; stage C over the whole record | **Epochs mixed in the stacks** at stations stage D could not scan (no event gains) | **Verified**; 21 tier-1 "good" stations with a shift of at least 5 % (`numbers.json` → `epochs`) | Up to tens of per cent (7D.G19B −36 %, 7D.J68A −30 %, 8J several) | Open: D3 |
| A3 | High (8J) | 8J hydrophones | **Network-wide drift** of \|Z/P\| (steps at about 25 stations, August 2021 – February 2022), following the event gains | **Verified** (synchronous steps; ρ = 0.71 with stage D gains) | Whole-record 8J compliance mixes sensitivities; the network-relative level test cannot see a common drift | Open: D4 |
| A4 | Med | `compy_pressure` sign choice; stage D second pass (D36) | **Polarity differs by DPG family**; the compliance phase gives the stage D sign (107/108) at no cost | **Verified** | None on \|η\|; saves the second stage D pass; fixes the BRUIT-FM sign column | Option `response_check.pressure_sign_from_phase`: D5 |
| A5 | Med | `compy/compliance.py` inclusive gate | **Selection bias measured:** median +0.24 %, above 1 % at 58 stations, about 5 % at the lowest acceptances | **Verified** (s01, 280 stations) | Compliance high by that much | Option `gate="exclusive"`: D6 |
| A6 | Med | YV.RR50, YV.RR36 (D40) | **The jumps are the seismometer**, not the DPG: Z steps +4 to +5 dB in both bands, P flat | **Verified** (relative levels, `which_sensor_yv.png`) | Stage E inputs are right (they use the stable epoch). The BRUIT-FM/atlas pressure gains of RR50 and RR36 after the jump are wrong if applied to P alone | Open: D7 (doc) |
| A7 | Med | Pressure atlas `src/s03b_tide.py` (branch `analysis/pressure-atlas`) | **NaN in the running median** makes `np.median` NaN, which excludes every hour: amplitude 0, "no tidal signal" | **Verified**: 7 APGs with `obs_M2_pa = 0` | Atlas tide check reports 40 of 47; 7 more APGs were never tested | Fixed in `response_check.tidal_harmonics`; atlas: D8 |
| A8 | Med | 7B tilt | **24-h line in the tilt angle** at 21 of 23 stations (S1), absent elsewhere | **Verified** (s06) | Possibly a daily instrument or processing cycle in the tilt correction | Open: test (§5.D) |
| A9 | Low | `scripts/stage_c_wide.py` (archive) | `Recorder` has no `min_median_coherence`, so `gate="exclusive"` raises | Code; s01 replicates the wrapper exactly | Blocks D6 for the wide pass | Open (fix if D6) |
| A10 | Low | `00_inventory/stationxml/XO.xml` (`build_stationxml.py`) | Station epochs do not cover the channel epochs: `select(time=…)` returns nothing at 36 XO stations | **Verified** | Any time-based selection fails silently | Open |
| A11 | Low | ObsPy `get_evalresp_response_for_frequencies` | Returns 0 at 22 µHz for some responses (XE); scaled from 0.1 Hz with the analog part instead | **Verified** | Tide checks | Handled in s03 |
| A12 | Low | `compy.anomaly.Anomaly.score` | Saturates at 40 (σ ≤ 8, effect ≤ 20 %), so hundreds of records tie and sort alphabetically | **Verified** | Ranking | Tie-break in s05; roadmap |
| A13 | Low | `examples/compy3_waveform_diagnostics.py` | On cleaned data the tilt azimuth splits into two bands 180° apart (tilt already removed, sign random) | **Verified** (RR52) | Misleading figure | **Fixed** (docstring) |

**Informational:**

- **Common mode:** monthly deviations are correlated across stations in 9 networks
  (permutation p ≤ 0.008).
- **7D, March 2013:** a synchronous +1.5 to +4.5 % step at 10 stations, 4–12 March 2013.
- **XE.CC07:** steps of −14 to −25 % in \|Z/P\| with no change in event gain, so the change
  is confined to the compliance band.
- **7D.G19B:** −48 % on 29 October 2012, the day after the M7.8 Haida Gwaii earthquake. It
  is the only such coincidence, and its phase is inconsistent anyway.

## 3. Diagnostic plots on real data

| Figure | Purpose | What it reveals here |
|---|---|---|
| [`phase_by_gauge.png`](figures/phase_by_gauge.png) | Residual phase by gauge type; APG vs DPG with the same seismometer | A1: the DPG, not the seismometer |
| [`dpg_tide_test.png`](figures/dpg_tide_test.png) | Tidal ratio predicted from the phase vs measured; within-network trend | A1, independently; unit-to-unit leak corners |
| [`gate_bias_measured.png`](figures/gate_bias_measured.png) | Exclusive vs inclusive gate vs acceptance, with #21's model range | A5, smaller than modelled |
| [`which_sensor_yv.png`](figures/which_sensor_yv.png) | Z and P levels minus the network, microseism and compliance bands | A6: RR50 and RR36 are the vertical; nothing in September at RR28 or RR52 |
| [`epoch_mixing.png`](figures/epoch_mixing.png) | Stack shift from mixed epochs, coloured by the sensor that stepped | A2, A3 |
| [`tilt_tides.png`](figures/tilt_tides.png) | Share of stations with tidal or 24-h lines in tilt, by network | A8, tidal currents |
| [`anomaly_ranking.png`](figures/anomaly_ranking.png) | Top 30 anomalies | Triage |
| [`YV_multistation.png`](figures/YV_multistation.png), [`XE_multistation.png`](figures/XE_multistation.png), `*_noise_vertical.png` | #21's network panels on the inverted stacks | Shape and level outliers |
| [`figures/waveform/`](figures/waveform) | #21's spectrograms, coherograms (Z–P, Z–H, P–H) and tilt, for RR52, XE.CC05 and the 7D APG/DPG pair J51A/J43A | RR52's vertical dies on 12 October 2013 (Z–P coherence lost); RR52's tilt azimuth shifts in April and June 2013 |

## 4. Anomaly and correlation detection on the archive

`s05_anomaly_scan.py` runs `compy.anomaly.scan_station` and `scan_network` on the restacked
results. Thanks to s01, these now carry per-window levels, so the step, tide and regime
detectors run on every station. Four cross-checks that need the data are added to each
step and station:

- **which sensor stepped:** s04 relative levels within 10 days;
- **stage D event gains** before and after the step;
- **the measured gate bias;**
- **the pressure-response model** with its tidal check.

**Results** ([`data/anomalies.md`](data/anomalies.md); 1,948 records at σ ≥ 2 over 287
stations and the networks):

| Detector | Records | Detector | Records |
|---|---|---|---|
| changepoint | 377 | tidal_acceptance | 266 |
| tidal_transfer_function | 215 | phase_residual | 207 |
| epoch_mixing | 144 | pressure_level_dependence | 141 |
| dpg_response | 116 | monthly_excursion | 115 |
| window_regime | 63 | gate_bias | 58 |
| gain_pull | 57 | poor_fit | 53 |
| monthly_drift | 46 | phase_inconsistent | 27 |
| outliers, synchronous_steps, common_mode, level outliers, lines, notches | 63 | | |

**Reading the table:**

- **Tier first:** tier-2 stations (not interpreted in D46) dominate the extreme effects
  (5F, 2D, 7D shelf), so filter on tier 1.
- **`gain_pull`:** marks stations whose stage E gain sits at a bound. The shallow 7D FS/FN
  stations pull to 0.60–0.78, and 2D.OBS03/06 to 1.40.
- **Phase:** `phase_inconsistent` (\|residual\| > 45°) marks 27 stations whose \|η\| should
  not be trusted.
- **`dpg_response` on hydrophones:** it also runs on hydrophones (8J, ZU and others), but
  their metadata corners sit at or above the compliance band (15–64 mHz). The implied
  biases there (up to +80 %) are large extrapolations of a one-factor model: read them
  as flags, not corrections.

## 5. Physical interpretation and discriminating tests

### A. Residual Z/P phase (#21 H1): resolved to the pressure gauge

| Hypothesis | Prediction | Test (run here) | Result |
|---|---|---|---|
| (a) Timing offset | Phase grows with f | — | Already rejected in #21 |
| (b) Pressure-gauge response ≠ metadata | Residual only with high-pass gauges; same seismometer gives none with an APG; tidal amplitude after deconvolution ≠ tide, as predicted by the phase | APG vs DPG (s02); tide vs phase (s03) | **Supported:** APG −0.6° (n = 33) vs DPG +6.0° (n = 168); 7D TC deep APG −0.2° vs DPG +5.4° (p = 0.02); tide ratio predicted (ρ = 0.48, n = 145) and within-network (ρ = −0.66 to −0.87) |
| (c) Seismometer corner ≠ nominal | Residual with any gauge, depending on seismometer model | APG stations with Trillium Compact | **Rejected:** no residual with an APG on the same seismometer. DPG residuals differ by seismometer (T240 +10.8°, CMG-3T +5.5°, TC +3.2°), but this is confounded with the DPG family, and the tide test involves no seismometer |
| (d) Poroelastic seafloor | Site-dependent, independent of gauge type | APG vs DPG at matched depth | **Rejected** (gauge type decides) |

**Model details:** `compy.response_check`, results in `data/dpg_tide.csv`.

- **Two pole-scaling models were compared:**
  - scaling only the lowest pole (the leak) fits the phase better: median rms 1.4° vs 2.2°
    at the 7D LDEO/WHOI DPGs;
  - it also predicts the tide better: ρ 0.48 vs 0.37.
- **Correction (research part, [RESEARCH.md](RESEARCH.md) §2):**
  - **The K1 statement was a bug.** The first version reported M2 and K1 ratios differing
    by 5–12 dB. `s03` had deconvolved K1 with the response at the M2 frequency. Fixed:
    M2 − K1 = −0.8 dB (median, n = 142).
  - **The ratio is flat.** Across O1, K1, M2 and S2 it has slope −0.11 (median), as a
    scaled leak corner predicts.
  - **The tide gives the corner.** The corner from the tide predicts the measured band
    phase (Spearman 0.82, median offset −1.8°).
  - **Where the +3.4 dB offset came from:** the phase-only corner. The phase constrains the
    corner weakly for LDEO/WHOI DPGs, whose metadata have a second pole at 34 mHz.

### B. Steps in \|Z/P\|: which sensor, and is it real?

| Hypothesis | Prediction | Test | Result |
|---|---|---|---|
| Sensitivity change of one sensor | The same step in the Rayleigh event gains (20–60 mHz); a level step in that channel relative to the network | s05 (event gains), s04 (levels) | **Supported in general:** ρ = 0.71 (n = 73), same sign at 85 % of steps ≥ 5 %. RR50 (30 July 2013) and RR36 (21 October 2013): **Z** +4 to +5 dB, P flat |
| Sea-state or noise regime | Step without an event-gain change; follows storms | Same | XE.CC07 and the 7D March 2013 cluster have no event-gain change: candidates |
| Common drift (whole network) | Synchronous steps; relative levels flat | Synchronous detector | **8J:** about 25 hydrophones from August 2021 to February 2022, event gains follow (UP36 −46 %, UP43 −45 %, UP46 −50 %) |
| YV September 2013 (#21 H2) | — | Same detector, independent processing | **Reproduced** (RR29 +6.5 %, RR40 +9.6 %, RR28 +7.6 %, RR52 +7.1 %, 10–19 September). Neither level moves by more than 1 dB, so the sensor is unresolved; event gains agree (+8 to +9 %) |

**Discriminating test for YV September and 8J:** the network-relative level cancels a common
change. Use a sensor-independent reference instead:

- M2 in the DPG counts (s03), month by month: a P change shifts it, a Z change does not;
- or teleseismic P-wave amplitudes on Z against neighbouring land stations.

### C. Tidal currents (#21 M4)

| Hypothesis | Test | Result |
|---|---|---|
| Current-driven tilt noise | M2 in the hourly tilt variance reduction and angle (s06) | **Supported:** variance-reduction M2 at 103 of 271 stations; angle M2 at 53 |
| Daily cycle (temperature, clock, processing) | S1 in the tilt angle | 7B: 21 of 23 stations at S1, 2 of 23 at M2 in the angle. Test: S1 phase common to all 7B stations in UTC (instrument or processing) or local solar time (thermal)? |

### D. Polarity (#21 M3, #22 P7)

| Hypothesis | Test | Result |
|---|---|---|
| Wiring or metadata convention differs by operator | Compliance phase polarity vs DPG metadata family; vs local-P; vs stage D sign | **Supported:** SIO DPG 0° (22/22 in 7D, all 7B, XE, ZU, 2D); LDEO/WHOI DPG 180° (38/38 in 7D, 21/21 XO, all 7K, X9, YO, ZZ). Local-P agrees at 61 of 64 channels; stage D sign at 107 of 108 |
| Low-SNR local-P estimates | — | Rejected (agreement 95 %) |

## 6. Implemented changes

| Commit | Change | Numbers changed? |
|---|---|---|
| Add `compy.response_check` | Corner fit, predicted tide ratio, compliance bias, NaN-safe tidal fit with spline drift, sign from phase; 13 tests | No (new module) |
| Add whole-archive audit scripts | `docs/audit_2026-11-data/scripts/s01–s07` | No |
| Add audit outputs | Figures, tables, `numbers.json` | No |
| Docs | This report, changelog, docs index, example docstring (A13) | No |

The archive run folder gets `13_audit/` (cache, restacks, tables) and decision D47 in
`DECISIONS.md`. Nothing in `06_compliance*` or `12_inversion_all_stations` was written.

## 7. Validation

- **Archived results reproduced bit-for-bit.** s01 recomputes every stage C result of the
  287 inverted stations from cached spectra:
  - 287 GJI-gate files and 370 wide-gate files (both gain labels, stable epochs included);
  - the path is `stacked_compliance_from_spectra` and `run_gain` of
    `pipeline/c_compliance.py`;
  - the maximum relative difference over compliance, uncertainty, coherence, phase and the
    per-window mean and scatter is **0.0**, with identical accepted windows;
  - 10 GJI epoch stacks have no archived counterpart.

  This also validates #21's streaming core and the replicated wide-gate wrapper on real
  data.
- **Tests:** the full suite passes with ObsPy, TiSKitPy and ArviZ (905 before this branch,
  plus 13 new in `tests/test_response_check.py`). `ruff check` and `ruff format --check`
  are clean.
- **Not validated:**
  - the corrected response itself (only its consistency with the tide);
  - the 8J and YV sensor attribution (§5.B).

## 8. Roadmap update (ComPy 3)

#21's ranking stands, with these changes (I = impact, E = effort, 1–5):

| Rank | Item | I | E | Change from #21 |
|---|---|---|---|---|
| 1 | **Pressure-response calibration per station:** joint fit of phase (compliance band), tides (M2/S2/K1/O1) and Rayleigh amplitude; deconvolve the fitted response in stage A; carry its uncertainty | 5 | 3 | Now targeted at the DPG (seismometer ruled out). Method validated by the tide; `response_check` has the parts |
| 2 | **Epoch segmentation for every station:** \|Z/P\| changepoints (stage C window levels, now saved by s01), confirmed by event gains or a tidal reference; one gain per epoch | 5 | 2 | 57 stations at ≥ 5 %, most without event gains |
| 3 | **Save per-window levels in stage C** (`<STA>_window_levels.npy`) | 4 | 1 | New: the anomaly scan needs them; s01 shows the format |
| 4 | **Automatic polarity:** `pressure_sign_from_phase` before stage D; record the sign | 3 | 1 | New, validated 107/108 |
| 5 | Systematic-error covariance | 4 | 2 | Same; the response term is now quantified (−6 to +3 %, p10–p90) |
| 6 | Exclusive gate | 3 | 1 | Lower than #21 thought (median 0.24 %); default for low-acceptance stations only? |
| 7 | Anomaly score: no cap on σ; rank by expected effect on the inverted stack | 2 | 1 | New (A12) |
| 8 | StationXML builder: station epochs that cover the channels (A10) | 2 | 1 | New |
| … | CLI and config, FDSN fetch, joint inversion, packaging, CI | | | As in #21 |

## 9. Decisions

Each is a yes/no question. "Results change" means published or archived numbers would move.

- **D1.** Accept that the residual Z/P phase comes from the DPG response, not the
  seismometer, and close #21's D1 on that basis? (Y/N) Impact: interpretation only.
- **D2.** Build the per-station response correction (roadmap 1) and rerun stage A → E for
  DPG networks with it? (Y/N) Impact: results change; about −2 % compliance at the band
  bottom (median), −3 to −7 % at YV, +9 % at RR40.
- **D3.** Rerun stage C/E on stable epochs for the 21 tier-1 "good" stations with an epoch
  shift of at least 5 %? (Y/N) Impact: results change at those stations; list in
  `numbers.json` → `epochs.tier1_good_abs_ge_5pct`.
- **D4.** Treat 8J whole-record compliance as provisional (downgrade 8J in the inversion
  report) until the drift is resolved? (Y/N) Impact: reporting.
- **D5.** Set the stage D pressure sign from the compliance phase
  (`pressure_sign_from_phase`) instead of a fixed −1 and a second pass? (Y/N) Impact: same
  gains at 107/108 stations; saves the second stage D pass.
- **D6.** Adopt the exclusive gate for new runs (at least below 30 % acceptance)? (Y/N)
  Impact: compliance −0.24 % median, up to −5 % at low acceptance; results change.
- **D7.** Reword D40 and the BRUIT-FM/atlas notes: RR50 and RR36 changed **seismometer**
  sensitivity, and their after-jump pressure gains are invalid for P alone? (Y/N) Impact:
  documentation; atlas gains of those two channels.
- **D8.** Fix the pressure atlas tide check (NaN-safe fit) and rerun it for the 7 APGs it
  lost? (Y/N) Impact: atlas numbers (40/47 → up to 47/54).
- **D9.** Save per-window levels in stage C by default (new file, no number changes)?
  (Y/N)
- **D10.** Merge this branch with #21 and #22 as ComPy 2.9.0 (new modules, defaults
  unchanged)? (Y/N)

## References

- Doran, A. K., et al. (2019). Calibration of differential pressure gauges through in
  situ testing. *Earth Space Sci.* 6. <https://doi.org/10.1029/2019EA000783>
- Hart-Davis, M. G., et al. (2021). EOT20: a global ocean tide model from multi-mission
  satellite altimetry. *Earth Syst. Sci. Data* 13, 3869–3884. doi:10.17882/79489
- Corela, C., et al. (2023). The effect of deep ocean currents on ocean-bottom
  seismometers records. *NHESS* 23, 1433–1451.
- Cox, C., Deaton, T., & Webb, S. (1984). A deep-sea differential pressure gauge. *J.
  Atmos. Oceanic Technol.* 1, 237–246.
- The references of [`../audit_2026-11/REPORT.md`](../audit_2026-11/REPORT.md) and
  [`../../pressure/docs/audit_2026-11/REPORT.md`](../../pressure/docs/audit_2026-11/REPORT.md).
