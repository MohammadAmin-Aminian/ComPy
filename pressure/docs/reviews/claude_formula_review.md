# Claude formula and signal-processing review

> Historical contributor review, imported selectively from `claude/formula-review`
> at `c767579` on 21 October 2026. Its proposed tasks are not blanket authorization
> to change scientific defaults. Current dispositions are in
> [the deployment review](deployment_readiness_20261021.md).
> Its historical figures and tables are in [formula_review/](formula_review/) and
> its optional plotting/audit scripts in `scripts/` (merged on 28 October 2026).
> They are not saved results and are not accepted calibrations.

| Item | Value |
|---|---|
| Reviewer | Claude (independent read-only review requested by the project owner) |
| Date | 2026-10-19 |
| Code reviewed | `main` at `9b43ce8`; `2fb1fc0` changes documentation only, so line numbers hold |
| Branch | `claude/formula-review` |
| Data access | Saved PPSDs, paired spectra and diagnostics only; no raw waveform, master, score or saved product was changed |
| Tests | 423 passed; `ruff check .` and `ruff format --check .` pass |
| Next owner | Codex finishes the tasks in §6. Scientific defaults change only with owner approval |

## 1. Verdict

The physics and spectral core is correct. I found:

- one software defect: wrapped-phase medians (§3, F1);
- two interpretation problems in saved results (F2, F3);
- two low-impact range or edge issues (F4, F5).

None of them changes a saved YV gain or YV PPSD.

## 2. Verified correct

Synthetic checks C1–C6 are in `docs/reviews/claude_formula_validation.py`.

| Item | Location | Check |
|---|---|---|
| P/a_up = ρ·tan(rH)/r, r = ω√(α⁻²−c⁻²), tanh branch, Taylor series | `rayleigh_calibration.py:43` | C1: independent ODE, 2.5e-8 |
| Same relation in the primary texts | — | Zha & Webb (2016) Eq. A5, read in full: z down, a = (1/ρ)∂P/∂z is upward. Thesis Eqs. 3.7–3.11. An et al. (2017) Eq. 5 is the incompressible tanh(kh)/kh limit |
| Gain direction: g multiplies pressure; Zha coefficient = 1/g | `fit_multiplier` | Zha Fig. A2 caption; thesis Fig. 3.13 caption |
| Microseism admittance with free-surface gravity, (1+e²)/(1+e·tan y) | `microseism_check.py:62` | Derived independently. C1 agrees. The incompressible limit gives ω² = gκ tanh κH |
| Structure-free bound L ≤ L(c=∞) below α/(4H); holds for mode mixtures under H1 | `microseism_check.py:97` | Derivation |
| Polarity: D/P = −1/(ρLω²), so correct metadata gives 180°; Rayleigh P·conj(a) gives 0° | `microseism_check.py:314`, `rayleigh_calibration.py:756` | Consistent |
| Local-P θ = ωhη/c₀; ρh·tanθ/θ; first transmission ρc₀/η | `local_p_physics.py` | C4: sign and scale. Plane-wave impedance |
| DPSS L2 norm (Kmax given) × √(2/fs); zero-coherence level; Maeda AIC | `local_p_physics.py:172` | Standard forms |
| H1 upper, H2 lower, noise-corrected gain between them | `local_p_check.py:596` | Derivation |
| Pole reflection (s+a)/(s−a) is all-pass; log-gain weights 2Kγ²/(1−γ²); high-pass phase lead atan(fc/f) | `calibration_experiment.py` | Bendat & Piersol; NumPy/evalresp convention |
| Partial coherence (Schur complement); fixed real projection commutes | `microseism_horizontal.py`, `microseism_fixed_projection.py` | Standard |
| PPSD pressure path: hydrophone handling, no ω², McNamara nfft and taper, dB-bin averaging (ObsPy 1.5.1 source) | `ppsd.py` | C6: white pressure noise level within −0.3 dB (log-χ² bias) |
| Brown digitisation pixel→dB/Hz and µPa²→Pa² (−120 dB) | `reference_data` | Recomputed rows |
| Godin form: surface S_η = B·tanh^(−3/2)(kH) with (ρg)²/cosh²(kH) | `infragravity_test.py:60` | Zabotin et al. (2015); Ardhuin et al. (2014). C2: 2000 Pa²/Hz lies within 0.05–0.29 dB of B(ρg)² |
| One-sided density 2/(fs·Σw²), cross-spectra averaged first, upward sign from SEED dip | `paired_spectra.py:81` | Standard |
| Noise-model histograms: union-edge pooling, per-station normalisation, CDF quantiles | `noise_model.py`, `plot_network_ppsd.py` | Code reading |

## 3. Findings

### F1 (medium, defect): linear median of wrapped phases

**Locations:**
- `microseism_diagnostic.py:177`;
- `microseism_check.py:373` and `:471` (`median_phase_deg`);
- `calibration_experiment.py:79`;
- `heldout_validation.py:69`.

The diagnostic value feeds `phase_family` (`microseism_compare.py:145`) and `scripts/compare_microseism_methods.py:132`.

**Evidence** (`formula_review/phase_baseline_46_stations.csv`):
- 15 of the 46 coherent baseline stations are off by more than 5°, up to 19°.
- F02W is reported at 160.6° against a circular median of 180.0°. C01W is −163.2° against −179.7°; B09W is −159.5° against −173.9°.
- No family label flipped. However, F02W's bins are split 51/49 across ±180°, and an even split would report about 0°.
- The "distance to family centre" values (up to 20°) in `microseism_baseline_comparison.csv`, and the phase column of `microseism_method_comparison_table.md`, are inflated.

**YV result** (`formula_review/YV_phase_all_methods_weeks.csv`; 308 rows, 41 stations, 4 completed weeks):
- 0 of 175 coherent rows (γ² ≥ 0.5) are off by more than 5°, because YV sits near 0°.
- 11 family flips occur, all in incoherent rows (γ² ≤ 0.04).
- RR47's opposite polarity (about 169°) is confirmed.

**Fix:** use a circular statistic. `local_p_physics.circular_mean_deg` and `microseism_compare.circular_mean_degrees` already exist.

### F2 (medium, interpretation): `phase_family` has no coherence gate and no "neither" class

`microseism_compare.phase_family` forces every station to 0° or 180°.

- RR31, RR53 and RR57 (γ² ≤ 0.02) are labelled near_0.
- RR22 and RR56 (about ±110°, where Claude's screen reports `unexplained_phase`) are labelled near_180.

**Fix:** add an "indeterminate" class for low coherence or a large distance from both centres. The owner chooses the thresholds.

### F3 (interpretation): the 0°/180° split is physically resolved

`microseism_baseline_comparison.md:24` calls the split unresolved. Below α/(4H), P = −ρLω²·D_up, so D/P = 180° when the metadata polarity is correct. YL/7D/XO therefore behave as expected, and YV is the inverted family, consistent with the provisional RR52 sign −1.

Neither P/Z method can tell whether the pressure or the vertical channel is the inverted one. An independent check of vertical polarity is needed, for example teleseismic P first motion against synthetics.

The labels `nominal` (0° = sign −1) and `reversed_nominal` (180° = metadata correct) read backwards. Add a docs alias; don't rename frozen rule IDs.

### F4 (medium for interpretation): the nearby test is not depth-neutral in the IG band

Godin's model (C3) predicts a level change from water depth alone:

| Depth | Depth difference | 15 mHz | 20 mHz |
|---|---|---|---|
| 4.5 km | 10 % | −3.5 dB | −6.3 dB |
| 3 km | 5 % | −1.1 dB | −2.1 dB |

The comparison takes residuals in dB without depth normalisation (`nearby_test.py:173`, `:311`).

**YV:** 0 compared bins are affected (`formula_review/YV_nearby_depth_bias.csv`). The four YV comparisons start at 0.0244 Hz, above f = √(g/2πH) ≈ 0.018 Hz. Other networks and shallow sites are not checked yet.

### F5 (low): two edge issues

- **Godin model below 2 mHz.** The model is compared from 1 mHz (`infragravity_test.py:80`). Ardhuin et al. (2014) note that atmospheric forcing dominates below about 2 mHz. YV has 0 affected bins (`formula_review/YV_godin_below_2mhz.csv`).
- **Linear interpolation across sub-sample tears.** `rayleigh_calibration.py:456/470` accepts a negative sub-sample tear between files and bridges it by linear interpolation. This attenuates the upper local-P bands, by up to about −46 % at 16 Hz for a half-sample tear at 50 Hz. It does not apply to YV.

### Still open from earlier audits

These are unchanged in code:
- SENT-010: taper before differentiating, `rayleigh_calibration.py:752`;
- SOC-002: INSU glitch removal;
- SOC-003: density convention;
- SOC-005: noise bias of the PSD-ratio estimator. Local-P already implements a noise-corrected ratio and an H1/H2 bracket that could be reused.

The thesis computes phase velocities with PREM and a 20-min window after a PhaseNet P pick. The code uses the RR52 layered model and an 1800 s window centred on the 40 mHz group arrival. Adding the water layer changes T by at most 0.22 % (C5).

## 4. Median PPSD with provisional gains (YV, display only)

**Command:** `scripts/plot_median_ppsd.py`.

**Figures:** `formula_review/YV_{pole_trial,baseline}_median_ppsd_{overview,per_station}.png`.

The calibrated curve is median + 20·log₁₀(g); nothing is applied to saved products. RR12 is excluded as a dead channel at −3076 dB.

- **Outside 0.03–0.07 Hz:** the gains reduce the station-to-station spread by about 1–2 dB (1 Hz: 2.3 → 1.1 dB). They also close the DPG–hydrophone gap in the microseism band (about +5 → 0 dB).
- **Inside the 0.03–0.07 Hz fit band:**
  - The gap opens from −0.9 to −5.4 dB (31-station pole-trial set).
  - The spread rises from 4.3 to 4.9 dB.
  - With the baseline set (16 stations), the gap goes from 0.0 to −4.0 dB and the spread from 2.3 to 3.2 dB.

This band is the pressure noise notch, so the medians there may be set by local or instrument noise. That is a hypothesis for Socrates or the owner, not a verified result.

- **Shapes worth review:**
  - RR31 has no IG or microseism structure.
  - RR17 and RR18 rise steeply at low frequency (responses with a positive pole).
  - RR53 and RR57 have low IG levels.

## 5. Gain distribution

**Command:** `scripts/plot_gain_histogram.py`.

**Figures:** `formula_review/gain_histogram_{pole_trial,baseline}.png`.

| Pole-trial group | n | Median g | Range |
|---|---:|---:|---|
| YV DEPAS hydrophones | 23 | 1.140 | 1.015–1.592 |
| YV INSU DPGs | 8 | 0.698 | 0.585–0.820 |
| YL | 22 | 0.602 | 0.526–0.793 |
| 7D | 4 | 0.905 | 0.780–1.080 |
| XO | 4 | 0.932 | 0.909–0.963 |

YV separates completely by sensor family. g is relative to the vertical, and the two pools also use different seismometers, so a metadata offset in either sensor would produce this.

The high hydrophone tail is RR44, RR20, RR48 and RR55 (1.31–1.59). RR55's hydrophone identity is already questioned (SOC-001).

## 6. Tasks for Codex

1. **F1.** Replace the linear phase medians with a circular statistic at all five locations. Add a regression test where phases straddle ±180°, e.g. [179, −179, 178, −178] gives about 180° and near_180. Regenerate only the derived comparison tables, in new folders, and keep the historical products.
2. **F2.** Add an `indeterminate` phase family gated by coherence and by distance from both centres. Ask the owner for the thresholds.
3. **F3.** Correct `microseism_baseline_comparison.md:24` and add a label glossary (`nominal` = sign −1).
4. **F4.** Flag or depth-normalise nearby residuals in the IG band. Extend `scripts/review_formula_findings.py` to the YL/7D/XO result roots to check whether any compared bin is affected.
5. **F5.** Choose the Godin lower bound (proposal: 2 mHz) with the owner. Refuse negative sub-sample tears in `native_resampling_samples`, or interpolate band-limited.
6. **Integration.** Integrate this branch:
   - Add a numbered README/research-log step.
   - Add small tests for `pressure_gain_db`, dead-channel exclusion and the sensor-family split used by the plotting scripts.
7. **Socrates or owner.** Explain why the gains widen the family gap inside the fit band but narrow it elsewhere (§4).

## 7. Historical reproduction on the source branch

The following commands are historical; they are not part of main's reviewed
deployment workflow. In particular, its plotting scripts overwrite
fixed output names and use station-level gain matching; they are not shipped as
deployment-ready tools. The synthetic `claude_formula_validation.py` is retained
on main as an attributed worksheet, not an assertion-based test suite.

```bash
PYTHONPATH=src python docs/reviews/claude_formula_validation.py                  # C1-C6
PYTHONPATH=src python scripts/review_formula_findings.py --output-dir docs/reviews/formula_review
python scripts/plot_median_ppsd.py --gain pole_trial    # also --gain baseline
python scripts/plot_gain_histogram.py --gain pole_trial
```

By default the plotting scripts write to the git-ignored `results/REVIEWS/`. The committed copies were generated from the main checkout's saved results.

## 8. Sources read

- **Full text, local:**
  - the owner's thesis, §3.1.6;
  - Zha & Webb (2016), Appendix A;
  - An et al. (2017), GRL 44, doi:10.1002/2017GL075137;
  - Stähler et al. (2016).
- **Open access:**
  - Ardhuin et al. (2014), Ocean Modelling 77 (Archimer 00181/29272/28639);
  - Rawat et al. (2014), GRL;
  - Zabotin et al. (2015), EPS, as quoted for B.
- **Not available:**
  - Godin et al. (2013): HTTP 403 at Wiley, IFREMER and WHOI;
  - Deng et al. (2022): the local file is An et al. (2017); the earlier Socrates audit read Deng in full.
