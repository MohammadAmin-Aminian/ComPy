# ComPy 3 audit: correctness, diagnostics, anomalies and roadmap (November 2026)

Branch `compy3-audit`, based on `main` at `cd41b4d` (ComPy 2.8.0). Nothing here changes a
published number by default; every proposed change that would is listed under
[Decisions](#9-decisions).

**How to read the evidence labels.** **Verified** means it was computed in this audit from
code or archived results, and the script or test that reproduces it is named. **Hypothesis**
means it is consistent with the evidence but not yet proven; each one comes with a test
that would settle it.

## 0. Scope and constraints

- **Read:** the whole repository: package, pipeline, examples, tests, docs, the archived
  YV and XE results, and the BRUIT-FM notes in `pressure/docs`. The audit builds on the
  October audit (`docs/AUDIT_2026-10.md`, items D1–D7, S1–S10) and the inversion notes
  (`docs/INVERSION.md`, V1–V8, A1–A4) and does not repeat them.
- **Environment:** no ObsPy, TiSKitPy or disba could be installed, and the raw waveforms
  are not in the repository.
- **What was run:** everything that does not need ObsPy. That includes 268 existing and
  new tests, and the full analysis of the archived products of 31 stations
  (8 YV, 23 XE): stacked spectra, monthly stacks, per-window levels (YV), accepted-window
  times and inversion records.
- **What was not run here:** the ObsPy-dependent tests run in CI on this PR. The waveform
  figures are implemented and tested on synthetic records. For the real deployment, run
  `examples/compy3_waveform_diagnostics.py` on the IPGP machine.

## 1. Executive summary

1. **The published compliance reproduces exactly (verified).** All 61 archived compliance
   files in `docs/results` are recomputed from their stored stacks with the current code,
   with a maximum relative difference of 2.2 × 10⁻¹⁶. The core estimator is correct:
   - gravity and free-air terms, with their signs;
   - the wave-attraction expression;
   - the dispersion solver;
   - the Bendat & Piersol error;
   - the S_pz convention.
2. **The main open science issue is the Z/P phase (verified).** In the compliance band,
   28 of the 29 usable stations of two networks with different operators (YV, INSU-IPGP;
   XE, Scripps) show a residual phase of +5 to +19° at the bottom of the band. It falls
   with frequency. All 21 usable XE stations have the same sign (binomial p = 10⁻⁶);
   RR40 is the one negative station.
   - **Ruled out:**
     - a timing offset (fit rms 2–7°, and the Rayleigh-band phase is 0 ± 3–8°);
     - a poroelastic seafloor (no correlation with sediment thickness).
   - **Fits (hypothesis):** a mismatch of the seismometer or DPG low-frequency corner
     (fit rms ≤ 1° at most stations).
   - **Why it matters:** if it is instrumental, a scalar Rayleigh gain fitted at 20–60 mHz
     cannot remove it. The phase alone implies a compliance bias of **about −1 % to −22 %
     at the bottom of the band** (positive at RR40), largest at the lowest frequencies,
     which carry the deep structure.
   - **Next step:** this is the highest-value item to resolve before ComPy 3 results.
3. **Steps inside stacked periods (verified).**
   - The changepoint detector finds a synchronous +6 to +10 % |Z/P| step on 10–18
     September 2013 at RR28, RR29, RR40 and RR52. Per-event Rayleigh gains show the same
     rise independently (+6 to +9 %, from only 1–3 events after the step), and no
     StationXML epoch changes.
   - RR28, RR29 and RR52 were inverted on the full period, so their stacks are high by
     +1.0, +0.9 and +0.6 %.
   - RR50's "stable" epoch still contains a +8 % step of 6 June 2013, so its stack is
     high by +2.6 %; RR40's epoch is high by +0.4 %.
4. **The coherence gate biases compliance high (verified by Monte Carlo).** Selecting
   windows on their own Z–P coherence biases the stacked H1 high, because windows are
   kept partly for noise that raised |S_pz|.
   - **Size:** +0.03 to +8 %, depending on acceptance and sea-state spread.
   - **Expected at YV**, which accepts 17–39 % of windows: about +1 to +4 % (model
     estimate). At most XE stations it is below 0.5 %.
   - **Fix (implemented, opt-in):** a per-bin "exclusive" gate that removes the bias.
5. **The calibration water column also biases compliance high (verified).** The Rayleigh
   calibration assumes a homogeneous 1028 kg/m³, 1500 m/s water column. The in-situ
   column mean at 3–5 km is about 1035–1039 kg/m³, so gains are low by 0.6–1.0 % and
   compliance high by the same amount, at every station.
6. **The error budget misses the systematics.** The inversion covariance holds the
   Bendat & Piersol error (0.1–0.3 %), the monthly scatter (0.5–0.8 %) and a 1 % floor.
   It does not hold the systematics in items 2–5, each 0.5–4 % (up to 20 % for item 2),
   mostly of one sign.
7. **New tools.** Three new modules and one figure module:
   - `compy.anomaly`: automatic, ranked anomaly detection; 120 flagged behaviours;
   - `compy.timefreq` and `compy.results`: array-level diagnostics that need no ObsPy;
   - `compy.plotting.diagnostics`: 13 diagnostic figures.

   The detector finds the documented sensitivity steps on its own (RR36, RR40, RR50) and
   the two failed XE stations (EC04, WC01), and diagnoses the latter as failed verticals.
   It also adds the findings above and two more:
   - strong tidal (M2) control of window acceptance at YV, up to 50 % modulation, which
     points to bottom currents;
   - a network-wide −1 %/yr drift of XE compliance (Spearman −0.86, p = 7 × 10⁻⁵).

## 2. Audit findings (severity-ranked)

**Severity scale:**

- **High:** can move published compliance by ≥ 1 % or invalidate a station.
- **Medium:** 0.3–1 %, or a robustness gap that could corrupt results silently.
- **Low:** code quality or documentation.

**Status:** Fixed = changed on this branch; Option = new behaviour, off by default;
Open = needs a decision or data.

| # | Sev. | Where | Finding | Evidence | Effect | Status |
|---|---|---|---|---|---|---|
| H1 | High | Stage A responses; `compy_pressure/rayleigh_calibration.py:80` (scalar `fit_multiplier`); `compy/compliance.py:158` | **Residual Z/P phase of +5–19° falling with frequency** at 28 of 29 usable stations (§5.A). A response mismatch at low frequency is the only family of hypotheses that fits. | **Verified** phase (`anomaly.phase_models`, `numbers.json`); cause is a **hypothesis** | Compliance bias at band bottom of about −1 % to −22 % (positive at RR40; station dependent), not removed by a scalar gain fitted at 20–60 mHz | Open: D1, D2 |
| H2 | High | `pipeline/c_compliance.py:190` (`--end`); YV inversion inputs | **Epochs mixed in the stacks:** steps on 2013-09-10…18 (RR28, RR29, RR40, RR52), on 2013-06-06 (RR50, inside its "stable" epoch) and on 2013-08-05 (RR36) | **Verified** on archived per-window levels; confirmed by per-event Rayleigh gains (`stability_steps.json`) | +0.6 to +1.0 % (RR28, RR29, RR52, full period); +2.6 % (RR50 epoch); +0.4 % (RR40 epoch) | Open: D3 |
| H3 | High (YV) / Low (XE) | `compy/compliance.py:130` (`SelectionCriteria.check`), `stacked_compliance` | **Selection bias of the coherence gate:** the same noise sets the per-window coherence (about 3.8 Welch averages) and \|S_pz\| | **Verified** by Monte Carlo (§4.3); station values **model-based** | +0.03 to +8 % in general; about +1 to +4 % at YV, < 0.5 % at most XE stations | Option `gate="exclusive"`: D4 |
| M1 | Med. | `compy_pressure/rayleigh_calibration.py:54, 1021–1022`; `pipeline/d_calibrate.py:60` | **Homogeneous water column (1028 kg/m³, 1500 m/s)** in the Rayleigh prediction; gain ∝ ρ | **Verified** (`water_column` in `numbers.json`) | Gain −0.6 to −1.0 % (−0.5 to −1.45 % for c = 1490–1510 m/s), so compliance +0.6 to +1.0 % at every station | Open: D5 |
| M2 | Med. | `compy/inversion/dram.py:543–563` (`ComplianceData.build`) | **The covariance omits the systematics** H1–M1. The B&P term is diagonal, although power estimates in adjacent Hann bins are correlated (ρ ≈ 0.44) | **Verified** (error budget figures) | χ²/N < 1 at most stations while systematic offsets of 1–4 % are not represented | Open: D6, roadmap |
| M3 | Med. | `compy/compliance.py:158` docstring; `pressure/docs/rayleigh_calibration_review.md` | **Polarity:** Z/P phase ≈ 0° at every station of both networks, where an elastic seafloor gives 180° for upward Z and positive compression. This agrees with BRUIT-FM's assumed pressure sign −1 and with Doran et al. (2019), whose DPGs output a negative pulse for a pressure increase | **Verified** | None on \|η\|. Breaks any phase-aware use (phase QC, ATaCR-style Z–P correction, joint inversion) | Open: D7 |
| M4 | Med. | Stage C acceptance | **Tidal control of window acceptance** (tidal lines at 7 of 8 YV stations, M2 up to 50 % modulation at RR34; S1 or O1 at CC11 and EE01) | **Verified** (empirical and tail p) | No measurable modulation of \|Z/P\| (≤ 1 %), so \|η\| is not biased. The weighting and the selection-bias regime vary with the tide | Info |
| M5 | Med. | XE monthly stacks | **Network-wide drift:** the median monthly deviation falls from +0.8 % (2018-04) to −0.7 % (2019-05); Spearman −0.86, p = 7 × 10⁻⁵; inter-station common mode p = 0.002 | **Verified**; cause is a **hypothesis** (§5.D) | About 1 %/yr. The monthly-subset term treats it as random scatter | Open: test |
| M6 | Med. | `pipeline/b_clean.py:215`; `compy/preprocessing.py:195–215` | **Silent fallback:** when TiSKitPy refuses the event spans, coherent-noise transfer functions are estimated *with* the events, and nothing is recorded unless both attempts fail | Code reading | Possible earthquake leakage into the transfer functions of affected 24-h windows | Open: D8 |
| M7 | Med. | `compy/preprocessing.py` (24-h DataCleaner windows) | **Cleaning-window edges** (hypothesis): FFT-domain transfer functions change at window edges, which could depress coherence of 35-min windows that straddle them. S1 clustering of acceptance at several XE stations is compatible | **Hypothesis**; test in §5.C | Fewer accepted windows near boundaries; possible local bias | Open: test |
| L1 | Low | `compy/compliance.py` (2.8.0 `stacked_compliance`) | The frequency axis came from the last window, so windows of different sampling rate were mixed silently; custom rejection reasons raised `KeyError` | Code reading, test | Silent corruption on mixed inputs | **Fixed** |
| L2 | Low | `compy/streams.py:4` | ObsPy imported at module load, so array-level compliance needed ObsPy | Verified | Testability | **Fixed** |
| L3 | Low | `compy/cleaning.py:170, 184–187` | `rotate_smoothly` median-filters over `smooth_hours` *samples*: skipped earthquake hours stretch the 24-h window. Coefficients are evaluated at integer knots but placed at float knots | Code reading | Negligible (smoothing semantics) | Open (doc) |
| L4 | Low | `compy/cleaning.py:53` | Fold index can equal `n_bins` through rounding | Code reading | Harmless | Open |
| L5 | Low | `compy/plotting/spectra.py:1008–1105` | Legacy `plot_spectrogram` titles the BDH panel "BHZ", indexes channels by position (`st[3]`) and uses `jet` | Code reading | Mislabelled legacy figure | Open (deprecate) |
| L6 | Low | `compy/compliance.py:733` | `gravitational_attraction(pw=1025)` vs 1028 elsewhere (affects only the returned `h_w`) | Code reading | Cosmetic | Open |
| L7 | Low | `tests/` | 72 tests need ObsPy, ArviZ or TiSKitPy, and no test ran the core estimator without the full stack | Verified | Testability | **Improved** (array path) |

**Verified correct (no change needed):**

- the gravity-correction form k(ω²|Z/P| + a_w/p)/(ω² + 2g/R) and its signs;
- the wave attraction a_w/p = 2πG/g · (1 + e^{−2kH})/2;
- the free-air constant 3.07 × 10⁻⁶ s⁻²;
- the bounded dispersion solver;
- the B&P amplitude-coherence error;
- `csd(p, z)` = E[P* Z];
- Welch effective averages;
- the gain algebra of the sampled-gain likelihood, η_meas = η_a + g(η − η_a).

The band search is bounded by f_c = √(g/2πH). That is a deep-water wavelength equal to H
(kH = 2π, pressure attenuated by 1/cosh 2π ≈ −49 dB) and is a safe upper limit; the
coherence then sets the band.

## 3. New diagnostic and insight plots

All figures live in `compy/plotting/diagnostics.py`. They take arrays or
`compy.results.StationProducts`, need no ObsPy, and use a fixed categorical palette, one
blue ramp for magnitude and blue–grey–red for signed deviations. The examples were made
by `examples/compy3_diagnostics.py`, except the waveform rows, which were made by
`examples/compy3_waveform_diagnostics.py` on real data or by the synthetic demo here.

| Figure (function) | Purpose | Anomaly it reveals | Example |
|---|---|---|---|
| Spectrograms (`plot_spectrogram`) | Z, P, H1, H2 power over the deployment | Storms, gaps, glitch lines, sensor changes | `figures/synthetic_spectrogram.png` |
| Coherograms (`plot_coherogram`) | Z–P, Z–H1/H2, P–H coherence vs time | When compliance or tilt dominates. P–H coherence means horizontals carry pressure signal, which coherent-noise removal would also take out of Z | `figures/synthetic_coherogram.png` |
| Tilt (`plot_tilt` + `timefreq.tilt_direction`) | Tilt azimuth and admittance vs time, folded on M2 | Fixed vs current-driven tilt, re-levelling | `figures/synthetic_tilt.png` (azimuth recovered: 72.5° for a true 70°) |
| Error budget (`plot_compliance_confidence`) | Compliance with nested ±1σ/±2σ bands: B&P, monthly SEM, floor | How small the formal error is next to the scatter | `figures/RR52_error_budget.png`, `RR50_error_budget.png` |
| Stability (`plot_compliance_stability`) | Month × frequency deviation map and band-mean time series with window counts | Steps, drifts, frequency-dependent changes (flat means gain; low-frequency only means tilt or DPG leak) | `figures/RR52_monthly_stability.png`, `RR50_…`, `CC07_…` |
| Noise models (`plot_noise_models`) | Stacked PSD of accepted windows vs Peterson NLNM/NHNM or Webb (1998) | Gain, units or failed channels (EC04, WC01 are 36–48 dB low) | `figures/noise_vertical.png`, `noise_pressure.png` |
| Network panels (`plot_multistation`) | Small multiples of η ± total σ and coherence | Stations whose shape or coherence departs | `figures/YV_multistation.png`, `XE_multistation.png` |
| Phase diagnostic (`plot_phase_diagnostic`) | Residual Z/P phase with response-model fits, and the implied bias range | H1 | `figures/phase_diagnostic.png` |
| Window levels (`plot_window_levels`) | Per-window \|Z/P\|, P and Z levels with detected steps | Which sensor changed: a Z-level step with a \|Z/P\| step means the seismometer | `figures/RR50_window_levels.png`, `RR52_…`, `RR28_…` |
| Tidal maps (`plot_periodicity`) | Station × constituent modulation, with significance | Tidal currents, diurnal processing artefacts | `figures/tidal_acceptance_map.png`, `YV_tidal_transfer_function.png` |
| Folded acceptance (`plot_acceptance_phase`) | Accepted windows folded on M2 or 24 h | Gate failures at a tidal phase | `figures/YV_acceptance_M2.png`, `YV_acceptance_S1.png` |
| Common mode (`plot_common_mode`) | Monthly deviations by station and their correlation | Network-wide instrument drift vs station-specific changes | `figures/YV_common_mode.png`, `XE_common_mode.png` |
| Ranking (`plot_anomaly_ranking`) | Ranked anomalies by category | Triage | `figures/anomaly_ranking.png` |

**Research basis:**

- ATaCR/OBStools produce daily QC, coherence, phase and admittance vs date, and tilt from
  the azimuth of maximum Z–H coherence (Janiszewski et al. 2019; Bell et al. 2015).
- Corela et al. (2023) show current-driven tilt noise modulated by tides.
- Doran et al. (2019) show DPG gains off by up to a factor of two, time constants off by
  more than 2×, and inverted polarity.
- Janiszewski et al. (2023) compile OBS noise properties.

`timefreq.tilt_direction` takes the principal direction of (S_1Z, S_2Z) instead of the
maximum coherence. With quiet horizontals the coherence is flat in azimuth, so the
maximum-coherence direction is noise-dominated; the principal direction is not (tested).

## 4. Anomaly and correlation detection system

### 4.1 Design (`compy/anomaly.py`)

**Detectors:**

- **Robust scale:** median/MAD z-scores everywhere.
- **Changepoints:** binary segmentation, accepting a split when the median shift exceeds
  5σ after deflating for lag-1 autocorrelation. The tests check that no changepoint is
  found in AR(1) noise.
- **Periodicity:** a sinusoid fit at M2, S2, N2, K1, O1, M4 and S1, compared with about
  250 non-tidal control periods (8–40 h). The p-values are empirical, plus a tail
  extension that scales an exponential law for amplitude² to the controls, so
  storm-driven red variability does not register as tidal. The same is done for event
  times with the Rayleigh resultant.
- **Level dependence and regimes:**
  - binned |Z/P| against pressure level;
  - k-means regimes, kept only when concentrated in time (KS ≥ 0.3), so the tail of one
    distribution is not reported as a cluster.
- **Network detectors:**
  - common mode of monthly deviations (permutation p);
  - synchronous steps (≥ 3 stations within 15 days);
  - PSD-level outliers in the compliance and microseism bands;
  - systematic phase sign.
- **Spectra and inversion:**
  - response-model fits to the residual phase, with the compliance bias they imply;
  - epoch mixing;
  - inversion gain pull and poor fit.

**Score:** min(σ, 8) × (1 + min(|effect %|, 20)/5). It is transparent and favours strong
evidence of a large effect on compliance.

### 4.2 Results

The ranked table is `anomalies.md` (or `.csv`), with 120 records at σ ≥ 2; the chart is
`figures/anomaly_ranking.png`.

- **Known issues found without supervision:**
  - the sensitivity steps of RR36 (2013-10-21, +38 %), RR40 and RR50 (2013-07-30, +58 %);
  - RR38's pressure-level dependence (−18.5 % per 10 dB) and its two regimes, consistent
    with its exclusion;
  - EC04 and WC01 (rejected for having fewer than 500 windows), now with a cause: verticals
    36–48 dB below the network in the microseism band, and Z/P phase −65°/−81°.
- **New:** the items of §1 (H1, H2, M4, M5) and CC07's −1.55 %/month drift (2019-05 is
  −22 %).

### 4.3 Coherence-gate selection bias (Monte Carlo, `anomaly.simulate_gate_bias`)

**Setup:**

- White pressure; vertical = pressure + incoherent noise.
- 2625-sample windows, nseg 1024, 9-bin band, gate on the band-median coherence ≥ 0.9.
- Log-normal window sea state with spread σ (natural log of amplitude).
- Three seeds; standard deviation 0.1–0.3 %; the exclusive-gate scatter at the
  0.8/σ = 0 case is ±0.65 %.

| Median coherence | σ | Acceptance | Inclusive bias | Exclusive bias |
|---|---|---|---|---|
| 0.80 | 0.0 | 0.15 | **+7.9 %** | +0.7 % |
| 0.80 | 0.5 | 0.32 | +0.8 % | −0.1 % |
| 0.85 | 0.0 | 0.40 | **+3.9 %** | +0.3 % |
| 0.85 | 0.5 | 0.46 | +0.5 % | −0.05 % |
| 0.90 | 0.0 | 0.80 | +1.1 % | +0.1 % |
| 0.90 | 0.5 | 0.63 | +0.3 % | 0.0 % |
| 0.95 | any | ≥ 0.71 | ≤ 0.1 % | ≤ 0.1 % |

**At each station's acceptance** (median coherence matched by bisection):

- **YV:** +1.4 to +1.9 % at σ = 0.5, and +1.6 to +4.1 % at the observed spread of
  accepted pressure levels (σ ≈ 0.28, a lower bound because selection narrows it).
- **XE:** +0.04 to +0.5 % at most stations.

The model is idealized, so these values are the order of magnitude to confirm with stage C
`COMPLIANCE_GATE = "exclusive"` on the IPGP machine.

## 5. Physical interpretation and discriminating tests

### A. Residual Z/P phase (H1)

The table lists the residual at the band edges, the rms of the timing fit and of the two
response-model fits, and the compliance bias at the band bottom allowed by acceptable fits
of each model (from `numbers.json`). XE stations show the same pattern, with
bottom-of-band residuals of +6 to +14°.

| Station | Band [mHz] | Residual [°] | Timing rms | DPG / seis. rms | DPG-model bias | Seis.-model bias |
|---|---|---|---|---|---|---|
| RR28 | 7.3–13.4 | +16.3 → +8.7 | 4.9 | 0.26 / 0.27 | −10 to −4 % | −14 to −7 % |
| RR29 | 7.3–12.2 | +19.3 → +13.9 | 4.8 | 0.14 / 0.40 | −19 to −12 % | −17 to −16 % |
| RR34 | 7.3–13.4 | +9.2 → +4.9 | 2.8 | 0.22 / 0.31 | −9 to −1 % | −10 to −4 % |
| RR36 | 7.3–15.9 | +8.2 → −1.9 | 3.8 | 2.4 / 1.5 | −9 to 0 % | +6 to +7 % |
| RR40 | 8.5–13.4 | −12.6 → −8.4 | 3.1 | 0.18 / 0.32 | +3 to +13 % | 0 to +6 % |
| RR50 | 7.3–13.4 | +15.2 → +10.9 | 4.0 | 0.15 / 0.64 | −18 to −12 % | −15 to −14 % |
| RR52 | 7.3–17.1 | +17.9 → +14.5 | 4.8 | 1.4 / 2.5 | −22 to −16 % | −20 to −10 % |
| CC11 | 4.9–14.7 | +11.9 → +0.5 | 4.9 | 1.6 / 1.2 | −10 to −1 % | −10 to +3 % |

| Hypothesis | Prediction | Test | Result |
|---|---|---|---|
| (a) Timing offset between Z and P | Phase = 360 f Δt, through zero, growing with f | Fit in the band; Rayleigh-band phase at 30–70 mHz | **Rejected:** rms 2–5° vs ≤ 1° for the response models; Rayleigh P/a phase rms 2.4–8.3° about 0 (BRUIT-FM `stations.csv`), where the fitted Δt of +0.5 to +4.5 s would give 10–110° |
| (b) DPG leak time constant ≠ metadata (first-order high-pass) | atan(f/f_true) − atan(f/f_meta), decaying above the corners | Fit (done); tidal-band amplitude vs a tide model gives τ directly (Doran et al. 2019) | **Fits** (rms ≤ 1°). Implied f_true < f_meta, so τ_true > τ_meta. Needs the M2 test on raw pressure |
| (c) Seismometer corner ≠ nominal 240 s | Second-order shape around 4 mHz | Fit (done); calibration pulses or step response; compare sensors of the same model | **Fits** (true period ≈ 180–220 s at most stations). Needs a calibration record |
| (d) Poroelastic or anelastic seafloor | Phase grows with sediment thickness (drained/undrained) and is site-dependent | Spearman of residual vs posterior sediment thickness | **Not supported:** ρ = −0.26 (YV, n = 8), −0.23 (XE, n = 21); wrong sign, not significant. Same pattern on two different crusts |
| (e) Coherent-noise removal | Quadrature leakage αc (α ≈ 0.3 horizontal/vertical, c = residual tilt) | Analytic: phase ≈ atan(αc) < 1° | **Rejected** as the main cause |

**Discriminating test still to run.** The two response hypotheses predict different
amplitude errors.

1. Measure the DPG τ from the M2 tidal line, comparing measured pressure with TPXO or FES
   tide.
2. Extend the Rayleigh calibration down to 10–20 mHz with M ≥ 7.5 events. Each acceptable
   fit predicts a P/a phase and amplitude there; the per-bin outputs of `phase_models`
   give the predictions.
3. If (b) or (c) holds, deconvolve the corrected response before stage C (Decision D2).

### B. Synchronous September 2013 steps (H2)

| Hypothesis | Test | Result |
|---|---|---|
| StationXML epoch change | One epoch per channel for the whole deployment (BRUIT-FM `stations.csv`) | Rejected |
| Sea-state change (non-instrumental) | A pressure-level step should come with it, and compliance is a ratio | Rejected: "no concurrent P/Z level step" at RR28, RR29, RR40, RR52 |
| DPG sensitivity change (e.g. battery voltage late in the deployment, common hardware) | Per-event Rayleigh gains before and after; P/Z ratio only | **Supported:** Rayleigh g rises by +7.8 % (RR28), +6.1 % (RR29), +6.3 % (RR40), +8.5 % (RR52) after the step, from 1–3 events each |
| Seismometer sensitivity change | Z-level step, as at RR50 on 2013-07-30 (+4.7 dB) | Not seen in September. RR50's July step is a vertical-channel change |

To separate DPG from seismometer in September, compare the microseism-band Z level of each
station with its neighbours before and after the step. Per-window microseism levels are
not archived; `examples/compy3_waveform_diagnostics.py` produces them.

### C. Tidal control of acceptance (M4)

| Hypothesis | Test | Result |
|---|---|---|
| Current-driven tilt noise lowers Z–P coherence at peak flow | Modulation at M2 alone (mean + tidal current), or M4 (purely tidal speed) | **Supported:** M2 dominates (RR34 50 %, M4 15 %), as expected for a mean bottom current plus tide (Corela et al. 2023) |
| DPG nonlinearity under tidal pressure | \|Z/P\| itself modulated at M2 | **Rejected:** \|Z/P\| modulation ≤ 1 %, not significant at M2 |
| Processing artefact at 24-h cleaning boundaries (M7) | Deficit concentrated at boundary times (S1 at XE) | Open: needs segment start times (`work/<STA>/clean/summary.json`). Fold acceptance on the distance to the nearest boundary |

### D. XE network drift (M5)

| Hypothesis | Test | Result |
|---|---|---|
| Common instrument drift (DPG batch ageing) | Monotonic trend; per-event Rayleigh gains vs time | Trend: Spearman −0.86. Gains vs time not archived (run stage D output) |
| Selection-bias modulation (acceptance changes) | Correlation with monthly acceptance | Partial: ρ = −0.53 (p = 0.05). But 2018-04 and 2019-05 both have low acceptance and opposite deviations, which **favours drift** |
| Seasonal infragravity directionality (2-D or 3-D effects) | Annual cycle rather than a trend; frequency dependence | 13 months cannot separate cycle from trend; deviations are flat in frequency, consistent with gain |

### E. RR38, CC07, EC04 and WC01

- **RR38:**
  - |Z/P| falls by 18.5 % per 10 dB of pressure;
  - there are two time regimes;
  - the phase residual is non-monotonic, and no response model fits it (rms 9–10°).

  This points to pressure-channel noise or a nonlinearity, or to model error, as already
  concluded. **Keep it excluded.**
- **CC07:** a −1.55 %/month drift. Test it with per-event Rayleigh gains vs time; if they
  drift too, the cause is instrumental, otherwise settling or coupling.
- **EC04 and WC01:** the verticals are 36–48 dB low in the microseism band and the Z/P
  phase is far from 0 or 180°, so the vertical channel failed or its response is wrong.

## 6. Implemented changes

| Commit | Change | Numbers changed? |
|---|---|---|
| `Import ObsPy lazily in compy.streams` | `compy.compliance` imports without ObsPy (L2) | No |
| `Stack compliance through a streaming core…` | `_ComplianceStack`, `stacked_compliance_from_spectra`, `gate="exclusive"` (H3), checks on mixed frequency axes and custom rejection reasons (L1), `COMPLIANCE_GATE` in `pipeline/c_compliance.py` (default `"inclusive"`) | **No by default** (bit-identical); yes if D4 is adopted |
| `Add array-level time-frequency helpers…` | `compy.timefreq`, `compy.noise_models`, `compy.results` | No |
| `Add an anomaly and correlation detection module` | `compy.anomaly` | No |
| `Add diagnostic figures and tests…` | `compy.plotting.diagnostics`, two test files (26 tests) | No |
| `Add ComPy 3 audit scripts and outputs` | `examples/compy3_diagnostics.py`, `examples/compy3_waveform_diagnostics.py`, this folder | No |

## 7. Validation (before/after)

- **Published results unchanged:**
  - `test_every_published_compliance_is_reproduced_from_its_stored_stack` recomputes all
    61 archived files: max relative difference 2.2 × 10⁻¹⁶, coherence identical;
  - `test_default_gate_reproduces_the_280_estimator_bit_for_bit` compares every output
    field of the refactored estimator with a frozen copy of the 2.8.0 loop using
    `np.array_equal`, with and without an output band;
  - `test_stacked_compliance_on_streams_matches_the_array_path` does the same through
    ObsPy streams (runs in CI).
- **New behaviour validated:**
  - the exclusive gate equals the inclusive one when every window passes, and removes the
    bias in the Monte Carlo;
  - detectors recover synthetic steps, tides and response mismatches, with false-positive
    controls;
  - figures render.
- **Run here:**
  - all 26 new tests pass (the ObsPy one is skipped);
  - existing tests: 219 passed before the changes and 242 after, because the lazy import
    lets 23 more run without ObsPy;
  - the remaining 49 failures and 8 modules that fail to import are all
    `ModuleNotFoundError` for ObsPy or ArviZ (checked one by one), so they run in CI;
  - `ruff check` and `ruff format --check` pass (ruff 0.16.9, as in CI).
- **Not validated here:** the waveform script on real MiniSEED (no data or ObsPy here),
  and a stage C rerun with `gate="exclusive"`.

## 8. ComPy 3 roadmap

Ranked by scientific impact over effort (I = impact, E = effort; 1 low, 5 high).

| Rank | Item | I | E | Why |
|---|---|---|---|---|
| 1 | **Frequency-dependent instrument calibration:** fit DPG τ (tides) and seismometer corner (Rayleigh 10–60 mHz, phase + amplitude) and deconvolve before stage C; carry the response uncertainty into the covariance | 5 | 3 | Resolves H1, the largest potential bias, at the lowest frequencies |
| 2 | **Automatic epoch segmentation:** changepoints on per-window \|Z/P\| plus per-event gains, then `ComplianceData.join` with one gain per epoch | 5 | 2 | Removes H2; the detector already exists |
| 3 | **Systematic-error covariance:** add terms for gate bias, water column, response and epoch uncertainty; correlated B&P across Hann bins; jackknife over days | 4 | 2 | Honest error bars; addresses M2 |
| 4 | **Unbiased window selection by default** (exclusive gate, or gating on pressure level) | 4 | 1 | H3; implemented as an option |
| 5 | **In-situ water column in the calibration:** depth-mean ρ (compressibility) and sound speed (Mackenzie with WOA profiles) | 3 | 1 | M1, at every station |
| 6 | **Automated QC report per station** (`compy qc`): this anomaly scan plus waveform coherograms and tilt, run after stage C | 3 | 2 | Catches failed channels and steps before inversion |
| 7 | **CLI and config files:** `compy run config.yaml` (stages A–G) replacing `exec` of a Python config; schema-checked YAML, run manifest with config hash | 3 | 3 | Reproducibility, safety |
| 8 | **FDSN fetch with caching** (EIDA/RESIF, IRIS) and StationXML sanity checks (unpaired poles, positive real poles, polarity) | 3 | 2 | New networks without an SDS archive; BRUIT-FM already found invalid YV responses |
| 9 | **Joint observables:** compliance + Rayleigh/Scholte dispersion or receiver functions, using the existing `Forward`, `disba` and gain sampling | 4 | 4 | Breaks the Vs–Vp/Vs and sediment trade-offs |
| 10 | **Packaging:** split `compy-core` (NumPy/SciPy only) from `compy[obspy]`; numpydoc site (Sphinx or MkDocs) with these figures as a gallery | 2 | 2 | Light install, docs |
| 11 | **CI:** a no-ObsPy job, a weekly full-pipeline smoke test on a 3-day public dataset, and coverage | 2 | 1 | Catches import and environment drift |
| 12 | **Retire legacy figure code** (`plotting/spectra.py` legacy layouts, L5) behind `compy.legacy` | 1 | 1 | Maintenance |

## 9. Decisions

Each line is a yes/no question. "Results change" means published numbers would move.

- **D1.** Treat the residual Z/P phase (H1) as an instrument-response problem to be
  resolved before new inversions? (Y/N) Impact: possible frequency-dependent compliance
  bias of about −1 to −22 % at the band bottom; ComPy 3 results wait for the calibration test
  (roadmap 1).
- **D2.** Implement frequency-dependent DPG/seismometer response fitting (τ from tides,
  corner from low-frequency Rayleigh) and apply it in stage A or C? (Y/N) Impact: results
  change, deep Vs most.
- **D3.** Re-stack RR28, RR29 and RR52 before 2013-09-10 and RR50 before 2013-06-06 (or
  join the epochs with separate gains) and re-invert? (Y/N) Impact: compliance −0.6 to
  −2.6 % at those stations; results change.
- **D4.** Make `gate="exclusive"` the default for new runs (keeping `"inclusive"` to
  reproduce the paper)? (Y/N) Impact: about −1 to −4 % compliance at YV, < 0.5 % at most
  XE stations; results change.
- **D5.** Use a depth-mean in-situ water density (and sound speed) in the Rayleigh
  calibration? (Y/N) Impact: gains +0.6 to +1.0 %, compliance −0.6 to −1.0 % at every
  station; results change.
- **D6.** Add systematic-error terms (gate, water column, response, epoch) to
  `ComplianceData` and drop or reduce the ad hoc 1 % floor accordingly? (Y/N) Impact:
  wider, better-founded posteriors.
- **D7.** Adopt an explicit pressure-polarity convention (sign recorded per station,
  checked from the Rayleigh or compliance phase) and report the compliance phase relative
  to 180°? (Y/N) Impact: none on |η|; enables phase QC.
- **D8.** Record the coherent-noise-removal fallback (event spans refused) in the stage B
  summary and warn? (Y/N) Impact: traceability only.
- **D9.** Add `compy.anomaly.scan_station` as an automatic QC step after stage C, failing
  the run on "phase_inconsistent" or level outliers? (Y/N) Impact: stations like EC04 and
  WC01 are stopped before inversion.
- **D10.** Merge this branch as ComPy 2.9.0 (new modules, defaults unchanged) before the
  ComPy 3 work? (Y/N)

## References

- Bell, S. W., Forsyth, D. W., & Ruan, Y. (2015). Removing noise from the vertical
  component records of ocean-bottom seismometers. *BSSA* 105, 300–313.
- Corela, C., et al. (2023). The effect of deep ocean currents on ocean-bottom
  seismometers records. *NHESS* 23, 1433–1451.
  <https://nhess.copernicus.org/articles/23/1433/2023/>
- Doran, A. K., Rapa, M., Laske, G., Babcock, J., & McPeak, S. (2019). Calibration of
  differential pressure gauges through in situ testing. *Earth Space Sci.* 6.
  <https://doi.org/10.1029/2019EA000783>
- Janiszewski, H. A., Gaherty, J. B., Abers, G. A., Gao, H., & Eilon, Z. C. (2019).
  Amphibious surface-wave phase-velocity measurements… (ATaCR). *GJI* 217, 1929–1948.
  <https://github.com/helenjanisz/ATaCR>
- Janiszewski, H. A., et al. (2023). Broad-band ocean bottom seismometer noise properties.
  *GJI* 233, 297–315.
- OBStools ATaCR documentation. <https://nfsi-canada.github.io/OBStools/atacr.html>
- Peterson, J. (1993). Observations and modeling of seismic background noise. USGS OFR
  93-322.
- Webb, S. C. (1998). Broadband seismology and noise under the ocean. *Rev. Geophys.* 36,
  105–142.
- The references of `docs/AUDIT_2026-10.md` and `docs/INVERSION.md` (Crawford et al.
  1991, 1998; Zha & Webb 2016; Bendat & Piersol; Aminian et al. 2025).
