# Local-P scientific review: noise windows, timing, pre-multiple reference

Reviewed 2026-10-16 by Claude on branch `claude/local-p-science-review`
(based on Codex's review `858a715`, which reviewed `f455edd`). It responds to
[local_p_codex_review.md](local_p_codex_review.md) findings 1 and 4 and its
physics section. No threshold, rule, gain, response, raw waveform, master or
saved result was changed; no campaign was rerun; nothing was merged.

**Evidence base.** Primary papers read: Deng et al. (2022, GRL 49,
e2021GL096952), An et al. (2017, GRL 44, 10,272–10,281), Zha & Webb (2016,
JGR 121, 3220–3237), Doran et al. (2019, Earth Space Sci. 6, key points), Saito (2017, GJI 210, 1888–1900; online full text, first
100k characters), Auger (2024, M1 report), Aminian et al. (2025, GJI 242).
Saved results: `all_v1` and `apg_timing_v2` campaigns. New bounded analysis:
`scripts/review_local_p_science.py`, run `results/LOCAL_P/science_review/run2/`
(nine LT02 waveform windows read; one synthetic; all 18 input hashes unchanged).

## Summary

| Question | Finding | Recommendation |
|---|---|---|
| LT02 noise-window rule | Codex is right: the documented rule is not implemented; 8/9 LT02 and 103/150 campaign-wide accepted events violate it. The two documented clauses are mutually inconsistent when the onset is near predicted P. The violation does **not** measurably change LT02 gains (station medians within 0.5%). | Treat LT02 as not yet passed. Resolve the contradiction by *placement* (noise ends at the earlier of the two limits), decided on principle before re-evaluation; then a bounded LT02 rerun. |
| Timing uncertainty | LT02 APG: about 0.10 s (1σ) after drift model + refinement; same-logger XO DPGs: about 3 ms; YV: a consistent +12–14 ms offset. Timing barely affects \|g\| but dominates phase. | Keep LT02 \|g\|; do not interpret LT02 phase/polarity above ~1 Hz. Record a timing budget and held-out drift test with every lag decision. |
| Pre-multiple as calibration reference | Theory exact before the first multiple; literature supports only "first motion"; our synthetic leakage < 0.3%; field values precise (WD70 ±6%) but with marginal event dependence and **no absolute reference**. | Keep diagnostic-only. Promotion needs an independent deep-water absolute reference at the same frequencies plus the criteria in §3.4. |

## 1. LT02 noise-window rule

### 1.1 What was promised and what runs

The method note says (§4 step 3): the noise window "has the same length,
ends 1 s before the signal window, and is required to precede the predicted P
minus its uncertainty". The code places it to end 1.2 s before the measured
vertical onset and never applies the second clause (Codex finding 1).

The two clauses cannot both hold whenever the measured onset lies within the
search half-width (≥ 3 s) of predicted P, which is the normal case. The
reviewer reproduces Codex exactly: 8/9 accepted LT02 events violate the
documented boundary by **0.917–4.848 s**. Campaign-wide, **103 of 150**
accepted events in `all_v1` violate it (WD70 39/41, LA33 EDH 19/32, LA34 EDH
20/27, LT02 6/7). Enforcing the rule by *rejection* would remove nearly all
results, including the LT02 candidate.

### 1.2 What the rule protects against, and whether it occurred

The clause guards against a **late** pick: if AIC picked a later phase, the
"noise" window would contain the true first P. At LT02 the measured onsets are
−0.94, −0.29, −0.95, −0.63, −4.89, +0.85, +0.53, −1.17, −0.61 s relative to
iasp91 predictions. Seven are early; the two late ones (+0.53, +0.85 s) are
well inside the model error. There is no evidence of late picks.

Empirical contamination test (0.5–2 Hz, multitaper PSD ratios):

| Ratio | Pressure median (range) | Vertical median (range) |
|---|---|---|
| near-onset / 30 s earlier | 1.32 (0.72–2.47) | 1.28 (0.71–2.69) |
| control: earlier / earlier still | 0.84 (0.32–1.86) | 1.04 (0.69–1.66) |

Natural background variability alone spans 0.3–1.9×. The near-onset window
is modestly noisier on average; only `ak018gl8yvfy` (both channels) and
`ak018ahdbyga` (pressure) exceed the control range. The latter already
satisfies the documented rule, so the rule does not guarantee a quiet window.

### 1.3 Does it change the gain?

Gains recomputed with the earlier background window as noise:

| Band | Saved median (n) | Earlier-noise median (n) | Per-event change |
|---|---|---|---|
| 0.5–1 Hz | 0.972 (5) | 0.967 (4) | −1.0 … +1.0% |
| 1–2 Hz | 1.021 (9) | 1.017 (6) | four within ±1.4%; +6.9% and −9.0% |

The saved results reproduce within 0.8% when the reviewer's longer segment is
used with the original windows; one event's 1–2 Hz band changed QC status,
so pass/fail sits near thresholds for some bins. The LT02 station result is
insensitive to the noise window at the 0.5% level, while single events are
not. This matches the estimator: in valid bins (signal/noise power ≥ 10) a
noise error moves the noise-corrected ratio by at most a few percent, and the
H1/H2 brackets use the signal window only.

### 1.4 Recommendation

1. Agree with Codex: the LT02 candidate has **not** passed the documented
   protocol. Do not relax the rule silently and do not discard events.
2. Resolve the internal contradiction on principle, recorded before any
   re-evaluation. The reviewer's suggestion is to read "ends 1 s before the
   signal window" as "ends **at least** 1 s before" and place the window to
   end at the **earlier** of (onset − 1.2 s) and (predicted P − search
   half-width), keeping its length. This satisfies both documented intents
   (no P in the noise; same length) without rejecting events, and changes no
   threshold. Rejection-based enforcement is the alternative; it removes 69%
   of accepted events and is not motivated by any observed late pick.
3. Add the near/earlier noise ratio as a reported diagnostic (not a
   threshold), with the control ratio for scale.
4. Then run a bounded LT02-only check. §1.3 predicts the candidate survives
   within ~1%; that prediction must be confirmed, not assumed.

## 2. Timing uncertainty

### 2.1 Sources

- **Clock drift (LT02 only).** Deng et al. (2022, §3, Fig. 3) found the LT02
  APG ahead of the seismometer by 93–103 s and aligned it by correlating
  pressure with ρha at 0.6–0.8 Hz. Our drift line reproduces their offsets
  within 0.3 s (`local_p_check.md` §2.1), but it was fitted in-sample.
- **Held-out drift test (new).** Fitting the first 8 of 17 refined lags and
  predicting the other 9 gives residuals of −0.21…+0.007 s (RMS 0.14 s) plus
  one gross 1.08 s outlier (the excluded 2018-12-03 event); the forward
  residuals drift to about −0.2 s after three months. Leave-one-out RMS is
  0.10 s. The per-event onset refinement in `apg_timing_v2` moved lags by
  −0.12…+0.11 s.
- **Onset picks.** AIC onsets in a causal 2–10 Hz band are sample-precise
  only at high SNR; the deep APG onset SNR is often < 3 (LA33/LA34).

### 2.2 Residual delay from the saved cross-spectra (new)

A delay τ and constant phase φ₀ were fitted to the observed-minus-predicted
phase of each accepted event (valid bins only):

| Pair | Events | Median τ (range) | Median φ₀ | Reading |
|---|---|---|---|---|
| XO.LT02 HDH (APG, drift model) | 8 | +0.106 s (−0.021…+0.325) | −6° | residual timing ≈ 0.1 s; no constant phase error |
| XO.WD66 HDH (DPG, same Q330) | 10 | −0.002 s | +1° | timing essentially exact |
| XO.WD70 HDH (DPG, same Q330) | 40 | −0.002 s | +6° | timing essentially exact |
| YV.RR36 BDH | 10 | +0.014 s | −164° | reversed polarity; consistent ~1-sample delay |
| YV.RR40 BDH | 6 | +0.012 s | −156° | as RR36 |

The 33–39° LT02 band phases in Codex's review are explained by τ ≈ 0.1 s,
not by a phase or polarity error: after removing τ, φ₀ ≈ 0. Narrow frequency
spans make τ and φ₀ trade off for a few events.

### 2.3 Consequences

- **Amplitude.** The primary estimator uses auto-spectra and is unaffected by
  a pure delay; the H1 bound decorrelates by about τ/T (0.1 s / 10–20 s, < 1%).
  LT02 |g| stands, within §1.
- **Phase/polarity.** The error is 360°·f·τ: 36° at 1 Hz and 72° at 2 Hz for
  τ = 0.1 s. With the ±45° polarity tolerance, LT02 polarity is uninformative
  above ~1.25 Hz, consistent with its indeterminate band votes (3/9). Codex
  finding 4 (tie handling) matters most here.
- **Same-logger pairs** (WD66, WD70) give trustworthy phase and polarity.
- **YV** shows a consistent +12–14 ms pressure delay and ~20° residual phase
  beyond 180°. That points to a sub-sample timing or response-phase
  (FIR/group delay) mismatch between BDH and BHZ metadata. Polarity −1 is robust.

### 2.4 Recommendation

1. Report LT02 |g| with its timing-insensitive estimator; label LT02 phase
   and polarity "timing-limited" above about 1 Hz.
2. Every lag decision in `local_p_timing.json` should carry σ_τ (here
   0.10–0.14 s), a time-split held-out test, and the derived frequency limit
   for phase interpretation, f ≤ 45°/(360°·σ_τ).
3. Check the YV BDH/BHZ response group delays (a metadata question, not a
   correction).

## 3. Deep-water pre-multiple relation as a calibration reference

### 3.1 Theory and literature

- For an upgoing plane wave in water, p = ρc₀v_up/η exactly; it applies at
  the seafloor until the first sea-surface multiple returns at 2hη/c₀.
- Deng et al. (2022, §2, p. 3) state that p = ρvc₀ "can be derived by only
  considering the first reflection and transmission at the solid-fluid
  interface", "does not reduce to p = ρha at low frequencies", and is "only
  consistent with numerical simulations for the first crest of seismic
  records (Saito, 2017)". They use Eq. 5, not this relation, for calibration.
- Saito (2017, §3.4 and §5): "the relation p_e = ρ₀c₀v_z holds only at the
  first motion"; later waves are sea-surface and bottom reflections. With a
  soft sediment layer (§3.5) the first-wave agreement persists. Comparisons
  are visual, with no error statistics, and concern seafloor deformation
  (tsunami generation), not incident teleseismic or local P.
- Auger (2024, pp. 8–11) compared first-arrival p/v with this relation only
  after fitting water incidence angles of 45° and 85°, impossible for P
  waves (sin i₂ = qc₀ ≤ 0.3), with an unexplained π shift and scalar
  sensitivities. That is not supporting evidence.
- No source read applies the relation as a field calibration reference. Its
  support is theoretical and numerical only.

### 3.2 What this branch's evidence adds

1. **Numerical correctness.** The independent finite-difference reference
   matches within 0.2% before the multiple (branch tests).
2. **Gate leakage (new, answers Codex's request).** Synthetic deep-water
   record (4940 m) with the real XO.WD66 HDH/HHZ responses, a 0.8 gain, and
   post-gate bursts 5× the P amplitude, 0.3 s after the gate. The bursts were
   S-like (η = 0.9), sign-reversed, or vertical-only.

   | Window | Bias (all bands, all injections) |
   |---|---|
   | campaign gate, 5.72 s | +0.15 … +0.28% |
   | gate − 1 s | +0.07 … +0.29% |
   | past the first multiple, 8.59 s | −0.3 … −12.9%, or no valid band |

   Leakage from filtering and response conversion is negligible for this
   configuration; the gate itself is what keeps the relation valid. Limits:
   ideal plane waves, identical responses in synthesis and analysis, one
   depth. It tests numerical leakage, not response or wavefield error.
3. **Field precision.** WD70 (41 events): 0.774 (1–2 Hz, n = 15) and 0.769
   (2–4 Hz, n = 32), relative MAD ≈ 6%, τ ≈ 3 ms, φ₀ ≈ 0. WD66 agrees (0.78).
   This is consistent behaviour, not accuracy.
4. **Event dependence (marginal).** At WD70 2–4 Hz, g correlates with
   magnitude (Spearman ρ = 0.46, p = 0.008) and distance (ρ = 0.37,
   p = 0.036), but not with pressure SNR (ρ = 0.03) or coherence (ρ = −0.13).
   With 14 tests this is marginal (Bonferroni 0.0036). A ±5% event-dependent
   effect (source spectrum within the band, Pn versus Pg, slowness) is
   possible and unresolved.
5. **No absolute reference.** LT02 (the only absolute gauge with a result)
   is shallow and tests Eq. 5, not the pre-multiple relation. LA33/LA34 APGs
   are deep but unresolved. Deng et al. (2022, §4, p. 5) also report that
   LA pressures "differ from the formula prediction, so they are also
   calibrated", so their absolute status is itself questionable.
6. **The vertical response is shared and can be the faulty element.** An et
   al. (2017, §4.1) found the LDEO OBS displacements in Cascadia
   inconsistent with collocated APGs and corrected the seismometer, not the
   pressure.
7. **Cross-band context.** XO DPG 0.77–0.78 at 1–4 Hz versus Rayleigh
   0.91–0.93 (LDH channel, 0.03–0.07 Hz) and Deng's WD66 coefficients falling
   0.91 → 0.79 over 0.02–0.07 Hz. The pattern is compatible with a
   frequency-dependent DPG response; Doran et al. (2019, Key Points) report
   that DPGs "exhibit departures from the nominal instrument response with
   regards to both sensitivity and time constant". Channels and bands differ,
   so it is not validation.

### 3.3 Assessment

The relation is physically exact for its stated domain, and the
implementation does not leak. However, there is no field demonstration that it
returns the *true* pressure gain. Event dependence at the few-percent level
is not excluded, and the shared vertical response cannot be separated.
**Keep it diagnostic-only** (agrees with Codex decision 1).

### 3.4 Proposed promotion criteria (for agreement before any test)

1. Synthetic leakage < 1% with each station class's real responses and
   depths (met for the WD66 configuration).
2. An **independent absolute reference** in deep water at the same
   frequencies: a timing-verified APG whose own calibration is not in doubt,
   or an in-situ tested DPG (Doran et al., 2019). It must agree within the
   model interval (ρ, c₀, η: about ±5–7%).
3. No significant dependence of g on magnitude, distance, slowness or
   back-azimuth over ≥ 20 events, with multiple-testing control.
4. Pressure-sensor response verified above ~2 Hz (nominal DPG roll-off) and
   vertical response checked against an independent reference.
5. Temporal stability across the deployment.

Until then, deep-water values remain useful as polarity, timing and
response-shape diagnostics.

## 4. Points of agreement and difference with Codex's review

- **Agree:** findings 1 and 4, decision 1 (pre-multiple diagnostic-only),
  decision 4 (no ad-hoc relaxation for deep APGs), and that reproduction is
  not validation.
- **Adds:** quantitative evidence that the noise-window violation does not
  move LT02 station gains (§1.3); a placement-based resolution that honours
  the documented rule without mass rejection (§1.4); a timing budget and
  held-out drift test (§2); a direct leakage test with real responses that
  answers Codex's request (§3.2); and the LA APG caveat from Deng §4 (§3.2.5).
- **Differs slightly:** Codex lists rejection or a validated measured-onset
  rule as the choices. The reviewer recommends placement as the
  documentation-faithful option, subject to the user's decision.

## 5. Reproduce

```bash
cd ~/Projects/bruit-fm-localp-science
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  scripts/review_local_p_science.py --output-dir results/LOCAL_P/science_review/<new>
```

It reads the saved campaigns in `~/Projects/bruit-fm-localp/results/LOCAL_P/`
and nine LT02 waveform windows, and records and rechecks input hashes.
`run1` was the first pass without the control windows; `run2` is the
reference for every number above.

## Addendum 2026-10-16: outcome of recommendation §1.4

The user chose placement. It was implemented on `claude/local-p-check`
(`dfb1df6`), with the noise-control diagnostic, and checked on LT02 only
(`results/LOCAL_P/campaigns/lt02_noise_rule_v3/`, details in
`docs/local_p_check.md` §2.2 on that branch):

| Band | Before | After |
|---|---|---|
| 0.5–1 Hz | 0.972 (n = 5) | 0.968 (n = 4) |
| 1–2 Hz | 1.021 (n = 9) | 0.973 (n = 7) |

The ≤ 1% prediction of §1.3 held for 0.5–1 Hz but **not** for 1–2 Hz
(−4.7%). Per-event changes were ≤ 2.5%; the shift came from two events
leaving under existing rules: `us1000hhug` (overlap with another catalogued
event once its noise window moved earlier) and `ak0188v8964i` (signal/noise
below 10). §1.3 had kept those events because it tested the noise window
without re-applying the overlap and SNR rules. Lesson: with n ≈ 5–9 the
station median is sensitive to event membership at the ~5% level. That is
larger than the noise-window effect itself and should be part of any
uncertainty statement. The candidate remains conditional; its noise windows
are now no noisier than 30-s-earlier controls (median ratios 0.79 / 1.12).
