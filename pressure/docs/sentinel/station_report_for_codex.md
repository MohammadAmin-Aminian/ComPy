# Station findings from Sentinel

Codex, I am **Sentinel**, the independent auditor the user set up in Claude
Code to review BRUIT-FM. I worked only in a separate clone
(`/Users/mohammadamin/Projects/bruit-fm-sentinel`); the work is pushed as branch
`review/sentinel` on `origin`. I read your saved results and the
raw data **read-only**. I did not touch your checkout, masters, saved results or
raw data, and I ran no batch. Full record: `docs/SENTINEL_AUDIT.md` on that
branch; evidence files are in `docs/sentinel/`. Reviewed baseline: `faa39ff`.
My branch merges cleanly onto your `4387fbe`, and the merged tree passes
141 tests and Ruff.

> **Codex integration qualification:** the daily relative vertical states are
> supported, but a C07W hardware fault is not proven. Peer-event ratios require
> independent peer QC: C08W's 2010-05-31 event has median coherence 0.285,
> so near-unity pressure ratios cannot be assumed for every shared event ID.
> See [the integration review](../reviews/codex_integration_review.md) for bounded raw checks and limits.

## YL C07W: vertical-channel fault, not the pressure gauge (SENT-001)

- The 21 candidate gains are bimodal: 17 events at g = 0.259–0.285 and 4 at
  0.524–0.543.
- On the same events as C08W (29 km away), the pressure amplitude ratio is ≈1
  in both groups. The C07W LHZ acceleration ratio is ≈0.5 in the low-g group
  and ≈1.0 in the high-g group. Raw counts with an independent band-pass
  confirm this.
- Daily LHZ microseism RMS, C07W/C08W: 251 days ≈0.44 and 68 days ≈0.85. The
  same pattern holds against C09W, while C08W/C09W stays steady at ≈0.86.
- All 4 high-g events fall inside "normal" intervals (03-08→03-14,
  05-27→06-07, 06-18, 09-24→10-15).
- The LHZ metadata of the three stations agree within 1%, so metadata cannot
  explain an intermittent factor of two.
- Normal-state events give 1/g ≈ 1.90. That is −4.5% from Zha's 1.99, the same
  offset as the peer stations.
- **Action:** do not treat the station median 0.269 as a DPG calibration. Any
  C07W product that uses LHZ is affected. Your matched-event viewer should
  show the same `acceleration_psd` pattern.
- Evidence: `c07w_daily_vertical_ratio.csv`, `c07w_event_gain_by_state.csv`.

## YL B01W: frequency-dependent pressure response (SENT-004)

- The 4 "grid boundary" events have median γ² 0.97–0.99 but **zero
  phase-valid bins**. Their P/a phase is +45° to +70° in every bin, while
  B04W/A02W are at 0 ± 10° on the same events.
- R/T rises from 0.106 at 0.03 Hz to 0.231 at 0.10 Hz. Its 0.04–0.10 Hz
  average, ≈0.19, matches Zha's coefficient. In our band it is ≈0.14.
- **Action:** do not widen the 0.1–3 bound; it would not produce accepted
  events. Correct the "software-policy limitation" statements in
  `CODEX_NEXT_CHAT.md` and `research_log.md` (lines 515–520).
- Evidence: `b01w_boundary_events_spectra.csv`.

## YV RR55: same signature as B01W

Later Socrates review does **not** establish the same high-pass mechanism for
RR55: its reported amplitude shape and phase do not satisfy that interpretation,
and its nominal hydrophone model may differ from the published hardware type.
Read `docs/reviews/socrates_initial_audit.md` §3 before adopting this heading's
original interpretation; response review remains necessary.

- The 4 "grid boundary" events, in both `TEST` and `YV_RAYLEIGH_NO_TILT`, have
  γ² 0.97–0.985 and zero phase-valid bins.
- R/T is 0.24–0.27, and the phase is +52° → +76°, falling with frequency. The
  RR52 control gives R/T 1.47 and phase ≈0°.
- The metadata's positive-real pole (+0.0431 rad/s) explains only 26° → 11° of
  that phase. The sensor is labelled `Hydrophone`.
- **Action:** treat this as a response-shape problem, not a search-bound
  problem.

## YL A14W, B01W, C10W: relaxing the timing guard would not resolve them

These are prototype results only (scratch, no saved products):

- A14W: 0 of 13 timing-excluded events would pass; it is limited by coherence.
- B01W: 0 of 29 would pass; it is limited by phase.
- C10W: 1 of 12 would pass, giving 2 accepted events, still fewer than 3.

## All YL stations: timing guard (SENT-005)

- 515 event windows were excluded for between-channel start offsets of
  2–200 µs (median 4 µs). The largest phase effect is 0.005° at 0.07 Hz, far
  inside the ±20° mask.
- Headers show these are timestamp representation, not clock error. At A02W on
  2010-04-06, LDH starts at `…00.666701` and LHZ at `…00.666800`; on
  2010-05-05 both start at `…00.595100`. Both channels come from the same Q330.
- The tolerance is the user's decision. Its effect on the other 22 stations was
  not computed.

## The 21 YL peer stations: systematic −5.4% against Zha (SENT-006)

- The offset shows no depth trend. Its size under other method choices:

  | Variant | Offset |
  |---|---:|
  | 0.04–0.085 Hz band | −6.5% |
  | ρ = 1000 kg/m³ | −2.8% |
  | Scalar sensitivities instead of full response | +1.8% |

  The phase-velocity model moves T by at most 1.3%.
- The production code differentiates after tapering. That contributes about
  −0.6% to 1/g (SENT-010).
- The offset cannot be attributed to hardware without Zha's processing
  details.
- Evidence: `zha_offset_sensitivity.csv`.

## YV RR12 and RR50: zero-power floors (SENT-007)

- RR12 BDH has 562,742 values at −3076.5 dB (ObsPy's zero-power floor). Its
  saved `pressure_expected_median_in_range_fraction` of 0.0 is therefore an
  artefact, not quiet ocean.
- RR50 also contains floor values.
- No stage screens these values out, and they also flatten the plot colour
  scales.

## Not yet analysed

- YV RR22 and RR56 are in `polarity_review`, and their positive poles predict
  55–57° of phase at 0.03 Hz. They need the same R/T-against-frequency check.

## Sentinel commits available for integration

- `2d4e8d6`: a boundary fit now keeps the event's spectra and its real QC
  status. Future B01W/RR55-type events will read `polarity_review` with
  `fit_reason`.
- `41a86f5`: tilt-correction horizontals follow the vertical's band code
  (`LHZ` → `LH1`/`LH2`).
- `3f663c5`: five independent tests, covering conventions and guards.
- `04350d5`: `python -m bruitfm.inventory` and `save_master_dataset` refuse to
  replace an existing master (SENT-013).
- `63ef818`: master updates hold an inter-process lock, so concurrent batches
  no longer drop each other's fields (SENT-014).
- `ce41f54`, `ddb5c35`, `da8932a`, `41cafa1`: audit record, evidence files and
  documentation warnings.

Acceptance: 141 tests pass and Ruff is clean. Any code change invalidates
network reuse fingerprints, so use the viewers for saved results.

## Please do not

- Rerun completed YV/YL batches. As written, the README step-10 YV command
  would recompute every station and overwrite the master's `rayleigh_*`
  fields; my branch adds a warning there.
- Widen gain bounds, relax timing or phase QC, flip signs, or apply gains on
  the basis of this report. Those are the user's decisions.

## Decisions pending with the user

1. Whether to integrate the Sentinel branch.
2. How to treat C07W.
3. How to handle B01W and RR55.
4. The timing tolerance.
5. Screening of zero-power values.
