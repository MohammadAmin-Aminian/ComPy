# Next task: RR52 microseism pressure–vertical pilot

**Current engineering checkpoint:** the six-check integration and conditional
accepted-station/noise-model workflow are implemented; read
[the current review workflow](six_check_review.md). The scientific acceptance
policy remains draft, all production scores are unchanged, and no conditional
microseism estimate was promoted. Review saved evidence; do not restart this
historical pilot or its completed calibration campaigns.

**Latest continuation, 16 October 2026:** the bounded microseism follow-ups and
reviewed local-P engineering phase are complete and integrated into main. Read
[the integration closeout](local_p_integration_closeout.md) and
`CODEX_NEXT_CHAT.md` for current status. The owner cancelled the separate long
all-methods run; its partial weeks are preserved, not adopted as a new frozen
review. Next is agreement on six-check acceptance/scoring, not a rerun of this
historical RR52 pilot. No absolute gain or automatic correction was accepted.

**Main integration, 15 October 2026:** Claude's pressure–vertical method is now
combined with the model-free diagnostic and a frozen frequency/week consistency
review. The original 35 best-band gain candidates are superseded. Read
[microseism_integration.md](reviews/microseism_integration.md) and
[microseism_frequency_table.md](microseism_frequency_table.md) for the current
decision and remaining validation. The brief below is historical context.

**Gain-band follow-up:** [the saved-spectrum audit](microseism_gain_band_audit.md)
checks 58 bands at all seven repeatable stations. Three RR40 weekly bands have
approximate phase and daily consistency is weaker. The frozen rule is unchanged.
An 89-row epoch-bound staging export leaves every score blank. The next bounded
task is a P/H1/H2/Z tilt investigation, not a repeat of calibration batches.

**Completed orientation/separate-training follow-up:** read
[the quiet-projection review](rr52_microseism_fixed_projection.md). A published
RR52 orientation reference is now recorded, without modifying XML. Saved prior
quiet weights were tested on separate saved June spectra; they do not consistently
improve target-band coherence. Do not promote an automatic correction. The
bounded task is complete with conditional/unresolved results. Local-P and final
status/scoring integration follow; do not restart this historical task brief.

**Completed follow-up, 14 October 2026:** the bounded saved-spectrum diagnostic
is implemented and documented in
[rr52_microseism_diagnostic.md](rr52_microseism_diagnostic.md). It establishes
that P/Z alone does not identify an independent pressure gain, preserves the
fixed Rayleigh multiplier as a display reference, and leaves all masters/scores
unchanged. The task brief below is retained as the pre-implementation scope.

**Campaign follow-up, 14 October 2026:** the user requested commands for all
stations having first-method estimates. The resulting 46-station manifest and
one-station runner are documented in
[microseism_baseline_run.md](microseism_baseline_run.md). RR52 keeps its saved
paired spectra; other commands create isolated products. No campaign station
job was run during preparation, and no gain is estimated or applied.

**Comparison follow-up, 14 October 2026:** the user subsequently completed all
46 stations with no failures. The exact-log, artifact-verified comparison is in
[microseism_baseline_comparison.md](microseism_baseline_comparison.md). It fits
no gain and finds an unresolved phase-family split: all YV stations are near 0°,
while all YL/7D/XO stations are near 180°.

Status: an independent implementation and RR52 result are in
[microseism_pressure_vertical_check.md](microseism_pressure_vertical_check.md)
(branch `claude/microseism-rr52`). The original brief below is unchanged.

User decision, 14 October 2026: finish the current calibration discussion and
continue the second pressure–vertical check in a fresh chat. Further Rayleigh
failure diagnosis remains open but is deferred, not the next execution task.
No microseism calibration implementation or new waveform run was performed
while preparing this handoff.

## Concrete first task

Establish from primary sources what theoretical pressure–vertical relationship
can be tested in the microseism band. Start with the existing RR52 pilot,
not another Rayleigh calibration or a new all-station batch. Explain the
proposed measurement, assumptions and bounded implementation before proceeding.

Determine whether the method identifies an independent pressure multiplier or
only a diagnostic/joint constraint on pressure gain, vertical response and
seafloor structure. Do not force a multiplier by matching an unconstrained
prediction. State what additional reference/model is needed if it is not
identifiable. Coherence alone is not a calibrated pressure gain.

Primary-source leads already recorded in the project:

- Ruan, Forsyth & Bell (2014), DOI `10.1002/2014JB011162`.
  Read the original method and relevant follow-up before treating it as a
  standalone gain recipe. A 7D2001/J44A reproduction is a previous lead;
  exact data availability has not been verified.
- OBSIC quality-assessment material referenced in the original WP2 brief.
  Identify the exact source; do not invent a ready numerical reference.

Confirm the usable sensor band, transfer direction, physical units, gravity/
water-column/structure assumptions, noise and tilt treatment. The user's earlier
request to make tilt optional concerned teleseismic pressure calibration;
microseism requirements must be assessed separately and stated explicitly.

## Reuse the saved RR52 foundation

Read `docs/paired_spectra.md`, `docs/rr52_paired_validation_summary.json` and
`src/bruitfm/paired_spectra.py` selectively. Saved production pilot:

```text
results/NETWORKS/TEST/YV/RR52/YV.RR52.00.BDH__20121019T090000__58afc2d6f005/paired_spectra/run_bfbef17f2c13/
```

It contains `comparison.json`, `spectra.npz`, `spectra.csv`, `windows.csv` and
`spectra_coherence_transfer.png`. Current master fields point to this run.
Confirm metadata/hashes before reuse; never assume a different/newest run is
scientifically equivalent merely from its timestamp.

- P/Z pair: `YV.RR52.00.BDH` / `YV.RR52.00.BHZ`.
- UTC interval: 2013-06-10 through 2013-06-17; 167/168 hourly pairs retained.
- Full complex nominal responses, upward displacement; no tilt correction.
- `S_dp = mean(D * conj(P))`, `H_dp = S_dp/S_pp`, units m/Pa.
- Individual auto/cross spectra and linear means are saved; coherence is computed
  after averaging complex cross-power, not averaged from single-periodogram
  coherences.
- Display range 0.01–0.5 Hz is not a verified sensor band. Shaded 0.1–0.2 Hz
  is a diagnostic band, not an adopted calibration interval.
- Median squared coherence 0.802142 and transfer amplitude
  4.985790411e-8 m/Pa in that diagnostic band are descriptive, not pass scores.

First reuse saved spectra. If new samples are scientifically necessary, use one
bounded pilot in a fresh output directory, preserve grid/gap/epoch/response
guards and hash/stat-check protected inputs. Keep the nominal and comparison
variants traceable. Do not update the production master or any scientific score.

## Fixed Rayleigh reference, not a new fit target

RR52's saved no-tilt gain is `g = 0.6787278140712512`, multiplying nominal pressure.
It is a fixed comparison reference; it has not been applied to original data.
Fourteen isolated held-out events give median prediction error +0.66414% under
the shared model. This is conditional consistency, not absolute accuracy.
Do not refit RR52, rerun its completed held-out study, or tune the microseism
method to reproduce this number. Any transformed diagnostic spectra must be
new, explicitly labelled products, never replacements of saved inputs.

RR52 is not usable for every event: for 15 October 2013 `usb000kdb4`, saved
mean coherence is only 0.098078; RR05 is 0.941675 and RR06 0.982187. RR06 still
has body-arrival review flags. Do not assume the RR52 control passed that event
or confuse individual failed events with rejection of its station gain.

## Desired next deliverable

A source-supported method note and a bounded RR52 diagnostic implementation
using the existing spectral foundation, with plots, provenance, relevant tests
and explicit limitations. If the reference does not support a unique pressure
gain, deliver the defensible diagnostic and identify the missing constraint.
Extend only after this pilot is scientifically interpretable. Local-P checks,
acceptance/scoring/classification and combined low/high ocean-bottom noise
models follow; they are not authorized as a new full batch by this handoff.

Keep quota low, use targeted reads, no subagents. Preserve raw data, XMLs,
response snapshots, masters, saved results and scores. Code/docs/curated results
go to GitHub; full data/products remain local under existing exclusions.
