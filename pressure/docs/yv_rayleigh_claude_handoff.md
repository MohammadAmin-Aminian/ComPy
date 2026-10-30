# YV Rayleigh unresolved-station diagnosis — Claude handoff

Date: 19 October 2026. Prepared from main `9b43ce8`.
Repository: <https://github.com/MohammadAmin-Aminian/bruit-fm>.

## Task and boundary

The owner wants to understand, calibration method by calibration method, why
not every available station yields an estimate. Start with **Rayleigh calibration
in the YV network only**: 43 pressure epochs, `YV.RRxx.00.BDH`, with their exact
vertical partners. This is a diagnostic/scientific-review task, not another
calibration implementation or an instruction to make every station pass.

Explain the actual bottleneck at each station and recommend the smallest
defensible next experiment. Distinguish a numerical estimate from an accepted
absolute calibration. Keep microseism, local-P, full-data coverage and ratings
outside this task. The owner is waiting for a new computer; no long run now.

Work on a **new branch from current origin/main**, suggested name
`claude/yv-rayleigh-diagnosis`, preferably in a separate worktree. Preserve every
existing contributor branch and worktree. Commit/push the report and any small
read-only audit helper/tests to that branch; do not merge it into main yourself.
Inspect Git and live processes before changing anything. Keep reads/output
targeted; no subagents.

## Essential distinction: 16 baseline estimates, 31 experimental estimates

The owner's previous table deliberately used the frozen **no-tilt nominal
baseline**: 16 YV station estimates, 25 unresolved, two reference-band exclusions.
It did not include experimental response alternatives in the baseline count.

The completed Socrates trial already produced **31 YV numerical estimates**
with `pole_absolute`: an in-memory alternative for positive-real pressure poles,
with the original fitting objective. This adds 15 stations; ten remain unresolved
and two retain their band exclusions. These 31 are **not physically certified**.
Do not repeat that trial or describe its 15 additions as newly discovered fixes.

The weighted-log objective is a separate experimental factor, not an approved
upgrade. Never pool it with pole reflection or choose an objective because it
approaches a paper's factor. Production defaults remain unchanged.

Older tilted YV reports use different gains and event totals: 48 candidate events
versus 49 in the no-tilt baseline. They are historical evidence, not replacements
for this baseline. Always state which workflow and estimator a number belongs to.

### Frozen station snapshot

Here `g` multiplies nominal-response pressure; `1/g` is the reciprocal sensitivity
multiplier. Positive magnitudes do not verify physical polarity. A dash is no
station-level estimate, not zero gain. Trial candidates are numerical QC passes,
not isolated or independently calibrated events. Values are rounded to six places.

| Station | Nominal baseline g | Pole-only trial g | Trial candidate events |
|---|---:|---:|---:|
| RR01 | — | 1.041380 | 4 |
| RR03 | 1.124056 | 1.142024 | 3 |
| RR05 | — | 1.120144 | 4 |
| RR06 | — | 1.139679 | 6 |
| RR07 | 1.079187 | 1.067288 | 3 |
| RR08 | 1.159651 | 1.160453 | 3 |
| RR09 | 1.193725 | 1.175480 | 3 |
| RR10 | — | 1.083792 | 6 |
| RR11 | 1.082658 | 1.053545 | 3 |
| RR12 | — | — | 0 |
| RR13 | 1.223614 | 1.201068 | 4 |
| RR14 | 1.095657 | 1.064746 | 3 |
| RR16 | — | 1.078309 | 8 |
| RR17 | — | — | 0 |
| RR18 | — | — | 0 |
| RR19 | — | 1.129611 | 4 |
| RR20 | — | 1.341531 | 7 |
| RR21 | — | 1.114118 | 6 |
| RR22 | — | — | 0 |
| RR25 | 1.162581 | 1.157092 | 4 |
| RR26 | — | 1.088703 | 4 |
| RR28 | 0.819513 | 0.819513 | 3 |
| RR29 | 0.729913 | 0.729913 | 3 |
| RR30 | — | 1.185682 | 8 |
| RR31 | — | — | 0 |
| RR34 | 0.585256 | 0.585256 | 3 |
| RR36 | 0.653573 | 0.653573 | 3 |
| RR38 | 0.636484 | 0.636484 | 3 |
| RR40 | 0.742586 | 0.742586 | 3 |
| RR41 | — | — | 0 |
| RR42 | — | — | 0 |
| RR43 | — | 1.140108 | 4 |
| RR44 | — | 1.307477 | 4 |
| RR45 | — | — | 2 |
| RR46 | — | 1.015093 | 7 |
| RR47 | — | — | 0 |
| RR48 | — | 1.386710 | 5 |
| RR50 | 0.716797 | 0.716797 | 3 |
| RR52 | 0.678728 | 0.678728 | 3 |
| RR53 | — | — | 1 |
| RR55 | — | 1.592099 | 3 |
| RR56 | — | — | 0 |
| RR57 | — | — | 0 |

### Remaining trial blockers: observations, not proven root causes

| Station(s) | Existing evidence to investigate |
|---|---|
| RR12 / RR42 | Nominal attempts were refused for missing/gapped/nonfinite/flat windows. RR12 also lacks synchronized nonflat paired segments; RR42 has limited shared epoch coverage. Separate pressure problems from vertical problems and selection outside coverage. |
| RR17 / RR22 | Nominal coherence and phase/arrival-review rejections; each also has a recorded 0.0001 s sample-grid mismatch. RR22 does not improve under the pole hypothesis. Do not infer reversed wiring or a hardware defect from the label. |
| RR18 / RR31 / RR47 | Predominantly coherence and/or phase/arrival-review screening; inspect saved frequency-bin and packet evidence before attributing the cause. |
| RR45 / RR53 | Pole-only trial has two / one candidates, below the existing three-event minimum. Older tilted runs also had rotation/training issues, but those are not the binding explanation for a no-tilt result. |
| RR56 | Coherence/phase review; an earlier audit records frequency-dependent phase. Test hypotheses against the actual complex response and event spectra. |
| RR41 / RR57 | Requested 0.03–0.07 Hz band fails the conservative water-column resonance preflight. RR57 is documented as borderline under a frequency-dependent model. Preserve the exclusions; a lower band would be a separate, justified protocol. |

RR55 is not among these ten unresolved trials, but its pressure hardware/response
identity remains uncertain despite three trial candidates. A passing mask is not
evidence that the assumed sensor response is physically correct.

## Read first; do not repeat completed work

1. `CODEX_NEXT_CHAT.md`, then this handoff. Its older task briefs are historical;
   this handoff is the owner's current Rayleigh-only delegation.
2. `docs/calibration_network_summary.md/.json` and
   `docs/socrates_all_station_comparison.md/.json/.csv`: current no-tilt baseline
   versus the completed 89-station experiment.
3. `docs/reviews/codex_integration_review.md` and relevant Sentinel/Socrates
   sections linked there. Initial audit claims have explicit Codex qualifications;
   do not repeat the unqualified original claims as established facts.
4. `docs/yv_heldout_results.md/.json` and
   `docs/yv_heldout_arrival_review.md/.json`: both studies and RR05 replacement
   are **already complete**. Read these before proposing more validation.
5. As needed: `docs/calibration_assumptions.md`,
   `docs/yv_rayleigh_saved_results_review.md`,
   `docs/yv_rayleigh_independent_validation.md`,
   `docs/yv_response_window_validation.md`,
   `docs/socrates_selected_station_trial.md`, `docs/rayleigh_window_review.md`.
   Their older counts/default descriptions do not supersede exact saved settings.

Completed checks include saved-spectrum arithmetic, selected independent raw CSD
checks, arrival/quiet reviews and response/window sensitivities. Do not spend this
task proving all of those again. Relevant implementation entry points:
`src/bruitfm/rayleigh_calibration.py`, `rayleigh_network.py`,
`calibration_experiment.py`, `rayleigh_validation.py`; scripts
`run_calibration_experiment.py`, `run_yv_heldout_validation.py`,
`review_yv_heldout_arrivals.py`, `review_saved_calibrations.py`.
Some existing review scripts also read raw windows: inspect before running.

### Already-completed held-out results

- RR52: 14 isolated usable events, median frozen-gain prediction error +0.664%.
  This is conditional consistency under shared response/model assumptions.
- RR06: one isolated usable event; RR01 and RR48: zero. RR48 had 38/50 events
  outside its selected epoch. Do not interpret absent validation as a wrong gain.
- Longer-epoch RR05 replacement: all 50 events attempted, zero QC-usable;
  40 coherence failures and ten input/window refusals. One example has median
  squared coherence above 0.95 but mean below 0.95; both were required.
- Three RR01/RR06 waveform reviews found stronger later packets. Those cases
  are now inspected development data, not untouched validation. Do not move
  their windows and call the resulting agreement independent validation.

## Exact saved inputs and provenance

Local project root: `/Users/mohammadamin/Projects/bruit-fm`.
Python: `/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python`,
with `PYTHONPATH=src:.`.

- `results/REVIEWS/SIX_CHECKS/final_20261019/review.json`:
  SHA-256 `241e9ab125a718fe92958b2129a267ebeff43076287142d85c0beb6d8cea1fd2`.
  Select YV exact rows and `checks.rayleigh.metrics`; follow explicit result paths.
- `docs/calibration_network_summary.json`: tracked baseline and saved result
  paths/hashes. No newest-directory selection.
- `results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json`: frozen no-tilt
  workflow. Do not substitute `results/NETWORKS/TEST/_workflow/rayleigh.json`,
  which represents the historical tilted workflow.
- `results/EXPERIMENTS/SOCRATES/all_20261013T195356Z_00cf3343/experiment.json`:
  SHA-256 `446b8e36b94c9952736b9a2c9b2353f0dfb0c95e05e8cc3719b9cea72fcc20cf`.
  Read `YV_RRxx/trial.json`, specifically `summaries.pole_absolute`, separately
  from the other variants. All 43 trial files exist locally, including exclusions.
- Held-out root: `results/VALIDATION/YV_HELDOUT/run_20261013T213126Z_87e24fe6/`.
- Completed RR05 root:
  `results/VALIDATION/YV_HELDOUT/RR05_replacement_20261013T215958Z_fee6a2d2/`.

Full products/raw data are intentionally not in Git. If Claude is remote, use
tracked curated evidence and state which local checks cannot be performed.
Do not fabricate fresh verification, download replacement data or rerun an
experiment just because a local product is absent. Request a bounded missing
artifact list. Never guess a channel/epoch/inventory join from station name alone.

## Requested workflow and deliverables

**First, read-only diagnosis.** Reconcile all 43 exact pressure epochs against
the baseline and trial. Produce a station table containing both gains/counts,
selection policy, response hypothesis, epoch/coverage, event refusal reasons,
frequency-bin coherence/phase evidence, physical/reference flags, existing
held-out evidence and one recommended next action. Keep exclusions and
24-hour selection skips separate from failed signal QC. Aggregate event counts
must not be presented as counts of failed stations. Flag missing evidence.

Separate: (a) observed screening/coverage limitation, (b) demonstrated software
defect, (c) supported response-convention concern, (d) unresolved physical
hypothesis. Low coherence is a measurement, not a complete causal diagnosis.
Pole reflection cannot improve coherence on unchanged spectra; identify changes
in phase-valid bins/events separately from changes in fit amplitude.

**Then, propose the smallest next protocol.** Prioritize unresolved mechanisms,
not the largest possible passing-station count. Use saved spectra first. Specify
controls (RR52 remains frozen at `g=0.6787278140712512`), one changed factor,
selection fixed before amplitude/error inspection, input hashes and disjoint
future validation data. A changed band/window is development, not retrospective
validation. Propose at most a small representative pilot before any network run.
Obtain owner approval before new raw-data experiments, production fixes or
scientific-default changes. A clear bug can be documented with a synthetic
regression; fixing production is not automatically authorized by this diagnosis.

Deliver on the new branch:

- `docs/yv_rayleigh_diagnosis.md`: evidence-backed station-by-station explanation,
  existing-work reconciliation, limitations and prioritized next actions.
- Compact machine-readable table/provenance with exact source identities/hashes;
  no copied waveform archive, NPZ campaign or full ignored results tree.
- If useful, a small saved-input-only audit helper, with synthetic tests and a
  safe command that cannot start calibration jobs or update masters. Test only
  new/changed code; no need to rerun every completed science check.
- A short return handoff with branch, commit, files, tests actually run,
  preservation scope, unresolved questions and any approval needed.

Stop at that review/protocol boundary. It is a successful outcome to establish
that a station cannot currently be calibrated with the available information.

## Primary-source leads for any new scientific claim

Read relevant original sections/equations and report access limitations; existing
audit summaries are leads, not substitutes for a newly claimed correction.

- Stähler et al. (2016), *Performance report of the RHUM-RUM ocean bottom
  seismometer network around La Réunion, western Indian Ocean*.
  DOI [10.5194/adgeo-41-43-2016](https://doi.org/10.5194/adgeo-41-43-2016):
  hardware pools, Appendix A response conventions and INSU periodic glitches.
- Zha & Webb (2016), *Crustal shear velocity structure in the Southern Lau Basin
  constrained by seafloor compliance*.
  DOI [10.1002/2015JB012688](https://doi.org/10.1002/2015JB012688):
  Appendix A reference and coefficient direction. Its YL coefficients are not
  additional YV calibration factors.
- Aminian et al. (2025), *Shallow crustal structures of the Indian ocean derived
  from compliance function analysis*.
  DOI [10.1093/gji/ggaf253](https://doi.org/10.1093/gji/ggaf253):
  eight matched YV pressure multipliers in Table 1, with differing preprocessing.
  Published agreement is not independent validation of every shared assumption.

## Non-negotiable preservation rules

Do not modify raw samples/headers, StationXMLs, response snapshots, production
masters (CSV/XLSX), historical spectra/results/workflows, manual scores or any
saved gain. Do not apply gains to original data, responses or PPSDs. Do not
reflect XML poles, flip signs per event, relax timing/phase/coherence/resonance
guards, change model/depth/density, or tune toward a published multiplier to
manufacture agreement. Do not assign ratings or adopt acceptance thresholds.

No completed calibration/held-out batch, `run_three_methods.sh`,
`run_all_calibrations.sh` or full-coverage execution is needed. Use fresh isolated
outputs for any new read-only report and verify the actual protected input scope
before/after; do not claim an archive-wide hash audit. Keep nominal and
experimental conclusions separate throughout.
