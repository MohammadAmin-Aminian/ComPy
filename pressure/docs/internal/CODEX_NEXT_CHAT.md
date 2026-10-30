# BRUIT-FM: compact continuation handoff

## Start here

Workspace: `/Users/mohammadamin/Projects/bruit-fm`, branch `main`.
Repository: https://github.com/MohammadAmin-Aminian/bruit-fm.
README progress is numbered through **68**.
The prior integration was merged
and pushed at the user’s request; frozen-gain studies and targeted waveform
review are complete. The engineering review
and optional scientific trial tools are integrated; inspect `git log -1` for
the latest checkpoint. **Production scientific defaults remain unchanged**:
pole reflection and weighted log fitting are available in the separate trial
runner, not automatically applied by normal calibration commands.

## Current ownership — transfer to Claude

The owner transferred ongoing project work to Claude on 21 October 2026. Read
[`docs/claude_project_handoff_20261021.md`](../claude_project_handoff_20261021.md)
first for the consolidated software/science state, frozen sources, gain medians,
safety boundaries and next task. Existing code is already pushed; this transfer
is documentation-only. The next bounded task is official response comparison
for RR22/RR56 and RR17/RR18, not another calibration campaign or implementation
of all historical formula-review suggestions.

## Latest task — pending-branch audit and deployment readiness

Owner request (21 October 2026): review/respond to/close all PRs and prepare
main for deployment before continuing YV response confirmation. GitHub CLI and
REST both showed **zero existing PRs**. The three unmerged contributor branches
were reviewed instead; see `docs/reviews/deployment_readiness_20261021.md` and
`docs/deployment.md`. Integration PR #1 is now merged/closed at main `7ac2e63`;
its exact post-merge CI run 39587156168 passed. No open PR remains at transfer.
The `codex/deployment-readiness` and all contributor branches are preserved.

Selective formula-review integration retains the attributed report/worksheet
and fixes wrapped-phase summary reporting (versioned circular L1 median) plus
negative fractional native-grid tears. No gain, per-bin QC mask, scientific
threshold, response correction, frozen rule ID or scoring policy is changed.
Historical tables stay unchanged; no science batch is rerun. Untested fixed-name
median/gain plotting tools remain on `claude/formula-review`, not in main's
deployment workflow. The old local-P branch tips are superseded by reviewed
v2 safeguards; keep their evidence, do not overwrite main with them.

Full local suite: **447 passed**, 117 warnings; Ruff/format/whitespace pass.
Wheel build and installed-package smoke pass for nine module help entry points;
CI now exercises packaging too. All 1,104 protected inputs match the preceding
review, including the owner's untracked `gain_histogram.png`. The software is
tested/installable on Linux/macOS Python 3.13, not scientifically accepted or
fully calibrated. Local archives/results are not shipped in Git/wheels and
new-machine paths/provenance need explicit reconciliation. No real deployment,
release publication or long-data campaign was launched.

After PR closeout, return to official response/calibration-record confirmation
for RR22/RR56 and RR17/RR18, then propose the small saved-spectrum pilot. No
automatic metadata repair, gain fitting or new threshold adoption follows.

## Previous task — reviewed Claude YV Rayleigh diagnosis

Owner handed control back on 21 October 2026. Claude's
`claude/yv-rayleigh-diagnosis` at `dbd4856` is reviewed on
`codex/yv-rayleigh-review`; preserve both branches/worktrees. Read
`docs/reviews/yv_rayleigh_diagnosis_integration.md`, then
`docs/yv_rayleigh_diagnosis.md`. Response/polarity hypotheses are not accepted
instrument corrections. Counts remain 16 nominal / 31 experimental YV estimates;
no new gain or production-default change. Ten trial unresolved plus two band
exclusions remain. Existing gains, spectra, masters, XMLs and scores are unchanged.

Four new regressions exposed audit-helper defects; the helper now checks the
completed trial's source hashes, exact epoch/P/Z binding, individual spectra,
optional microseism provenance and before/after input integrity. Scientific
wording now separates screening failures from physical causes and post hoc
H4 from calibration. Claude's original JSON/CSV evidence remains unchanged.
Targeted suite: **59 passed**, 24 warnings; Ruff passes. Real saved-input replay:
43 stations, 583 NPZ spectra, 875 hashed inputs, every prior station metric
identical. Output: `results/REVIEWS/YV_RAYLEIGH_DIAGNOSIS/codex_review_20261021/`.
All 1,104 protected inputs match before/after, including the owner's untracked
`gain_histogram.png`; do not delete, add or overwrite that user file.

Next: official response/calibration-record confirmation for RR22/RR56 and
RR17/RR18, then a proposed bounded saved-spectrum pilot for owner approval.
Do not silently repair metadata or fit new gains. The separate
`claude/formula-review` at `c767579` was inspected but is **not merged**; its
wrapped-phase, phase-family and timing proposals need their own review and do
not authorize adopting thresholds or changing production defaults. No completed
calibration/held-out batch or long full-data run was restarted.

## Previous task — YV Rayleigh diagnostic handoff to Claude

Owner decision (19 October 2026): investigate why each calibration method does
not yield estimates at all stations, starting with **YV Rayleigh**, delegated to
Claude. Read `docs/yv_rayleigh_claude_handoff.md` for the full scope, 43-station
baseline/trial snapshot, exact saved sources, already-completed validation and
requested deliverables. This is diagnosis/protocol work, not another batch.

The prior user-facing gain table used 16 frozen nominal no-tilt YV estimates;
the completed pole-only experiment has 31 provisional numerical estimates.
Do not confuse either with accepted absolute calibration or forget the finished
RR01/RR06/RR48/RR52 held-out study and RR05 replacement. Ten ordinary trial
unresolved stations plus RR41/RR57 band exclusions remain. Keep all old branches,
gains, production defaults, raw data, responses, masters, results and scores.

Claude should branch from current origin/main (suggested
`claude/yv-rayleigh-diagnosis`), inspect saved evidence first, push the report and
any read-only helper/tests to that branch, and propose a bounded protocol for
owner approval. No raw-data experiment, production fix, automatic merge,
threshold change, rating or completed campaign rerun is authorized by this
handoff. Main changes in this handoff are documentation only.

## Previous task — pressure-reference correction and RR52 display

Read `docs/pressure_reference_review.md`. Eight core references/relevant
sections were reviewed; this is not an exhaustive literature or absolute
calibration validation. Handout values attributed to Webb remain preserved
but are excluded from new numerical comparisons, not promoted to a global
envelope. The original Brown CSV remains unchanged; a new 37-row transcription
uses actual solid-curve support about 0.01139-78.10 Hz (axes are 0.01-100 Hz).
New plots distinguish environmental regions, nominal sensor response and
hypothetical full-band gains. Local-P shifts remain legacy band-centre
illustrations because saved PSDs are already smoothed. Gains are unchanged.
Owner's subsequent display decision: show everything supplied by Wayne on the
plots. All five reference frequencies (ten low/high markers) are restored by
default as "Wayne-supplied (Webb)", including points outside the numerical
comparison band. Numerical exclusion and provenance caveats remain unchanged.
Focused display/pressure/gain regression: 37 passed, 5 ObsPy warnings.

Latest RR52 command:
`PYTHONPATH=src:. /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python scripts/plot_rr40_gain_comparison.py --station RR52 --show`
This opens the new four-panel comparison and gain-only dB plot. It reads only
saved evidence, creates isolated figures/manifest and verifies source hashes.
No raw/calibration batch, master/score/result rewrite or response correction.
Historical reference fractions (including RR52 0.8421) remain historical and
must not be presented as verified Webb acceptance. Do not rerun old workflows
because a new code/reference fingerprint differs.
Verification: **423 tests passed** (110 ObsPy warnings); Ruff check/format pass.
Fresh reviews/plans flag historical reference protocols without replacing the
saved statistics. Of 353 historical protected inputs, only the authorized
quality-review code changed; the other 352 data/evidence/master inputs match.

## Previous task — reviewed full-data coverage workflow

Codex reviewed Claude's `87ff816`; see
`docs/reviews/full_data_coverage_integration.md` for integration corrections,
verification, the bounded real saved-PPSD smoke and preservation checks.
No long data run was launched. The owner is waiting for a new computer for
about two weeks; do not start the 198-task campaign now. Continue descriptive
comparisons of existing evidence, not ratings or calibration acceptance.
This runner fills missing PPSDs and three spectral checks, **not** missing
calibration campaigns. New-machine paths/provenance need reconciliation.

Owner decision (19 Oct 2026): **no ratings, thresholds, scores or approved
noise models yet.** First inventory all local data, reuse valid products and
run only missing PPSDs and applicable checks. This supersedes the
acceptance-policy review below as the immediate next task. Branch
`claude/full-data-coverage` adds `python -m bruitfm.coverage` with
`plan`/`run --execute`/`report` and `scripts/run_full_coverage.sh`. Read
`docs/full_data_coverage.md` and `docs/full_data_coverage_handoff.md`.
Dry run: 89 saved pressure PPSDs reused, 15 new pressure PPSDs plus 61 × 3
spectral checks planned, 19 blocked. Calibration checks are reuse/deferred
only. Claude executed no task; Codex's one-task smoke is recorded above.
The owner runs the later full `execute`; do not rerun
completed campaigns.

## Previous task — six-check integration and conditional noise-model workflow

The owner asked to implement the remaining integration after the pressure PPSD
plots. Read `docs/six_check_review.md` and `docs/six_check_phase_report.json`.
The completed checkpoint is developed on `codex/six-check-integration` and
integrated/pushed to main after verification. Contributor branches are preserved.

New `python -m bruitfm.quality_review` and
`bash scripts/run_quality_review.sh` read explicit saved sources, never newest
folders. Real review: `results/REVIEWS/SIX_CHECKS/final_20261019/` — 89 stations,
118 exact channel/epoch/source identities, 89 saved PPSD paths, 54 unjoined
local-P rows. Baseline 46 Rayleigh estimates, seven microseism conditional
estimates, one current LT02 conditional and 63 legacy local-P rows retained.
Two stations have numerical values from three methods; 14 more from two.
This is not accepted calibration. Full scientific scores in production remain
unchanged, and old manual values are retained in the derived review.

`docs/six_check_policy_v1.json` is **draft**, not an adopted scientific policy.
The classifier uses explicit human verdicts, never calculation completion or
an average of six scores. Conditional/diagnostic outcomes do not auto-pass.
Noise approval additionally requires a fingerprint-bound reviewed decision,
independent-reference record, explicit nominal-response acceptance and approved
frequency/time support. No real approvals were fabricated: accepted-channel
export is empty, and the real noise run is a selection plan only.

`python -m bruitfm.noise_model` has a tested synthetic approved path: equal
station weighting, per-frequency coverage, complete smoothing-bin/time-window
support, no duplicate overlapping station records, no gain application, protected
sources and optional all/network/channel/network-channel grouping. Actual model
generation refuses an empty approved pool. Its 5th/95th percentile outputs are
descriptive, not certified global noise models. Environmental/seasonal grouping,
uncertainty and field acceptance remain future work.

The plotting code (`plot_network_ppsd.py`, docs and 20 focused tests), previously
uncommitted, is included in this checkpoint. See `docs/network_ppsd_overlay.md`.
User's preferred plot: pressure only, probability heat map, PQLX, no means or
medians. Occurrence cutoff is display-only; counts are preserved. Full-range
plots retain extreme levels. One-octave stored smoothing cannot be undone by
plotting. Never rerun raw PPSDs to improve the image without an explicit task.

Verification: **388 passed, 78 existing warnings, 62.74 s**; Ruff check/format
and shell syntax pass. Real wrapper + selection plan passed; all 353 protected
review inputs and 356 noise-plan inputs retain their hashes, including production
CSV/XLSX. No real waveform/calibration batch, master write or gain application
was performed. Synthetic approvals exercise code only. Full provenance and the
saved review SHA are in the curated phase report.

Next concrete action: scientifically review the draft acceptance policy and
the YV evidence table with the owner, including independent reference criteria.
Do not mark any candidate passed, apply gains, write master scores or call this
the completed WP2/global noise model. No completed campaign needs to be rerun
to review the table. Safe reproduction command: `bash scripts/run_quality_review.sh`.

## Prior completed task — reviewed local-P phase closeout

The owner initially kept main pinned for the long all-methods batch, then
requested cancellation and final testing/documentation/push before leaving.
The exact batch process group was interrupted with SIGINT and its exit verified
before main integration. Completed and partial products are preserved; nothing
was restarted. The reviewed `codex/local-p-finish` branch is integrated into
main; Claude's branches/worktrees remain available. See
`docs/local_p_integration_closeout.md` for interrupted-run provenance,
Claude-review disposition and next steps.
Read `docs/local_p_phase_report.md` and its curated JSON first. The engineering
phase is complete; absolute calibration and scoring are NOT accepted.

V2 enforces the already documented noise-window guard by moving noise earlier,
without changing the signal window or thresholds. Bounded LT02/LA33/LA34 HDH
validation: all three completed; LT02 retains seven events (was nine), gains
0.967805 at 0.5-1 Hz (n=4) and 0.972951 at 1-2 Hz (n=7). All accepted events
pass the guard. LA33/LA34 remain unresolved. Raw data/responses unchanged.
No all-station rerun was launched by Codex. New `review_local_p_phase.py`
composes the explicit prior 120-pair user campaign with the three v2 rows.
63 other executed pairs retain explicit legacy noise-guard flags; do not
present them as v2-validated. Numerical deep-water results remain diagnostic.

Hardening: unique full-key/epoch exports and explicit hash-bound 7D inventory
mapping; per-pressure-epoch evidence; waveform-root/stat/dependency resume;
parent failure isolation/campaign lock; exclusive protected helper outputs;
tie-safe, qualified polarity; no fallback to legacy microseism gains.
`docs/local_p_table_staging.csv` and `local_p_station_comparison_v2.*` are
staging/reporting artifacts, not production master updates. The restricted
local_p_* updater and read-only-default importer were verified on a private
copy: 16,884 rows and 140 original columns unchanged, 66 exact diagnostics
patched, 54 unresolved/not-run joins skipped. Scores/check statuses preserved.
All 144 phase-review input hashes matched. Combined suite: 332 passed,
78 warnings, Ruff check/format and shell syntax clean (see final report).

Final repeat verification: **332 passed, 78 warnings, 62.50 s**; Ruff check and
format clean; shell syntax and the no-waveform preflight pass. Fresh saved-result
composition `results/LOCAL_P/reviews/main_closeout_20261016/` matches the prior
phase-report SHA exactly; all 144 protected inputs remain unchanged. Main's
catalogues are byte-identical copies of the six saved catalogue/sidecar files,
not a refetch. No production master update or score assignment was performed.

The engineering workflow is now implemented above; scientific policy adoption
is still pending. Do not restart completed
campaigns, treat the interrupted -1 week as complete, or rewrite historical
results. **Late corroboration already arrived:** `438c4ca` on
`claude/local-p-check` and scientific review `7ebd872` reproduce the exact
LT02 0.967805/0.972951 medians and seven events with the placement fix. Those
late code additions are not wholesale merged; this branch keeps conservative
floor rounding and all tested export/resume/safety fixes. See the phase-report
addendum before asking Claude to duplicate that work. Field validation, independent vertical/reference constraints
and six-check scoring/classification policy are distinct future decisions.

## Prior completed task — integrated experimental microseism workflow

**Orientation/separate-training review completed:** read
`docs/rr52_microseism_fixed_projection.md/.json` and
`docs/rr52_orientation_reference.json`. Scholz et al. (2017), Table 1, gives
RR52's standardized BH1 circular mean as 29 degrees (P-pol 6 measurements,
reported 2STD 10 degrees; R-pol 27, 17 degrees). The PDF table was visually
checked. Native INSU vs standardized handedness mapping remains unverified;
do not write 29/119 degrees into XML or treat horizontal yaw as physical tilt.
New read-only runner: `python -m scripts.review_microseism_fixed_projection`.
It reuses the 46 saved hourly P/H1/H2/Z coefficient windows and saved historical
quiet weights, not raw data. Primary rule: latest completed prior two-hour
quiet fit, regardless of event QC or target amplitudes. This selects 18 April
2013; four earlier fits are sensitivities and July is excluded. Training/test
data are disjoint, but the June test intervals were previously inspected, so
this is not preregistered blind validation.
Result: low-band Z power falls, but highest-band power increases ~8-9% and
coherence decreases on both days. **Do not adopt an automatic correction.**
All gains, frozen selections, responses, masters and scores remain unchanged;
2,545 protected files match. Full review:
`results/REVIEWS/MICROSEISM_FIXED_PROJECTION/run_20261016T100814Z_584abbcf/`.
The bounded microseism validation task is complete with conditional/unresolved
outcomes, NOT accepted absolute calibration. Do not extend rotations or rerun
completed batches simply to force agreement. Next: separate local-P review and
six-check status/scoring integration, preserving experimental microseism labels.
A new physical geometry/tilt or independent reference study needs an explicit
protocol; this review does not justify a broader RR40/7D correction campaign.
Verification: 276 tests pass (78 existing warnings); Ruff check/format pass.

**RR52 horizontal-coupling pilot completed:** read
`docs/rr52_microseism_horizontal_pilot.md/.json`. The new bounded runner is
`python -m scripts.pilot_microseism_horizontal`. Two 24-hour target days
(10/14 June 2013) retain 46/48 windows; both 23:00 hours are refused for a
0.0024-second cross-file native-grid change. Every retained P/Z spectrum
reproduces its saved counterpart. All 2,498 protected inputs match.
In the lowest gain band, 67.5%/96.9% of Z power is statistically predictable
from H, while P/Z coherence is weak. These are NOT literal tilt fractions:
physical Rayleigh motion and finite-sample projection bias remain possible.
Both BH1/BH2 azimuths are recorded as 0 degrees, so do not infer geographic
axes or a physical rotation. No cleaned waveform, correction, gain or score
was produced. Do not rerun these days; hourly four-channel coefficients are
saved under `results/REVIEWS/MICROSEISM_HORIZONTAL/run_20261015T155523Z_d3564085/`.
The subsequent review above completes the bounded follow-up; its original next
brief was to verify orientation/deployment information or review a relative-
axis approach, design quiet-trained signal-preserving diagnostics with held-out
distortion checks, then RR40/7D controls as data permit. A blind full batch or
in-sample H subtraction would not solve the identifiability problem. Final
absolute response/model constraints and scoring policy are still open.
Verification: 262 tests pass (78 ObsPy/MiniSEED warnings); Ruff check/format pass.

**Saved gain-band audit completed:** read `docs/microseism_gain_band_audit.md`
and its curated JSON before interpreting the seven conditional estimates.
All seven stations, 22 passing weeks and 58 eligible gain bands were audited:
55 bands meet the legacy 10-degree phase screen; three low-frequency RR40
bands are approximate (10.46–11.89 degrees). Daily cases are weaker: five RR52
day/band cases exceed 20 degrees; sparse 7D days include single-hour records.
These are descriptive flags, not a new rule or retrospective reselection.
Complex responses reproduce exactly from saved XML; absolute calibration,
physical polarity and tilt are not thereby verified. Cached water/dry model
choice changes 7D band medians by <0.6%, not enough by itself to explain the
4–6% Rayleigh offset. No physical cause is established.
New runner: `python -m scripts.audit_microseism_gain_bands --review PATH`.
It checks 2,466 protected inputs and writes a fresh read-only review.
`docs/microseism_table_staging.csv` has 89 full epoch keys and blank scores:
seven conditional/unvalidated, 35 insufficient weeks, 42 no passing weeks,
five pairing-unresolved. It has NOT been merged into production masters.
The first bounded RR52 P/H1/H2/Z pilot is now complete above. Its original brief
was RR52 base week, including 10 and flagged 14 June 2013, and RR40 base control.
Header-only
preflight confirms the four local YV channels; the configured 7D inputs have
only P/Z. Freeze the procedure before inspecting new amplitudes, avoid removing
physical Rayleigh motion, and preserve every v1 selection and source product.
Final scoring and independent vertical/model constraints remain open.
The original local-P development remains preserved in `claude/local-p-check`;
its reviewed implementation is now integrated as described above. This older
microseism handoff is not a request to merge later contributor changes implicitly.

**Independent weeks completed:** read `docs/microseism_independent_weeks.md`.
The requested eight new YV/7D shifts all finished using main `cc323b9` and the
unchanged frozen rule. 232 station-weeks attempted, 188 paired/checked, zero
model-free diagnostic errors; 40 unresolved pairing and four without a base
interval. The current 24-run review gives seven conditional repeatable
estimates: stable RR40, RR52, YL.C05W/C09W, 7D.J45A/J54A; moderate 7D.J67A.
Five are outside the YL design set. 35 stations lack three passing weeks and
47 have none. YV medians are within 1% of Rayleigh; the two matched 7D values
remain 4-6% lower. No accepted absolute gains or scores are claimed.
All 5,274 protected files match. New review:
`results/REVIEWS/MICROSEISM_INTEGRATION/run_20261015T141224Z_2a1d41d2/`.
Do not rerun the completed weeks or Rayleigh batches. The gain-band phase
review is now complete above; next is coherent-tilt/vertical-response work, then
the 7D model/band offsets and C07W time states. The current P/Z-only spectra
cannot exclude coherent tilt. Preserve the frozen rule and original inputs.

**Superseded by the 15 October integration:** read
`docs/reviews/microseism_integration.md` and `docs/microseism_frequency_table.md`
first for the current microseism decision. Claude's branch at `b603948` is
preserved and merged with its history; Codex implemented its frozen
frequency/week rule, retaining the earlier model-free diagnostic and legacy
best-band results. The reviewed 16 saved runs currently yield two conditional
repeatable YL development-set estimates (C05W/C09W); 25 stations have fewer
than three passing weeks and 62 have none. No gain is accepted or applied.
The old 35 best-band candidates are superseded, and C07W's single-week rejection
of its baseline is retracted. Source hashes were checked; no science batch was
rerun during integration.
Verification: 236 tests and Ruff checks pass; 4,121 protected files still match
their recorded hashes. See the review for the exact scope and qualifications.

The initially requested `bash scripts/run_microseism_validation_next.sh` is
now complete; do not start it again. Existing YL/XO weeks and completed shifts
were reused. Recheck live status when continuing. Do not retune the frozen rule
to match the publications. Gain-band phase is audited; coherent tilt, vertical
response and model systematics remain open. `run_three_methods.sh` now defaults
to main's integrated code and saved Rayleigh results, with an explicit
`--rerun-rayleigh` override. Keep long jobs with the owner; no subagents.

The sections below retain historical completion and source context.

## Previous completed task — 46-station microseism comparison

The user completed every job in
`results/DIAGNOSTICS/MICROSEISM_BASELINE_LOGS/run_20261014T170400Z` with no
failures. Read `docs/microseism_baseline_comparison.md` and its 46-row CSV/JSON
plus figure. The new `bruitfm.microseism_compare` aggregator uses only explicit
paths from that run, verifies 462 protected artifacts/paths and reproduces all
saved spectral aggregates. It reads no new waveform and changes no response,
master, prior result or score.

Weekly median coherence² is at least 0.8 for 37/46 stations; 30/46 contain a
Ruan-reference >0.95 bin; 149/275 station-days across 36 stations meet the
Saikia daily >0.8 reference. These are diagnostics, not scores. Coverage ranges
from 24–168 hourly segments and weeks are not synchronized. RR52 stays at
0.802142 median coherence, 29/361 reference bins and 6/7 selected days.

All 16 YV phases are near 0° while all 30 YL/7D/XO phases are near 180°. A
positive scalar gain cannot create this split; audit channel/response and signed
processing conventions before pooling modal admittances. Fixed-gain display
reduces within-network amplitude dispersion for YV/YL/XO but raises it for 7D;
this is not independent validation because structure, depth, vertical response
and wavefield differ. No microseism gain was estimated, so paper differences
remain those of the first method: eight YV +0.35…+6.08% (median +3.20%); 21 YL
peers −7.04…−3.37% (median −5.43%). Ruan needs a modal reproduction; Deng's XO
body-wave bands remain non-equivalent. Complete suite: 186 tests; Ruff passes.

The next defensible task is a read-only phase/sign audit of exact P/Z response
and processing conventions by network. Only after resolving that should one
assemble independently verified water depths, vertical amplitudes and shallow
structure for a bounded modal reproduction. Do not infer gains from the current
cross-station amplitude spread.

## Prior task completed — RR52 microseism diagnostic

Read `docs/rr52_microseism_diagnostic.md` and its curated JSON/figure. The
master-linked saved spectra were reused without a waveform/calibration rerun.
Primary sources establish `H_nominal=(gP/gZ)*eta(f; water, structure)`, so the
P/Z pilot does not uniquely identify pressure gain. RR52's Rayleigh
`g=0.6787278140712512` remains fixed and was displayed only through `H/g`; it
was not fitted or applied. No master or score changed.

New analyzer: `src/bruitfm/microseism_diagnostic.py`. It verifies the five saved
artifacts and recorded inputs, reproduces the weekly spectra, and re-averages
the stored 167 hourly pairs by day. Median coherence² remains 0.802142; 29/361
bins exceed the Ruan 0.95 reference and six/seven daily band means exceed the
later Saikia 0.8 reference. Neither is adopted as a score because the pilot is
one week, uncorrected for tilt, and has no event screen or modal inversion.
Full suite: 176 passed; Ruff check/format passed.

## Completed campaign implementation context

Read `docs/microseism_baseline_run.md`. The scope was the 46 saved
first-method (`previous_g`) estimates only: YV 16, YL 22, 7D 4, XO 4. The
manifest is `src/bruitfm/reference_data/microseism_baseline_stations.csv` and
the runner is `python -m bruitfm.microseism_campaign --station NET.STA`.
`--list` prints exact P/Z IDs and frozen gains. RR52 reuses its designated saved
paired run; other stations write fresh products under
`results/DIAGNOSTICS/MICROSEISM_PAIRED/` and derived diagnostics under
`results/DIAGNOSTICS/MICROSEISM_BASELINE/`. No long campaign station was run
while preparing the code; the user subsequently completed the campaign above.

The user requested a one-command launcher. Use
`bash scripts/run_microseism_baseline_all.sh`; it attempts all 46 sequentially,
continues after station failures, and writes per-station logs plus `summary.tsv`
under a fresh `results/DIAGNOSTICS/MICROSEISM_BASELINE_LOGS/run_<UTC>/` folder.

Most pairs have equal native rates. XO LA33/LA34 explicitly reuse the audited
125/100 Hz to 1 Hz working-copy resampling; WD66/WD70 use its guarded 1 Hz path.
The original data/headers, response snapshots, masters, scores and prior results
are never changed. All 46 exact response pairs reproduce at 0.10/0.15/0.20 Hz.
The implementation checkpoint had 181 tests; the comparison checkpoint has
186. The combined table uses only the printed explicit diagnostic paths and
does not silently select newest runs. YV publication values remain in `g`, YL
in `1/g`, and no numeric matched 7D/XO paper values were invented. Microseism
P/Z still does not identify an independent gain without vertical, water-depth
and coupled-modal constraints.

The next defensible microseism extension needs independent RR52 water depth and
shallow structure, verified vertical amplitude response, and a station-specific
coupled modal curve. A new bounded P/H/Z product is needed to assess tilt. Do
not fit an unconstrained model to force the Rayleigh gain, rerun calibration
batches, or update scores. Rayleigh troubleshooting and 7D/XO benchmarking stay
deferred unless the user reorders them. Use targeted reads, low quota and no
subagents; the user prefers to run long jobs.

Latest agreed counts: 89 distinct stations; 46 baseline estimates (51.7%);
61 optional-trial estimates (68.5%), all provisional. YV 31/43, YL 22/25,
7D 4/15, XO 4/6; 26 unresolved plus two reference exclusions. The 7D/XO gains
exist; the missing part is a verified directly comparable publication benchmark.
Eight matched YV INSU gains differ +0.35–6.08%, median +3.20%; 21 YL peer
reciprocals differ −7.04…−3.37%, median −5.43%. C07W's normal-state diagnostic
is −4.45% versus its mixed-state +86.76%. These are publication discrepancies,
not an absolute accuracy percentage. No newly estimated DEPAS gain is certified.

The saved 15 October 2013 `usb000kdb4` comparison has mean coherence² RR05
0.941675, RR06 0.982187 and RR52 0.098078. RR52 is not a valid control for that
particular event; RR06 passes coherence but retains body-arrival review. This
does not undo RR52's other 14 isolated held-out successes. No new fit/run was
performed for this read-only comparison.

Current checkpoint verification: no calibration or paired-processing job was
live. All five designated RR52 artifacts exist; their hashes are recorded in
the curated diagnostic. Its 16 raw-source size/mtime records and inventory/two
response-snapshot hashes match, and the saved aggregate reproduces from all 167
window spectra. Git whitespace, Ruff check/format and all 176 tests pass. Raw
data and the complete fresh diagnostic remain local/excluded from Git; only its
curated JSON/figure and method/code/tests are tracked.

## Latest completed held-out study

Latest follow-up: read `docs/yv_heldout_arrival_review.md/.json` and
`docs/yv_heldout_rr05_replacement.csv`. Three native RR01/RR06 overviews show
stronger calibration-band packets later than predicted body arrivals; this
supports prospective window research, not retrospective isolation acceptance.
Original QC labels stay unchanged. The display tool is
`scripts/review_yv_heldout_arrivals.py`; Ruff and its three real-window checks pass.
Four candidates were screened using headers before new errors; RR05 had 42
covered epoch-valid windows versus RR48's nine. Completed follow-up:
`results/VALIDATION/YV_HELDOUT/RR05_replacement_20261013T215958Z_fee6a2d2`:
g frozen at 1.12014351665, all 50 events attempted, zero QC/isolated usable,
40 coherence failures, ten input refusals. Do not rerun either completed study.
Review root: `results/VALIDATION/YV_ARRIVAL_REVIEW/review_20261013T215958Z_fee6a2d2`.
Sources/raw stats still match, including the original four-station study.
Deferred Rayleigh work: inspect signal/noise and packet placement; any changed
windows on reviewed cases are development, not fresh validation. Keep gains
provisional, production defaults unchanged and scores blank.

The user completed `results/VALIDATION/YV_HELDOUT/run_20261013T213126Z_87e24fe6`.
Read `docs/yv_heldout_results.md`, `.csv` and `.json`. Frozen gains were never
refitted/applied. RR52 has 14 usable isolated events, median error +0.664%; this
supports event consistency under the shared reference, not absolute acceptance.
RR06 has one isolated event (-1.041%); RR01/RR48 have zero. RR01 has one and
RR06 three QC-usable events in total when body-arrival review cases are included;
these do not establish isolated validation. RR48's epoch ended 2013-06-10 and
38/50 attempts lie outside it; nine other attempts fail coherence. Do not label
that gain wrong from absence of usable tests. RR01/RR06 primarily fail coherence.
All 71 protected metadata/code/product hashes and 2,208 indexed raw-file stats
still match. No scores or source data changed. No rerun is needed. The targeted
arrival review and longer-epoch replacement test above are now complete;
retain the new DEPAS gains as provisional.

## Held-out runner and previous task context

The user authorized validation of RR01/RR06/RR48 with RR52 as control and asked
for code to run. New utility: `scripts/run_yv_heldout_validation.py`.

```bash
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  scripts/run_yv_heldout_validation.py --open-plots
```

It freezes pole-only/original-fit gains from the completed all-station trial;
excludes all 42 previously recorded YV events plus the native development probe;
and evaluates all **50 unused mainshocks** under Mw≥6/depth/24 h rules. No gain
refitting, source gain application, XML/master/score writes or production-default
change. Primary summaries flag/exclude other packets and predicted P/S arrivals;
event scatter/bootstrap remains conditional, not absolute calibration confidence.
Sensor/serial metadata are recorded; physical hardware certification is pending.
Native raw checks passed; full software suite **173 passed**, style checks passed.
The complete four-station study is finished. For future studies, read the
latest completed root `validation.json` under `results/VALIDATION/YV_HELDOUT/`,
then station `validation.json`, CSVs and response audits. Do not substitute the
`YV_HELDOUT_CODE_CHECKS` development folder for the formal study. Instructions:
`docs/yv_heldout_validation.md`. README numbered progress reaches **49**;
step 46 prepares validation, not completed physical acceptance.

## Previous completed integration checkpoint

- **173 tests passed**; Ruff check/format and runner syntax passed. No repeated
  science batch is needed for the prior integration checkpoint. Check GitHub CI for the final
  merge if needed; do not repeat completed tests without a new code change.
- The user completed all **89 station trials**. Root:
  `results/EXPERIMENTS/SOCRATES/all_20261013T195356Z_00cf3343`.
  Read `docs/socrates_all_station_comparison.md` and its 89-row CSV/JSON.
- Preserved baseline: **46 provisional station gains**. Pole/state trial:
  **61 estimates**, including C07W's three-event normal-state subset; 26 ordinary
  unresolved stations and two reference-band exclusions. YV changes 16 → 31.
- All 2,432 protected source hashes stayed unchanged. No gain, XML, master,
  waveform, previous result or score correction was applied. Full products
  remain local under existing exclusions; curated comparisons are tracked.
- Keep the original fitter as default: the weighted objective worsens YL peer
  median published discrepancy from −5.43% to −6.72%; eight matched YV INSU
  comparisons change from +3.20% to +3.50%. Improved numerical QC does not
  certify a physical calibration.
- Newly estimated YV stations: RR01/RR05/RR06/RR10/RR16/RR19/RR20/RR21/RR26/
  RR30/RR43/RR44/RR46/RR48/RR55. RR55's response/hardware identity remains
  unverified. RR22 remains unresolved. Review held-out phase, arrival isolation,
  absolute vertical response and uncertainty before accepting these gains.
- C07W's normal-state subset gives 1/g=1.90135 versus Zha 1.99, instead of its
  mixed-state 3.71649. Its hardware mechanism is unconfirmed. B01W has no
  accepted scalar calibration; high-pass fits remain diagnostic.
- INSU glitch removal, noise subtraction, local model/density changes and timing
  relaxation are not implemented by this merge.
- Roles/history: `docs/contributions.md`. User: scientific direction, downloads,
  real runs and decisions. Sentinel: independent engineering audit/fixes.
  Socrates: independent scientific review/proposals. Codex: development,
  experiment implementation, verification, reporting and integration.

The all-station command below is historical reproduction, not a request to rerun:

```bash
PYTHONPATH=src conda run --no-capture-output -n bruitfm python \
  scripts/run_calibration_experiment.py --network all --all-stations --open-plots
```

The runner reads all stored spectra and recoverable legacy boundary windows;
keeps input/selection/overlap guards; includes XO APGs using their saved rate;
and retains RR41/RR57 reference exclusions. Each invocation uses a fresh output
folder. Instructions: `docs/socrates_selected_station_trial.md`.

Read this note first; consult `CODEX_HANDOFF.md` only as needed. Preserve the
separate review worktree's untracked products. The earlier sections below retain
baseline/history details; do not mistake them for requests to restart batches.

## Latest independent review

- **155 tests pass**; Ruff and runner syntax pass. All 416 saved accepted-event
  fits reproduce exactly across 89 stations; 46 provisional station gains remain.
  Eighteen bounded raw event probes reproduce candidate gains within 8.9e-10
  and retain rejection diagnoses. No full YV/catalogue batch was rerun.
- Sentinel fixes preserve spectra/QC at gain boundaries, use the vertical band
  for tilt horizontals, refuse master overwrite and lock cooperating writers.
  Codex also locks fresh export creation; its new regression failed before fix.
- C07W daily neighbor ratios support two relative vertical-amplitude states;
  hardware cause is unconfirmed. Do not accept its mixed-state median as DPG gain.
- B01W has no phase-valid bins in its four coherent boundary events. Widening
  the scalar range alone will not resolve it. RR55's same phase rejection does
  not establish the same high-pass mechanism; check sensor identity separately.
- Aminian (2025) Table 1 provides eight verified YV pressure multipliers.
  Socrates compares the original tilt baseline: +0.3–6.1%, median +3.0%.
  This is method-lineage reproducibility, not independent absolute validation.
- Real-pole reflection changes phase but not amplitude at fixed bins; recomputed
  QC can nevertheless change station gains. No XML/pole correction was applied.
  IID day-bootstrap intervals omit serial dependence; they are descriptive.
- New review products remain in the separate integration worktree
  `/Users/mohammadamin/Projects/bruit-fm-sentinel-review/results/REVIEWS/`.
  Preserve that worktree's untracked products. Curated reports are tracked.

## Completed work and current results

- New 7D/XO runners and current status: `docs/7d_xo_runs.md`. All writes use
  isolated `_master` and result roots. Original exports remain unchanged.
- 7D completed 30 PPSDs / 15 calibration jobs: J30A/J38A/J45A/J54A provisional
  gains; 11 unresolved. No verified published station-matched gains located.
- **Latest XO all-gauge run is complete**: LA33 APG g=.9342752751 / 7 events;
  LA34 APG g=.9626801225 / 8; WD66 DPG g=.9094978847 / 16; WD70 DPG
  g=.9298885005 / 21. LT02/WD48 have zero accepted events. Root:
  `XO_ALL_GAUGES_NO_TILT`. LT02 has 40 coherence rejections; WD48 has 29 missing
  windows and seven coherence rejections. Timing/coverage review remains.
- Original `XO_ALL_NO_TILT` EDH/DPG run is preserved: two DPG estimates / four
  unresolved. LA33/LA34/LT02 used EDH hydrophones; this was not an APG test.
  WD66/70 gains change +.199244%/−.007922% in the latest run, retaining all earlier
  accepted events and adding six/five. This does not establish published agreement.
- User completed `bash scripts/run_xo_apg.sh`: **XO_APG_NO_TILT** has two APG
  gains (LA33/LA34, seven/eight events), LT02 unresolved. These 2 Hz gains differ
  from the latest 1 Hz APG results by less than .00003%. All runs are preserved.
  Explicit 2 Hz anti-alias working-copy resampling, native UTC gridding, 0.1 ppm
  tolerance for 125.000007629 Hz header encoding; no raw/header writes. Independent
  exact response snapshots remove the need for PPSD computation. Default unequal-
  rate refusal remains. Logger clock offsets (especially LT02) are not corrected.
  Read `docs/7d_xo_runs.md` for command, verified method and limitations.
- Deng et al. 2022 DOI 10.1029/2021GL096952 is a matching XO paper: body-wave
  frequency-dependent calibration; WD48 about .95, WD70 .91–1.71. Direction/exact
  band values unverified, supplement blocked. Do not claim scalar gain agreement.
- All downloaded XO command (already completed by the user):
  `bash scripts/run_xo_all_gauges.sh`, separate `XO_ALL_GAUGES_NO_TILT` root,
  no PPSD, 1 Hz working copies, all eligible events, no gain/master update.
  Hardware-based APG-first / DPG fallback selects LA33/LA34/LT02 HDH/HHZ and
  WD48/WD66/WD70 LDH/LHZ. WD66/70 HDH is DPG, not APG. No EDH selected.
  Added 1 Hz synthetic validation; that implementation checkpoint had **143 passed**.
- User asked where Ruan (2014) recorded 7D gain factors. They estimated per-
  station scale factors, but no verified original numeric factor table located
  in publisher/author PDF, related studies or indexed public repository searches.
  Do not claim they did not compute them or that no public record can exist.
  Bell/Ruan/Forsyth 2015 DOI 10.1785/0220150101 corrects the 2014 STARTING profile
  c=.029 (not .29); do not confuse with fitted c or pressure gain. Source notes:
  `docs/7d_xo_runs.md`. No gain/model/raw changes during this search.
- Combined latest workflows: **four networks, 89 distinct stations, 46 provisional
  gains, 41 unresolved and two reference-band exclusions**. Repeated runs count
  once. The reporting check rederives saved medians, reciprocals and counts,
  checks workflow/result hashes after reading, and confirms blank scientific
  scores in the original and both new XO masters. No batch/refit or gain
  application. This is not a fresh comprehensive raw-archive audit.

- YV: all 172 channel PPSD/base products and pressure diagnostics completed.
  Original Rayleigh run: 43 stations, 16 provisional candidates, 25 unresolved,
  two conservative reference-band exclusions. The user's separate no-tilt run
  has the same station counts. **Do not rerun the completed YV batch.**
- YL: downloaded daily MiniSEED names and colocated StationXML are supported.
  Signed negative pressure sensitivity is accepted without changing its polarity.
  First-five products remain preserved; all-station preparation reuses complete
  PPSDs only when inventory hashes and source metadata match.
- The user completed **all 25 YL stations**: **22 provisional gains**, three
  unresolved (**A14W/B01W/C10W**). The workflow is
  `results/NETWORKS/YL_ALL_NO_TILT/_workflow/rayleigh.json`.
- Calibration tilt is optional: explicit `--tilt-method none` or network
  `--rayleigh-tilt-method none` uses P/Z only. Historic default remains `tiskitpy`.
- Event sampling defaults to **all eligible events**. Three events is the current
  minimum for a station estimate, not a cap. Explicit `bounded` retains old limits.
- Optional `--earth-model-file` / `--rayleigh-earth-model-file` generates
  fundamental Rayleigh phase/group curves using disba 0.7.0 and Numba 0.68.0.
  **No verified numeric Zha/Lau model has been obtained.** Bundled RR52 layers
  must not be presented as Zha's model.
- **No calibration gain was applied. All six scientific scores remain blank.**
  YL calibration did not update masters. Raw data, XMLs, masters, previous results
  and the original five-station comparison were preserved.

## Short Zha comparison

Read `docs/yl_all_zha_comparison.md`; full values are in its CSV/JSON sidecars.
`scripts/compare_yl_zha.py` audits saved fits without refitting or changing masters.

Our `g` multiplies pressure. Comparing **1/g** with Zha & Webb (2016) Table 1
aligns correction direction, but does not prove identical fitting objectives or
nominal response versions. **21 stations are 3.4–7.0% lower (median −5.4%).**
**C07W is +86.8%: 3.716 versus 1.99.** Its already accepted M>=7 subset still gives
3.708. Individual event gains range .259–.543, with interleaved dates; do not
assume a simple gain step or sensor failure. Subsequent independent evidence
supports intermittent relative vertical-amplitude states; the hardware mechanism
is not established. Listed nominal P/Z sensitivities alone do not explain it.

A14W has no accepted events, chiefly coherence/timing/window rejection. C10W has
only one. B01W has no accepted events and four search-boundary rejections; the
published coefficient .19 implies our pressure factor **5.26**, beyond the
current search upper bound **3**. Subsequent raw checks find zero phase-valid
bins in those four events: phase is binding, so widening the bound alone will not
resolve them. Their original saved labels/results remain unchanged.

Current YL run: no tilt, 0.03–0.07 Hz, borrowed historical RR52 dispersion,
extra pressure sign +1, all eligible events under existing selection/QC.
The paper uses 0.04–0.1 Hz, regional Lau dispersion and a large-event set.
Model, response, events, windows and averaging still need reconciliation.
Do not tune gains to the publication or claim final calibration acceptance.

## Deferred response/Rayleigh work — historical priorities

Sentinel and Socrates completed their independent reviews; Codex's integration
verdict and limits are in `docs/reviews/codex_integration_review.md`.
User previously requested code for a C07W/C08W/C09W saved-event
comparison. New standalone tool: `scripts/compare_yl_matched_events.py --show`;
commands and limits: `docs/yl_matched_event_review.md`. Default includes all shared
saved spectral events, preserving rejection labels. The M>=7 runtime check found
six shared events and one all-three candidate; 105 protected hashes stayed intact.
No raw samples, new fits or master writes; next inspect original arrival links.

1. Investigate C07W using matched saved event spectra, response conventions and
   event windows; reconcile the published method before selecting a new run.
2. Investigate B01W's frequency-dependent phase/response mismatch; bounds alone do not resolve it.
3. Review A14W/C10W usable windows and timing; never silently relax safeguards.
4. Review LT02 logger timing and WD48 coverage. Obtain Ruan's fitted scale factors
   and Deng coefficient direction/band values before numerical comparisons.
   Validate local dispersion/depth and held-out uncertainty.
5. YV response/window validation remains unfinished. Consult
   `docs/yv_response_window_validation.md`,
   `docs/yv_rayleigh_independent_validation.md` and
   `docs/teleseismic_dpg_calibration_report.md` only as needed.
   Positive-real poles, absolute Z calibration and mode isolation are unresolved.

Later work includes held-out validation/uncertainty, microseism/local-P checks,
station classification, other networks and accepted-spectrum ocean-bottom noise
models. There is no measured whole-project completion percentage.

## Verification and live state

Latest implementation checks: **167 tests passed**, Ruff check/format and Bash
runner syntax passed. Git and all eight saved network workflows were inspected;
no BRUIT-FM preparation/calibration process was live at the reporting checkpoint.
The APG changes add analytic gain/phase/alias checks, gap/overlap/epoch rejection,
explicit rate-roundoff handling and calibration without PPSD arrays. The previous
read-only audit checked 222 PPSD channels / 980,630 segments,
2,219 saved Rayleigh spectra across 116 workflow station runs, other diagnostics
and every CSV/Excel value. No arithmetic/provenance mismatch was found. All 50
YL P/Z products cover their full available intervals. All 93,784 baseline source/
production file stats stayed unchanged; current master and protected product
hashes match. Historical YV master hash changes are the authorized YL append.

Five software fixes: preparation no longer mistakes bounded complete PPSDs for
full-archive products; pressure selection shares case/alias handling; network
preflight honors explicit depth; cache dependencies include Lat/Lon; auxiliary
H1/raw/median-curve fit failures preserve valid primary results with recorded
diagnostic reasons. Primary gain bounds/QC/model/sign/timing policy did not change.
No saved fit or workflow was edited and no science batch was rerun. Code changes
invalidate cache fingerprints on future computation; use viewers for old results.
Follow-up: 38 channels have 684,714 saved values outside the histogram range.
RR12 pressure/all four RR50 channels include numerical zero-power floors near
−3076.53 dB; these are not ultra-low measured noise. Histogram bins collect
under/overflow at their edges. Preserve originals and review screening before
accepted noise models. This did not alter raw-event calibration gains.
Reinspect Git and workflows before acting. Earlier RR52 GUI processes are
historical observations; do not assume old PIDs remain current or stop a viewer.

## Constraints and paths

Keep quota low: targeted reads, concise output, **no subagents**, no repeated full
history scans or full batches. The user prefers to run long computations.
Preserve raw waveforms, StationXML, saved spectra/results, masters, scores and
Excel styles. Keep alternatives in separate folders. Commit/push completed work.
Code, docs and curated comparison artifacts go to GitHub; raw data, masters and
full results remain local and excluded. Do not force-add them.

Python: `/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python`.
Data: `/Users/mohammadamin/Desktop/Data_Bruit_FM/Data/{YV,YL,7D,XO}`.
YL XML: `Data/YL/<station>/inventory.xml` beside downloaded waveforms.
Masters: `exports/master_inventory.csv` and `.xlsx`.
All-YL preparation/run command: `bash scripts/run_yl_all.sh` (already completed).
Commands/settings: `docs/yl_all_station_run.md`.
Model schema/provenance: `docs/rayleigh_earth_model_options.md`.
Original project scope: `docs/project_context.md`.

View a saved YL plot without calibration computation:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network view \
  --network YL --station A02W --view-kind rayleigh \
  --view-workflow results/NETWORKS/YL_ALL_NO_TILT/_workflow/rayleigh.json
```
