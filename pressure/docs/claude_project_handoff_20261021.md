# BRUIT-FM: project transfer from Codex to Claude

> **Status update, 22 October 2026:** the immediate task below (response
> confirmation) and the YV reproduction are done. Read
> [project_status_20261022.md](project_status_20261022.md) first for the current
> state, the open owner decisions and the new-computer runbook.

21 October 2026. The owner is transferring ongoing project work to Claude.
This is the current continuation guide; older task briefs are historical when
they conflict with it. This transfer changes documentation only, not science.

## Repository and verified software checkpoint

- Repository: <https://github.com/MohammadAmin-Aminian/bruit-fm>.
- Main workspace: `/Users/mohammadamin/Projects/bruit-fm`.
- Reviewed software checkpoint: `7ac2e639e84268f1ef2317412d5e1ef4537874af`.
  This handoff is committed after that checkpoint; use the latest main.
- [PR #1](https://github.com/MohammadAmin-Aminian/bruit-fm/pull/1) is merged and
  closed. A fresh remote check found no open PRs and no local branch commits
  absent from all origin refs. Contributor branches/worktrees are preserved.
- [Post-merge main CI](https://github.com/MohammadAmin-Aminian/bruit-fm/actions/runs/39587156168)
  passed. Full local suite at the software checkpoint: **447 passed**, 117
  warnings; Ruff/format/whitespace checks passed. Wheel installation outside the
  checkout and nine safe module-help entry points passed locally and in CI.
- This is installable/tested research software, **not** completed all-data
  processing, accepted absolute calibration, adopted ratings or a global noise
  model. No service/release was published.
- No matching calibration/coverage batch process was found during transfer.
  Recheck live processes before future changes or execution; historical output
  folders and missing summaries alone do not establish live/completed status.

## Read in this order

1. This handoff, then `CODEX_NEXT_CHAT.md` for detailed recent chronology.
2. `docs/reviews/yv_rayleigh_diagnosis_integration.md`, then
   `docs/yv_rayleigh_diagnosis.md` (read the review's qualifications first).
3. `docs/reviews/deployment_readiness_20261021.md` and `docs/deployment.md`.
4. As needed: `docs/socrates_all_station_comparison.md`,
   `docs/yv_heldout_results.md`, `docs/yv_heldout_arrival_review.md`.
5. Before other methods: `docs/microseism_next_chat.md`,
   `docs/local_p_integration_closeout.md`, `docs/six_check_review.md`.
6. Before the future long campaign: `docs/full_data_coverage.md` and
   `docs/reviews/full_data_coverage_integration.md`.

Do not reread the whole repository or restart old briefs. `CODEX_HANDOFF.md`
is additional historical context, not the latest task instruction.

## What is implemented and reviewed

- PPSD production, response handling, protected-input checks, plotting and
  network execution; expected-pressure, nearby and infragravity diagnostics.
- Rayleigh pressure/vertical calibration, saved complex spectra, isolated
  opt-in response/fitting experiments, and frozen-gain held-out evaluation.
  Sentinel engineering changes and the Socrates review/trial dispositions are
  already integrated; do not blindly remerge older branches.
- Model-free microseism P/Z diagnostics plus Claude's forward-model check,
  combined with frozen frequency/week consistency criteria. P/Z alone does
  not identify an independent pressure gain without vertical/model constraints.
- Reviewed local-earthquake P-wave v2 check: response/epoch identity, strict
  pre-multiple timing/noise placement, tied-polarity and output guards. Most
  old results remain explicitly legacy-unvalidated; deep-water gains remain
  diagnostic, not automatically accepted broadband calibrations.
- Exact-identity six-check evidence review, draft categorical policy, explicit
  approval records and approval-gated station-balanced noise-model machinery.
  No real approvals were fabricated; accepted-channel export is empty.
- Resumable full-coverage planner/executor: reuse compatible saved products,
  compute only missing PPSDs and spectral diagnostics, isolate outputs, keep
  provenance and original masters. It does not fill missing calibration gains.
- Read-only YV Rayleigh diagnosis with exact P/Z epoch/source binding, completed
  experiment source-hash checks, per-spectrum integrity and before/after guards.
  Replay covered 43 stations and 583 saved event NPZs, reproducing prior metrics.
- Latest software fixes: versioned circular L1 phase summaries, null rather
  than fabricated phase for directionless inputs, and refusal of negative
  fractional native-grid tears. Per-bin QC, gain objectives and frozen gains
  were not changed. A new code fingerprint is not a reason to rerun old jobs.

## Frozen numerical availability, not calibration acceptance

Baseline pressure multipliers from `docs/calibration_network_summary.json`:

| Network | Stations with a baseline estimate / available | Median g |
|---|---:|---:|
| YV | 16 / 43 | 0.949350 |
| YL | 22 / 25 | 0.602219 |
| 7D | 4 / 15 | 0.904638 |
| XO | 4 / 6 | 0.932082 |
| All networks | 46 / 89 | 0.730075 |

These medians give equal weight to each non-null station estimate, not each
event. `g` multiplies nominal pressure; a sensitivity multiplier is `1/g`.
Do not invert a median and call it the median reciprocal without calculating it.

The completed experimental pressure-pole trial has **31 / 43 YV estimates**,
median **1.088702813782616**, ten ordinary unresolved stations and two reference
band exclusions. The nominal baseline has 25 ordinary unresolved and two
exclusions. No physical response correction is accepted by these counts.

Across networks, the historical trial summary's **61 / 89** explicitly includes
C07W's separate three-event normal-state diagnostic. Strict non-null
`summaries.pole_absolute.median_pressure_multiplier` entries number **60 / 89**;
C07W's mixed-state scalar is null there. Do not conflate these selections.

Frozen microseism consistency review: seven conditional estimates overall
(two YV, two YL, three 7D). Current local-P v2 has one conditional station, LT02;
63 legacy executed pairs remain flagged. Execution/availability is not pass/fail.

For YV's original baseline comparison, RR40/RR52 have numbers from all three
methods; RR11/RR36/RR38/RR50 from two. Rayleigh and microseism are not independent
when they share the same nominal Z response and modelling assumptions. Local-P
gains must remain frequency dependent, not collapsed into a single scalar.

## Immediate next task: response confirmation, not another implementation

Start with **RR22, RR56, RR17 and RR18**. The metadata comparison has **not yet
been performed**. Explain the bounded plan briefly to the owner before acting.

1. Inspect exact local P/Z response epochs and saved source hashes. Retrieve
   official public response metadata only within the owner's approved scope;
   save fresh downloads in a new review folder, never over `YV.xml`. Record
   provider, exact query, retrieval UTC, response bytes/hashes and failures.
2. Compare roots, conjugacy, units, stage gain, sensitivity, normalization,
   sample rate and epoch. Check Stähler et al. (2016), Appendix A/Table A1
   against the exact stations, not just a generic sensor family.
3. Distinguish an archive agreeing with local XML from independent physical
   calibration evidence. Manufacturer/instrument-specific records can still
   be needed. If unavailable, document the limitation; do not invent a pole.
4. Deliver a concise station-by-station evidence table and provenance report,
   with any reusable comparison helper/tests. Propose a bounded saved-spectrum
   pilot before fitting gains or changing production responses/defaults.

No raw waveform or full-network rerun is needed for metadata comparison.
Do not contact institutes or send messages without the owner's authorization.

### Existing leads and limits

- RR22/RR56: local vertical XML has duplicated unpaired complex poles.
  Published Table A1 prints conjugate pairs. H2 restores the pairs after the
  pressure-pole hypothesis, reducing inspected phase residuals; it does not
  verify station-specific hardware/normalization/absolute sensitivity.
- RR17/RR18: unusually large positive pressure poles. H4 (reflect and divide
  that pole by ten) improved inspected phases but was selected **post hoc**.
  DEPAS/AWI instrument-specific resistance/capacitance records are needed to
  establish a physical pole. Do not promote H4 as a measured correction.
- RR47: relative P/Z inversion hypothesis; no evidence yet identifies which
  physical channel is inverted. Independent polarity evidence is needed.
- RR53: the same generic vertical repair does not resolve the phase; hardware
  and response identity need confirmation. RR55 also has hardware ambiguity
  despite a numerical trial gain.
- RR12: 4,769 / 4,778 saved pressure PPSD window medians are below -3000 dB.
  This is an archive/export/source issue to investigate, not proven failure of
  the deployed hydrophone. RR31 has published artefacts but usable intervals
  may exist. RR42/RR45 have coverage/vertical limitations. RR41/RR57 are
  reference-band exclusions, not ordinary gain-fitting failures.

After the metadata review, seek owner agreement on fixed hypotheses, saved
events, unchanged QC and fresh output paths. **RR52 stays unchanged** at
`g = 0.6787278140712512`. Do not tune other methods to it. Keep previously
untouched events for independent validation; inspected cases are development,
not blind tests. Existing RR52 held-out evidence is 14 isolated usable events,
median prediction error +0.664%, conditional on the shared model.

## Exact saved sources and local-only data

Paths below are relative to the main workspace unless absolute:

- Nominal YV workflow:
  `results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json`.
  Do not substitute the older tilted `TEST/_workflow/rayleigh.json`.
- Completed trial:
  `results/EXPERIMENTS/SOCRATES/all_20261013T195356Z_00cf3343/experiment.json`;
  SHA-256 `446b8e36b94c9952736b9a2c9b2353f0dfb0c95e05e8cc3719b9cea72fcc20cf`.
  Individual `YV_RRxx/trial.json` files retain variant-specific outcomes.
- Claude's original diagnosis:
  `docs/yv_rayleigh_diagnosis/evidence.json` and `stations.csv`, unchanged.
  Reviewed local replay:
  `results/REVIEWS/YV_RAYLEIGH_DIAGNOSIS/codex_review_20261021/`.
- Frozen combined review:
  `results/REVIEWS/SIX_CHECKS/final_20261019/review.json`;
  SHA-256 `241e9ab125a718fe92958b2129a267ebeff43076287142d85c0beb6d8cea1fd2`.
- Completed held-out studies:
  `results/VALIDATION/YV_HELDOUT/run_20261013T213126Z_87e24fe6/` and
  `results/VALIDATION/YV_HELDOUT/RR05_replacement_20261013T215958Z_fee6a2d2/`.
- Additional microseism sources remain in the preserved
  `/Users/mohammadamin/Projects/bruit-fm-claude/results/` worktree. Use explicit
  pointers in the frozen review; never pick the newest directory by glob.
- Raw archive: `/Users/mohammadamin/Desktop/Data_Bruit_FM/Data/`;
  original inventories: `/Users/mohammadamin/Desktop/Data_Bruit_FM/INVENTORIES/`.
- Primary paper: Stähler et al. (2016), *Performance report of the RHUM-RUM
  ocean bottom seismometer network around La Réunion, western Indian Ocean*,
  <https://doi.org/10.5194/adgeo-41-43-2016>. Relevant Table 3/section 3.1 and
  Appendix A/Table A1 were checked. Local PDF:
  `results/LITERATURE/PRESSURE_REFERENCES_20261019/Staehler_2016_RHUM_RUM.pdf`.

Git contains code, documentation and curated evidence, **not** the full private
archive, masters, cached catalogues or saved campaigns. On another computer,
missing local files are a missing-input condition, not authorization to refit
or regenerate everything. Reconcile paths/provenance explicitly; do not edit
historical JSON or copy completion flags blindly to make validation pass.

## What remains later

- Complete justified response/hardware identification and a scoped pilot, then
  genuinely independent validation and interpretation of remaining failures.
- On the new machine, reconcile archive/source paths, inspect a coverage dry
  run, then have the owner deliberately execute only missing work. The saved
  pressure-scope plan has 89 reused pressure PPSDs, 15 new pressure PPSDs plus
  61 each of three spectral checks (**198 tasks**), and 19 blocked epochs.
  These are tasks/epochs, not newly calibrated stations. One isolated saved
  YL.A02W expected-pressure smoke already completed; do not repeat it needlessly.
- Only after coverage/distributions and scientific agreement: adopt evidence
  requirements/ratings, record real approvals, and build a scoped noise model.
  The owner explicitly deferred ratings and thresholds; do not turn draft
  policy into adopted policy or fill scores to manufacture a finished table.

The separate October 16 all-methods batch was interrupted: base, -12, -6 and
-2 week microseism stages finished, -1 week partial, later weeks not started;
campaign-level combination/local-P stages were not reached. Preserve it, do not
call it complete, and do not automatically resume or promote it.

## Non-negotiable workflow constraints and useful commands

- Preserve original waveforms/headers, XML, response snapshots, masters CSV/XLSX,
  saved spectra/results, fixed gains, scores, branches and worktrees. Never
  apply calibration gains to these originals. New experiments need isolated,
  explicitly labelled outputs and before/after integrity checks.
- The untracked `gain_histogram.png` at the main root is the owner's personal
  file: do not delete, overwrite, add or clean it. Stage named task files only.
- Keep quota low, read targeted files, explain the next concrete task briefly,
  and use **no subagents**. No heavy campaign is authorized by this handoff.
- Keep the original scientific fitter/defaults. Weighted fitting, pole repair,
  new phase-family thresholds, lower Godin bands, depth normalization and
  polarity flips are not silently adopted fixes.
- Pressure PPSD figures: frequency-domain probability heatmaps, ObsPy PQLX,
  no mean/median overlays by default; show all Wayne-supplied markers, but do
  not treat unverified handout values as numerical acceptance. Brown is an
  in-water reference, not seafloor instrument bounds. Read
  `docs/pressure_reference_review.md` before reference changes. A display-only
  gain shift is `20*log10(g)` dB, not `10*log10(g)`.
- Preserve contributor attribution. The formula branch was selectively
  integrated; untested fixed-output plotting helpers stayed on its original
  branch. Old local-P branches must not replace main's stronger reviewed v2.

Current Python:
`/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python`.
Use `PYTHONPATH=src:.` for checkout scripts. Supported deployment is Python
3.13 on macOS/Linux; wheel smoke is not a locked dependency or field test.

Safe software verification, when relevant to actual changes:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src:. /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python -m ruff check .
PYTHONPATH=src:. /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python -m ruff format --check .
PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python -m pytest -q
```

Create a new branch from current `origin/main`, e.g.
`claude/yv-response-confirmation`. Read any existing changes before switching;
do not reset/clean another worktree. Push focused code/docs/tests, report tests
actually run and unresolved evidence, and open a PR. Do not force-push, delete
contributor branches or merge scientific changes without the owner's decision.
