# Codex review of Claude's full-data coverage workflow

Reviewed contributor tip: `87ff816`, branch `claude/full-data-coverage`, based
on main `db9b045`. Contributor history/branch is preserved. Review date:
19 October 2026. This accepts engineering integration, not calibration gains,
station ratings, scientific acceptance or a validated pressure-noise model.

## Scope and disposition

Retain the archive planner, explicit product/evidence sources, isolated outputs,
campaign locking, per-task resume/logs, failure isolation and pressure heatmaps.
Existing estimators and scientific defaults are unchanged. Sentinel/Socrates
audit qualifications remain applicable; see `codex_integration_review.md` and
the limitations table in `../full_data_coverage.md`.

Additional review fixes:

1. Require exact PPSD source-epoch boundaries for automatic reuse. An overlapping
   epoch with an equivalent response does not authorize mixing its PSD windows
   into a different epoch. Report `review_required`, preserving the old product.
2. A terminal empty/insufficient-data task is not reported as covered input time.
   Coverage tables explicitly distinguish processed input spans from continuous
   valid PSD support. Three regression cases cover epoch mismatch and the two
   empty-result states.
3. Header-cache writes use unique temporary names so concurrent planners do not
   race over a shared `.tmp` file. This does not authorize concurrent execution;
   the campaign execution lock remains mandatory.

## Real saved-data smoke and preservation

Only one missing expected-pressure comparison was executed, for
`YL.A02W..LDH`, on an isolated verified copy of its existing nominal-response
PPSD. No real PPSD, Rayleigh, microseism or local-P batch was started.

Local campaign: `results/COVERAGE/codex_review_20261019/` in the review worktree.
First run: `runs/20261019T093803Z_ac046156_1d9ba7/run.json`.
Result: `comparison_only`; one task completed and 197 tasks excluded by the
explicit kind/seed filter. All **1,056 hashed protected files** and **97,519
waveform size/mtime records** remained unchanged. These checks are not a
byte-level hash audit of all raw waveforms. The archive scan reused Claude's
copied header cache; it did not recompute the saved PSDs.

Repeat invocation: `runs/20261019T093847Z_ac046156_eb5a13/run.json` reported
one `already_complete` and 197 excluded tasks, with protected inputs unchanged.
The saved-data comparison was not executed a second time.

The real plan remains **198 missing tasks**: 15 pressure PPSDs and 61 each of
expected-pressure, nearby and infragravity comparisons. There are 19 blocked
pressure epochs. These count channel epochs/tasks, not new stations. Historical
calibration evidence is reused/deferred; this runner does **not** automatically
fill missing Rayleigh/microseism/local-P gains or update legacy results to v2.

## Verification

- Claude's original 22 coverage tests pass independently.
- Reviewed full suite: **413 passed, 110 warnings, 105.92 seconds**. The
  warnings are ObsPy synthetic-unit/encoding warnings, not failed assertions.
- Ruff check and formatting checks pass; shell runner syntax passes.
- Three added regression cases pass for mismatched source epochs and both
  terminal empty-result statuses. Resume, output integrity, campaign locking,
  protected source preservation and worker-process execution are exercised by
  the inherited synthetic tests.
- Real saved-PPSD execution/resume is verified above. New real raw-data PPSDs
  and new real merge/nearby tasks are **not** claimed tested by that smoke.

## Owner's next step

No long computation is needed now. Use existing data/results for descriptive
comparisons while waiting for the new computer. Later, configure the actual
archive/source paths, run a dry run, inspect blocked prerequisites, then execute
only planned missing work. No command should be described as running every
calibration method automatically.

From the main workspace after integration:

```bash
bash scripts/run_full_coverage.sh plan
bash scripts/run_full_coverage.sh dry-run
```

Later, on the properly configured computation machine:

```bash
bash scripts/run_full_coverage.sh execute --workers 4
bash scripts/run_full_coverage.sh report
```

Reusing the same campaign on the same machine resumes it. Migration to a new
machine still needs path/provenance reconciliation; absolute file paths and
mtime-bound source records are not a portable resume guarantee. Do not blindly
copy task states or present the archived dry-run as a current machine's plan.
