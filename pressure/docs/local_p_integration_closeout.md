# Local-P final integration and continuation handoff

16 October 2026. Owner request: cancel the long batch, finish testing, review
Claude's work, document and push. This updates the earlier
[phase report](local_p_phase_report.md), which deliberately kept main pinned.
Engineering integration is complete; absolute scientific acceptance is not.

## Code and peer-review disposition

Main integrates `codex/local-p-finish`, including Claude's initial implementation
and Codex's review/hardening. The original `claude/local-p-check`,
`claude/local-p-science-review`, `codex/local-p-review` and integration branches
are preserved. No branch was deleted or force-pushed. Remote fetch confirmed
the latest Claude tips were still `438c4ca` and `7ebd872`.

The [peer science review, including its outcome addendum](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/7ebd8726a759c5fb3071352b8dfb88f79522edd7/docs/reviews/local_p_science_review.md)
and [late implementation changes](https://github.com/MohammadAmin-Aminian/bruit-fm/compare/c828e56...438c4ca)
were read. Disposition:

| Contribution or claim | Final decision |
|---|---|
| Earlier equal-length noise placement | Retained in reviewed v2; conservative floor-to-grid boundary rather than rounding |
| LT02 independent implementation | Saved results match exactly: 0.967805 (0.5–1 Hz, n=4), 0.972951 (1–2 Hz, n=7), seven accepted events |
| Earlier-noise control diagnostic | Useful bounded sensitivity evidence; remains on Claude's branch, not a new acceptance rule or required extra waveform read |
| Initial prediction of less than 1% median change | Superseded by full-QC outcome: −4.72% at 1–2 Hz when two events leave; retain membership sensitivity |
| Timing study | Retain timing-limited interpretation; nominal vertical response and refined timing do not provide an independent reference |
| One response-aware leakage synthetic | Evidence for that configuration only, not proof for every response, depth or field wavefield |
| Deep-water absolute calibration | Not accepted; numerical pre-multiple results remain diagnostic-only |
| Restricted export/resume/output/master handling | Keep the tested integration fixes; no wholesale replacement by the late peer branch |

These are development checks on the same data/model, not independent field
calibration. The original peer report remains preserved at its cited commit;
its preliminary summary must be read with its final addendum and this disposition.
Its proposed promotion thresholds are research proposals, not adopted policy.

## Interrupted all-methods run

Session in the main checkout:
`results/EXPERIMENTS/ALL_CALIBRATIONS/run_20261016T144245Z/`.
Nested experiment:
`results/EXPERIMENTS/THREE_METHODS/run_20261016T144248Z_2515/`.
Its saved microseism settings identify source commit `ace6ee3`.

The validated terminal job's process group was interrupted with SIGINT, and
exit of the wrapper/children was verified before moving main. No output was
deleted, rewritten, resumed or promoted. The old wrapper did not record an
explicit cancellation marker; absence of its final session summary alone is
not a reliable completion test. This note records the observed interruption.

Under the nested experiment's `microseism/` directory:

| Interval | Folder | Observed state |
|---|---|---|
| Base | `run_20261016T144432Z_536542f9` | Summary and method comparison present; stage 1 finished |
| −12 weeks | `run_20261016T145753Z_90e7e1bb` | Summary/comparison and wrapper completion line present |
| −6 weeks | `run_20261016T151040Z_c36ae13d` | Summary/comparison and wrapper completion line present |
| −2 weeks | `run_20261016T152342Z_601f267d` | Summary/comparison and wrapper completion line present |
| −1 week | `run_20261016T153718Z_834e7d85` | Partial; no run summary or method comparison at interruption |
| +1, +2, +6, +12 weeks | None | Not started in this session |

The base summary has 89 station rows (84 checked, five unresolved pairing).
"Checked" is execution status, not calibration acceptance. The campaign-level
frozen combination and local-P stages were not reached. Do not label this a
completed all-methods validation or replace the existing frozen microseism
review with it. The earlier completed campaigns remain the comparison sources.
Do not rerun completed weeks; any future continuation needs explicit scope and
must exclude or separately handle the partial −1-week directory.

## Final verification

Verified on code `7dde4d4`; subsequent integration edits are documentation only.

- Combined tests: **332 passed, 78 existing warnings, 62.50 s**.
- Ruff check passed; format check clean for 158 files; three runner scripts
  pass `bash -n`.
- The all-methods `--check` preflight passed without waveform execution:
  66 local-P eligible pairs, 39 sample-rate exclusions, 15 response exclusions.
- Fresh read-only saved-campaign composition:
  `bruit-fm-localp-finish/results/LOCAL_P/reviews/main_closeout_20261016/`.
  All 144 protected input hashes match. Its `phase_report.json` SHA-256 is
  `798eda5c2a5dcc67701476ed856cb95142755e29967436d10e1ddb93800ffe4f`,
  exactly the earlier phase report. The import plan is read-only: 66 exact
  diagnostic joins, 54 skipped; no production write.
- Production master SHA-256 values still match the curated phase report:
  CSV `d929f43ad006decd638ac64b9bb9388f5c0b9d82e9402c82a33b5c9b3e8c45f8`;
  XLSX `7041a9d0142d62764a591d3af62e56cc3acd5ea6a7b8ee4701d3ea470c3335e5`.
- To make main's local preflight usable, the three saved catalogue CSVs and
  their JSON provenance sidecars were copied byte-for-byte into its previously
  absent `results/LOCAL_P/catalogs/`. No network refetch or source modification;
  all six source/destination hashes match. Results/catalogues remain Git-ignored.

No new real-data calibration batch was run. The prior bounded v2 APG validation
and private-master-copy test remain valid evidence, not newly repeated tests.
No original waveform, response, saved result, master, scientific score or fixed
Rayleigh reference was changed. The staging still covers 120 pairs/89 stations,
with only three v2 pairs; 63 legacy executed pairs remain explicitly flagged.

## Safe commands and next task

These commands run tests/preflight only, not a data campaign:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m pytest -q
bash scripts/run_all_calibrations.sh --check
```

The preflight's "nothing was run" means no waveform analysis; it does perform
imports and archive/metadata planning. Omit `--check` only when intentionally
launching a new campaign. Defaults reuse saved Rayleigh and skip new shifted
weeks, but still compute a new base microseism and local-P campaign. The runner
does not process every archive day/event. Do not launch it simply to verify the
merge. Historical results live in their original worktrees and remain local.

Next concrete project task: propose and agree a versioned six-check
acceptance/scoring policy, including unperformed/not-applicable/conditional/
diagnostic/unresolved distinctions and minimum evidence. Then implement its
dry-run table/classification export before any production score write, followed
by the accepted-station/noise-model pipeline. This closeout does not adopt
thresholds, classify instruments or promote candidates. Claude's useful next
research is independent sensor/reference/timing evidence, not another duplicate
method or a relaxation to force agreement.
