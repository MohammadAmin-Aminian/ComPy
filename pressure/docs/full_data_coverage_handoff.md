# Handoff: full-data coverage workflow (for Codex review)

**Reviewed continuation:** read
[Codex's integration review](reviews/full_data_coverage_integration.md) for
additional fixes, the bounded saved-data smoke, verification and the owner's
decision to postpone the long run until the new computer arrives. The account
below preserves Claude's original contributor checkpoint.

Branch `claude/full-data-coverage`, based on main `db9b045` (still the latest
origin/main when work started). Worktree:
`/Users/mohammadamin/Projects/bruit-fm-coverage`. Not merged into main.

Method, decisions and current numbers: [full_data_coverage.md](full_data_coverage.md).
Curated real dry-run outputs: [coverage/](coverage/).

## Changes

| File | Purpose |
|---|---|
| `src/bruitfm/coverage.py` | Archive/header inventory, StationXML epochs, response identity, saved-product matching, evidence joins, deterministic task plan, CLI |
| `src/bruitfm/coverage_run.py` | Dry run / `--execute`: campaign lock, task state and logs, resume, failure isolation, symlink input view, hash-verified product views, PPSD/merge/check execution, protected-input verification |
| `src/bruitfm/coverage_report.py` | Coverage tables; pressure probability heat maps (window-pooled and station-balanced); curated export |
| `scripts/run_full_coverage.sh` | `plan`, `dry-run`, `execute`, `report` wrapper |
| `docs/coverage_sources.json` | Explicit sources and PPSD parameters (equal to the `network ppsd` defaults) |
| `tests/test_coverage.py` | 22 tests (below) |
| `docs/full_data_coverage.md`, `docs/coverage/` | Method, results, limitations; curated tables and 14 detail heat maps |

Existing modules are unchanged. The workflow calls the following existing
code rather than re-implementing it:

- `network.ppsd_job` and `ppsd.process_epoch`: PPSDs.
- `network.diagnostic_job` → `run_pressure_test`, `run_nearby_test` and
  `run_infragravity_test`: the three spectral checks.
- `quality_review` loaders: historical calibration evidence.
- `plot_network_ppsd.accumulate_histogram`, `spectrum` and `add_references`:
  the heat maps.

## Commands

Read-only inventory and dry run:

```bash
cd /Users/mohammadamin/Projects/bruit-fm-coverage
bash scripts/run_full_coverage.sh plan
```

Execute all planned missing work (owner runs this):

```bash
cd /Users/mohammadamin/Projects/bruit-fm-coverage
bash scripts/run_full_coverage.sh execute --workers 4
```

Afterwards, write the report (tables and heat maps including the new
products):

```bash
bash scripts/run_full_coverage.sh report --curated-dir docs/coverage
```

## Verification

- `tests/test_coverage.py`: 22 passed. The tests use synthetic 1 Hz SDS and
  daily archives, a StationXML, a partial saved product with a saved check
  result, registry and evidence files, and run the real PPSD and estimator
  code. They cover:
  - interval arithmetic and the filename layouts;
  - the `00`/blank location alias, and refusal of other header mismatches;
  - reuse of a complete product and its hash-bound saved result;
  - supplements planned only for uncovered intervals, then a merge that
    keeps the base windows;
  - explicit all-channel scope;
  - response or parameter differences → review, without recomputation;
  - epoch identity splitting coverage, and data outside epochs reported;
  - calibration evidence reused or deferred, never executed;
  - legacy local-P flagged, and `not_applicable` taken from saved rows;
  - missing configured sources refused, and a campaign inside protected
    inputs refused;
  - a dry run writing no task state;
  - an interrupted run resuming;
  - a repeated run executing nothing, with estimators patched to fail if
    called;
  - failure isolation with dependent blocking;
  - the nearby pool as a soft dependency, reported as `stale` and rerun
    when the pool changes;
  - the campaign lock refusing concurrent execution;
  - checks on reused products writing only to campaign views;
  - protected trees byte-identical;
  - heat maps keeping zero-power floors, with each station-balanced column
    summing to one;
  - heat maps split by declared sensor family; the keyword rule is checked
    on the real YV/YL/7D/XO labels;
  - PPSD tasks executed by two worker processes (the `--workers` path).
- Full suite: **410 passed** (388 existing + 22 new), 104 warnings (78
  existing; 26 ObsPy unit-mapping/encoding warnings inside existing
  estimators on synthetic responses). Ruff check/format and `bash -n` pass.

## What was executed, and what was not

**Executed (read-only on inputs, writing only to this worktree's
`results/COVERAGE/`):**

- `plan`, `dry-run` and `report` on the real archive: about 97.5k
  first/last-record header probes, plus full-header verification of the
  files of the 15 channels scheduled for computation;
- heat maps from the 89 saved pressure products.

**Not executed:**

- no PPSD or check task (the campaign has 198 not-started tasks). A bounded
  real smoke run, and then the full run, were attempted from Claude's
  session; the session's permission classifier refused them. The owner runs
  `execute` (the long part is about 9.1 × 10⁹ samples of PPSD input, mostly
  the four XO HDH channels);
- no merge, Rayleigh, microseism or local-P run;
- no waveform or StationXML fetch;
- no master, score, response or saved-product write;
- no gain application, threshold, rating, approval or noise model;
- nothing merged into main.

## Points for Codex review

1. **Code fingerprints.** Saved products are reused on provenance (epoch,
   response amplitude/phase, parameters, source files, covered interval),
   not on `network.code_digest`. Code/ObsPy differences are reported only.
2. **Primary sources.** For 7D the FDSN snapshots saved by the completed run
   are primary because `INVENTORIES/7D.xml` has no responses. Local-P evidence
   keyed to `INVENTORIES/7D.xml` is joined as
   `same_epoch_response_unverifiable` and flagged.
3. **Supplement and merge.** Supplement PPSDs use the base product's
   parameters and start at uncovered intervals. The merge is ObsPy
   `add_npz`, which omits already-covered windows. The real plan currently
   needs no supplement (J43A's is blocked by its 0.15 ppm rate), so this
   path is covered only by tests.
4. **Views.** The existing estimators resolve symlinks and write beside the
   PPSD, so reused products are cloned (`cp -c`, SHA-256 verified) into
   `views/`.
5. **Header probe.** The fast probe gives file spans. Channels scheduled for
   computation get a full-header verification of identity and rate.
6. **Nearby.** The reference pool is soft: failed members are excluded and
   listed, and the result reruns when the pool signature changes. The report
   labels such results `stale`. Saved YV nearby results are reused. The
   larger pool cannot change their candidate sets: there is no cross-network
   pair within 100 km, and the networks are ≥ 1,707 km apart.
7. **Real-data checks of the execution path.** The checks, views, symlink
   input view and merge have run only on synthetic data. The first real
   `execute` is therefore their first real-data exercise. Review its
   `runs/<stamp>/run.json`, which records protected-input verification, and
   any `failed` task logs.

## Decisions left to the owner

1. Obtain full 7D responses for horizontals, `BXH` and `LDI`? A fetch is not
   authorized by this task.
2. Allow the existing 0.1 ppm rate-encoding tolerance in PPSDs (XO LA33 HDH,
   7D J43A)?
3. Run `--ppsd-scope all-channels` (YL/XO horizontals, 7D BHZ, XO HHZ)?
4. Should alternate pressure streams (7D BDH, XO HDH) get calibration runs,
   under an explicitly frozen protocol?
5. After the run: examine distributions, coverage, method agreement and
   limitations before any acceptance rule (the owner's stated next step).

Resolved without recomputation: the YV nearby results do not need a rerun
against the four-network pool (see point 6 above).
