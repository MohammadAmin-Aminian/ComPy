# Full-data coverage: inventory, missing-work plan and resumable execution

Branch `claude/full-data-coverage` (from main `db9b045`). The owner decided
to finish data coverage **before** any rating, threshold, station score or
approved noise model. Nothing here adopts one. The draft six-check policy
stays draft, and no gain is applied anywhere.

## What it does

`python -m bruitfm.coverage` has three subcommands.

| Command | Writes | Purpose |
|---|---|---|
| `plan` | a new `plans/<stamp>/` folder and a header cache in the campaign | Read-only inventory and dry run |
| `run` | nothing (only the header cache) | Lists the tasks that would run and those already complete |
| `run --execute` | campaign `tasks/`, `products/`, `views/`, `_inputs/`, `runs/` | Runs only the missing tasks |
| `report` | a new `reports/<stamp>/` folder; optional `--curated-dir` | Coverage tables and pressure heat maps |

The campaign folder defaults to `results/COVERAGE/<campaign>/` in the checkout
that runs the module. `--source-root` names the checkout that holds the saved
`results/` and `exports/`, i.e. the main workspace. All configured sources must
exist, otherwise the command refuses to run. This prevents a wrong source root
from silently planning to recompute everything. The campaign folder may not
lie inside the data, StationXML, registry or saved-product folders.

### Sources (explicit, never "newest")

`docs/coverage_sources.json` lists every input:

- data directory;
- StationXML per network, in priority order;
- product registries (the three masters with PPSD pointers);
- product roots;
- historical evidence: the six-check review's masters, the Rayleigh summary,
  and the microseism and local-P staging tables.

Within a network, the epochs of the first file that defines a SEED ID are
primary. A later file is never added as a competing epoch. It is listed with
an explicit amplitude/phase response comparison.

`INVENTORIES/7D.xml` has no response stages. The full-response FDSN snapshots
saved by the completed 7D runs (P/Z only) are therefore listed first for 7D.

### Inventory

**Archive.** Every recognized daily file is classified by layout:

- SDS: `NET.STA.LOC.CHA.D.YYYY.DDD`;
- SDS with `.mseed`: the XO archive;
- daily download: `NET.STA.LOC.CHA.YYYYMMDD.mseed`.

The identity comes from the MiniSEED headers: the first and last record of
each file, cached by path, size and mtime. `--full-headers` reads every
record instead.

- The only filename/header difference accepted is the `00`/blank location
  alias already authorized by the 7D/XO preparation scripts. Any other
  mismatch, or an unreadable header, excludes the file and reports it.
- Several copies of one stream (SDS plus `waveforms/` downloads) are never
  merged where they overlap. The preferred layout is kept, and other copies
  only fill time that the preferred layout lacks.

**Epochs.** For each StationXML epoch the inventory records:

- whether the response is usable, using the existing
  `ppsd.response_handling`;
- the channel class from the declared units: pressure, vertical, horizontal
  or other;
- the available intervals from the file spans clipped to the epoch.

Data outside every epoch is reported as unmatched.

**Prerequisites for new computation:**

- a usable response;
- header sample rate equal to the StationXML rate within the existing PPSD
  code's `math.isclose` tolerance (rel 1e-9).

Every header of a channel scheduled for computation is read in full before a
task is created.

### Completeness of saved products

A saved PPSD is found only through a registry pointer or a configured root.
The match to an epoch requires the same SEED ID and an overlapping epoch.
It is **compatible** when all of the following hold:

- its snapshot response equals the primary response, in amplitude and phase;
- its channel epoch boundaries exactly match the primary epoch; an overlapping
  epoch alone requires review and explicit PSD-window selection, not reuse;
- its length, overlap, minimum frequency, dB bins and gap handling equal the
  plan parameters;
- its status is complete or partial;
- no recorded source file was modified.

The covered interval is the product's requested interval. Chunks that failed
and an interrupted tail are subtracted, and data in files the product never
read are added as uncovered, unless they are another copy of a file it did
read.

`product_covered_days` describes processed input spans, **not** continuous
valid PSD-window support or independent samples; gaps and rejected windows can
remain inside those spans. The table discloses this as `coverage_basis`.
A terminal `no_data`/`insufficient_data` task is retained for resume bookkeeping,
but it does not make its interval covered.

Some differences are reported but do not cause recomputation:

- code fingerprints (the old `network` workflow invalidated everything on
  any code change);
- the ObsPy version;
- being built from another inventory file with an identical response.

Each epoch receives one PPSD decision.

| Decision | Meaning |
|---|---|
| `reuse_existing` | One compatible registered product covers all available data |
| `supplement` | Compatible product, but data are uncovered. A supplement PPSD per uncovered interval (≥ one window), then a merge task (`PPSD.add_npz`, base windows preserved) |
| `reuse_existing_partial` | As above, but the supplement is blocked by a prerequisite |
| `compute` | No saved product: one new PPSD over the available interval |
| `blocked` | A prerequisite is missing; no task |
| `review_required` | Several, unregistered or incompatible products. Preserved, not recomputed |
| `not_planned_out_of_scope` | Non-pressure channel without a product; `--ppsd-scope all-channels` plans it |

### Checks

The six checks remain distinct from the three calibration approaches.

- **Spectral checks** (expected pressure, nearby, infragravity) use the
  existing estimators through `network.diagnostic_job`.
  - A saved result is reused when its registry row is the reused product's
    exact key and its `source_npz_sha256` equals the product hash.
  - Otherwise a task is planned. Supplemented epochs are checked on the
    merged product.
  - Nearby uses all available pressure products of the plan as its reference
    pool. Saved products are read in place.
  - The pool is a soft dependency. Failed members are excluded and recorded,
    and each result stores a pool signature. If the pool changes later, that
    nearby result reruns.
- **Calibration checks** (Rayleigh, microseism, local-P) are never executed
  here. Historical evidence is joined with the six-check review's exact-key
  loaders (`quality_review.load_masters`, `attach_rayleigh`,
  `attach_staging`). The join records one of these relations:
  - exact key;
  - same epoch / same response;
  - same epoch with the response unverifiable: the evidence key points to
    `INVENTORIES/7D.xml`, which has no responses;
  - response differs: not used.

  Legacy local-P rows keep their `legacy_unvalidated` state, and the reason
  names the missing v2 noise-window guard. An unjoined staging row that says
  `not_applicable*` gives `not_applicable`; any other missing evidence gives
  `deferred`. When another pressure stream of the same station has evidence,
  it is named and **not transferred**: metadata do not prove one physical
  sensor.

### Execution safety

- **Lock.** A non-waiting `fcntl` campaign lock refuses a second concurrent
  execution.
- **State and logs.** Each task has a deterministic ID, and per-task
  `state.json` and `log.txt`.
- **No repeated work.** A task is complete only while its recorded outputs
  still exist with the recorded SHA-256. A repeated invocation therefore
  runs nothing. A nearby result additionally needs its pool signature.
- **Interruptions.** Ctrl-C marks running tasks as interrupted. PPSD
  checkpoints are kept, and `process_epoch` resumes them on the next run.
- **Failure isolation.** A failure is recorded, and only its dependents are
  blocked.
- **Read-only inputs.** Waveforms are read through header-named symlinks
  (`_inputs/waveforms`), as in the 7D/XO preparation.
- **No writes into saved folders.** The estimators write outputs beside the
  PPSD they read, and resolve symlinks. Checks on reused products therefore
  run on hash-verified copy-on-write clones (`views/`).
- **No master writes.** No master is created or patched. Estimator updates
  stay in the task state.
- **Protected-input check.** Before and after execution, SHA-256 is computed
  for config, StationXML, registries (CSV+XLSX), evidence, saved products,
  snapshots, provenance and reused results, and size/mtime for every waveform
  file. The run report states whether anything changed, and the run fails if
  it did.

## Commands for the owner

Read-only inventory and dry run (writes only the campaign plan/cache):

```bash
cd /Users/mohammadamin/Projects/bruit-fm-coverage
bash scripts/run_full_coverage.sh plan
```

Execute all planned missing work (pressure scope, four PPSD workers):

```bash
cd /Users/mohammadamin/Projects/bruit-fm-coverage
bash scripts/run_full_coverage.sh execute --workers 4
```

Rerun the same command after any interruption; completed tasks are skipped.
The report's `tasks_status.csv` marks each task as follows:

- `complete`: still valid;
- `stale`: recorded as complete, but its outputs or (for nearby) its
  reference pool changed since; it reruns on the next `execute`;
- `failed`, `partial`, `interrupted`, `blocked_by_dependency` or
  `not_started`.
`dry-run` lists what would run. `report --curated-dir DIR` writes the coverage
tables and heat maps. Bounded variants: `--only-kind ppsd`,
`--only-seed XO.LA34..HDH`. Broader PPSD execution is explicit:
`--ppsd-scope all-channels` (vertical/horizontal channels with usable
responses).

## Current inventory (real dry run, 19 October 2026)

Curated outputs are in [`docs/coverage/`](coverage/); the full plan, task
list and report stay local in `results/COVERAGE/full_data_coverage_v1/`.
Plan `ac046156fadf1625`. No task has been executed.

**Archive:**

- 97,519 recognized waveform files (SDS 55,945; daily 33,496; SDS+`.mseed`
  8,078), 454 GB.
- 566 header streams at 89 stations (YV 43, YL 25, 7D 15, XO 6).
- 566 channel epochs with data: 123 pressure, 126 vertical, 252 horizontal,
  65 other.
- Excluded files:
  - 661 XO `RKO` state-of-health files whose headers carry `AN.R`/`SO.R`;
  - 16 7D J43A download copies with mixed sample rates;
  - 1,152 overlapping duplicate copies (7D/XO `waveforms/`).
- YV RR53: 11.2 h on each of four channels lies outside every inventory epoch.
- 6,316 StationXML epochs have no local data (mostly other 7D/XO stations).

**Pressure PPSDs** (123 epochs):

| Decision | Epochs | Detail |
|---|---:|---|
| Reuse saved product | 88 | YV 43 BDH, YL 25 LDH, 7D 3 BDH + 11 LDH, XO 3 EDH + 3 LDH |
| Reuse, supplement blocked | 1 | 7D J43A BDH: 26.46 uncovered days whose headers say 50.0000076 Hz (0.15 ppm) |
| Compute | 15 | 7D BDH ×11 (40 Hz; stations calibrated on LDH); XO LA34/LT02 HDH (APG, 125 Hz); WD66/WD70 HDH (100 Hz) |
| Blocked | 19 | 7D BXH ×15 and LDI ×3: no full response locally; XO LA33 HDH: headers 125.0000076 Hz (0.061 ppm) |

The 15 new PPSDs cover 1,190 channel-days, about 9.1 × 10⁹ samples. The four
XO HDH channels account for about 7.7 × 10⁹ of these (125/100 Hz, 163–225
days each), similar in size to the completed XO EDH products. Of the existing
vertical and horizontal products, 88 and 86 are reused.

**Nearby candidates (planned targets).** These counts are descriptive. They
count candidates meeting the estimator's distance and depth limits only;
common days and spectra still decide the result.

| Network | Targets with ≥ 2 candidates |
|---|---:|
| YL | 13 of 25 |
| 7D | 24 of 26 |
| XO | 1 of 10 |
| YV | 5 of 43 (saved results) |

Expect many `insufficient_reference_stations` / `no_suitable_neighbors`
outcomes. These describe the deployment geometry, not instrument failures. Not planned in pressure scope:
13 vertical and 66 horizontal epochs with usable responses, including all
50 YL horizontals and 16 XO horizontals. 7D horizontal/other channels are
blocked: there is no local response metadata.

**Checks per pressure epoch:**

| Check | Reuse saved | Planned | Blocked | Other |
|---|---:|---:|---:|---|
| Expected pressure | 43 (YV) | 61 | 19 | — |
| Nearby | 43 (YV) | 61 | 19 | — |
| Infragravity | 43 (YV) | 61 | 19 | — |
| Rayleigh | 89 | 0 | — | 34 deferred: all alternate streams at stations with evidence on another stream |
| Microseism | 89 | 0 | — | 34 deferred: same |
| Local-P | 66 | 0 | — | 39 not applicable (1 Hz streams, saved record); 18 deferred (BXH/LDI) |

The 66 local-P rows are 63 legacy-unvalidated, 2 unresolved and
1 conditional (LT02). This reproduces the six-check review: 89 PPSDs,
46/7/1/63 calibration evidence states, and 54 unjoined local-P rows.

**Heat maps of the saved products** (pressure only, ObsPy PQLX, occurrence
per 1-dB bin, no mean/median/percentile):

| Group (declared sensor family) | Stations | Saved windows | Zero-power floor values | Above saved histogram |
|---|---:|---:|---:|---:|
| YV BDH, hydrophone | 34 | 140,602 | 562,742 | 3,758 |
| YV BDH, DPG | 9 | 40,973 | 28,556 | 88 |
| YL LDH, DPG | 25 | 127,447 | 0 | 3,660 |
| 7D BDH, DPG | 4 | 3,854 | 20,296 | 0 |
| 7D LDH, DPG | 11 | 6,580 | 0 | 17 |
| XO EDH, hydrophone | 3 | 10,692 | 0 | 474 |
| XO LDH, DPG | 3 | 5,727 | 0 | 0 |

- **Grouping.** Figures are drawn per network, pressure channel code and
  declared sensor family. Two streams of one sensor are never pooled, and
  one channel code covering different instruments is split.
  - The family comes from a keyword rule on the StationXML sensor label:
    APG/Parosci → APG; differential pressure/DPG/Cox-Webb → DPG;
    hydrophone/HTI → hydrophone.
  - YV BDH is 34 hydrophones and 9 DPGs. After the run, XO HDH will be
    2 APGs (LA34, LT02) and 2 DPGs (WD66, WD70).
  - The labels are declared metadata, not hardware checks. The exact labels
    of each group are in `heatmap_manifest.json`.
- **Weighting.** Window-pooled and station-balanced figures are separate
  files, and each title states its weighting and denominator.
  - Window-pooled: each window counts once; the denominator is the number
    of windows covering that frequency.
  - Station-balanced: each station is normalized per frequency, then
    weighted equally; the denominator is the number of stations covering
    that frequency.
- **Extreme values.** Values are kept. ObsPy's numerical zero-power floor
  (−3076.53 dB) is counted separately from measured low levels, and values
  outside the saved −250…100 dB histogram range are counted. The detail view
  is −100…120 dB. The full-range figure and the NPZ counts keep every level.
- **No band censoring.** The usable sensor band is not censored: no band is
  invented.

After the owner's run, `report` adds the 15 new and merged products.

## Limitations and review status (not changed by this work)

| Topic | Status | Effect here |
|---|---|---|
| YV DEPAS positive-real RC pole (SOC-001) | Open. Pole reflection exists only as an optional trial; not adopted | PPSDs use \|R\|², so pole sign leaves them unchanged; it matters for calibration phase QC |
| INSU hourly levelling glitch (SOC-002) | Not removed | Affects INSU vertical/horizontal PPSDs and Rayleigh windows. No glitch subtraction applied |
| Usable bandwidth / self-noise (SOC-007) | No censoring adopted | Heat maps show the full saved band (0.001 Hz–Nyquist); low-frequency levels may be instrument noise |
| Zero-power floors (SENT-007) | Open | Counted separately, never dropped |
| Timing: 1 µs inter-channel guard (SENT-005), LT02 drift, LA33/LA34 | Guards unchanged; LT02 drift model only for local-P | Spectral levels are insensitive; calibration checks are not rerun |
| C07W intermittent vertical-amplitude state (SENT-001) | Supported; cause unconfirmed | C07W vertical PPSD mixes both states; no state selection |
| B01W/RR55 frequency-dependent response (SENT-004, SOC-006) | B01W supported; RR55 open | Scalar-gain evidence is reused as historical, not reinterpreted |
| Calibration identifiability | Microseism P/Z does not identify an independent gain; Rayleigh and microseism share the nominal vertical response and models; local-P deep-water results are diagnostic | Agreement between methods is not independent validation |
| Network default pressure sign −1 (SENT-012) | Open (YL/7D/XO runs used +1 explicitly) | One reason new calibration runs are deferred, not automated |
| Horizontal correction, relaxed QC, weighted fit | Not adopted | Not used |
| 7D responses | `INVENTORIES/7D.xml` lacks response stages | 7D horizontals, BXH, LDI and state-of-health channels are blocked until full responses are obtained |
| Rate encodings (XO LA33 HDH 0.061 ppm; 7D J43A 0.15 ppm) | The 0.1 ppm tolerance was authorized only for calibration | The existing PPSD code refuses these files; widening the tolerance is the owner's decision |
| Saved YV nearby results | Computed with the earlier pool | Reused. The larger pool cannot change them: networks are ≥ 1,707 km apart (YV ≥ 11,060 km from any other network), and no pressure product has a cross-network candidate within the estimator's 100 km / 10 % depth limits |
| Header probe | Spans come from the first/last record | Internal gaps appear as skipped PSD windows, not in planned intervals (`--full-headers` gives exact intervals) |
