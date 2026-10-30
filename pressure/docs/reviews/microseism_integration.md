# Reviewed integration of the microseism methods

15 October 2026. Source main: `ed5cb16`. Claude handoff:
`origin/claude/microseism-all` at `b603948`. The owner authorized this integration
and preservation of Claude's branch. The merge retains its commit history.

## Decision and scope

Integrate the two implementations as an experimental research workflow:

* Claude's homogeneous-water theory, modal calculations, H1/H2 estimates,
  phase hypotheses, paired-spectrum runners and historical comparisons;
* Codex's model-free coherence/phase/daily diagnostics on the same spectra;
* Claude's frozen frequency-consistency proposal, implemented by Codex with
  input-hash, station/model/epoch and non-overlap guards.

Rayleigh calibration defaults, response XMLs, waveforms, masters, prior numerical
results and scientific scores are unchanged. No gain is applied or accepted.
The original 35 best-band candidates remain available as historical diagnostics;
they are superseded for calibration interpretation.

## Why the additional rule is necessary

Claude's repeated YL weeks identified a selection artifact: a best-conditioned
sub-band can sit where the historical RR52 mode-0 model is unsuitable. C01W's
result changes by about a factor of two when different sub-bands are selected,
although the model-free diagnostic amplitude is approximately stable. This
supports the model/band-selection explanation, not a hardware certification.

C07W differs across time and cannot be resolved from the one favorable week.
The earlier claim that its microseism result rejects Rayleigh g=0.269 is
retracted. Existing evidence for relative vertical-amplitude states remains
relevant; matching dates and spectra is needed before assigning a pressure gain.

## Frozen rule and interpretation

The exact numeric rule is recorded as `claude-frequency-v1-20261015` in
`src/bruitfm/microseism_validation.py`. It comes from the handoff at `b603948`.
YL is the development set. Thresholds are not adjusted to improve agreement with
Rayleigh, the publications, or the later XO results.

1. A sub-band must end at or below 0.10 Hz, contain finite positive H1/H2/bound
   values, and have conditioning (H2 vertical-incidence bound / model H1) <=1.5.
2. A week needs an exact explained phase label and at least two eligible bands.
   Their midpoint spread, max/min minus one, must be <=10%.
3. The week value is their median midpoint; the descriptive envelope is min H1
   to max H2. No reference gain enters the rule.
4. At least three independent passing weeks are required. Their midpoint range
   divided by median must be <=10%; <=5% is labeled stable.

These are conditional repeatability criteria. They do not establish absolute
calibration. The v1 proposal retains the original whole-band phase screen and
adds no explicit coherence threshold. Coherence is recorded per eligible band;
phase in the actual gain bands and coherent tilt remain separate open checks.
The vertical response, water parameters and historical RR52 dispersion model
are still assumed. H1/H2 bounds rely on additive uncorrelated channel noise;
pressure-correlated tilt or mode mixing can bias both estimates.

## Review of the saved data

The counts in this section are the initial integration snapshot at `cc323b9`.
The subsequent [completed independent-week evaluation](../microseism_independent_weeks.md)
updates the current curated tables to 24 runs and seven conditional repeatable
estimates. No numeric threshold was changed.

The curated [frequency/week table](../microseism_frequency_table.md) and
[JSON provenance](../microseism_weeks.json) use 16 completed run folders:
four base runs, four YL shifts, four out-of-coverage 7D shifts, and four XO
shifts. Four repeated YL jobs are deduplicated. Each selected run and source
artifact is explicitly identified and hash-checked; the recorded review checks
1,269 distinct files. No waveform processing or calibration batch is rerun.

| Frozen-rule status | Stations |
|---|---:|
| At least three passing weeks, stable | 2 |
| At least one passing week, fewer than three | 25 |
| No passing weeks | 62 |

The two repeatable development-set results are:

| Station | Passing / checked weeks | Conditional g | Weekly spread | Rayleigh baseline g |
|---|---:|---:|---:|---:|
| YL.C05W | 3 / 5 | 0.57455 | 4.35% | 0.58466 |
| YL.C09W | 3 / 5 | 0.58433 | 1.95% | 0.57530 |

Neither is independent validation because YL informed the rule. None of the
five-week XO stations currently has three passing weeks. YV has only base-week
coverage in this review; the original 7D shifts fall outside its short archive.
Absence of a passing week does not identify an instrument fault.

Both legacy and new calculations are saved. The old
[week table](../microseism_weeks_table.md) is labeled legacy. Fresh derived
folders live under `results/REVIEWS/MICROSEISM_INTEGRATION/`; earlier experiment
and review products are not overwritten. Numeric historical comparison CSV/JSON
files are retained unchanged.

## Software corrections

* The combiner requires completed comparisons, deduplicates shifts, records
  excluded/mismatched runs, refuses changed input hashes or differing
  water/model/channel epochs, and prevents overlapping intervals from counting
  as independent weeks.
  Repeated shifts must have identical scientific comparison values; differing
  repeats require review instead of silently choosing the newest.
* Cached model-free diagnostics are reused only after checking source binding,
  multiplier and every recorded artifact hash.
* New Deng comparisons require the entire sub-band to fall within a published
  band. WD66's original partial overlap is not an exact match. Frequency-dependent
  published values without a matching band are omitted from gain plots.
* Shift runners inherit the base duration and refuse another network/duration.
  The weekly wrapper reuses completed shifts instead of recomputing them.
* The terminal-safe three-method launcher uses the integrated main code and
  reuses saved Rayleigh results by default. A new Rayleigh replay requires
  `--rerun-rayleigh`; `--check` performs only preflight checks.
* Source branches and mathematical reference calculations are retained. Claude
  receives implementation/proposal credit; Codex receives integration, guard and
  independent review credit. A self-review is not described as independent.

Verification: **236 tests pass** (77 existing ObsPy unit/encoding warnings),
Ruff checks and formatting pass, and Bash launcher syntax passes. The complete
suite includes numerical water-ODE and independent dispersion checks plus
regressions for the frozen rule, the high-frequency selection artifact,
cross-week provenance, changed models, duplicate/overlapping windows, stale
diagnostics and preservation of prior review output directories.
The final preservation check verified **4,121 distinct protected files** across
the completed three-method run and the new saved-week review: zero missing or
changed hashes. This is the recorded-input scope, not an archive-wide byte audit.

## Next concrete task

**Completed follow-up:** the eight requested YV/7D shifts below have finished;
see [the outcome](../microseism_independent_weeks.md). Do not rerun them. The
next scientific work is a bounded phase/tilt/vertical-response review of the
passing subset and the remaining 7D model/band discrepancy.

Evaluate the unchanged rule on independent YV weeks and coverage-appropriate 7D
weeks. The owner runs the longer jobs with:

```bash
bash scripts/run_microseism_validation_next.sh
```

This uses 7D shifts -2, -1, +1, +2 weeks and YV shifts -12, -6, +6, +12 weeks,
reuses completed shifts, and writes a new combined review. Existing YL and XO
weeks are read as saved. Do not rerun their completed jobs or a calibration batch.

Then assess per-band phase and tilt on the passing subset and investigate
C07W's time states using matched existing event dates/responses. If the frozen
rule fails, report that result before proposing another rule; preserve v1 and
clearly separate development data from subsequent tests. Microseism acceptance
and noise-model/scoring work remain downstream of this evidence.
