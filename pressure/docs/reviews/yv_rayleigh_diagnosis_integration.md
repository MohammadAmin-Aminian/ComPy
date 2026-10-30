# Review and integration of Claude's YV Rayleigh diagnosis

Date: 21 October 2026. Source branch: `claude/yv-rayleigh-diagnosis`,
commit `dbd4856`, based on main `2fb1fc0`. Review branch:
`codex/yv-rayleigh-review`. The contributor branch/worktree is preserved.

## Decision and scope

Integrate the **diagnostic report, saved evidence and read-only audit helper**
after the corrections below. This accepts no physical response correction,
pressure gain, polarity assignment, calibration score or production-default
change. No new gain was fitted and no calibration/held-out batch was run.

The separate `claude/formula-review` at `c767579` was inspected for context but
is **not merged here**. Its wrapped-phase, phase-family, nearby-depth and timing
findings require their own review. Its instructions to change production code
or adopt new thresholds are not automatically owner authorization.

## Findings retained, with scientific qualifications

- The frozen nominal workflow has **16** YV estimates, 25 ordinary unresolved
  stations and two band exclusions; the completed pole-only trial has **31**,
  ten ordinary unresolved and two exclusions. These are numerical availability
  counts, not accepted absolute calibrations.
- Unpaired complex vertical poles are confirmed in the original local YV XML
  at RR22/RR23/RR24/RR33/RR39/RR42/RR53/RR56. Four are in the 43-station local
  pressure review. [Stähler et al. (2016)](https://doi.org/10.5194/adgeo-41-43-2016),
  Appendix A/Table A1, prints conjugate pairs; that table was visually checked.
  Restoring pairs is a supported **hypothesis**, not confirmation of the actual
  station-specific hardware, normalization or absolute sensitivity.
- RR17/RR18's large pressure poles and H4 phase reduction reproduce. Dividing
  the pole by ten is post hoc development, not a manufacturer calibration.
  RR47's H3 result supports opposite relative P/Z polarity under the assumed
  forward relation; it does not identify which physical channel is inverted.
- Table 3 and section 3.1 of the same primary paper support the reported failed
  verticals, RR45 noise and RR31 artefact history. They do not prove that all
  archived RR31 records are irrecoverable. RR12's stored pressure PPSD has
  **4,769/4,778** median window levels below -3000 dB, confirmed directly from
  its saved NPZ; this is a source/archive question, not an automatic failure
  diagnosis for the original hydrophone.
- Immediate software screening blockers must not be confused with complete
  physical causes. Structure/mode mixtures, response errors and window placement
  remain possible contributors. RR57's pointwise resonance ratio 0.898 and its
  conservative maximum-phase-speed preflight exclusion are different calculations.

The primary PDF was already local under
`results/LITERATURE/PRESSURE_REFERENCES_20261019/Staehler_2016_RHUM_RUM.pdf`.
Relevant Table 3/section 3.1 and Appendix A pages were checked; the publisher
article page was also verified. No new station metadata or waveform was fetched.

## Helper corrections

Four newly added regressions failed on Claude's unmodified helper: a foreign
network snapshot was accepted, a trial could point at an unrelated baseline,
source changes during the run were not detected, and a positive complex pressure
pole was silently converted to a real one. These are audit-helper defects,
not changes to the production Rayleigh calibration.

The integrated helper now:

- verifies workflow/comparison/response/spectral hashes against the completed
  trial's recorded source manifest, with the explicit expected trial-manifest
  hash available on the CLI;
- binds baseline pressure full epoch/inventory identity, trial baseline pointer,
  network/station, exact snapshot channel epochs and vertical partner;
- verifies the fixed no-tilt sign -1, band and phase/coherence settings before
  using those diagnostic conventions;
- hashes every individual saved event NPZ, plus the optional microseism summary,
  check JSON, paired-source JSON and bins, binding both full P/Z epochs;
- rechecks all consumed input hashes before writing a fresh exclusive output;
- refuses unsupported complex-positive pole reflection, unknown event outcomes
  and invalid/nonfinite/mismatched saved spectral arrays;
- includes full pressure/vertical IDs, epoch, inventory source and trial state
  in fresh CSV exports. Existing saved CSV/JSON files are not regenerated.

It still reads raw **directory names only** when `--data-root` is supplied;
file-name day counts are not validated waveform coverage. It does not correct
metadata, choose new events, fit gains or start raw-data processing.

## Verification actually performed

Targeted suite: **59 passed, 24 existing/ObsPy warnings**, 14.09 seconds:

```bash
PYTHONPATH=src:. /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python \
  -m pytest -q tests/test_yv_rayleigh_diagnosis.py \
  tests/test_rayleigh_calibration.py tests/test_calibration_experiment.py \
  tests/test_calibration_experiment_runner.py tests/test_rayleigh_validation.py
```

There are 16 diagnosis tests (seven original plus nine added), including changed
source, network/epoch mismatch and optional microseism epoch-binding refusals.
Ruff check/format pass on the two changed Python files. No claim is made that
the entire historical 423-test suite was rerun during this integration.

One saved-input-only replay completed for **43 stations and 583 saved event
spectra**, tracking **875 individual input hashes**. Every metric already in
Claude's 43 station records matches exactly. Only identity/provenance fields
were added. No raw waveform was opened and no calibration gain was refitted.

New local output, never overwriting Claude's evidence:

`results/REVIEWS/YV_RAYLEIGH_DIAGNOSIS/codex_review_20261021/`

Its `evidence.json` SHA-256:
`953e40d5facf2e42c52dc92a38fae60699a1b2c12d126bcfa1ce12f2322e986c`.
Use the explicit CLI in the source report and a **different new output folder**
if another saved-input audit is justified; no completed science batch is needed.

Original Claude evidence remains unchanged:

| Artifact | SHA-256 |
|---|---|
| `docs/yv_rayleigh_diagnosis/evidence.json` | `ebd4f95f1a6f43a4e1fe1a70528c15153f6db02be2744dfbfc2e972bcbb6f854` |
| `docs/yv_rayleigh_diagnosis/stations.csv` | `5c740eb661953182192e33c72ee38e10cfacc58014faeac6c1bbdd252010d783` |

Before/after the audit, **1,104 unique protected files** retained their hashes:
the union of Claude's 210 recorded inputs, the baseline's event NPZs, the
historical six-check protected-input list and the owner's untracked
`gain_histogram.png`. This includes the saved PPSDs, responses, inventories,
masters and evidence in that bounded scope. Its sorted-path/hash aggregate is
`7536e90035e12ec2767967dd775a554c75c69ccdc32ec0d150365d21a8dc9646`.
This is not an archive-wide raw-byte audit; no raw waveform byte hashes are claimed.

## Next concrete task and approval boundary

First confirm official response metadata/calibration records for RR22/RR56 and
the unusual RR17/RR18 pressure poles, keeping downloaded metadata separate from
the original inventory. Official archive agreement alone would not certify the
hardware response. RR47 still needs independent polarity evidence and RR53
needs hardware/response identification.

Then present a bounded saved-spectrum response-hypothesis pilot for owner
approval, with frozen events, unchanged QC and RR52 control. Any experimental
gains remain development results; no XML/master update follows automatically.
Keep untouched validation events for later testing. Do not rerun the completed
RR01/RR06/RR48/RR52 or RR05 held-out studies, launch network campaigns, adopt
new thresholds, apply gains or mark stations calibrated from this review.
