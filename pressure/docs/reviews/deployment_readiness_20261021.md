# Pending-branch review and deployment-readiness checkpoint

Date: 21 October 2026. Starting main: `6380338`.
Review branch: `codex/deployment-readiness`.

## Repository/PR audit and decisions

The owner requested review, responses and proper closure of all pull requests,
and a deployment-ready main. Both `gh pr list --state all` and the GitHub REST
`pulls?state=all` endpoint returned **zero existing pull requests**. There were
therefore no old PR conversations to answer or requests to close. Three remote
branches remained outside main's ancestry; these were reviewed as the pending
work instead. [Integration PR #1](https://github.com/MohammadAmin-Aminian/bruit-fm/pull/1)
records this review and was merged after green checks at main `7ac2e63`.
[Exact post-merge CI](https://github.com/MohammadAmin-Aminian/bruit-fm/actions/runs/39587156168)
also passed. Review and closeout comments were posted; no open PR remains at
the project transfer. No contributor branch was deleted or force-pushed.

| Branch / tip | Disposition | Reason |
|---|---|---|
| `claude/formula-review` / `c767579` | Selectively integrate the attributed report and synthetic worksheet; implement corrected phase reporting and the native-boundary refusal with tests | Two reproducible software defects. Do not import untested fixed-output plotting scripts or silently adopt scientific-policy proposals. |
| `claude/local-p-check` / `438c4ca` | Preserve; retain main's already-reviewed v2 implementation, not a wholesale merge | Direct replacement would remove protocol identification, conservative floor rounding, exact epoch summaries, tied-polarity handling and stronger output guards. The useful placement change and matching LT02 outcomes were already integrated. |
| `claude/local-p-science-review` / `7ebd872` | Preserve its attributed report/addendum; retain the existing closeout disposition | The independent LT02 result is already documented, including the -4.7% membership-sensitive 1–2 Hz change. Its raw-window review helper is not needed for deployment and is not rerun. |

The older local-P decisions are fully documented in
[local_p_integration_closeout.md](../local_p_integration_closeout.md). Keep the
late outcome addendum with the initial science review; the original less-than-1%
prediction is not a validated conclusion. Extra noise-control reads and proposed
promotion thresholds remain research options, not deployed defaults.

## Formula-review dispositions

The attributed [Claude review](claude_formula_review.md) is historical evidence,
not an instruction to implement every suggestion. Its optional figures/CSVs
remain accessible at the preserved source commit. Main's pressure probability
heatmap defaults, PQLX palette and absence of median overlays are unchanged.

| Finding | Reviewed action |
|---|---|
| F1: linear median of wrapped phases | Correct five producer locations with a shared, versioned circular L1 median. A summary straddling +/-180 degrees no longer reports zero. Empty/directionless phases are null, not fabricated polarity. |
| F2: forced family interpretation | Keep legacy nearest-centre labels descriptive. Null/nonfinite phases yield `indeterminate`; no new coherence/distance thresholds are selected. A physically meaningful coherence-gated family policy remains for the owner. |
| F3: 0/180 interpretation | Interpret the expected family only under stated positive-pressure/upward-displacement response conventions and the appropriate real forward relation. Relative inversion does not identify the physical pressure versus vertical channel. Do not flip responses or rename frozen rule IDs. |
| F4: nearby depth sensitivity | Retain as an interpretation warning; no model-dependent depth normalization or new station thresholds are silently applied. |
| F5: lower Godin band | No change to the existing frequency support; a new lower bound is a scientific protocol decision. The newer pressure-reference review records additional primary-source access beyond the historical Claude access notes. |
| F5: negative fractional native-grid tear | Refuse the boundary instead of interpolating across it. Preserve valid aligned duplicate overlaps, declared near-rate tolerances and anti-alias guards. |

The shared phase statistic minimizes shortest-angle L1 distance. Prefix sums
avoid an all-pairs matrix; the circular-mean-centred median is used only when it
is a global minimizer. The numerical directionlessness check is machine epsilon,
not an adopted scientific-coherence threshold. Fresh outputs record the summary
method; diagnostic/model-free outputs also record the phase-helper code hash.
Legacy summaries remain labelled `legacy_linear_median` when reused.

Only phase **reporting** changes. Per-bin phase masks, coherence criteria,
gain-fitting objectives, frozen gain values, experiment selections and acceptance
rules are not changed. The native-boundary guard can refuse future affected
reads; no completed run is repeated merely to replace historical decisions.

## Regression and packaging verification

Four new regressions failed before the fixes: trial, held-out and model-free
phase summaries around +/-180 degrees, and a negative half-sample native-boundary
tear. They pass after the corrections. Additional tests compare circular L1
results with brute-force minima on 40 seeded samples, wrapping/permutation
invariance, null/antipodal handling and preserved per-bin QC/frozen gains.

- Full local suite: **447 passed, 117 warnings**, 106.67 seconds. Warnings are
  existing/ObsPy pressure-unit and encoding warnings; they are not suppressed.
- Ruff check and format check pass; `git diff --check` passes.
- The project wheel builds successfully. A disposable environment installs that
  wheel without altering the base environment, runs outside the checkout with
  isolated Python, verifies packaged reference JSON/CSV and successfully invokes
  **nine safe module `--help` entry points**.
- CI now repeats wheel build and isolated-install smoke after style/tests.
  The smoke reuses installed scientific dependencies; it does not claim a fully
  dependency-locked environment or field calibration validation.

See [deployment instructions](../deployment.md). A checkout with its source
scripts/configuration and a package-only installation have different contents;
the wheel does not contain local archives or private saved results. Python 3.13
on Linux/macOS is the supported target; POSIX file locks remain in use.

## Data and history preservation

No real calibration, raw PPSD, all-data campaign or saved-result regeneration was
run. Raw directory contents, inventories, response snapshots, masters, fixed
gains, saved scores and original scientific evidence are unchanged. The existing
untracked `gain_histogram.png` belongs to the owner and remains untouched.

The same **1,104 protected scientific/evidence files** used in the preceding
YV integration retained their hashes after this review and full test run.
The sorted-path/hash aggregate remains
`7536e90035e12ec2767967dd775a554c75c69ccdc32ec0d150365d21a8dc9646`.
This bounded check includes saved spectra/PPSDs/response/master evidence; it is
not an archive-wide raw-byte hash audit. Intentional source-code/doc changes
are reviewed in Git, not represented as unchanged historical code hashes.

## What “ready to deploy” does and does not mean

The software checkpoint is installable and tested, with supported CLI entry
points and a conservative scientific configuration. It is **not** a completed
all-station calibration/noise-model campaign or acceptance of conditional gains.
New-machine path/provenance reconciliation, official YV response confirmation,
independent calibration evidence, ratings and the long missing-data campaign
remain separate work. No service is published and no data campaign is launched
by merging this PR.
