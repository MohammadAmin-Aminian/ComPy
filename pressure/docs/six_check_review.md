# Six-check integration and approved-station noise workflow

## Scope and scientific boundary

The six calculation methods exist. This integration adds one read-only evidence
review, a versioned categorical scoring/classification engine, an explicit
approval record and a station-balanced pressure noise-model builder. It also
retains the new PPSD plotting tools. It does not invent numerical pass thresholds
or turn agreement among methods sharing a nominal vertical response into an
independent calibration. The default policy is **draft**, not scientifically
adopted. No real station approval was created during development.

Current saved data cover 89 stations in YV, YL, 7D and XO. There are 118 full
channel/epoch/inventory-source identities in the combined review, not 118
independently calibrated instruments. Keeping source identities separate avoids
copying a pressure multiplier across channels, responses or epochs.

## One safe command

```bash
cd /Users/mohammadamin/Projects/bruit-fm
bash scripts/run_quality_review.sh
```

Reads the explicitly configured saved sources in
`docs/quality_review_sources.json`. No newest-directory guessing. It writes a
fresh folder below `results/REVIEWS/SIX_CHECKS/` containing:

- `station_summary.csv`: 89-station availability/state overview. It does not
  average gains or infer whole-station validity from a single channel.
- `station_channel_review.csv`: 118 full-key rows with six source states,
  source outcome labels, reviewer verdicts and blank/unassigned review scores.
- `review.json`: exact source hashes, original numerical/band evidence,
  limitations, old manual-score values, source/epoch identities and all joins.
- `accepted_channels.csv`: selected exact rows, currently header-only.
- `decisions_template.json`: blank `decisions` plus fingerprint-bound
  `pending_review` records. These are prompts, not fabricated approvals.
- `noise_plan/`: a selection audit, not an estimated noise model.

Original master CSV/XLSX files, raw data, responses, saved spectra, reference
gains and scientific score/status fields are never written. Existing manual
scores are retained as source evidence but are not silently adopted or erased.
Derived CSV files are review exports; open them in Excel as text for SEED IDs
and source paths. No generated XLSX or production-master import is performed.

## Source bindings and current results

The key is `(Seed_id, channel_start_time, channel_end_time, inventory_file)`.
An open-ended epoch end may be blank; start and source must be present.

| Evidence | Binding |
|---|---|
| Expected pressure, nearby, Godin | Exact master-row pointer; saved result SEED and available source-PPSD hash checked |
| Rayleigh | Frozen baseline summary, verified result hash and explicit pressure epoch key inside the saved result |
| Microseism | Frozen epoch-bound staging export, never an old best-band fallback |
| Local-P | Reviewed staging export with `master_epoch_match=exact`; unresolved/not-run joins remain outside exact evidence |
| Saved PPSD | Exact master-row path and content hash; downstream modelling rechecks SEED, pressure units, bins and time windows |

Duplicate full keys in an in-scope source master are refused. Conflicting PPSD
or diagnostic pointers across masters are flagged and block noise approval.
Out-of-scope inventory duplicates do not contaminate this four-network review.
Inventory files and master XLSX siblings are also hash-protected when present.
No fuzzy station-only, channel-only or matching-date join is used.

The real-data dry run reproduces:

- 46 provisional Rayleigh estimates, seven conditional microseism estimates.
- One current conditional local-P candidate (LT02), two unresolved v2 APGs and
  63 explicitly legacy-unvalidated executed pairs.
- 54 local-P not-run/unresolved joins recorded separately, not guessed onto
  pressure epochs and not labelled as failed instruments.
- Three-method numbers at two stations, two-method numbers at another 14.
  These are availability counts, not validated agreement.
- 43 expected-pressure comparisons, four exploratory nearby comparisons and
  19 usable saved-PPSD Godin comparisons. Raw Godin pilot fields remain in the
  original masters; they are not substituted as an automatic alternative to
  the saved-PPSD Godin test.

Status counts in the exact table count source epochs. Do not add them and call
the total a station count. Processing errors, missing events/references,
inapplicable frequency bands, unresolved physics and failed instruments are
different concepts. Unknown scientific outcome labels fail closed as unresolved.

## Draft classification/scoring policy

`docs/six_check_policy_v1.json` implements the categories recorded in the
project brief, without claiming their scientific adequacy has been agreed:

| Explicit reviewer verdicts | Proposed classification |
|---|---|
| None performed (only not tested/not applicable) | `not_assessed` |
| At least one performed; every performed check passes | `valid_under_reviewed_checks` |
| At least two performed; every performed check fails | `invalid_under_reviewed_checks` |
| Partial, unresolved, mixed pass/fail, or only one failure | `unresolved` |

Only explicit review verdicts enter this rule, never numerical execution
outcomes. Pass score is 1, fail 0, partial strictly between 0 and 1. Unperformed,
inapplicable and unresolved scores remain blank. Scores are not averaged into
a station rating. No `nonworking` or `rescued` diagnosis is inferred from low
counts or disagreement; those require independent diagnosis/correction work.

Even `valid_under_reviewed_checks` does **not** automatically select a channel
for a noise model. Approval is a separate frequency/time-scoped scientific
decision. Every approval needs an independent-reference record; two methods
using the same nominal vertical response are not such a reference.

Adoption requires an explicitly reviewed policy copy with `status=adopted`,
`approved_by`, timezone-aware `approved_utc` and a scientific `rationale`.
The fixed rule/version and safety fields are validated. Do not switch this flag
just to obtain curves. Quantitative acceptance criteria for the six methods,
minimum usable calibration coverage, handling of the shared vertical/model
uncertainty and independent reference standards remain scientific decisions.
This task does not pretend to settle them.

## Approval record schema

After scientific agreement, use a separate source-config copy pointing to the
adopted policy, build a fresh review, then copy **only genuinely reviewed**
records from `pending_review` to `decisions`. Keep all original files.

Each record requires:

```json
{
  "row_id": "exact identifier from the review",
  "evidence_fingerprint": "exact fingerprint from that review",
  "reviewer": "actual reviewer",
  "reviewed_utc": "timezone-aware ISO timestamp",
  "rationale": "scientific basis and limitations",
  "checks": {
    "expected_pressure": {"verdict": "unresolved", "score": null, "rationale": "reason"},
    "nearby": {"verdict": "unresolved", "score": null, "rationale": "reason"},
    "infragravity": {"verdict": "unresolved", "score": null, "rationale": "reason"},
    "rayleigh": {"verdict": "unresolved", "score": null, "rationale": "reason"},
    "microseism": {"verdict": "unresolved", "score": null, "rationale": "reason"},
    "local_p": {"verdict": "unresolved", "score": null, "rationale": "reason"}
  },
  "noise_model_approved": false,
  "nominal_response_accepted": false,
  "independent_reference": {"citation": "", "description": ""},
  "approved_frequency_band_hz": [],
  "approved_time_range_utc": []
}
```

The template intentionally cannot authorize anything as written. Scoring absent
evidence or a processing failure is refused. Stale fingerprints, duplicate or
unknown identifiers, missing check rationales, boolean/invalid scores and naive
timestamps are refused. All six decisions must be explicit. Human review is a
trust boundary, not a digitally authenticated measurement or proof of accuracy.

For noise approval, classification must be valid, nominal response explicitly
accepted, an exact PPSD available, no source conflict/missing-inventory issue,
and independent reference, finite frequency interval and increasing time
interval provided. The model builder requires full PSD windows and full
smoothing-bin support inside these intervals. A gain is never automatically
applied. An instrument requiring response correction is not rescued by this
nominal-response-only builder.

## Conditional downstream model builder

Safe plan, valid even while no channel is approved:

```bash
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  -m bruitfm.noise_model --review PATH_TO_REVIEW/review.json \
  --output-dir results/REVIEWS/SIX_CHECKS/new_noise_plan --plan
```

Only after an adopted policy and real approvals, omit `--plan` and choose a new
output directory. It then writes pressure low/high percentile CSVs, PNGs and
station-balanced probability NPZs. Defaults are descriptive 5th/95th percentiles
at 1-dB histogram resolution, configurable by `--percentiles`; these choices
are not newly adopted global standards. `--group-by` supports `network_channel`
(default), `network`, `channel`, or `all` for an explicitly scoped combined pool.

Epoch windows are pooled within each station; distributions are normalized
before equal weighting across available stations at each frequency. Long
records do not dominate the model. Different native frequency grids use the
documented piecewise-constant display-bin edge union, without inventing PSD
resolution. Coverage is saved per frequency. Overlapping approved records of
one station within a group are rejected to avoid double counting, including
collocated channels when using `all`. Disjoint epochs can be combined.
The original overlapping windows are correlated; their number is not a count
of independent observations, and these percentile curves are not confidence
intervals.

No approved channels means an explicit refusal to generate noise-model curves,
not empty curves labelled success. Review acceptance flags are reconstructed
from the original decisions and fingerprints; manually toggling a Boolean in
`review.json` is refused. Every protected input hash is rechecked downstream.
The model output remains labelled descriptive, not a validated global model.
Depth/coastal/seafloor/seasonal grouping, representativeness, uncertainty and
independent field validation require additional approved scientific protocols
and metadata. They are not completed by this engineering checkpoint.

## Verification and continuation

Synthetic tests exercise the approved path without approving any real station,
including weighting, coverage, full-bin/time support, unchanged inputs and
tamper/stale-source rejection. The real run is read-only and deliberately has
zero accepted channels. See the curated `six_check_phase_report.json` for final
test counts and protected-input results. Saved integration reviews remain local
under `results/REVIEWS/SIX_CHECKS/`; source spectra are not copied into Git.

Next: review this draft policy and the 89-station table scientifically, starting
with YV and the best-supported reference cases. Record actual approval evidence
before any accepted-station noise-model run. Do not rerun completed calibration
campaigns, relax thresholds to manufacture passes, or populate production scores
merely because the engineering workflow is ready.
