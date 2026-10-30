# Local-P phase closeout: reviewed v2

**Final integration update:** the owner subsequently requested cancellation of
the long batch and final testing/documentation/push. After verified interruption,
the reviewed branch is integrated into main. See the
[final closeout](local_p_integration_closeout.md) for current execution status
and Claude-review decisions. The pinned-main statements below describe the
earlier checkpoint; scientific results and limitations are unchanged.

16 October 2026. Engineering implementation and bounded validation are complete
on `codex/local-p-finish`. **This is not acceptance of an absolute calibration.**
At the owner's explicit request, the separate long all-methods batch continues
undisturbed, main stays at `ace6ee3`, and this integration branch is pushed
separately. Claude's original branches and worktrees are preserved.

## Result

The missing noise-window safeguard is fixed. Rather than discarding an event
solely because the old noise endpoint was too late, v2 moves the equal-length
noise window earlier to meet both existing conditions: at least one second
before the signal and no later than predicted P minus its search uncertainty.
The signal window, timing decisions, physical model and numerical QC thresholds
are unchanged. The expanded noise/signal span is checked for other catalogue
arrivals before analysis. No reference gain was used to select the fix.

Only three HDH APG pairs were rerun, at clean code `4c70941`:

| Channel | Outcome | Accepted v1 → v2 | Band | v1 g | v2 g | v2 events |
|---|---|---:|---|---:|---:|---:|
| XO.LT02..HDH | conditional_gain_candidate | 9 → 7 | 0.5-1 Hz | 0.972350 | 0.967805 | 4 |
| XO.LT02..HDH | same run | 9 → 7 | 1-2 Hz | 1.021197 | 0.972951 | 7 |
| XO.LA33..HDH | unresolved | 0 → 0 | none | — | — | 0 |
| XO.LA34..HDH | unresolved | 0 → 0 | none | — | — | 0 |

All seven accepted LT02 events pass the recorded noise guard. Relative changes
from v1 are -0.47% and -4.72% by band. With the earlier noise span,
`us1000hhug` is excluded for an overlapping catalogue event and `ak0188v8964i`
has no valid band. These are existing acceptance rules, not new thresholds.
LT02's band MADs are 0.0223/0.0429; model-corner envelopes are 0.890-1.067 and
0.879-1.115. These envelopes are not statistical confidence intervals. Residual
phase means remain about 25.9/32.3 degrees. Per-event timing refinement and the
shared nominal vertical response prevent independent phase/absolute validation.

## Review findings closed in code

| Finding | Resolution |
|---|---|
| Missing documented noise guard | Earlier conservative noise endpoint; persisted endpoint/earliest-P/guard fields; regression tests; bounded real-data check |
| Pooled gains copied to other epochs | Per-pressure-epoch summaries; no numerical export to an epoch without evidence; ambiguous pooled attribution suppressed |
| Weak join labelled exact | Full master key resolution using verified source inventory, unique epoch metadata and explicit hash-bound aliases; unrelated inventory matches/duplicate keys refused |
| Fetched 7D response versus master source | [Explicit mapping](local_p_inventory_map.json), checked source/target hashes and channel metadata; all 15 BDH epochs resolve; no response XML changed |
| Stale resume after changed data | Fingerprint includes resolved waveform root, matching file paths/size/mtime-ns, relevant code and numpy/scipy/ObsPy versions; pre/post dependency comparison; added/missing files invalidate reuse |
| Parent failure aborts campaign | Recorded failed row, next pair continues; corrupt summaries are archived, not reused; previous pair logs archived; nonblocking campaign writer lock |
| Polarity ties reported positive | Ties return indeterminate; counts and classified fraction retained; dominant classified sign explicitly not a confidence statement |
| Legacy microseism gains reintroduced | Only the frozen epoch-bound staging data are used; no best-band/weekly fallback; channel/source/epoch identity retained |
| Unsafe helper outputs | Resolved LOCAL_P-only destinations; raw/exports roots refused; exclusive writes; existing sidecars refused; unique catalogue staging directory |
| Missing restricted import path | `update_local_p_results` plus read-only-default `apply_local_p_results.py`; full keys and unique targets required; scientific scores and all six check-status fields forbidden |

The original independent [review](reviews/local_p_codex_review.md) is preserved
as a dated review of `f455edd`, not rewritten to conceal the original findings.
Development changes were not blind preregistration. Original hardware/response
diagnoses remain hypotheses; a phase match is not an independently known sign.

### Late independent implementation cross-check

Before final handoff, Claude pushed `dfb1df6`/`438c4ca` and the separate
[scientific peer review at `7ebd872`](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/7ebd8726a759c5fb3071352b8dfb88f79522edd7/docs/reviews/local_p_science_review.md).
Its saved `lt02_noise_rule_v3` summary was read, not rerun: it independently
reproduces **exactly** the two v2 band medians above, n=4/n=7 and seven
accepted events, with unchanged inputs. This is implementation corroboration
on the same data/model, not an independent absolute calibration reference.

The peer review also recommends placement rather than blanket rejection and
retaining the deep-water relation as diagnostic. Its preliminary prediction of
less than 1% station-median change failed for 1-2 Hz when the full QC rules were
reapplied; the addendum acknowledges the same two event removals and -4.7%
change recorded here. Event-membership sensitivity belongs in interpretation.
Its limited response-aware leakage synthetic is useful evidence for that
configuration, not proof of negligible leakage for all instruments/depths.

This branch retains its tested hardening and conservative floor-to-grid noise
boundary (Claude rounds the boundary). Claude's optional extra earlier-noise
control reads remain separate, not a new acceptance rule or prerequisite for
this engineering closeout. No late wholesale merge overwrites the reviewed
export/resume/safety changes; both contributor branches remain available.

## Staged integration and coverage

The [120-row epoch-bound staging export](local_p_table_staging.csv) composes
the explicit user campaign with only the three v2 replacements. The
[comparison table](local_p_station_comparison_v2.md) and
[full comparison CSV](local_p_station_comparison_v2.csv) preserve bandwise gains,
event counts and cross-channel caveats. They do not overwrite the v1 campaign.

| Outcome | Pairs |
|---|---:|
| Conditional LT02 candidate | 1 |
| Deep-water diagnostic-only | 17 |
| Unresolved | 22 |
| No suitable selected catalogue events | 26 |
| Insufficient sample rate | 39 |
| Missing full pressure response | 15 |

These are **120 channel pairs across 89 stations**, not 120 calibrated stations.
Of 66 executed pairs, three have v2 results and **63 remain explicitly legacy
noise-guard-unverified**. The 17 deep-water numerical results were not rerun
or promoted. The published references and Rayleigh baseline remain fixed.
Microseism's seven repeatable conditional estimates and frozen rule are unchanged.
Do not interpret this engineering closeout as full-archive v2 execution.

Production masters were **not updated**. A private copy of the actual master
was patched using the explicit importer: 66 unique exact diagnostic updates,
54 unresolved/not-run joins skipped. All **16,884 rows and 140 original CSV
columns** matched the source after patching; only local_p_* columns were added.
Existing score and check-status cells are preserved, not set to zero or pass.
The updater uses the established backup/transaction machinery and rejects
non-local-P edits or duplicate/unmatched full keys. Not-run channel-only rows
are not guessed onto individual master epochs.

## Verification and preserved inputs

- Full combined suite: **332 passed, 78 existing warnings**, 66.60 s at the
  first closeout check; final verification is recorded in the curated JSON.
- Ruff check/format and shell syntax pass.
- Bounded APG campaign: three completed, no failures; inventory/catalogue
  guards and waveform/dependency guards true; seven of seven accepted events
  satisfy the noise guard. No new all-station waveform batch was launched.
- All three completed APG pairs resume without executing a child waveform job.
- Fresh full legacy export plus explicit three-pair supersession works; no
  legacy-protocol conditional candidate is labelled as a current candidate.
- Private-master-copy updater check passes, including existing field equality.
- **144 protected phase inputs** retain their hashes: master CSV/XLSX,
  inventories, catalogues, source reports/events/saved spectra and reference
  staging data. The private verification copy is not a production master.
- Original production master SHA-256: CSV
  `d929f43ad006decd638ac64b9bb9388f5c0b9d82e9402c82a33b5c9b3e8c45f8`;
  XLSX `7041a9d0142d62764a591d3af62e56cc3acd5ea6a7b8ee4701d3ea470c3335e5`.

Waveform resume uses file identity/stat manifests, not a full byte hash of the
large raw archive. This detects normal additions/removals/changes; it is not
proof against a deliberately modified file with restored size and timestamps.
Original raw data, responses, saved campaigns, scores and frozen selections
were not modified. Generated staging/report files are new artifacts.

## Reproduce without a full campaign

Use the integration worktree and the `bruitfm` environment. Do not run another
all-station batch to repeat this closeout. Choose a new output name on repeat.

```bash
cd /Users/mohammadamin/Projects/bruit-fm-localp-finish
export PYTHONPATH="$PWD/src:$PWD"
conda run --no-capture-output -n bruitfm python -m pytest -q
conda run --no-capture-output -n bruitfm ruff check .

# Read-only campaign planning; no output folder or waveform samples.
conda run --no-capture-output -n bruitfm python scripts/run_local_p_campaign.py \
  --campaign-dir results/LOCAL_P/campaigns/plan_only --dry-run

# Read-only import plan. Omit --apply for production masters.
conda run --no-capture-output -n bruitfm python scripts/apply_local_p_results.py \
  --export docs/local_p_table_staging.csv \
  --master-dir /Users/mohammadamin/Projects/bruit-fm/exports
```

Saved bounded run:
`results/LOCAL_P/campaigns/validation_v2_apg/`.
Saved composed review:
`results/LOCAL_P/reviews/phase_v2_final_20261016/phase_report.json`.
Private updater verification:
`results/LOCAL_P/reviews/master_copy_v2/`.
These are git-ignored; [the curated JSON](local_p_phase_report.json) retains
key results and provenance in Git.

Recomposition is explicit and read-only for its inputs:

```bash
conda run --no-capture-output -n bruitfm python scripts/review_local_p_phase.py \
  --base-campaign /Users/mohammadamin/Projects/bruit-fm-localp-review/results/LOCAL_P/campaigns/user_all_20261016T120124Z_17172 \
  --supersede results/LOCAL_P/campaigns/validation_v2_apg \
  --output-dir results/LOCAL_P/reviews/phase_v2_repeat
```

For a future explicitly requested full campaign,
`bash scripts/run_all_calibrations.sh --check` is preflight only. Its defaults
reuse Rayleigh and do not launch additional shifted weeks; `--rerun-rayleigh`
and `--new-weeks` are explicit opt-ins. A new base microseism/local-P campaign
is still a real computation. The runner covers configured stations/selected
event-week windows, not every day or event in the archive; local-P retains the
80-event cap. Do not start it while the existing owner batch is running.

## Remaining science and next stage

1. Keep the pre-multiple relation diagnostic-only. Ideal-plane-wave synthetic
   tests do not establish real OBS coupling; filtering before gating can
   admit later-wave leakage. A separately designed field/response/gate study
   is required before treating deep-water numbers as calibration references.
2. Keep local-P gains bandwise, with no response fitting, cross-band scalar
   extrapolation or gain application to originals. Independent vertical scale,
   physical polarity and local propagation remain unresolved assumptions.
3. Do not bypass the deep-APG pressure-onset rule to rescue LA33/LA34.
   A response-appropriate detection/timing protocol needs independent evidence.
4. Leave scores and six-check classification undecided. The next project stage
   is a documented acceptance/scoring policy and accepted-station/noise-model
   workflow, not automatically declaring every computed estimate calibrated.
5. Completed in the final integration: the owner's batch was interrupted, live
   state checked, and this reviewed branch integrated into main. No branch
   deletion or force-push is needed. See the final closeout before any new run.

Claude's independent cited scientific review has arrived and was reviewed;
the addendum above records its limits and corroboration. Do not duplicate it.
Further research should seek independent reference/timing evidence, not relax
selection rules or treat implementation agreement as physical validation.
