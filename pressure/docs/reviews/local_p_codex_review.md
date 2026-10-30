# Local-P implementation review

Reviewed 2026-10-16: Claude `f455edd`, based on `2b4fef8`.
Review branch: `codex/local-p-review`. Current main: `ace6ee3`.

**Recommendation: retain this implementation, but do not merge it unchanged or
promote its gains to accepted calibration.** Numerical reproduction passes;
the noise-window contract and export attribution need resolution first.
This review changes no implementation, scientific settings, gains or scores.
Claude's worktree and main remain untouched. No full campaign was rerun.

## Verification

Commands ran in `/Users/mohammadamin/Projects/bruit-fm-localp-review`, using
the `bruitfm` environment and `PYTHONPATH=src:.`.

| Check | Result |
| --- | --- |
| Catalogue copies versus JSON SHA-256 | PASS: XO `101995d9...`, 7D `b6237a9b...`, YV `a05a5376...` |
| Full branch test suite | PASS: 271 tests, 77 warnings, 62.13 s |
| Ruff check / format check | PASS: 138 files already formatted |
| Full dry-run plan | PASS: 66 eligible, 39 insufficient sample rate, 15 missing full pressure response; no output folder created; planning code does not read waveform samples |
| Timing plan | PASS: LT02 HDH `drift_model`; LA33/LA34 HDH `none` |
| Regenerated 120-row CSV | PASS: byte-identical to `docs/local_p_campaign_v1.csv` |
| Three-pair APG rerun | PASS: scientific summary fields match `apg_timing_v2`; all three `events.csv` files byte-identical |
| Fresh APG export | PASS: three rows, all `master_epoch_match=exact`; no score field assigned; protected-input hashes unchanged |
| Saved campaigns, read only | PASS: 66 completed pairs in `all_v1`, zero failed; all 69 pair summaries across both campaigns report unchanged inputs; current inventory/catalogue hashes still match their recorded hashes |
| Saved all-station export | PASS for this dataset: 120 rows; all 15 high-rate 7D BDH rows marked exact; see general attribution defect below |

The final combined table has 1 conditional candidate, 17 diagnostic-only,
22 unresolved, 26 no-suitable-events, 39 sample-rate exclusions and 15
missing-response exclusions. These are channel-pair outcomes, not scores.

| APG channel | Outcome | Accepted | Band | Multiplier | Events in band |
| --- | --- | --- | --- | --- | --- |
| XO.LT02..HDH | conditional_gain_candidate | 9 | 0.5-1 Hz | 0.9723496861 | 5 |
| XO.LT02..HDH | same run | 9 | 1-2 Hz | 1.0211971547 | 9 |
| XO.LA33..HDH | unresolved | 0 | none | none | 0 |
| XO.LA34..HDH | unresolved | 0 | none | none | 0 |

LT02's exploratory scalar is 0.9992287421, not an adopted correction.
Its band phase means are 33.31 and 38.79 degrees, not zero residual phase.
Expected rerun values match; differences are timestamps, paths, code provenance
and path-dependent campaign fingerprints. The handoff's n=7 / scalar=0.998
paragraph is stale relative to its final n=9 drift-model results. Its open
question saying the drift model is unimplemented is also stale.

Verification artifacts, git-ignored, are under `results/LOCAL_P/`:
`verify_table.{csv,md}`, `campaigns/verify_apg/`, and
`campaigns/verify_apg_export/`. `campaigns/verify/` does not exist.

## Findings, in priority order

Severity: P1 = resolve before authoritative results/master integration;
P2 = correctness or safety issue to fix before routine production reuse.
Line numbers refer to the reviewed Claude commit, not current main.

1. **P1: documented noise-window safeguard is absent.**
   `src/bruitfm/local_p_check.py:964` computes the noise window and checks
   array bounds, but never requires its end to precede predicted P minus
   search uncertainty, as promised by `docs/local_p_method.md:158`.
   Recomputing that boundary with `search_half_width` from saved event records
   finds **8/9 accepted LT02 events** outside the stated condition, by
   0.917-4.848 s. Example: `us1000g6g6`, 4.848 s.
   This proves a code/documentation mismatch, **not observed P contamination**:
   the noise still precedes the measured onset. Decide explicitly whether the
   conservative theoretical boundary or a validated measured-onset rule is
   intended, add its regression tests, then perform only a bounded LT02 check.
   Do not silently relax the rule or automatically discard these events during
   integration. The present candidate has not passed the documented protocol.

2. **P1: export labels a weaker join as an exact full-key match and can spread
   gains across untested epochs.**
   `scripts/export_local_p_results.py:136` matches SEED ID and coverage epoch
   dates, then copies the same aggregate numbers to every match at line 181.
   `inventory_file` is copied from the master rather than independently
   resolved as part of the result's identity. Accepted-event epochs are not
   used to allocate gains. An in-memory mocked export with two coverage
   epochs, acceptance only in the first, and a second inventory for that first
   epoch produced three `exact` rows, all with the same gain. Thus an untested
   epoch and an unrelated inventory both received numbers. Existing campaigns
   have one coverage epoch per pair, so this is **not evidence that their
   present epoch values changed**. Require explicit full-key resolution,
   unique/verified inventory mapping and epoch-specific evidence; suppress
   numerical promotion when attribution is ambiguous. Fetched response paths
   may legitimately differ from master inventory paths: verify that mapping,
   do not merely demand path equality.

3. **P2: resume fingerprint omits waveform identity/state.**
   `scripts/run_local_p_campaign.py:233` hashes code, options, inventory and
   catalogue, but not waveform root, files or their state. Changing only
   `waveform_dir` to a completely different archive produced the same
   fingerprint (`21cb6d974baf8d53` in the probe). Line 258 then permits reuse
   without checking the previously read waveforms. Persist/check an input
   manifest and archive root, and detect added coverage as well as modified or
   missing files. Record relevant dependency versions for reproducibility.

4. **P2: polarity summary resolves a tie as +1 and hides indeterminate votes.**
   `src/bruitfm/local_p_check.py:1122` drops all zero classifications and uses
   `max` without a tie check. Synthetic votes `[+1,-1]` return +1;
   `[+1,0,0]` also return +1. Actual RR52 has one negative and one
   indeterminate band vote; RR40 has 9 negative and 12 indeterminate votes,
   yet each receives -1. Preserve counts, return indeterminate for ties,
   distinguish a dominant *classified* sign from a reliable overall polarity,
   and agree an event-level confidence policy before scoring. Do not invent a
   new confidence threshold within this review.

5. **P2: comparison generator reintroduces superseded microseism numbers.**
   `scripts/summarize_local_p_campaign.py:33` falls back to legacy weekly
   medians when no frozen-rule candidate exists. CSV basis labels help, but
   the Markdown still presents these in the same gain column. Main explicitly
   supersedes the original 35 best-band candidates. Use main's current frozen
   staging/provenance; leave non-reportable gains blank in an authoritative
   comparison, or move legacy values to an explicitly historical column.
   Retain channel/epoch/band identities, not just station names. Local-P Hz
   bands cannot validate a Rayleigh 0.03-0.07 Hz scalar by percentage agreement.

6. **P2: helper output protection is weaker than the station runner's.**
   `scripts/fetch_local_catalog.py:56` checks only CSV existence; lines 71/85
   can overwrite an existing sibling `.tmp` or `.json`. Its output is not
   confined to LOCAL_P. The summarizer (`:115`, `:150`) also overwrites explicit
   output paths, and the drift fitter accepts arbitrary fresh output paths.
   These helpers were **not used to overwrite any existing artifact** here.
   Before treating the whole workflow as protected, validate every destination,
   refuse existing sidecars, and use unique atomic temporary files. The
   campaign runner intentionally updates its own reports/logs on resume;
   that is not append-only archival behavior.

Additional operational gap: child-process failures normally become failed
rows and allow continuation, but exceptions in fingerprinting, prior-summary
JSON parsing or directory moves can escape `run_pair` and stop the campaign
(`scripts/run_local_p_campaign.py:411`). Add a per-pair parent exception
boundary and a regression test; preserve failed attempts/logs.

## Physics and identifiability

The PDF review included visual checking of Deng's Eq. 5 on page 3 and reading
Zha & Webb's Appendix A, pp. 3234-3236. The implementation's
`p/a_up = rho*h*tan(theta)/theta`, `theta=2*pi*f*h*eta/c`, agrees with those
relations. Compression-positive pressure, upward Z, NumPy derivative `+i*w`,
the pressure-multiplier convention and separate polarity handling are
consistent. Before a surface reflection, the plane-wave impedance is
`p/v_up=rho*c/eta`; the delay `2*h*eta/c` is a fixed-horizontal-slowness phase
delay, not a geometric ray round-trip length divided by c. Deng explicitly
limits the simpler velocity relation to the early crest; it is not blanket
validation of local-P coda calibration. See [Deng et al. (2022), sections 2-4](https://doi.org/10.1029/2021GL096952).

The low-frequency limit is also consistent with [An et al. (2017), section 2](https://doi.org/10.1002/2017GL075137).
Gain conventions must be inverted when comparing with the recorded/predicted
factor in [Zha & Webb (2016), Appendix A](https://doi.org/10.1002/2015JB012688).
Deng's supporting information was not independently checked; this review
does not claim new per-station published benchmarks from it.

The independent finite-difference tests are useful: both relations, two
half-space properties, nonzero slowness, and synthetic gain/sign recovery
passed. They test an ideal layered wavefield, not real OBS coupling or an
independent absolute sensor reference. The eta bound assumes vP >= 5 km/s;
that is conditional, not universal for slow sediments or scattered arrivals.

The identical-window forward model avoids directly dividing a smoothed
observed spectrum by a rapidly varying pointwise theoretical curve. However,
it filters the **whole segment before gating** (`local_p_check.py:548`, `:568`).
Zero-phase filtering and response conversion have temporal support outside
the selected window. A geometrically pre-multiple end time alone does not
prove absence of later-wave leakage, especially with non-flat real responses.
Require controlled post-gate multiple/S-wave injections, realistic responses
and gate/padding sensitivities before promoting this regime. No claim that a
specific saved diagnostic is thereby biased has been established here.

Four tapers provide actual spectral averaging; zero-padding does not create
independent observations. The minimum three retained FFT bins are not three
independent estimates. H1/H2 noise brackets and model-corner envelopes are
diagnostics, not calibrated statistical confidence intervals. Station
frequency assessment compares band medians against event MAD, with potentially
different events in each band: it does not prove a frequency-independent
instrument response.

Timing is explicitly conditional: LT02's drift comes from the earlier
campaign's timing measurements, and the rerun still refines lag per event.
This is neither held-out timing validation nor independent residual-phase
validation. Reproducing results does not make the remaining phase offset zero.
The shared nominal vertical response also prevents Rayleigh/microseism/local-P
agreement from independently proving absolute pressure gain.

## Decisions and integration path

1. **Pre-multiple reference:** keep diagnostic-only. An ideal wave equation
   and successful synthetics are necessary but insufficient for field-gain
   acceptance. First validate gate/filter leakage and synchronized absolute
   gauges under a separately agreed, bounded protocol.
2. **Frequency dependence:** retain per-band results and uncertainties;
   do not extrapolate high-frequency multipliers to Rayleigh/microseism bands,
   fit a new response or apply a scalar to original data.
3. **Scoring:** leave the local-P calibration score blank. Separate execution,
   applicability and scientific outcome. Do not translate unresolved into
   failure/zero, or conditional into pass/one. Agree the computed-but-unscored
   status mapping with the six-check schema before updating anything.
4. **Deep APG onset rule:** do not remove it ad hoc to rescue LA33/LA34.
   A separate diagnostic protocol could use independently supported timing
   plus a response-appropriate pressure detection band, retaining S exclusion,
   conservative gating and explicit lost-independence flags. Phase fitted to
   the same waveform is not an independent timing check. Keep current results.

After approval: resolve findings 1-6 and the failure-isolation gap with focused
tests, preserving all old outputs; reconcile current-main documentation and
microseism provenance; test the combined branch. The 271 branch tests versus
276 on current main reflect divergent additions, not evidence of a regression
or a verified combined test count.

Only then implement a restricted `local_p_*` master updater, modelled on the
existing restricted updater, using all four master keys and refusing
duplicates, ambiguous mappings, cross-epoch copying and non-local-P writes.
Carry execution/scientific status, regimes, bandwise evidence and provenance,
plus timing/diagnostic flags. Keep
`calibration_check_pressure_vertical_local_p_waves_score` blank and preserve
all existing other-method values/statuses/scores. Protect/have before-after
hashes for masters, responses and saved inputs. This review grants no approval
to merge, update the master or launch an all-data rerun.
