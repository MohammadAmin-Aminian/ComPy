# Handoff: local-P check (Claude, branch `claude/local-p-check`)

**Archived contributor handoff.** Codex subsequently reviewed and hardened
this implementation. Read [the current phase report](local_p_phase_report.md)
first for v2 status and numbers. The n=7/scalar=0.998 pilot paragraph and the
open question saying the drift model is unimplemented below are historical
inconsistencies, not current instructions. Hardware diagnoses remain hypotheses.

Base: `main` at `2b4fef8`. Worktree used: `~/Projects/bruit-fm-localp`.
Not merged into `main`. Codex's checkout and `claude/microseism-all` were not
touched.

## What exists

| File | Role |
|---|---|
| `docs/local_p_method.md` | Method note: sources/pages, conventions, two regimes, applicability to local P, estimator, pre-registered settings, post-pilot revision, open questions |
| `docs/local_p_check.md` | Data inventory, catalogue provenance, pilot results, exact commands, output fields |
| `src/bruitfm/local_p_physics.py` | Pure relations (Deng Eq. 5, pre-multiple ρc₀/η, η bounds, poles, round trip, regime classifier), multitaper cross-spectra, AIC |
| `src/bruitfm/local_p_reference.py` | Independent 1-D finite-difference water/half-space simulation (validation only) |
| `src/bruitfm/local_p_check.py` | One-pair pipeline and CLI (`python -m bruitfm.local_p_check`) |
| `scripts/fetch_local_catalog.py` | USGS ComCat CSV into new files with provenance sidecar |
| `scripts/run_local_p_campaign.py` | All-station plan, `--dry-run`, resumable run, per-pair logs, summary |
| `scripts/export_local_p_results.py` | Read-only export keyed to master epochs, and cross-check table |
| `scripts/summarize_local_p_campaign.py` | Station table of a campaign beside Rayleigh/microseism values |
| `docs/local_p_campaign_v1.md`, `.csv` | Full-campaign report and 120-row table |
| `scripts/fit_local_p_clock_drift.py` | Robust linear clock-drift fit from refined lags (timing only) |
| `docs/local_p_timing.json` | Per-channel timing decisions with evidence (LT02 drift model; LA33/LA34 none) |
| `tests/test_local_p_{physics,check,campaign}.py` | 35 tests |

Reused without modification: `archive.WaveformFile`,
`rayleigh_calibration.load_aligned/load_resampled/catalog_events`,
`ppsd.response_handling/atomic_json`. Production Rayleigh/microseism code is
unchanged. README step 57 and the Calibration 3 row of `project_context.md`
were updated.

**Tests:** full suite 271 passed, Ruff check/format clean (the local-P files
add 35 tests). Coverage: unit/slowness conversions, sign and Fourier
conventions, incompressible limit, resonance guards, consistency with
`water_ratio`, an independent finite-difference simulation of both relations
(R-invariance included), synthetic gain/polarity recovery through the full
estimator in both regimes, window/regime rules, sample-quantum, lag sign,
onset/clipping, outcome rules, catalogue parsing, header-based file indexing,
output protection, launcher rows and export flags.

## What was run (bounded)

- Catalogue fetches (3 files, `results/LOCAL_P/catalogs/`, git-ignored).
- Development runs, `results/LOCAL_P/dev/` (LT02 EDH/HDH, WD66; ≤ 25 events).
- Pilot campaign `results/LOCAL_P/pilot/XO_LT02_WD66/` at commit `fadfb06`,
  plus export `results/LOCAL_P/pilot/XO_LT02_WD66_export/`. Results in
  `docs/local_p_check.md` §2. All input hashes unchanged.
- The user ran the all-station campaign `results/LOCAL_P/campaigns/all_v1/`
  (66 pairs completed, 0 failed, inputs unchanged) and its export.
- APG timing follow-up at commit `ff45286`: rerun of the three XO APG pairs
  only, `results/LOCAL_P/campaigns/apg_timing_v2/` (`docs/local_p_check.md` §2.1).

## Key findings to carry forward

1. Only XO.LT02 (152 m) is shallow enough for the steady-state regime with
   local P. Every other eligible station is 2.6–5 km deep, where only the
   time-gated pre-multiple relation applies (diagnostic). The local-P check
   will therefore produce at most one conditional gain candidate in this
   archive unless the pre-multiple relation is accepted as a calibration
   reference (scientific decision; see below).
2. LT02 APG: g = 0.972 (0.5–1 Hz, n = 5) and 1.021 (1–2 Hz, n = 7), consistent
   (scalar 0.998), model interval ±10–17%.
3. LT02 APG clock drift −0.483 s/day (−5.6 ppm); it predicts Deng's four
   teleseismic offsets within 0.3 s. Refit as a drift model (16 events,
   RMS 0.09 s), it is now the LT02 timing mode (9 events accepted, g 0.97/1.02).
   LA33/LA34 APGs do **not** drift (synchronized within ~0.02 s at LA33); they
   stay unresolved because their 2–10 Hz onset SNR is below 3 (APG noise floor
   above ~2 Hz). Timing decisions live in `docs/local_p_timing.json`
   (`docs/local_p_check.md` §2.1).
4. LT02 EDH hydrophone is dead for all of 2018.
5. WD66 DPG (pre-multiple, diagnostic): ≈0.78 at 1–4 Hz, rising above 4 Hz
   where the nominal response rolls off; polarity +1 vs nominal response.

## Master-table integration requirements

- Keys: `Seed_id`, `channel_start_time`, `channel_end_time`, `inventory_file`
  from the export (`master_epoch_match = exact`). Epoch strings match the
  master's `str(start or "")` form. All 15 7D BDH epochs of the fetched
  response file match master keys exactly. 7D BXH has no full response and
  stays `no_full_pressure_response`.
- New prefix `local_p_` for results (as `rayleigh_*`): status, execution
  status, vertical SEED ID, regimes, accepted events, per-band JSON, frequency
  assessment, scalar-if-consistent, polarity, depth, result/plot files,
  review flags, error. A restricted updater like
  `inventory.update_rayleigh_results` with prefix `local_p` is needed. The
  export already provides the rows; this branch does not write any master.
- `calibration_check_pressure_vertical_local_p_waves_score` must stay blank.
  Its `_status` should stay `not_tested` until a scoring policy is agreed.
  Proposed mapping when one is: `not_applicable_sample_rate` and
  `no_full_pressure_response` → `not_applicable`; `failed` → keep the error;
  everything else → a computed-but-unscored status. Outcome names are not
  scores.

## Open questions (scientific)

1. Accept or reject the strictly gated pre-multiple relation (Eq. B) as a
   calibration reference for deep stations. It is exact for an upgoing plane
   wave before the first multiple (finite-difference check: 0.2%), but it
   has not been validated on real data in the cited literature.
2. How to treat frequency-dependent DPG multipliers above ~2 Hz, where nominal
   responses roll off, versus the 0.03–0.07 Hz Rayleigh gains.
3. APG timing: per-event lag removes timing/phase as an independent check. A
   station drift model (fitted line) could rescue the 35 LT02 events whose
   coarse lag was off by > 3 s; this is not implemented.
4. Vertical absolute response is shared with the Rayleigh and microseism
   checks; agreement among them is not independent proof.
5. Scoring: how a conditional candidate from n = 5–9 events, with a model
   interval of ±10–17%, should map onto 0–1.
6. Deep APGs (LA33/LA34): should channels with recorded synchronized timing
   accept events whose 2–10 Hz pressure onset is unmeasurable, using the
   per-bin phase slope as timing QC instead? They would give an absolute-gauge
   check of the shared vertical response (pre-multiple, diagnostic). Not
   implemented: it changes an acceptance rule.

## Implementation limits

AIC onsets in a fixed 2–10 Hz band; the coarse envelope lag is good to only a
few seconds; TauP `iasp91` is continental (q is bounded by η ≥ 0.954 anyway);
80 events maximum per pair (deterministic priority); catalogue completeness
and location errors are as published by ComCat; ComCat revisions make
refetched catalogues differ (hashes recorded).

## For Codex: how to verify (read-only unless stated)

Work in a checkout of `claude/local-p-check` (the code is not on `main`). The
result folders below are git-ignored and live in
`~/Projects/bruit-fm-localp/results/LOCAL_P/`; read them in place.

1. `python -m pytest -q` (expect 271 passed) and `ruff check . && ruff format --check .`.
2. Plan without reading waveforms:
   `python scripts/run_local_p_campaign.py --campaign-dir results/LOCAL_P/campaigns/verify --dry-run`
   (expect 66 eligible, 39 `not_applicable_sample_rate`, 15 `no_full_pressure_response`;
   LT02 HDH `timing=drift_model`).
3. Reproduce the report table from the saved campaigns into a new path with
   `scripts/summarize_local_p_campaign.py` (`docs/local_p_campaign_v1.md` §5) and
   diff it against `docs/local_p_campaign_v1.csv`.
4. Bounded rerun (about 40 s) into a fresh folder, then compare the three
   `summary.json` files with `results/LOCAL_P/campaigns/apg_timing_v2/`:
   `python scripts/run_local_p_campaign.py --campaign-dir results/LOCAL_P/campaigns/verify_apg --networks XO --stations LT02 LA33 LA34 --pressure-ids XO.LT02..HDH XO.LA33..HDH XO.LA34..HDH`
5. Review points: sign/Fourier conventions (`local_p_physics.py`, method note
   §2–3), the forward model (`local_p_check.forward_model`), window/regime
   rules (`choose_window`), timing modes (`process_event`), output protection
   (`assert_fresh_output`) and that nothing writes outside `results/LOCAL_P/`.
