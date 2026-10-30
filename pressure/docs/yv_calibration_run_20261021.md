# YV all-methods calibration run, 21 October 2026

Session `results/EXPERIMENTS/YV_ALL_METHODS/run_20261021T122108Z` (main
checkout, local only), started by the owner with `scripts/run_yv_all_methods.sh`
at 12:21 UTC and finished at 13:04 UTC. This report records what ran, what it
shows and what it does not show. No master, response, raw waveform, saved
result, score or rating changed. No gain was applied. **RR52's frozen reference
stays g = 0.6787278140712512.**

## Update 22 October: clean reproduction (confirmed)

* The owner's uncommitted edit was saved as
  `~/Projects/bruit-fm-uncommitted-rayleigh-guard-removal-20261021.patch`
  (SHA-256 `9bcd9007…`; re-applies cleanly). The committed file was then
  restored. The edit is preserved, not discarded.
* Session `results/EXPERIMENTS/YV_ALL_METHODS/run_20261022T082359Z`: code
  `d572bdf`, no uncommitted changes, `--rayleigh-events bounded`, 08:24–08:40
  UTC.
* **Rayleigh:** **RR52 = 0.6787278140712512 exactly**.
  * **43 of 43 stations identical** to the frozen workflow in status, g and
    candidate count; same 16 / 25 / 2 split.
  * Every result records the clean code hash `89cfb7b8…`.
  * Only event-label difference: RR55, where 4 non-candidate events changed
    from "unresolved, grid boundary" to "polarity_review". No gain or station
    status changed.
* **Local P:** 26 / 9 / 8 again; 41 of 43 pairs identical to the frozen v2
  campaign. **The RR11 and RR50 event flips also occur on clean code**, so the
  uncommitted edit did not cause them. They follow from reviewed code changes
  made after 16 October.
* **Conclusion:** the latest committed package reproduces the frozen YV
  Rayleigh reference exactly. The 21 October microseism results (identical to
  14 October) and local-P results (identical to today's clean run) stand
  despite the edit. Only their recorded provenance is affected.
* The provenance fix (code hash at import) was merged on 22 October as PR #5
  (`90adc26`).

## Bottom line (21 October run)

* **Rayleigh.** The run used the package's current default of evaluating
  **all** eligible events, not the frozen run's bounded selection. It is
  therefore a new, larger-sample result, not a reproduction.
  * The same 16 stations get an estimate as in the frozen run (16 with an
    estimate / 25 unresolved / 2 reference-band exclusions, identical sets).
  * Accepted events: **48 → 442**.
  * Per-station gain change: median **+0.01%**, range −3.74% to +4.02%, with
    **12 of 16 within ±2%**.
  * RR52: **0.6750** from 40 events, against the frozen 0.6787 from 3 (−0.54%).
  * More events rescued **none** of the 25 unresolved stations, so their
    blockers are systematic, as the diagnosis said.
* **Microseism.** 41 stations checked; 2 unresolved pairing (RR12, RR42).
  * All key numbers are **identical to the 14 October run** at all 43 stations.
  * The frozen Rayleigh g lies inside the best microseism bracket at 7 of 16
    stations; the pole-trial g at 14 of 31.
* **Local P.** 26 `no_suitable_events`, 9 `unresolved`, 8 `diagnostic_only`;
  no station value.
  * 41 of 43 pairs are identical, event by event, to the frozen v2 campaign of
    16 October.
  * Two events changed outcome; the cause is not established (see below).
* **Code provenance problem.** An uncommitted edit to
  `src/bruitfm/rayleigh_calibration.py` appeared in the main checkout at
  **12:33:58 UTC, during the run**. The microseism and local-P stages
  **executed the edited code**.

## 1. Code actually executed (confirmed)

* Committed code: main `d572bdf`. `rayleigh_calibration.py` SHA-256 starts
  `89cfb7b8…` (committed) and `e4c09ff4…` (edited).
* The edit removes the reviewed refusal of fractional native-grid boundary
  tears (two lines in `native_resampling_samples`). The deployment review lists
  this refusal as one of the latest software fixes. The edit is not committed
  and was not made by Claude.
* **Rayleigh:** the process started at 12:21, before the edit. Python had
  already loaded the module, so it executed the committed code. But each result
  records `code_sha256` by reading the file **at write time**:
  * 15 results recorded `89cfb7b8…`;
  * 26 recorded `e4c09ff4…` (written after 12:33:58);
  * 2 recorded none (reference-band exclusions).
  This recorded hash does not prove what code ran.
* **Microseism** (started 12:52) and **local P** (started about 13:00): new
  processes that executed the edited code. Both import
  `bruitfm.rayleigh_calibration`.

## 2. Rayleigh: all events versus the frozen bounded selection

* Cause of the difference (confirmed in code and data):
  * The frozen results were written on 2026-10-09 at 08:49 UTC with the old
    bounded pilot: stop after 3 candidates or 8 attempts per magnitude level.
    The recorded outcomes are `target_reached` or `attempt_limit_reached`.
  * Commit `a19bc3b` (2026-10-09 11:48 UTC) made `all` the default.
  * The runner passed `--rayleigh-event-sampling all`. It now defaults to
    `bounded` for reproduction.
* Same 16 stations with an estimate. The table shows candidate / attempted
  events.

| Station | Frozen g (bounded) | Events | All-events g | Events | Change |
|---|---:|---:|---:|---:|---:|
| RR03 | 1.1241 | 3/6 | 1.1664 | 5/71 | +3.76% |
| RR07 | 1.0792 | 3/4 | 1.1021 | 7/67 | +2.12% |
| RR08 | 1.1597 | 3/4 | 1.1802 | 12/68 | +1.77% |
| RR09 | 1.1937 | 3/4 | 1.1843 | 10/68 | −0.79% |
| RR11 | 1.0827 | 3/4 | 1.0643 | 7/68 | −1.70% |
| RR13 | 1.2236 | 3/6 | 1.2236 | 5/63 | +0.00% |
| RR14 | 1.0957 | 3/4 | 1.0547 | 12/68 | −3.74% |
| RR25 | 1.1626 | 3/6 | 1.1788 | 9/71 | +1.39% |
| RR28 | 0.8195 | 3/4 | 0.8113 | 52/72 | −1.01% |
| RR29 | 0.7299 | 3/4 | 0.7238 | 46/72 | −0.84% |
| RR34 | 0.5853 | 3/4 | 0.5738 | 22/34 | −1.95% |
| RR36 | 0.6536 | 3/4 | 0.6539 | 45/75 | +0.05% |
| RR38 | 0.6365 | 3/4 | 0.6366 | 63/75 | +0.02% |
| RR40 | 0.7426 | 3/4 | 0.7724 | 61/76 | +4.02% |
| RR50 | 0.7168 | 3/4 | 0.7210 | 46/76 | +0.59% |
| **RR52** | **0.6787** | 3/4 | **0.6750** | 40/76 | **−0.54%** |

* Network median of the 16: **0.9494 frozen → 0.9330 all-events**.
* INSU DPG stations (RR28–RR52) keep many more events (22–63) than DEPAS
  stations (5–12), which still lose most events to coherence.
* **Interpretation.** The frozen values agree with a nine-times larger event
  sample within about ±4%. This supports their stability under event selection;
  it is not an independent absolute calibration, because both use the same
  model and vertical response. Which policy becomes the reference is the
  owner's decision. The frozen RR52 value is unchanged.

## 3. Microseism (week of 2013-06-10)

* 41 `checked`; 2 `unresolved_pairing`: **RR12** (flat archive) and **RR42**
  (failed vertical).
* Compared with run `run_20261014T175553Z_30002f23` (preserved worktree): **43
  of 43 stations identical** in status, week, retained windows, coherence, the
  h1/h2 gains and the bounds. Only the newer phase-summary label differs.
  The edited guard therefore did not change any microseism number here.
* Frozen Rayleigh g inside the best microseism bracket: **7 of 16**.
  * Inside: RR08, RR09, RR13, RR14, RR29, RR40, RR52.
  * Outside: RR03, RR07, RR11, RR25, RR28, RR34, RR36, RR38, RR50.
* Pole-trial g inside the bracket: **14 of 31**.
* RR52: best sub-band coherence² 0.986 with g 0.671–0.681, which brackets the
  frozen 0.6787. Phase −1.0°.

## 4. Local P (reviewed v2 protocol)

* 43 eligible pairs, all run: 26 `no_suitable_events`, 9 `unresolved`,
  8 `diagnostic_only`.
  * Unresolved: RR09, RR12, RR14, RR41, RR43, RR44, RR46, RR47, RR48.
  * Diagnostic only: RR10, RR11, RR36, RR38, RR40, RR50, RR52, RR53.
* Compared with the frozen campaign `user_all_20261016T120124Z_17172`: **41 of
  43 pairs identical** event by event. Two events changed, in opposite
  directions, both with reason `no_valid_band`:
  * RR11 `us2013pjc9`: rejected → accepted.
  * RR50 `usb000ja9l`: accepted → rejected.
* Cause **not established**. It could be the removed guard or other reviewed
  code changes since 16 October. Only a rerun on clean code can separate them.
* RR52: 1 of 22 events accepted; g 0.56 (4–8 Hz), 0.44 (8–16 Hz). Diagnostic
  only; vertical onset not measured.

## 5. RR52 across methods

| Method | g | Basis | Status |
|---|---:|---|---|
| Rayleigh, frozen | **0.6787278140712512** | 3 events | Reference; unchanged |
| Rayleigh, pole trial | 0.6787278140712512 | 3 events | Identity (no positive pole) |
| Rayleigh, all events (today) | 0.6750 | 40 events | −0.54%; development comparison |
| Rayleigh, held-out | +0.66% median prediction error | 14 events | Earlier validation |
| Microseism (today) | 0.671–0.681 best sub-band | 167 windows | Diagnostic; brackets the reference |
| Local P (today) | 0.56 / 0.44 (4–8 / 8–16 Hz) | 1 of 22 events | Diagnostic only |
| Published (Aminian 2025) | 0.67 | — | −1.3% |

## 6. Why 25 stations stay unresolved (from the diagnosis; unchanged today)

* **Coherence** (coh² ≥ 0.95) removes most events everywhere.
* **Positive hydrophone pole in the metadata**, contradicting Stähler Eq. A6:
  15 stations recover in the pole trial (RR01, RR05, RR06, RR10, RR16, RR19,
  RR20, RR21, RR26, RR30, RR43, RR44, RR46, RR48, RR55).
* **Metadata defects** (pilot of 21 October; development values only):
  RR22 and RR56 (unpaired vertical poles), RR17 and RR18 (pole about 10×).
* **Hardware or data:** RR12 (flat −3077 dB archive), RR31 (DPG artefacts),
  RR42 (failed vertical), RR45 (noisy vertical, too few events), RR53
  (unexplained).
* **Polarity:** RR47.
* RR41 and RR57 are reference-band exclusions, counted separately.

## 7. Runner changes in this branch (software only)

* `--rayleigh-events bounded|all`, default **bounded**.
  * In bounded mode the run stops unless RR52 returns exactly
    0.6787278140712512.
  * The output folder name records the mode.
* The runner **refuses uncommitted code** in `src/` or `scripts/` unless
  `--allow-dirty` is given; the code state is recorded either way.
* It **stops if the code changes during the run**. It checks a hash of HEAD
  plus the uncommitted diff before every stage.
* `--check` against the current main checkout correctly refuses because of the
  uncommitted edit, and passes with `--allow-dirty`.

## 8. Open items and decisions (owner)

1. **The uncommitted edit** in the main checkout. Keep it, or restore the
   committed file? It removes a reviewed safety guard. It is not pushed and
   was not changed by Claude.
2. ~~Bounded Rayleigh reproduction on clean code~~. **Done 22 October:** RR52
   exact; 43/43 stations identical.
3. ~~Local-P rerun on clean code~~. **Done 22 October:** the RR11/RR50 changes
   persist on clean code, so they are not caused by the edit.
4. **Which Rayleigh event policy is the reference:** frozen bounded (3 events
   per station) or all events (9× more). The owner decides; nothing is adopted
   here.
5. ~~Provenance fix~~. **Merged 22 October as PR #5:** the Rayleigh code
   hash is now taken at import.
6. PRs #3 and #4 were merged on 22 October. The DEPAS/RESIF drafts in
   `docs/yv_response_confirmation.md` are still unsent. Validation of the 21
   October pilot gains waits for the new computer. The current overall status
   is in [project_status_20261022.md](project_status_20261022.md).

## Sources (local; not in Git)

| Item | Path |
|---|---|
| Session logs | `results/EXPERIMENTS/YV_ALL_METHODS/run_20261021T122108Z/` |
| Rayleigh all-events workflow | `results/NETWORKS/YV_RAYLEIGH_NO_TILT_20261021T122108Z/_workflow/rayleigh.json` |
| Frozen Rayleigh workflow | `results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json` |
| Microseism today / 14 Oct | `results/EXPERIMENTS/MICROSEISM_YV/run_20261021T125204Z_b2a2ad63/`; `~/Projects/bruit-fm-claude/results/EXPERIMENTS/MICROSEISM_YV/run_20261014T175553Z_30002f23/` |
| Local P today / frozen | `results/LOCAL_P/campaigns/yv_20261021T122108Z/`; `~/Projects/bruit-fm-localp-review/results/LOCAL_P/campaigns/user_all_20261016T120124Z_17172/` |

All comparisons above were read-only and wrote nothing into `results/`.
