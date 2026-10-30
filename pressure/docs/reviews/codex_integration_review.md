# Codex integration review of Sentinel and Socrates

Decision: **approve the corrected engineering changes for merge into main**.
This decision accepts the software fixes and records the scientific audits;
it does not accept every audit interpretation or certify the provisional gains.

Reviewed main: `66a3d6a`. Reviewed independent branch:
`review/sentinel` at `6f074c9d5a6f82f21b523f35f9b58a13b12cf990`.
The independent work branched from `faa39ff`; this integration also retains
main's subsequent 7D/XO/APG implementation and four-network reporting.
Review date: 12 October 2026. No conflicts in the trial merge.

## Code verdict and correction

- Decide event QC before attempting a scalar fit. When fitting hits the search
  boundary, preserve spectra, plots, phase/coherence/overlap status and a separate
  `fit_reason`. A QC-valid boundary failure remains unresolved. No additional
  event becomes accepted through this change.
- Select tilt horizontals from the paired vertical's band/instrument prefix:
  `LHZ` uses `LH1/LH2`, `BHZ` uses `BH1/BH2`. No-tilt runs are unaffected;
  horizontal N/E support remains outside this change.
- Refuse fresh exports over either existing master file. Serialize inventory
  append and restricted result updates with the shared master lock.
- **Additional defect found and fixed by Codex:** the public fresh-export path
  initially bypassed that lock, leaving a race between simultaneous exporters
  or an exporter and inventory synchronization. Fresh creation now holds the
  same lock across checking existence and writing. The already-locked append
  path uses an internal creation helper to avoid acquiring the lock twice.
  The regression failed before this fix and passes afterwards; another writer's
  manually scored CSV remains unchanged and is not paired with a fabricated XLSX.

The fitter, gain bounds, acceptance thresholds, response polarity, timing guard,
reference formula and default tilt behavior have not changed. Locking uses
`fcntl`, supported by the project's macOS and Linux environments. Cooperating
writers serialize; an external editor is not controlled by this lock.

## Evaluation on existing networks

The integrated checkout passed **155 pytest tests** (77 existing ObsPy/MiniSEED
warnings), Ruff check/format, and runner Bash syntax. Socrates' read-only
V1–V11 calculations completed and reproduced its reported numerical examples;
these printed diagnostics are not eleven additional assertion-based tests.

`scripts/review_saved_calibrations.py` replays saved event fits and reconstructs
selected raw windows. These are consistency/regression checks using the reviewed
production estimator, not an independent absolute calibration measurement.

| Latest preserved workflow | Stations | Station estimates | Accepted event fits replayed |
|---|---:|---:|---:|
| YV no tilt | 43 | 16 | 49 |
| YL all stations, no tilt | 25 | 22 | 297 |
| 7D all stations, no tilt | 15 | 4 | 18 |
| XO all gauges, no tilt | 6 | 4 | 52 |
| **Total** | **89** | **46** | **416** |

All 416 saved accepted-event gains reproduce exactly; saved station medians,
reciprocals and counts agree. Event totals include accepted events at stations
that have too few events for a station estimate. The original tilt YV workflow
has 48 accepted events and is a different preserved baseline.

Raw probes covered **18 event windows at 17 stations**:
YV RR03/RR52/RR55; YL A02W/C07W/B01W/C08W;
7D J30A/J38A/J45A/J54A; XO LA33/LA34/WD66/WD70/LT02/WD48.
C07W has two probes, spanning its smallest and largest accepted gains.
Fourteen candidate event gains reproduce with maximum absolute difference
`8.8124e-10`, including the XO APGs. RR55 and B01W have zero phase-valid bins;
LT02 remains incoherent, and WD48's selected window is missing.

A separate B01W end-to-end one-event run saved spectra and returned
`polarity_review`, zero candidates, and the boundary `fit_reason`. This confirms
the new diagnostic behavior on real data without overwriting the original run.

Two additional raw daily microseism checks corroborate C07W's relative vertical
amplitude state: C07W/C08W is 0.901 on 2010-03-14 and 0.431 on 2010-08-04;
C08W/C09W remains 0.835 and 0.846. These are raw-count RMS diagnostics in
0.12–0.3 Hz, not pressure-calibration factors.

**Preservation:** 557 checked source/product/master/response/model SHA-256 hashes
and 34 raw-file size/mtime records remained unchanged during the replay.
The B01W run separately rechecked its original workflow/result/master hashes;
the daily checks hashed six raw files before and after. These scopes overlap;
they are not an archive-wide byte audit. No full catalogue batch was rerun, no
gain or response correction was applied, and no original score/result was edited.
New products are isolated in the review worktree's
`results/REVIEWS/SENTINEL_INTEGRATION/`; curated records are in
[codex_integration_validation.json](codex_integration_validation.json).

## Scientific conclusions and qualifications

1. **C07W:** the two relative vertical-amplitude states are supported by daily
   neighbor checks and the bimodal event gains. Its saved median is unsuitable
   as an accepted pressure calibration without separating these states.
   A hardware wiring/feedback fault is a hypothesis, not an established cause.
   Origin-day labels are coarse: the low-g event on 2010-04-17 occurred after
   a late-day state transition. The normal-state reciprocal near 1.90 is a
   diagnostic subset, not a newly accepted station estimate.
   Peer event comparisons must also screen the peer's own coherence and window:
   C08W's saved 2010-05-31 event has median coherence only 0.285. The claim
   that pressure amplitudes agree near unity on every matched event is too broad;
   it cannot be applied to that event. Shared event IDs alone are insufficient.
2. **B01W:** the phase failure is binding. Widening scalar bounds alone would
   not make the coherent boundary events pass. The single-pole high-pass model
   is a useful amplitude/phase diagnostic; it does not establish the exact
   hardware cause or independently validate the full published-band coefficient.
3. **RR55:** zero phase-valid bins are reproduced, but Socrates does not confirm
   the same single-pole mechanism as B01W. Resolve its actual sensor/response
   identity separately. Do not label it repaired or faulty from this audit alone.
4. **YV poles:** the published causal RC response has a negative pole, whereas
   several nominal XML responses have positive-real poles. An in-memory RR03
   diagnostic confirms that reflecting a real pole leaves response amplitudes
   unchanged at fixed bins but changes phase. Re-running phase QC could change
   accepted bins/events and therefore station gains. No inventory was modified;
   pole reflection is not approved as a production correction by this merge.
5. **INSU glitches and YV published comparison:** Aminian et al. (2025) Table 1
   supplies eight pressure multipliers, in our multiplication direction.
   Socrates V11 compares the **original tilt YV baseline**, not the latest no-tilt
   run: deviations are +0.3% to +6.1%, median +3.0%. This is a method-lineage
   reproducibility comparison. Their glitch treatment/event selection differs;
   the positive offset does not prove that missing glitch removal caused it.
6. **Timing and uncertainty:** microsecond header offsets are small in phase,
   but timing guards were not relaxed. Actual logger clock errors are distinct.
   The C07W IID day-bootstrap omits serial dependence within multi-day states;
   its intervals are descriptive, not calibration confidence intervals.
   Model/response/window systematics remain beyond event-to-event scatter.

Primary sources checked for the YV response/glitch and gain conventions:
[Stähler et al. (2016), Appendix A](https://adgeo.copernicus.org/articles/41/43/2016/adgeo-41-43-2016.pdf)
and [Aminian et al. (2025), Table 1](https://www.ipgp.fr/~stutz/2025-aminian-crawford-stutzmann-al-GJI.pdf).
The Sentinel and Socrates reports remain attributable to their authors;
Codex qualifications are explicitly marked where their initial statements were
too strong. Historical saved numerical snapshots are preserved.

## Next concrete scientific task

Design a separate response-convention experiment for selected YV pressure epochs:
verify hardware identity and the causal pole against the primary response equation,
then compare phase/QC on a few existing event windows using an explicit alternative
response in a new output root. Report changes in accepted bins/events before any
gain comparison. Do not silently replace XMLs or rerun the full YV batch.
Handle RR55's sensor identity and C07W's vertical-state selection separately.
LT02 timing, WD48 coverage, local dispersion/depth, held-out uncertainty and
published 7D/XO coefficient matching remain open.
