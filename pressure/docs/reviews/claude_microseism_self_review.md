# Self-review of the Claude microseism check (Sentinel + Socrates roles)

Branch `claude/microseism-all`, 15 October 2026. Requested by the user before the
repeated-week runs.

I reviewed my own work in both team roles:

* **Sentinel:** software (bugs, reliability, tests).
* **Socrates:** science (equations, conventions, claims).

This is a self-review, not the team's independent review. Findings are
corrected in this branch unless marked open.

**Scope.** I reviewed the following:

* `src/bruitfm/microseism_check.py`;
* `scripts/run_microseism_check.py`, `compare_microseism_methods.py`,
  `summarize_yv_microseism.py` and the wrapper scripts;
* tests;
* `docs/microseism_all_networks.md`, `yv_microseism_check.md` and
  `microseism_pressure_vertical_check.md`;
* the claims made in chat.

## Sentinel (software)

| ID | Severity | Finding | Resolution |
| --- | --- | --- | --- |
| SC-1 | High | Without `PYTHONPATH`, pytest in this worktree imported Codex's installed `bruitfm` from `~/Projects/bruit-fm/src`. New modules failed to import, and **shared modules (e.g. `paired_spectra`) would silently be tested from the wrong checkout**. All earlier test runs set `PYTHONPATH=src:.`, so their results stand. | `pyproject.toml`: `pythonpath = ["src", "."]`. Verified: tests pass without `PYTHONPATH`. |
| SC-2 | Medium | Scripts have the same hazard. During this review a regeneration ran without `PYTHONPATH` and failed on import. A `\| tail` hid the failure behind exit code 0. | The wrappers already export `PYTHONPATH`. Manual commands in the docs include it. Reruns now check exit status. |
| SC-3 | Medium | `compare_microseism_methods.py` aborted the whole comparison if one Codex diagnostic failed or one Codex campaign file was missing. | Per-station isolation (`codex_error`); missing campaign files are skipped. |
| SC-4 | Medium | Codex's `run_diagnostic` refuses an existing folder (correctly), so a compare could not be rerun. Isolation turned this into 84 silent `codex_error` rows on the first regeneration. | A saved diagnostic is reused only if its `source_comparison` and multiplier match; otherwise an error is raised. Tested. After the fix: 0 errors, 84 Codex phases, 46/46 identical to Codex's campaign. |
| SC-5 | Low | `--combine` overwrote phase-label counts instead of summing them, and did not detect a station supplied twice. | Counts are summed; duplicates raise. Tested. |
| SC-6 | Low | The wrapper left temporary logs. There was also dead code in `channel_poles`. | `trap` cleanup; dead code removed. |
| SC-7 | Low | Two pole sources: the YV summary reads `yv_rayleigh_independent_validation.json`, while the comparison reads response snapshots. Labels can differ. | Documented in `yv_microseism_check.md` (superseded note). |
| SC-8 | Medium | The new week-shift code and the weekly combiner were untested. | Tests for `shifted_start`, `week_runs`, `aggregate`, frequency matching, combine and Codex reuse. |

## Socrates (science)

| ID | Severity | Finding | Resolution |
| --- | --- | --- | --- |
| SS-1 | High | Frequency-dependent Deng et al. (2022) values were compared without matching frequency. "WD70 inside 0.91–1.71" tests nothing. WD48's sub-band (0.115–0.171 Hz) lies outside Deng's ~0.02–0.07 Hz, so the "+15% lead" was not like for like. | A published value is used only when the sub-band overlaps exactly one Deng band. Otherwise it is marked † with no percentage. WD66 is matched to 0.92 (+1.2…+2.7%). |
| SS-2 | High | The reported g bracket (H1–H2, mode 0) is a noise bracket only. It hides Earth-model dependence of up to the conditioning ratio (≤ 1.5). | A structure-free maximum (vertical incidence, H2) is added to the tables. Docs list the dependence per station: +1–3% (most YL) to +19–36% (RR28/29/36, A10W, C03W, J29A/J30A, LA33/34). |
| SS-3 | Medium | "The YL 5% offset points to a convention difference with the publication, not to either method." This is too strong. Both of our methods share the nominal responses, master depths and ρw. | Reworded: the offset comes from a shared input or from a convention in Zha & Webb; this cannot be resolved here. **Open.** |
| SS-4 | Medium | "RR50: Earth-model limit, not evidence against the Rayleigh gain." The Rayleigh reference uses the same model. | Reworded as a hypothesis in both docs. |
| SS-5 | Medium | The YV "published" values (Aminian et al. 2025) use the same Rayleigh P/a method family as column 2, so they are not an independent check of it. | Stated in both docs and in the source label. |
| SS-6 | Low | "69 explained" mixed fits ≤ 10° with approximate fits (10–20°). | Split: 65 within 10°, 4 approximate. Wording about polarity is tightened: relative P/Z polarity agrees; which channel carries the inversion is not established. |
| SS-7 | Medium | "RR38, RR40, RR52 reproduce across weeks." The weeks overlapped by 4–7 days. | Corrected earlier (`1f817d8`). Repeated independent weeks are now the next step. |
| SS-8 | Low | ρw = 1028 kg/m³ is a surface value. A 4–5 km column averages about 1035–1040 kg/m³ in situ, a sensitivity of about +1% in g (RR52 table). This is shared with the Rayleigh reference, so it cancels in the comparison but not in the absolute value. | Noted; no change. |
| SS-9 | Medium | One week per station gives no uncertainty estimate. Some brackets change between nearly independent weeks (RR29, RR34, RR36). | Repeated-week runner and combiner added. The stability thresholds were fixed before running: stable ≤ 5%, moderate ≤ 10%. |

| SS-10 | High | *(found with the repeated weeks)* The best-sub-band rule can select a 0.1–0.2 Hz sub-band where the mode-0 model is wrong. The true admittance can lie below mode 0, outside the noise + structure-free bracket. This produced a ×2 week-to-week jump at YL.C01W, while Codex's model-free amplitude stayed flat. 63 of 100 informative station-weeks fail a sub-band consistency check. | 35-station gain set marked superseded. Proposed rule (post hoc, to be applied blind): coherent sub-bands ≤ 0.1 Hz must agree within 10%, and repeated weeks must be stable. **Open.** |
| SS-11 | High | "C07W: microseism rejects the baseline 0.269" was based on one week. Three other weeks give 0.27–0.31, agreeing with the baseline. | Retracted in `microseism_all_networks.md`. C07W is a lead for time-varying response. |

Checked and found sound (no change):

* the D/P sign mapping (upward acceleration raises bottom pressure, so correct
  polarity gives D/P at 180°);
* H1 biased low by pressure noise and H2 biased high by vertical noise;
* infragravity leakage, which is negligible above 1.5 f_IG;
* the gravity term, which cancels at a Lagrangian gauge in the long-wave limit
  (RR52 sensitivity −1%);
* the selection of the best sub-band, which uses conditioning, not reference
  agreement;
* thresholds fixed before the YL/7D/XO runs;
* Deng's Cp direction (multiplies pressure, the g direction);
* Zha & Webb inversion (1/g).

## Open

1. **SS-3.** Find the YL conventions in Zha & Webb (2016): the nominal DPG
   sensitivity, ρw and depth source. The paper was not among the supplied PDFs.
2. **SS-4 / SS-2.** Run local Earth models for RR50, 7D and XO, and verify
   master depths for the stations whose Earth-model dependence is above 10%.
3. **SS-9.** Run the repeated weeks
   ([microseism_all_networks.md](../microseism_all_networks.md#repeated-weeks-empirical-uncertainty)).
