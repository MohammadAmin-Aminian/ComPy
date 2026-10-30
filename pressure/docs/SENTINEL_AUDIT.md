# Sentinel audit record

> **Codex integration qualification:** the daily relative vertical states are
> supported, but a C07W hardware fault is not proven. Peer-event ratios require
> independent peer QC: C08W's 2010-05-31 event has median coherence 0.285,
> so near-unity pressure ratios cannot be assumed for every shared event ID.
> See [the integration review](reviews/codex_integration_review.md) for bounded raw checks and limits.


Independent review of BRUIT-FM. Software fixes live only on branch
`review/sentinel` (pushed to `origin`; same history as the local working branch
`sentinel/audit-2026-10`) until the user decides on integration. Nothing here
applies a calibration gain, edits a master, response, waveform or saved result,
or reruns a YV/YL batch.

**Roles (from 2026-10-12).** Codex develops and integrates. Sentinel reviews
software architecture, bugs, performance, reliability and testing. Socrates
reviews scientific correctness. The scientific findings below (SENT-001, -004
to -007, -010, -012) are handed to Socrates for independent scientific review;
see `docs/sentinel/HANDOFF_SOCRATES.md`.

## 1. Baseline and environment

| Item | Value |
|---|---|
| Reviewed baseline | `faa39ff` (local `main` = `origin/main`, clean tree, last five GitHub Checks passed) |
| Sentinel checkout | `/Users/mohammadamin/Projects/bruit-fm-sentinel`, separate clone, branch `sentinel/audit-2026-10` |
| Environment | `.venv` in the Sentinel clone, built from the conda `bruitfm` Python 3.13.15 and pinned to the shared env's versions (ObsPy 1.5.1, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.6, TiSKitPy 2.3.1, disba 0.7.0, Numba 0.68.0, pytest 9.1.1, Ruff 0.16.9). The shared env is untouched; it imports `bruitfm` from Codex's checkout, so it was not used for tests. |
| Baseline checks | 130 passed; Ruff check/format and `bash -n scripts/run_yl_all.sh` passed. Tests ran with `HOME` redirected to a scratch folder; production master stats were unchanged afterwards. |
| Upstream during review | `main` advanced to `4387fbe` (Codex: read-only C07W/C08W/C09W matched-event viewer; no core changes). `git merge-tree` of this branch onto it is conflict-free; the merged tree passes 138 tests and Ruff check/format (68 files). |
| Current branch checks | 141 passed (about 43 s); Ruff check/format passed. Trial merge of `63ef818` onto `origin/main` `4387fbe`: conflict-free, merged tree passes 141 tests and Ruff check/format (69 files). |
| Read-only access | Local `Data/YL`, `Data/YV` StationXML/MiniSEED and Codex's `results/`, `exports/` (read, never written). |
| Not accessible | Zha & Webb (2016) full text (Wiley 403) and Zha's thesis (bot protection). Table 1 values were taken from the repository transcription. |

## 2. Review coverage

"Full" means every line was read.

| Area | Depth |
|---|---|
| All 15 modules in `src/bruitfm/` | Full (second pass completed `network.py`, `pressure_test.py`, `infragravity_test.py`, `godin_pilot.py`, `nearby_test.py`, `plot_ppsd.py`) |
| All 6 scripts in `scripts/` | Full |
| All 15 test modules | Full; each checked for whether its expected values are independent of the code under test |
| Package data (`__init__.py`, all 5 `reference_data/` files) | Full. All 17 Brown rows recomputed from their stored pixel coordinates and axis mapping (max error 5e-11 dB, including the −120 dB conversion). RR52 layer and dispersion files match their provenance SHA-256. disba 0.7.0 regenerates the table from the layers exactly (phase) and to 1.5e-11 m/s (group). The Webb points and the Brown source raster cannot be checked against primary sources here. |
| Method and report docs (21 files under `docs/`, plus `README.md`, `CODEX_NEXT_CHAT.md`) | Full, checked against current code |
| `docs/research_log.md`, `CODEX_HANDOFF.md`, `SESSION_NOTES.md` | Historical narratives: searched for runnable commands and current-behaviour claims, and final sections read; not re-read line by line |
| JSON/CSV sidecars | Read where a finding depended on them (`project_audit.json`, `yl_all_zha_comparison.csv`) |

### Test-suite independence

Most tests use independent expected values: hand-built FFTs and Parseval,
SciPy periodograms, analytic high-pass corners, `brentq` dispersion roots,
the acoustic ODE, a known 2× amplitude (6.02 dB), and a known delay or complex
gain. Gaps found and closed in `3f663c5` (all pass on current code):

1. The Rayleigh synthetic fixture builds pressure as `2jπf·V·T`, the same
   derivative the code uses, so a flipped convention would cancel. A new test
   uses an analytic time-domain derivative at leakage-free bins. Production
   `event_spectra` returns ratio exact to 1e-9 and phase 0 to 1e-6°.
2. Paired spectra were only tested with displacement (`M`) sensors. A
   velocity-sensor test now checks ObsPy `DISP` integration (exact).
3. The seismic PPSD path was only checked by its unit string. It is now checked
   numerically as acceleration density (ω²/|R_vel|²), to 1e-4 dB.
4. The CSV/Excel row-disagreement guard had no test. It refuses, leaves both
   files byte-identical and makes no backup.
5. Conflicting overlapping records are shown to become missing samples.

The prior audit's saved-product script recomputes coherence, ratios and gains
from saved arrays with the production formulas. That is arithmetic consistency,
not independent estimation, and it does not recompute `expected_pa`.

### Prior-audit fixes re-checked

| Prior fix | Sentinel verdict |
|---|---|
| 1. Short PPSD pilot reuse | Correct. The full-end computation mirrors `process_epoch`; tests cover both bounds. |
| 2. Pressure-unit aliases | Correct for selection. No remaining exact-`"PA"` matches in `src/`. |
| 3. Explicit-depth preflight | Correct. |
| 4. Coordinate fingerprint | Correct. |
| 5. Auxiliary fit failures | Correct for H1, raw and median-curve fits. The primary fit of non-accepted events had the same failure mode; see SENT-002. |

### Independent numerical checks

| Kind | Check | Result |
|---|---|---|
| Derivation | Water-column reference: z down, P(0)=0, ∂P/∂z(H)=ρ·a_up gives ρ·tan(rH)/r. The tan/tanh series 1 ± x²/3 + 2x⁴/15 is correct. | Agrees with `water_ratio` |
| Implementation | Separate estimator (`scripts/sentinel_yl_review.py`: SciPy Welch CSD, cross-spectral response division, own water formula) on 16 saved YL candidate events at C07W/B04W/A02W/C08W | Gains reproduced within 0.07% |
| Implementation | ObsPy 1.5.1 PPSD averages **dB** within period bins and clamps zero power to `dtiny`, i.e. 10·log10(2.2e-308) = −3076.5 dB | Confirms the Godin dB bin-mean and the audit's floor interpretation |
| Shares responses | Rayleigh phase velocity 4000 → 3500 m/s changes T by −0.06% to −1.3% (H = 1655–2822 m, 0.03–0.07 Hz) | The borrowed Earth model is a minor term in the gains |
| Field consistency | Cross-station comparisons for C07W, B01W and RR55 (below) | Not absolute calibration |
| Convention | Analytic-derivative synthetic through production `event_spectra` (new test) | Ratio exact, phase 0°: iω sign, response division and upward Z agree |
| Metadata algebra | Positive-real pole phase 2·atan(p/ω): RR03 28.2°, RR22 55.4° at 0.03 Hz | Reproduces `yv_response_window_validation.md` |
| Timing | MiniSEED headers, A02W: 2010-04-06 LDH `…00.666701`, LHZ `…00.666800`; 2010-05-05 both `…00.595100` | Same-digitizer offsets at the 100 µs header quantum (SENT-005) |

### Software-engineering pass

| Area | What was checked | Outcome |
|---|---|---|
| Static analysis | Ruff with bug-prone families (B, BLE, S, PLW, PLE, RUF, ARG, PERF, TRY, SIM, C4, DTZ, PTH) run from the command line; project config unchanged | 336 hits, mostly style (221 `TRY003` long messages, 31 `dict()` calls). Each bug-prone hit was read: the 10 broad `except Exception` handlers all record a `failed` status per channel/event (failure isolation, intended); the PPSD chunk handler still advances its cursor; `B023` in `validate_yv_rayleigh.py:286` is a false positive (closure used in the same iteration); `assert` appears only in scripts. One unused loop variable in Sentinel's own script was fixed. |
| Type checking, coverage | mypy and coverage are not installed; not installed by Sentinel | Not measured. Recommend adding `pytest-cov` to the `dev` extra. |
| Write paths | Every write site in `src/` | JSON state uses `atomic_json` (temp + rename); PPSD checkpoints commit NPZ before manifest; masters are staged, backed up and restored on failure. Diagnostic JSON/CSV/NPZ use direct writes (reproducible, see Minor observations). Gaps: SENT-013, SENT-014. |
| Concurrency | Process pool, cancellation, locks | Workers never write masters (parent only). Rayleigh runs hold an flock per workflow folder. Master files had no inter-process lock (SENT-014). |
| Performance | Test-suite timings; one master commit timed on a scratch copy of the production master (16,884 rows × 140 columns) | Suite 43 s (slowest test 9.5 s). One master commit 8.4 s, of which 2.7 s is reading the Excel file, plus an 11 MB backup per commit (SENT-015). |
| Architecture | Function sizes (AST) | 151 functions; `run_calibration` is 683 lines, and nine more exceed 180 lines (`run_godin_pilot`, `run_nearby_test`, `run_paired_spectra`, `process_epoch`, `run_stage`, `comparison_plots`, `network.main`, `run_infragravity_test`, `_run_network`). Each mixes I/O, estimation, QC and plotting, which is why SENT-002 hid in one branch. Recommend extracting pure estimation/QC functions when they are next changed, not as a standalone refactor. |

## 3. Findings

Severity and confidence are assessed separately.

### SENT-001: C07W's outlier comes from an intermittent vertical-channel amplitude state, not the pressure gauge

- **Category / severity / confidence:** Scientific concern / High / Strong evidence.
- **Baseline:** `faa39ff`; saved `YL_ALL_NO_TILT` workflow; raw `Data/YL`.
- **Location:** saved C07W run `run_621bf2f89b89`; affects `docs/yl_all_zha_comparison.md`.
- **Problem:** C07W's 21 candidate gains are bimodal, not scattered.
  - 17 events have g = 0.259–0.285 (1/g ≈ 3.7).
  - 4 events have g = 0.524–0.543 (1/g ≈ 1.90).
- **Evidence:**
  1. Saved spectra, C07W against C08W (29 km away) on the same events: pressure amplitude ratio ≈ 1.0 in both groups. Vertical acceleration ratio ≈ 0.5 for the low-g group and ≈ 1.0 for the high-g group.
  2. Raw counts, band-passed 0.03–0.07 Hz by independent code, reproduce this: LHZ ratio 0.43–0.56 against 0.80–1.05.
  3. The daily microseism-band (0.12–0.3 Hz) ratio of C07W LHZ to C08W LHZ is bimodal over 325 days:
     - 251 days near 0.44 and 68 near 0.85 (state ratio 1.91); 6 other days (deployment start 2009-11-25/26/27, and three transitional days near 0.62–0.64).
     - C08W/C09W is steady at 0.86; pressure C07W/C08W is steady near 1.1.
  4. Every high-g candidate falls in a "normal" episode (03-08→03-14, 05-27→06-07, 06-18, 09-24→10-15). The 2010-04-17 23:15 event (g = 0.275) is in the reduced state by its event-window ratio (0.506); the state switched late that day.
  5. Metadata cannot explain an intermittent factor: C07W/C08W/C09W LHZ sensitivities agree to 1%, and the DPG responses are identical.
- **Interpretation:** for most of the deployment, C07W's vertical records about half the ground motion its metadata imply.
  - The normal-state events give 1/g ≈ 1.90, which is 4.5% below Zha's 1.99. That is the same offset as the other 21 stations.
  - The cause (sensor or feedback fault, mass centring, digitizer setting) is not established.
- **Impact:** the existing C07W station median (g = 0.269) measures vertical-channel behaviour, not DPG sensitivity. The same contamination would affect any C07W product that uses LHZ, including paired and compliance products.
- **Correction:** none applied. Recommendations:
  - Do not use the C07W station median as a pressure calibration.
  - Consult the OBSIP/WHOI deployment logs.
  - Consider a cross-station vertical-consistency diagnostic before accepting gains. That is a new QC policy and needs user approval.
- **Verification:** `docs/sentinel/c07w_daily_vertical_ratio.csv`, `c07w_event_gain_by_state.csv`, `summary.json`.
- **Status:** Open.

### SENT-002: a gain-grid boundary discarded the event's spectra and its QC outcome

- **Category / severity / confidence:** Confirmed defect / Medium / Confirmed.
- **Location:** `run_calibration`, `src/bruitfm/rayleigh_calibration.py` (event fit block).
- **Problem:** for a coherent event, `fit_multiplier` raising "Grid minimum lies on a boundary" escaped to the outer handler. As a result:
  - The spectra NPZ/CSV, comparison and waveform plots and grid CSV were never written.
  - An event that had already failed the phase guard (`polarity_review`) was relabelled `unresolved` with the boundary reason. The binding QC reason was hidden.
- **Evidence:**
  - New test `test_boundary_fit_keeps_spectra_and_actual_qc_outcome`. At baseline both polarities gave `unresolved` / boundary, with only `arrival_review` saved.
  - Real data: all four B01W "boundary" events have median γ² of 0.97–0.99 but 0 phase-valid bins (see SENT-004).
- **Impact:** existing saved runs are unchanged; they keep the old labels. This led to the documented reading that B01W is limited by the 0.1–3 search bound.
  - Every saved boundary-labelled event in every workflow has this pattern: 12 records in total. These are YV RR55 ×4 in both `TEST` and `YV_RAYLEIGH_NO_TILT`, and YL B01W ×4. All have median γ² 0.97–0.99 and **0 phase-valid bins**, so all should read `polarity_review`.
- **Correction:** commit `2d4e8d6`.
  - The QC status is decided before the fit.
  - A fit failure records `fit_reason` and keeps spectra and plots.
  - Accepted-but-boundary events stay `unresolved`.
  - No acceptance rule, bound or gain changed; no new candidate can result.
- **Verification:** new test (2 cases); full suite 132 → 133 passed after SENT-003; Ruff passed.
- **Status:** Fixed on Sentinel branch.

### SENT-003: tilt correction hard-coded `.BH1`/`.BH2` horizontals

- **Category / severity / confidence:** Confirmed defect / Low / Confirmed.
- **Location:** `select_calibration_channels`, `src/bruitfm/rayleigh_calibration.py`.
- **Problem:** with `tiskitpy` or `local_eigen` tilt correction, any station whose horizontals use another band code (YL: LH1/LH2; HH*) failed with "Expected one overlapping .BH1 epoch". The vertical and pressure are found by units, but the horizontals were found by channel letters.
- **Evidence:** new test `test_tilt_horizontals_follow_the_paired_vertical_band` (LDH/LH1/LH2/LHZ) failed at baseline with exactly that error. After the fix TiSKitPy recovers g = 0.67.
- **Impact:** no saved result (the YL runs used `none`). Future YL tilt runs would still need horizontal PPSD snapshots.
- **Correction:** commit `41a86f5`.
  - Horizontal codes are derived from the paired vertical (BHZ → BH1/BH2 unchanged for YV; LHZ → LH1/LH2).
  - N/E pairing is still not implemented.
  - `calibration_assumptions.md` was updated.
- **Verification:** 133 passed; Ruff passed.
- **Status:** Fixed on Sentinel branch.

### SENT-004: B01W has a frequency-dependent pressure response; the search bound is not the binding limit

- **Category / severity / confidence:** Scientific concern / High / Strong evidence.
- **Evidence:** the four boundary events (2010-02-27 M8.8, 06-12, 07-18, 08-13) were recomputed with the independent estimator:
  - Coherence is high (median γ² 0.97–0.99).
  - P/a phase is **+45° to +70° in every bin** (0.03–0.085 Hz). Same-event neighbours B04W and A02W are at 0 ± 10°.
  - B01W's R/T rises monotonically from 0.106 (0.03 Hz) to 0.231 (0.10 Hz). B04W is flat (1.28–1.47).
  - Averaged over Zha's 0.04–0.10 Hz band, R/T ≈ 0.19, which reproduces the published 0.19. Over our band it is ≈ 0.14.
  - The pattern is consistent with an extra high-pass of roughly 0.07 Hz corner. That is a hypothesis; it would fit, for example, a DPG leak time constant far shorter than nominal.
  - In a scalar-gain world, expanding the bound would yield g ≈ 7 with every bin failing phase.
  - Relaxing the timing guard adds one more coherent B01W event, which also fails phase.
- **Impact:** corrects `yl_all_zha_comparison.md` ("a software-policy limitation"). A scalar gain is not an adequate model for B01W. The published scalar is a band average of a strongly frequency-dependent ratio.
- **YV analogue, RR55** (the YV saved review's four "gain-grid boundary" events), independent estimator, no tilt, saved sign −1:
  - R/T = 0.24–0.27 (pressure about 4× low, so g≈4 lies beyond the bound).
  - Phase is a coherent +52° to +76°, falling with frequency. Control RR52 on the same code has R/T≈1.47 and phase ≈0°.
  - RR55's one positive-real metadata pole (+0.0431 rad/s) predicts only 26°→11° (2·atan(p/ω)), about a third of the observed lead.
  - RR55's inventory label is `Hydrophone`, which has a separate review flag.
  - Evidence: scratch output `rr55_frequency_dependence.txt`, reproducible with the method in `scripts/sentinel_yl_review.py`.
- **Correction:** none. Recommendations:
  - Do not widen the gain bound for B01W.
  - Treat it as a response-shape problem (frequency-dependent fit or flag).
- **Verification:** `docs/sentinel/b01w_boundary_events_spectra.csv`.
- **Status:** Open.

### SENT-005: the 1 µs timing tolerance excludes 515 YL event windows over µs offsets

- **Category / severity / confidence:** Scientific concern (policy) / Medium / Confirmed (counts).
- **Location:** `load_aligned` (`max(starts) - min(starts) > 1e-6`).
- **Evidence:**
  - 515 event windows across all 25 YL stations were excluded for between-channel start spreads of 2 µs to 200 µs (median 4 µs).
  - The largest implied phase error at 0.07 Hz is **0.005°**, against a 20° phase guard.
  - LDH and LHZ come from the same Q330, and 100 µs is the MiniSEED 2 BTIME quantum. These look like timestamp representation, not clock error.
- **Header evidence:** at A02W, 2010-04-06, the first samples are LDH `…00.666701`, LHZ `…00.666800` (the documented −0.000099 s). On 2010-05-05 both are `…00.595100`. An intermittent same-digitizer offset at the 100 µs quantum points to timestamp representation; Q330 per-channel timing documentation would settle it.
- **Prototype (scratch only):** relaxing the tolerance would not resolve A14W (0/13 pass), B01W (0/29) or C10W (1/12 pass, giving 2 < 3 events). Its effect on the other 22 stations was not computed; that is a network-wide calculation that needs your approval.
- **Correction:** none. A defensible option is a tolerance tied to the phase error at band maximum, e.g. 360·f_max·Δt ≤ 0.1° (Δt ≤ 4 ms at 0.07 Hz) and less than one sample. That is a user decision.
- **Status:** Open.

### SENT-006: the systematic −5.4% offset is within methodological sensitivity

- **Category / severity / confidence:** Scientific concern / Low / Strong evidence for the sensitivities; hypothesis for the cause.
- **Evidence:** recomputed from saved spectra (`docs/sentinel/zha_offset_sensitivity.csv`), median offset against published values:

  | Variant | Median offset |
  |---|---:|
  | As run | −5.4% |
  | 0.04–0.085 Hz band (dispersion table ends at 0.085 Hz) | −6.5% |
  | ρ = 1000 kg/m³ | −2.8% |
  | Scalar nominal sensitivities instead of full response | +1.8% |

  The offset shows no depth trend.
- **Interpretation:** the band does not explain the offset. Density and response treatment span it. Without Zha's processing details the offset cannot be attributed to hardware.
- **Status:** Open (needs the primary source's parameters).

### SENT-007: zero-power floor windows enter pressure diagnostics; RR12's "0% in range" is an artefact

- **Category / severity / confidence:** Scientific concern (missing guard) / Medium / Confirmed.
- **Location:** `validate_pressure` accepts any finite dB. Percentiles in `pressure_test`, `infragravity_test` and `nearby_test` therefore include −3076.5 dB floors.
- **Evidence:** YV.RR12.00.BDH has 562,742 floor values over 4,778 windows, and its saved `pressure_expected_median_in_range_fraction` is 0.0. RR50 also has floors (28,644 out-of-range values).
- **Impact:** existing diagnostics are `comparison_only` with blank scores, but RR12's outcome must not be read as "below expected noise".
- **Correction:** recommend counting and excluding windows at the floor before percentiles, recorded per channel. This changes window selection and needs user approval.
- **Status:** Open.

### SENT-008: documentation statements contradicted by SENT-001 and SENT-004

- **Category / severity / confidence:** Documentation / Medium / Confirmed.
- **Location:** `docs/yl_all_zha_comparison.md`, `CODEX_NEXT_CHAT.md`, README step 30.
- **Problem:**
  - C07W: "do not assume a simple gain step". The data show a two-state vertical amplitude, not a pressure step.
  - B01W: "software-policy limitation". The bound is not binding; phase is.
- **Correction:** a labelled "Sentinel follow-up" section was appended to `yl_all_zha_comparison.md`, with existing text and numbers unchanged. `CODEX_NEXT_CHAT.md`, README step 30 and `research_log.md` (lines 515–520) still carry the earlier reading; they belong to Codex's handoff.
- **Status:** Partly fixed on Sentinel branch.

### SENT-009: stale run commands could overwrite the completed YV Rayleigh workflow and master fields

- **Category / severity / confidence:** Documentation / Medium / Confirmed.
- **Location:** README step 10; `docs/yv_rayleigh_calibration.md` ("Run all available YV stations"); `network.py:57-67` against README and `network_workflow.md`; sampling descriptions in `rayleigh_paper_figure.md`, `teleseismic_dpg_calibration_report.md` §5 and `rayleigh_calibration_review.md`.
- **Problem:**
  1. The YV command is presented as "run or resume" against the default `results/NETWORKS/TEST`.
     - The reuse stamp hashes every package `.py` (`rayleigh_network.py:65`), and the code has changed since that batch.
     - So it would recompute all 43 stations under the new `all` sampling. `--rayleigh-max-events 8` is ignored without `--rayleigh-event-sampling bounded`.
     - It would replace the historical workflow entries and overwrite the master's `rayleigh_*` fields; master update is the default.
  2. `network all` also runs `paired` (seven raw days per pressure epoch, with master writes). README and `network_workflow.md` described eight stages.
  3. Three documents describe "stop at three candidates, eight attempts per level" as current behaviour. Since `a19bc3b` that requires `bounded`.
- **Correction:** warning and clarification text only, on the Sentinel branch. No code changed.
- **Status:** Fixed on Sentinel branch (documentation).

### SENT-010: acceleration is formed after tapering, a small estimator bias

- **Category / severity / confidence:** Scientific concern (estimator) / Low / Strong evidence.
- **Location:** `event_spectra` uses A = iω·FFT(w·v). The windowed transform of true acceleration is FFT(w·a) = iω·FFT(w·v) − FFT(w′·v).
- **Evidence:**
  - Synthetic broadband case: per-bin ratio errors up to 6.6%. They vanish for bin-centred sinusoids, which is why the new convention test uses those.
  - Real data, 48 YL candidate events at six stations: the production estimator was compared with FFT(w·a), using a from an untapered spectral derivative over a 1.5 h window.
  - Per-bin differences reach ±10%. Event-mean effect on g: median **+0.65%** (range −1.8% to +2.4%).
  - So 1/g reads about 0.6% low, the same direction as the −5.4% Zha offset but about a tenth of it.
- **Correction:** none. Differentiating, or response-correcting, over the full window before segmenting would remove it. That changes gains and is a method decision.
- **Status:** Open.

### SENT-011: missing tests for conventions and guards

- **Category / severity / confidence:** Missing test / Medium / Confirmed.
- **Correction:** five tests added in `3f663c5` (see "Test-suite independence"); all pass on unchanged code.
- **Status:** Fixed on Sentinel branch.

### SENT-012: the network CLI applies YV's provisional pressure sign to any network by default

- **Category / severity / confidence:** Scientific concern (default) / Low / Confirmed.
- **Location:** `network.py:666` (`--rayleigh-pressure-sign` default −1) and `rayleigh_network.py:150`. The standalone runner refuses to run outside RR52 without an explicit sign (`rayleigh_calibration.py:878-880`).
- **Impact:** a new network run without the flag uses YV's assumed −1. With signed metadata such as YL's, phase QC would then reject events (`polarity_review`); that is visible, but costs a batch.
- **Correction:** recommend requiring an explicit sign outside RR52 in the network CLI too. This is a CLI behaviour change, so it is your decision.
- **Status:** Open.

### SENT-013: the full inventory export overwrote an existing master without a backup

- **Category / severity / confidence:** Confirmed defect (data preservation) / High / Confirmed.
- **Location:** `save_master_dataset` and `python -m bruitfm.inventory`, `src/bruitfm/inventory.py`.
- **Problem:** the module CLI defaults `--output-dir` to the repository's `exports/`, the production master location. `save_master_dataset` wrote `master_inventory.csv/.xlsx` unconditionally, so one run replaced the master with a freshly built table. All processing results, manual scores and extra sheets were lost, with no backup (`exports/` is not in Git). `load_master_dataset`'s error message even suggests "rebuild it with inventory.py".
- **Evidence:** new test `test_full_export_refuses_to_replace_an_existing_master` (both files present; CSV only) failed at baseline: the scored master was replaced.
- **Correction:** commit `04350d5`. Both entry points raise `FileExistsError` when either file exists and name `network inventory` (which appends and keeps every cell). New exports into new folders, and `synchronize_inventory`'s first creation, are unchanged. README note added.
- **Verification:** 2 new cases; 140 passed; Ruff passed.
- **Status:** Fixed on Sentinel branch.

### SENT-014: concurrent batches could silently drop each other's master fields

- **Category / severity / confidence:** Confirmed defect (reliability) / Medium / Confirmed.
- **Location:** `_update_master_results` and `synchronize_inventory`, `src/bruitfm/inventory.py`.
- **Problem:** every master update reads both files, patches them and replaces them. Nothing serialized writers across processes; the only lock is per Rayleigh workflow folder, and `ensure_master_idle`/`assert_workflows_finished` are, as their docstrings say, not locks. Two batches sharing one master (YL and YV in two terminals, or a diagnostic stage beside a Rayleigh run) could both read, and the later commit erased the earlier one's fields. The docs do not warn against concurrent runs.
- **Evidence:** two real processes on the previous code, the first slowed inside its commit: the second process's committed `rayleigh_status` column was absent afterwards. Same scenario with the fix: the second process printed "Waiting for another master update to finish", and both `ppsd_status` and `rayleigh_status` survived. New test `test_master_update_waits_for_another_writer` holds the lock and checks the update waits (fails at baseline).
- **Correction:** commit `63ef818`. An exclusive `fcntl.flock` on `.master_inventory.lock` in the master folder is held from read to replace. A waiting writer prints a message, then proceeds. No field, selection or value changes.
- **Limits:** the lock is advisory and Unix-only (as is the existing Rayleigh lock); a manual edit in Excel is not covered.
- **Verification:** 141 passed; Ruff passed.
- **Status:** Fixed on Sentinel branch.

### SENT-015: per-flush master commits are slow and their backups grow without bound

- **Category / severity / confidence:** Performance and storage / Medium / Confirmed (measured).
- **Location:** `run_stage` in `network.py` (`flush()` after every 4 completed epochs) with `_commit_master_files`.
- **Evidence:** on a scratch copy of the production master (16,884 rows × 140 columns; 9.9 MB CSV, 2.0 MB Excel), one commit took 8.4 s (2.7 s to read the workbook) and wrote an 11 MB backup. The production `exports/backups/` already holds 168 backups (1.7 GB). A 4,000-epoch batch would spend about 2.3 h in commits and add about 11 GB of backups. Nothing prunes backups.
- **Impact:** batch wall time and disk use; no correctness effect. With SENT-014, a commit holds the master lock for those seconds.
- **Correction:** none applied; this trades crash-safety granularity against cost. Options for the user: flush on a time interval (e.g. every 10 minutes) instead of every 4 epochs; keep one pre-batch backup plus the latest; or compress backups. Deleting existing backups is the user's decision only.
- **Status:** Open.

### Minor observations (no action proposed)

- TauP distance uses ellipsoidal metres / 111.195 km per degree, an error below 0.5% of epicentral distance. It affects review-only P/S times.
- The cluster rule excludes events within 24 h of any *attempted* event, including attempts that later failed for data gaps. This is conservative.
- The overlap test considers only other events' fundamental-mode windows (±15 min), not their long coda or R2/R3 arrivals. This is already flagged as a heuristic.
- The RR52 layer model has no water layer. Per the sensitivity row above, its effect on T is small.
- `compare_yl_zha.py:149` rewrites the tracked `docs/yl_all_zha_comparison.{json,csv}` on every run, in addition to its new run folder.
- Pressure, infragravity and nearby outputs are written to fixed folders beside the NPZ, so a forced rerun replaces them. They are reproducible from the NPZ; Godin pilot folders are keyed by interval.
- The plots and percentiles of RR12/RR50 are dominated by −3076.5 dB floors (`values.min()` sets colour limits); see SENT-007.
- Godin's 2000 Pa²/Hz normalization and Webb's bounds are restated in tests, not independently sourced. Primary texts were not accessible here.

## 4. Commits on `review/sentinel`

Based on `faa39ff`; in order:

| Commit | Content |
|---|---|
| `2d4e8d6` | SENT-002 fix and regression test |
| `41a86f5` | SENT-003 fix, regression test and assumptions-table line |
| `ce41f54` | `scripts/sentinel_yl_review.py`, `docs/sentinel/*` evidence, first audit record |
| `3f663c5` | SENT-011: five independent convention/guard tests (138 pass) |
| `ddb5c35` | SENT-008/009 documentation, completed coverage, SENT-009–012 |
| `da8932a` | Upstream `4387fbe` and the clean trial merge recorded |
| `41cafa1` | Reference-data verification recorded |
| `e142d81` | `docs/sentinel/station_report_for_codex.md` |
| `04350d5` | SENT-013 fix, regression test, README note |
| `63ef818` | SENT-014 fix and regression test (141 pass) |
| (this commit) | Software-engineering pass, SENT-013–015, Socrates handoff |

**Integration caveat:** any change to `src/bruitfm/*.py` changes the code digests in network and Rayleigh fingerprints. The next network invocation therefore recomputes rather than reuses saved stages. Use the viewers for existing results.

## 5. Codex-ready tasks

1. **Integrate SENT-002/003/013/014.**
   - Cherry-pick `2d4e8d6`, `41a86f5`, `04350d5` and `63ef818` onto `main`.
   - Also cherry-pick `3f663c5` (tests only) and the documentation commits.
   - Acceptance: 141 tests pass, Ruff check/format pass, no saved result edited, README test count updated.
   - After integration, expect `.master_inventory.lock` (empty) beside the masters; it is harmless and should stay out of Git like the rest of `exports/`.
2. **Correct the C07W/B01W narrative (SENT-008).**
   - In `yl_all_zha_comparison.md` and the continuation note, state that:
     - C07W's vertical LHZ shows a two-state amplitude (≈0.44 vs ≈0.85 of C08W), and the 4 normal-state events give 1/g ≈ 1.90.
     - B01W's four boundary events fail phase in every bin, and R/T rises from 0.11 to 0.23 over 0.03–0.10 Hz.
   - Link `docs/SENTINEL_AUDIT.md`. The appended section in `yl_all_zha_comparison.md` already does this; update `CODEX_NEXT_CHAT.md` and README step 30 to match.
   - Acceptance: no numbers changed in the existing tables; the new statements cite `docs/sentinel/*`.
3. **Decisions for the user, not Codex:**
   - Timing tolerance (SENT-005).
   - Zero-power screening (SENT-007).
   - A vertical-consistency QC diagnostic (SENT-001).
   - Treatment of frequency-dependent responses (SENT-004, B01W and RR55).
   - Estimator ordering (SENT-010) and the network default sign (SENT-012).
   - Master flush cadence and backup retention (SENT-015).

## 6. Next targeted review

- A network-wide check for SENT-001-type vertical-state problems across YL and YV, using the daily-ratio method (read-only, but network-wide; needs approval).
- The SENT-005 timing-tolerance effect on the 22 resolved YL stations (network-wide prototype; needs approval).
- RR55 and other YV `polarity_review` stations (RR22, RR56): R/T and phase against frequency, as for B01W/RR55.
- Primary sources still unavailable here: Zha & Webb (2016) Appendix A parameters, Godin (2013) Eq. 6 normalization, and the Webb bounds.
