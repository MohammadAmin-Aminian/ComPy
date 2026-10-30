# BRUIT-FM project status, 22 October 2026

This is the current continuation guide after Claude took over from Codex on 21
October. It summarises what was done, with the numbers, what the owner still
has to decide, and how to start on the new PPSD computer. Detailed evidence is
in the linked reports. Nothing here adopts a calibration, threshold, rating or
response correction.

## Rules followed throughout (unchanged)

* No gain applied to raw data, responses, saved PPSDs or masters. No XML repair.
* Frozen results kept. **RR52 stays g = 0.6787278140712512.**
* New work goes into new folders, and every input is hashed before and after.
* Confirmed facts are kept separate from hypotheses. Inspected events count as
  development, not validation.
* Scientific changes, outside contacts and long campaigns need the owner's
  approval. The owner's untracked `gain_histogram.png` was never touched.

## What was merged (main `90adc26`)

| PR | Content |
|---|---|
| [#2](https://github.com/MohammadAmin-Aminian/bruit-fm/pull/2) | Official response confirmation for RR17, RR18, RR22, RR56; saved-spectrum response pilot |
| [#3](https://github.com/MohammadAmin-Aminian/bruit-fm/pull/3) | Plots: network median PSD with limits; pooled PPSD by sensor type; four-panel station PPSD |
| [#4](https://github.com/MohammadAmin-Aminian/bruit-fm/pull/4) | YV all-methods runner (reproducible, code-safe) and the 21–22 October run report |
| [#5](https://github.com/MohammadAmin-Aminian/bruit-fm/pull/5) | Rayleigh code hash recorded at import, not at write time (provenance only) |

CI (Ruff, full test suite and wheel smoke test) passed on every PR. The last
local full suite, at #5, gave 469 passed.

## Results, with numbers

### 1. Official response metadata: RR17, RR18, RR22, RR56 ([report](yv_response_confirmation.md))

* Current RESIF StationXML is **identical to the saved responses in all 8 P/Z
  channels**: epochs, gains, A0, poles/zeros, decimation and FIR. Only the
  letter case of the unit names differs.
* RESIF is the **only data centre**, so archive agreement is not independent
  confirmation. The metadata also carries no sensor model or serial number.
* **RR56:** the recorded vertical A0 matches the conjugate-paired roots
  (+0.097%), not the stored duplicated roots (−3.41%). This supports H2.
* **RR22:** A0 is closer to the paired roots (+0.29%) than to the stored ones
  (+1.51%). Weaker support.
* **RR17/RR18:** A0 is consistent with the stored 0.747/0.766 rad/s poles to
  within 3×10⁻⁶. **H4 stays post hoc.**
* **Stähler et al. (2016)** prints S_ADC = 3.62×10⁵ counts/V. That disagrees
  with its own printed overall sensitivity and with the XML's 373380 counts/V,
  a 3.1% difference.

### 2. Saved-spectrum response pilot (development only)

| Station | Hypothesis | g | Events | Phase RMS before → after |
|---|---|---:|---:|---|
| RR56 | H2: restore pole pairs | 1.166 | 4 | 115° → 9° |
| RR18 | H4: pole ÷10 (post hoc) | 1.221 | 4 | 53° → 4° |
| RR17 | H4: pole ÷10 (post hoc) | 0.989 | 10 | 56° → 4° |
| RR22 | H2: restore pole pairs | 1.419 | 5 | 110° → 21°; 2 of 7 still fail |

* All controls passed: RR52 exactly 0.6787278140712512; RR16 identity; RR10
  negative control (phase 2.7° → 17.3°).
* Caveat from the negative control: RR10 still passes the 20° screen under the
  wrong H4, and its g moves by 6%. Passing the phase screen alone does not
  validate a pole.

### 3. YV all-methods calibration on the latest package ([report](yv_calibration_run_20261021.md))

* **Rayleigh, clean reproduction (bounded events, 22 October): RR52 =
  0.6787278140712512 exactly.**
  * **43 of 43 stations identical** to the frozen workflow: 16 estimates,
    25 unresolved, 2 reference-band exclusions.
  * The only difference is labels: at RR55, 4 non-candidate events changed
    their rejection label.
* **Rayleigh, all events (21 October):** the same 16 estimates.
  * Accepted events: 48 → 442.
  * Median gain change **+0.01%**, range −3.74% to +4.02%; 12 of 16 within ±2%.
  * Network median: 0.9494 → 0.9330.
  * RR52: 0.6750 from 40 events (−0.54%).
  * More events rescued none of the 25 unresolved stations.
* **Microseism:** 41 checked; 2 unresolved pairing (RR12, RR42).
  * **43 of 43 identical** to the 14 October run.
  * The frozen Rayleigh g lies inside the best microseism bracket at 7 of 16
    stations; the pole-trial g at 14 of 31.
* **Local P:** 26 no suitable events, 9 unresolved, 8 diagnostic only; no
  station value.
  * 41 of 43 pairs identical to the frozen 16 October campaign.
  * Two event flips (RR11 `us2013pjc9`, RR50 `usb000ja9l`) also occur on clean
    code. They come from reviewed changes after 16 October, not from the
    owner's edit.

### 4. RR52 across methods

| Method | g | Basis |
|---|---:|---|
| Rayleigh, frozen and reproduced | **0.6787278140712512** | 3 events |
| Rayleigh, all events | 0.6750 | 40 events |
| Rayleigh, held-out | +0.66% median prediction error | 14 events |
| Microseism | 0.671–0.681 (best sub-band) | 167 windows; coh² 0.986 |
| Local P | 0.56 / 0.44 (4–8 / 8–16 Hz) | 1 of 22 events; diagnostic |
| Published (Aminian 2025) | 0.67 | — |

### 5. Why 25 YV stations have no Rayleigh estimate

* **Wrong pole sign in the metadata:** 15 stations recover with the pole
  reflection trial.
* **Metadata defects:** 4 stations (RR17, RR18, RR22, RR56) recover in the
  pilot, as development values only.
* **Hardware or data:** RR12 (flat −3077 dB archive), RR31 (DPG artefacts),
  RR42 (failed vertical), RR45 (noisy vertical), RR53 (unexplained).
* **Polarity:** RR47.
* Coherence (coh² ≥ 0.95) removes most events at every station.
* RR41 and RR57 are reference-band exclusions, counted separately.

### 6. YV pressure PPSD by sensor type (saved PPSDs, no gain)

| Sensor type | Channels | Windows | Source of grouping |
|---|---:|---:|---|
| INSU differential pressure gauge (DPG) | 9 | 40,973 | StationXML sensor field |
| Geomar HTI-01-PCA hydrophone | 4 | 13,671 | Stähler et al. (2016), App. A2 |
| DEPAS HTI hydrophone | 30 | 126,931 | All others; HTI-01 vs HTI-04 is not resolvable per station |
| **Total** | **43** | **181,575** | |

* RR12's saved median is −3077 dB at every frequency. It is shown as off scale.
* One DPG has a distinctly low curve, about 0 to −50 dB, under the main band.
* 14 YV pressure channels have no saved PPSD yet.

## Owner decisions still open

1. **The uncommitted edit removing the fractional native-grid guard.** It is
   saved as `~/Projects/bruit-fm-uncommitted-rayleigh-guard-removal-20261021.patch`
   and is not in the checkout.
   * Recommendation: **discard it** unless there is a reviewed reason. It
     reverses a reviewed safety fix, and the clean rerun showed no YV result
     depends on it.
2. **Rayleigh reference policy.**
   * Recommendation: **keep the frozen bounded values as the reference**. They
     reproduce exactly, and the held-out validation used them.
   * Report the all-event values as robustness evidence, within ±4%.
   * Revisit at the ratings stage.
3. **Outside requests (unsent drafts in the response report).** DEPAS/AWI for
   RR17/RR18 R, C_hyd and S_d plus the GEOLON ADC gain; RESIF for the unpaired
   roots.
   * Recommendation: send the DEPAS request. Only DEPAS can settle RR17/RR18
     and the 3% ADC question.
4. **Not adopted, waiting:** pilot gains for RR17, RR18, RR22 and RR56;
   polarity at RR47; hardware identity at RR53 and RR55.

## Runbook for the new PPSD computer

Saved results record **absolute paths, file sizes and modification times**.
The coverage planner reuses a saved PPSD only when every source waveform still
has the same path, size and `mtime_ns`. A changed time blocks reuse, and a
changed path makes the files look uncovered, so the planner would recompute
almost everything. The move must preserve all three.

0. **Old computer, before anything is deleted.** Baseline, read-only, about
   3 s:

   ```bash
   cd ~/Projects/bruit-fm
   /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python scripts/check_migration.py
   ```

   On 22 October this gave **READY**: 264 saved PPSDs (YV 172, YL 50, 7D 30,
   XO 12), **65,572 source waveforms unchanged**, all frozen hashes and RR52
   matching, and all 574 paths in the frozen review present.
1. **Same user and same folders.** Keep the account name `mohammadamin` and
   copy to identical paths:
   * `~/Desktop/Data_Bruit_FM/` (raw `Data/`, `INVENTORIES/`, `RR52/`);
   * `~/Projects/bruit-fm/`, including `exports/`, `results/` and `.git`;
   * the sibling folders that the frozen review points into:
     `~/Projects/bruit-fm-claude/`, `~/Projects/bruit-fm-localp-review/` and
     `~/Projects/bruit-fm-localp-finish/`. Simplest: copy all of `~/Projects/`.
2. **Copy, preserving timestamps.** Use Migration Assistant, or `rsync -a`;
   `-a` keeps modification times. Never overwrite existing files.
   * **Do not re-download the waveforms** that saved products were made from.
     New downloads have new times and would force recomputation.
   * Extra data, such as the missing YV October–December 2012 days or a fresh
     RR12, goes into a **separate folder** and is treated as new input later.
3. **Environment.** The shell wrappers look for a conda environment named
   `bruitfm`, or the Homebrew Miniforge path, so recreate it with the same name:

   ```bash
   brew install miniforge
   conda create -n bruitfm python=3.13 && conda activate bruitfm
   cd ~/Projects/bruit-fm && pip install -e '.[dev,rayleigh,dispersion]'
   python -m pytest -q                      # expect all tests to pass
   conda info --base                        # expect /opt/homebrew/Caskroom/miniforge/base
   ```

   The coverage and YV wrappers default to
   `/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python`, which is
   exactly where Homebrew's Miniforge puts this environment on an Apple-silicon
   Mac. If `conda info --base` differs, run
   `export BRUITFM_PYTHON="$(conda run -n bruitfm which python)"` before the
   wrappers.

4. **Check the copy** before planning anything. It must say **READY** with
   the same counts as step 0:

   ```bash
   python scripts/check_migration.py
   ```

   If it reports `modified` or `missing` source files, fix the copy (paths or
   timestamps). Do not plan or execute until it passes.
5. **Plan only**, which writes just a plan folder and a header cache:

   ```bash
   bash scripts/run_full_coverage.sh plan
   bash scripts/run_full_coverage.sh dry-run
   ```

   Compare with the 19 October plan: 89 reused pressure PPSDs, 15 new pressure
   PPSDs and 61 × 3 spectral checks (198 tasks), with 19 blocked epochs. Explain
   any difference before executing.
6. **Execute only missing work** after reviewing the dry run. It is resumable;
   rerun the same command after an interruption.

   ```bash
   bash scripts/run_full_coverage.sh execute --workers 4
   bash scripts/run_full_coverage.sh report --curated-dir <new folder>
   ```

7. **Confirm the YV calibration still reproduces on the new machine.** The
   runner refuses uncommitted code and stops unless RR52 is exact.

   ```bash
   bash scripts/run_yv_all_methods.sh --check
   bash scripts/run_yv_all_methods.sh --only rayleigh
   ```

8. **Then, in order:**
   * remake the sensor-type and median PPSD plots from the complete data;
   * decide thresholds and ratings from the distributions;
   * validate the pilot gains on untouched M ≥ 6 events with RR52 as the
     control;
   * only after that, build the scoped noise model.

## Plot commands (from `~/Projects/bruit-fm`)

```bash
PY=/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python
PYTHONPATH=src:. $PY scripts/plot_network_median_ppsd.py --network YV --channel BDH --show
PYTHONPATH=src:. $PY -m bruitfm.plot_network_ppsd --network YV --channels BDH --probability --show
PYTHONPATH=src:. $PY scripts/plot_yv_ppsd_by_sensor.py --output-dir <new folder>
PYTHONPATH=src:. $PY scripts/plot_station_ppsd_panels.py --station RR52 --out <new file>.png
```

Every plot command refuses to overwrite an existing output.
