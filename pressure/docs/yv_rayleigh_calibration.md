# YV network Rayleigh calibration comparisons

The explicit `rayleigh` network stage applies the working method to available
pressure-channel epochs, with optional tilt correction, physical-unit complex
response correction, pressure/upward-acceleration spectra, coherence/phase QC
and analytically checked gain fitting. It uses the original figure's band and
layout. The current default evaluates all eligible events. The completed historical
YV batches used a bounded magnitude fallback seeking three candidate events;
those saved products remain unchanged.

## Comparison without tilt correction

Current sampling defaults to `--rayleigh-event-sampling all`: the minimum event
count governs availability of the station estimate and does not stop sampling.
`--rayleigh-max-events` applies only to explicit `bounded` mode. Old saved runs
retain their recorded limited sampling policy.

`--rayleigh-tilt-method none` uses pressure and the original response-corrected
vertical component. It skips horizontal-channel requirements, quiet training and
rotation. Existing event-window and spectral safeguards remain. The historical
default is still `tiskitpy`; select `none` explicitly for this comparison.

Run the available YV pressure epochs with a separate workflow and no master writes:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m bruitfm.network rayleigh \
  --network YV \
  --results-dir results/NETWORKS/YV_RAYLEIGH_NO_TILT \
  --rayleigh-no-master-update \
  --rayleigh-tilt-method none \
  --rayleigh-catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json \
  --rayleigh-dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign -1 --rayleigh-min-events 3 --rayleigh-max-events 8
```

View the new RR52 plot, rather than the old master-linked tilt-corrected plot:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network view \
  --network YV --station RR52 --view-kind rayleigh \
  --view-workflow results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json
```

Replace `RR52` with another YV station to view its saved comparison. Repeating the
same network command resumes/reuses results when the recorded dependencies match.
Output artifacts use new per-station `run_*` folders beside the saved pressure
PPSD; `--results-dir` selects the separate workflow. Existing workflow results,
masters, scores, responses, raw records and PPSDs are retained.

For a first single-station comparison that shows the plot immediately and never
writes the master:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.rayleigh_calibration \
  --network YV --station RR52 --tilt-method none \
  --catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json \
  --dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --frequency-min 0.03 --frequency-max 0.07 --pressure-sign -1 \
  --min-events 3 --max-events 8 --show
```

To restrict that first check to Alaska, add `--event-id ak0138esnzr --max-events 1`.
One event can produce an event-level comparison but cannot establish a station
gain under the existing three-event rule. No real-data run without tilt correction
has yet been performed as part of this code change.

## Where the earlier workflow stands

Before this stage, YV had 43 stations / 172 completed channel PPSDs and base plots.
All 43 pressure channels had expected-pressure, nearby and Godin comparisons.
Nearby outcomes were four exploratory comparisons, 36 without suitable neighbors,
two lacking reference stations and one lacking common days. Saved Godin had 19
comparisons and 24 without eligible bins. Raw Godin pilots had 39 comparisons,
two without eligible bins and two processing/coverage failures (RR12/RR42).
The continuous pressure–vertical stage had one RR52 interior-week comparison.
Those results and their original statuses are preserved.

## Run all available YV stations

> **Historical command; the YV batch below is complete.** Reuse fingerprints
> include the package code, so after code changes this command recomputes every
> station with the current all-event sampling (`--rayleigh-max-events` needs
> `--rayleigh-event-sampling bounded`), replacing the saved workflow entries and
> the master's `rayleigh_*` fields. View saved results instead; run new
> comparisons with a new `--results-dir` and `--rayleigh-no-master-update`.

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m bruitfm.network rayleigh \
  --network YV \
  --rayleigh-catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json \
  --rayleigh-dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign -1 --rayleigh-min-events 3 --rayleigh-max-events 8
```

Use `--station RR52` to select one station, or repeat `--station` for a subset.
The explicit stage does not rerun PPSDs or existing diagnostics. It is kept out
of `network all`, because it requires declared catalogue/model inputs and
provisional scientific assumptions. Other networks need appropriate inputs and
metadata; this does not establish validation of the YV assumptions elsewhere.

The stage processes stations serially, records each completed epoch immediately,
and appends only `rayleigh_*` master fields in one parent write with CSV/Excel
backups. No scientific score, sensitivity, source waveform or response is changed.
Existing extra Excel sheets, metadata and manual results are retained.

Rerun the same command to resume an interrupted run with unchanged code and
options. Reuse checks source file
size/mtime, all partner epochs/statuses/responses, depth/sensor information,
catalogue/model/review inputs and scientific code. A missing saved event product
also invalidates reuse. `--force` makes new computations in separate run folders.
No previous station run folder is overwritten. A file lock prevents competing
Rayleigh runs; an unfinished other workflow blocks the master-writing batch.

## Files and statuses

- Live checkpoint: `results/NETWORKS/TEST/_workflow/rayleigh.json`.
- Station results: the pressure PPSD epoch's `rayleigh_calibration/run_*/` folder,
  including `comparison.json`, event spectra/CSV, arrival/phase/waveform plots,
  reviewed-window template and the combined `paper_figure5_comparison.png`.
- Metadata/reference exceptions have a saved audit under
  `_workflow/rayleigh_unresolved/`, without a fabricated gain.
- Master fields link the station audit, Figure 5 plot, candidate count, candidate
  pressure/signed/sensitivity multipliers and review flags. Gains stay blank when
  fewer than the selected number of events pass. A diagnostic aggregate-curve
  fit can still exist with fewer events; it is not a station-level calibration.

`candidate_comparison_under_assumed_sign` means enough events pass the automatic
pilot rules; physical/window/uncertainty flags remain. `unresolved` means the
run lacks enough candidates. `unresolved_metadata` identifies missing/ambiguous
channels or responses. `reference_band_unusable` identifies an unsuitable water
reference for the requested band/depth. `failed` is a software/input exception;
none of these statuses declares instrument failure or changes scores.

All results record the [assumptions and remaining checks](calibration_assumptions.md).
Using the shared historical RR52 model and assumed network pressure sign −1 is
explicitly provisional. Depth comes from each station's inventory, not RR52's
depth. Generic hydrophone hardware labels receive a separate review flag.

## View without recomputing

The agent's batch was initially stopped at the user's request after **7/43**.
The user then completed the network run: **43/43**, with **16 candidate
comparisons, 25 unresolved and two unsuitable reference bands (RR41/RR57)**.
The final master update ran; all pre-existing CSV columns match the pre-run
snapshot exactly. Scores and input responses/data were preserved. Candidate
gains are still provisional and no correction was applied.
[Completed audit](yv_rayleigh_completed_summary.json) ·
[Historical partial audit](yv_rayleigh_validation_summary.json) ·
[Reviewed RR03 example](figures/yv/RR03_rayleigh_calibration.png).

The [initial saved-results review](yv_rayleigh_saved_results_review.md) checks
gain arithmetic, candidate plots, event rejection reasons and depth sensitivity.
It also identifies RR57's conservative reference preflight flag. Full packet,
response and uncertainty validation remains pending; no correction was applied.

The [subsequent raw/window validation](yv_rayleigh_independent_validation.md)
adds seven independent raw spectral checks, all 48 candidate arrival reviews and
quiet/model sensitivity. Future calibration reads reject fractional native
record boundaries before merging. Preserve yesterday's run as historical evidence;
do not rerun it simply to replace its statuses with the new guard's decisions.

The master links are now available:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network view \
  --network YV --station RR52 --view-kind rayleigh
```

Repeat `--station`, for example `--station RR03 --station RR52`, to view a subset.
The viewer requires a station or SEED-ID selection to avoid opening the entire
network's plot windows accidentally.
Per-event plots remain in each linked station folder; an unresolved metadata
audit may have no plot. Full data/products stay local. Code, documentation,
curated plots and a compact network summary are tracked in Git.

To view an example immediately, without finishing or recomputing the network:

```bash
open docs/figures/yv/RR03_rayleigh_calibration.png
```

For numerical results after completion, open `exports/master_inventory.xlsx`
and inspect the pressure rows' `rayleigh_*` columns. During a run, the live
checkpoint contains each completed station's update, result path and plot path;
master links are appended only at batch completion. Ctrl-C retains completed
checkpoints and releases the run lock; rerun the same command to continue.
