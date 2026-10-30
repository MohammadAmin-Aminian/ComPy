# All available YL stations, without tilt correction

Wait for the waveform download to finish before preparing/running this sequence.
The current XML coverage includes 25 station folders and 50 LDH/LHZ channels,
with deployment bounds 2009-11-20 through 2010-12-03. Preparation has cached
568 USGS events in a new all-station catalogue, retaining the five-station
catalogue and original waveforms/XMLs.

`--all-stations` discovers original station XMLs and appends their metadata to
the existing CSV/Excel master. It prepares a complete P/Z seed list plus a PPSD
execution selector. Existing complete PPSDs are skipped only when their saved
inventory hash and raw file paths/sizes/mtimes still match; newly downloaded
inputs must be processed. The execution selector is regenerated from live state.
Reuse also requires the full available interval: a complete short pilot is not
treated as a full-deployment product. Headers establish the available end, clipped
to the response epoch. The completed 50 YL P/Z products passed this check in the
[correctness audit](project_audit.md).
New/changed PPSD products use a fresh results root, preserving the original pilot
products. Existing manual master fields and scores remain.

The complete sequence is packaged as a runner that stops on command errors:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
bash scripts/run_yl_all.sh
```

It does not automatically open plot windows. Use the viewer after completion.

Run in order, stopping if a command fails:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python scripts/prepare_yl_test.py --all-stations

YL_ALL_ARGS=()
while IFS= read -r seed; do
  YL_ALL_ARGS+=(--seed-id "$seed")
done < results/RAYLEIGH_REVIEW/yl_all_pressure_vertical_seeds.txt

conda run --no-capture-output -n bruitfm python -m bruitfm.network check \
  --network YL "${YL_ALL_ARGS[@]}" \
  --results-dir results/NETWORKS/YL_ALL_NO_TILT

YL_PPSD_ARGS=()
while IFS= read -r seed; do
  YL_PPSD_ARGS+=(--seed-id "$seed")
done < results/RAYLEIGH_REVIEW/yl_all_ppsd_needed_seeds.txt

if (( ${#YL_PPSD_ARGS[@]} )); then
  conda run --no-capture-output -n bruitfm python -m bruitfm.network ppsd \
    --network YL "${YL_PPSD_ARGS[@]}" --workers 4 \
    --results-dir results/NETWORKS/YL_ALL_PPSD
fi

conda run --no-capture-output -n bruitfm python -m bruitfm.network rayleigh \
  --network YL "${YL_ALL_ARGS[@]}" \
  --results-dir results/NETWORKS/YL_ALL_NO_TILT \
  --rayleigh-no-master-update --rayleigh-tilt-method none \
  --rayleigh-event-sampling all \
  --rayleigh-catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_YL_ALL_20091120_20101204.json \
  --rayleigh-dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign 1

conda run --no-capture-output -n bruitfm python -m bruitfm.network view \
  --network YL --station A02W --view-kind rayleigh \
  --view-workflow results/NETWORKS/YL_ALL_NO_TILT/_workflow/rayleigh.json
```

Change A02W in the viewer to another available station. Horizontal channels and
the other spectral diagnostic stages are not needed for this comparison.
The sampling minimum remains separate from the number attempted: all eligible
events are evaluated under the current magnitude/depth, epoch, 24-hour separation,
timing and QC rules. No gain is applied and calibration does not update the master.

## What is aligned with Zha, and what is still different?

The physical pressure/upward-acceleration equation is based on Zha & Webb.
The code supports all-event sampling and optional disba layer-model generation,
with hashes and solver provenance. Signed input responses remain unchanged.

This command deliberately retains the existing provisional model and 0.03–0.07 Hz
band for comparison with the pilot. It is not an exact published-method
reproduction. Verified numeric Lau model inputs, the paper's 0.04–0.1 Hz band,
specific large-event selection, Rayleigh windows, coherence convention and
averaging/fitting details still need reconciliation. Existing timing and phase
guards are not silently relaxed. See [model options and comparison](rayleigh_earth_model_options.md).
