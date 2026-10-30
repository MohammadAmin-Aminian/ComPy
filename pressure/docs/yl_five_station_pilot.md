# YL five-station pilot without tilt correction

Latest user correction: **use all eligible earthquakes**, not the first three
successes. The current default is `--rayleigh-event-sampling all`; it does not
stop at the existing three-event station-estimate minimum or the old eight-attempt
budget. It exhausts the eligible catalogue under the existing magnitude/depth,
epoch, 24-hour separation and event QC rules. The old exploratory behavior is
available only through explicit `--rayleigh-event-sampling bounded`. Earlier
saved YL products used the old limited policy; they have not been rerun.

Scope: A02W, A06W, A08W, A10W and A12W. Original XMLs are under
`Data/YL/<station>/inventory.xml`; waveforms use downloaded daily MiniSEED names
with `--` for blank location. Discovery now accepts those names and validates
their SEED IDs against actual MiniSEED headers, without renaming/copying records.
These stations have LDH/LHZ at 1 Hz and deployment metadata spanning late 2009
through late 2010. The original signed pressure response has negative sensitivity;
the initial test uses `+1` as the additional pressure multiplier, preserving the
metadata sign rather than importing the YV `-1` assumption. Hardware confirmation
remains pending.

The preparation script registers only the five stations, appending metadata
while preserving existing master fields, and fetches/caches a deployment-specific
USGS magnitude >=5.5 catalogue. The initial snapshot contains 566 events; its
date limits come from these inventories. Existing catalogue files are retained.
It does not start PPSD or calibration computations.

The first user run registered all five stations and completed their vertical
PPSDs. All five pressure jobs failed on BRUIT-FM's positive-only sensitivity
guard; the following calibration stage therefore recorded five
`unresolved_metadata` entries. This software restriction is now corrected:
finite nonzero signed sensitivities are accepted without changing the original
XML, stage gains or complex response. Zero/nonfinite values remain rejected.
Synthetic tests show PSD invariance and complex-response sign preservation.
A separate six-hour A02W pressure pilot completed with three segments and
unchanged sensitivity −753.235; its raw/XML/master hashes were preserved.
Pilot output: `results/VALIDATION/YL_signed_response/run_1bbfc0ec3565/`.

To resume after that failure, retry only pressure channels (omit `--force` and
`--rebuild`), then rerun the separate calibration command below:

```bash
YL_PRESSURE=()
for station in A02W A06W A08W A10W A12W; do
  YL_PRESSURE+=(--seed-id "YL.${station}..LDH")
done
conda run --no-capture-output -n bruitfm python -m bruitfm.network ppsd \
  --network YL "${YL_PRESSURE[@]}" --workers 4 \
  --results-dir results/NETWORKS/YL_FIVE_NO_TILT
```

The completed vertical PPSDs do not need to be recalculated. The original timing
guards still apply after the signed-response blocker is removed.

The initial Rayleigh comparison explicitly borrows the existing RR52 dispersion
table as a provisional model, with a 0.03–0.07 Hz band. It is a pipeline test;
it is not an exact reproduction of Zha & Webb (2016), whose Appendix A uses a
local average Rayleigh model and 0.04–0.1 Hz. Their Table 1 includes all five
stations. Model, event policy and gain convention must be reconciled before
comparing numbers directly with the paper.
[Original paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688).

The [Earth-model option](rayleigh_earth_model_options.md) can now generate curves
from a supplied layer CSV using disba. Verified numeric Zha/Lau inputs are still
needed before selecting that model. The same note compares gain conventions:
the original pilot's reciprocal coefficients are about 5–6% below the published
ones, rather than differing by a factor of two.

Header-only checks on 2010-04-06 found pressure–vertical native-grid residuals:
A02W −0.000099 s, A06W −0.000002 s, A08W +0.000002 s, A10W −0.000002 s,
A12W −0.0001 s. These exceed the current 1e-6-second cross-channel guard.
They do not by themselves prove meaningful physical phase error, but the initial
run can reject such event windows for timing review. No timestamp rounding,
interpolation or tolerance relaxation is introduced by this preparation.

Run these in order from the repository root, stopping if a command fails:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python scripts/prepare_yl_test.py

YL_STATIONS=(--station A02W --station A06W --station A08W --station A10W --station A12W)
YL_SEEDS=()
for station in A02W A06W A08W A10W A12W; do
  YL_SEEDS+=(--seed-id "YL.${station}..LDH" --seed-id "YL.${station}..LHZ")
done

conda run --no-capture-output -n bruitfm python -m bruitfm.network check \
  --network YL "${YL_STATIONS[@]}" "${YL_SEEDS[@]}" \
  --results-dir results/NETWORKS/YL_FIVE_NO_TILT

conda run --no-capture-output -n bruitfm python -m bruitfm.network ppsd \
  --network YL "${YL_STATIONS[@]}" "${YL_SEEDS[@]}" --workers 4 \
  --results-dir results/NETWORKS/YL_FIVE_NO_TILT

conda run --no-capture-output -n bruitfm python -m bruitfm.network rayleigh \
  --network YL "${YL_STATIONS[@]}" "${YL_SEEDS[@]}" \
  --results-dir results/NETWORKS/YL_FIVE_NO_TILT \
  --rayleigh-no-master-update --rayleigh-tilt-method none \
  --rayleigh-catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_YL_2009_2010.json \
  --rayleigh-dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign 1 --rayleigh-event-sampling all

conda run --no-capture-output -n bruitfm python -m bruitfm.network view \
  --network YL --station A02W --view-kind rayleigh \
  --view-workflow results/NETWORKS/YL_FIVE_NO_TILT/_workflow/rayleigh.json
```

This computes PPSDs for ten pressure/vertical channel epochs; horizontal PPSDs
and other diagnostics are unnecessary for this test. Registration/PPSD append
their metadata/status fields. Calibration does not update the master or apply
gains. Change A02W in the viewer to another selected station. If a station has
only timing rejections, its combined plot can have no fitted curve; read its
comparison JSON before concluding that the pressure instrument failed.
