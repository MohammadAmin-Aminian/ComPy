#!/usr/bin/env bash
# First 7D calibration test; all writes stay in a separate result/master root.
set -eo pipefail
cd "$(dirname "$0")/.."
PRESSURE_PY=(conda run --no-capture-output -n compy python)
RUN_ROOT=results/NETWORKS/7D_J67A_NO_TILT
"${PRESSURE_PY[@]}" scripts/prepare_7d_run.py \
  --station J67A --low-rate-only --fetch-responses --results-dir "$RUN_ROOT"
ARGS=(--network 7D --station J67A \
  --data-dir "$RUN_ROOT/_inputs/waveforms" \
  --master-dir "$RUN_ROOT/_master" --results-dir "$RUN_ROOT")
"${PRESSURE_PY[@]}" -m compy_pressure.network ppsd "${ARGS[@]}" --workers 2
"${PRESSURE_PY[@]}" -m compy_pressure.network rayleigh "${ARGS[@]}" \
  --rayleigh-catalog "$RUN_ROOT/_inputs/catalog_usgs_m5p5.json" \
  --rayleigh-dispersion-file ../compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign 1 --rayleigh-tilt-method none \
  --rayleigh-event-sampling all --rayleigh-no-master-update
"${PRESSURE_PY[@]}" -m compy_pressure.network view "${ARGS[@]}" \
  --view-kind rayleigh --view-workflow "$RUN_ROOT/_workflow/rayleigh.json"
