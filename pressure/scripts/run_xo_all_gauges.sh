#!/usr/bin/env bash
# All downloaded XO stations: prefer APG, otherwise DPG; no hydrophone fallback.
set -eo pipefail
cd "$(dirname "$0")/.."
BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
RUN_ROOT="${BRUIT_XO_RUN_ROOT:-results/NETWORKS/XO_ALL_GAUGES_NO_TILT}"
"${BRUIT_PY[@]}" scripts/prepare_xo_run.py \
  --pressure-kind gauges --event-resample-rate 1 --native-rate-tolerance-ppm 0.1 \
  --fetch-responses --exclude-incompatible-files --results-dir "$RUN_ROOT"
ARGS=(--network XO --data-dir "$RUN_ROOT/_inputs/waveforms" \
  --master-dir "$RUN_ROOT/_master" --results-dir "$RUN_ROOT")
"${BRUIT_PY[@]}" -m compy_pressure.network rayleigh "${ARGS[@]}" \
  --rayleigh-catalog "$RUN_ROOT/_inputs/catalog_usgs_m5p5.json" \
  --rayleigh-dispersion-file ../compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign 1 --rayleigh-tilt-method none \
  --rayleigh-event-resample-rate 1 --rayleigh-native-rate-tolerance-ppm 0.1 \
  --rayleigh-event-sampling all --rayleigh-no-master-update
STATIONS=()
while IFS= read -r station; do
  STATIONS+=(--station "$station")
done < <("${BRUIT_PY[@]}" -c 'import json,sys; print("\n".join(s["station"] for s in json.load(open(sys.argv[1]))["stations"]))' "$RUN_ROOT/preparation.json")
"${BRUIT_PY[@]}" -m compy_pressure.network view "${ARGS[@]}" "${STATIONS[@]}" \
  --view-kind rayleigh --view-workflow "$RUN_ROOT/_workflow/rayleigh.json"
