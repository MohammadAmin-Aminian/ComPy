#!/usr/bin/env bash
# Run after downloading YL; stop on errors and preserve the five-station products.
set -eo pipefail
cd "$(dirname "$0")/.."

BRUIT_PY=(conda run --no-capture-output -n bruitfm python)
"${BRUIT_PY[@]}" scripts/prepare_yl_test.py --all-stations
YL_CATALOG=$(< results/RAYLEIGH_REVIEW/yl_all_catalog_path.txt)

YL_ALL_ARGS=()
while IFS= read -r seed; do
  YL_ALL_ARGS+=(--seed-id "$seed")
done < results/RAYLEIGH_REVIEW/yl_all_pressure_vertical_seeds.txt
if (( ${#YL_ALL_ARGS[@]} == 0 )); then
  exit 1
fi

"${BRUIT_PY[@]}" -m compy_pressure.network check \
  --network YL "${YL_ALL_ARGS[@]}" \
  --results-dir results/NETWORKS/YL_ALL_NO_TILT

YL_PPSD_ARGS=()
while IFS= read -r seed; do
  YL_PPSD_ARGS+=(--seed-id "$seed")
done < results/RAYLEIGH_REVIEW/yl_all_ppsd_needed_seeds.txt
if (( ${#YL_PPSD_ARGS[@]} )); then
  "${BRUIT_PY[@]}" -m compy_pressure.network ppsd \
    --network YL "${YL_PPSD_ARGS[@]}" --workers 4 \
    --results-dir results/NETWORKS/YL_ALL_PPSD
fi

"${BRUIT_PY[@]}" -m compy_pressure.network rayleigh \
  --network YL "${YL_ALL_ARGS[@]}" \
  --results-dir results/NETWORKS/YL_ALL_NO_TILT \
  --rayleigh-no-master-update --rayleigh-tilt-method none \
  --rayleigh-event-sampling all \
  --rayleigh-catalog "$YL_CATALOG" \
  --rayleigh-dispersion-file ../compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign 1
