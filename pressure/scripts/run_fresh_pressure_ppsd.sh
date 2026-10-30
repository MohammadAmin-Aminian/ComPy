#!/usr/bin/env bash
# Fresh pressure-PPSD run from a data drive into a new Desktop folder.
#   bash scripts/run_fresh_pressure_ppsd.sh /Volumes/<drive>/Data_Bruit_FM [stage] [options]
# The data folder must contain INVENTORIES/ (StationXML) and Data/ or SDS/ (waveforms).
# Stages: inventory, check, ppsd, plots, or all (default: all, in that order).
# Extra options go to python -m compy_pressure.network (e.g. --network YV --workers 8).
# Raw data are only read. Outputs, by default under ~/Desktop/BRUIT-FM/fresh_run_<date>:
#   master/   master_inventory.csv and .xlsx
#   results/  <network>/<station>/<epoch>/ppsd.npz ... and <station>/plots/
#   logs/     one log per stage
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run with bash; do not source this script.\n' >&2
  return 2
fi
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="${BRUITFM_PYTHON:-/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python}"
SOURCE="${1:-}"
WAVEFORMS=""
for name in Data SDS; do
  if [[ -n "$SOURCE" && -d "$SOURCE/$name" ]]; then
    WAVEFORMS="$SOURCE/$name"
    break
  fi
done
if [[ -z "$WAVEFORMS" || ! -d "$SOURCE/INVENTORIES" ]]; then
  printf 'Usage: bash scripts/run_fresh_pressure_ppsd.sh DATA_FOLDER [stage] [options]\n' >&2
  printf 'DATA_FOLDER must contain INVENTORIES/ and Data/ or SDS/.\n' >&2
  exit 2
fi
shift
STAGE="all"
if [[ $# -gt 0 && "$1" != --* ]]; then
  STAGE="$1"
  shift
fi
OUT="${BRUITFM_FRESH_DIR:-$HOME/Desktop/BRUIT-FM/fresh_run_$(date +%Y-%m-%d)}"
mkdir -p "$OUT/master" "$OUT/results" "$OUT/logs"
cd "$ROOT"
export PYTHONPATH="$ROOT/.."
# One numerical thread per PPSD worker process.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
COMMON=(--data-dir "$WAVEFORMS" --inventory-dir "$SOURCE/INVENTORIES"
  --master-dir "$OUT/master" --results-dir "$OUT/results")

run() {
  local stage="$1"
  shift
  local log
  log="$OUT/logs/$(date +%Y%m%dT%H%M%S)_$stage.log"
  printf '== %s -> %s\n' "$stage" "$log"
  {
    printf '# git %s\n' "$(git -C "$ROOT" rev-parse HEAD)"
    printf '# %s -m compy_pressure.network %s' "$PYTHON" "$stage"
    printf ' %q' "${COMMON[@]}" "$@"
    printf '\n'
  } >"$log"
  # caffeinate keeps the Mac awake while a stage runs.
  caffeinate -i "$PYTHON" -m compy_pressure.network "$stage" "${COMMON[@]}" "$@" 2>&1 | tee -a "$log"
}

# check lists waveform channels without StationXML and exits nonzero; the
# later stages skip those channels (recorded in _workflow/uninventoried.json).
check() {
  run check "$@" || printf 'check reported channels without StationXML; continuing\n'
}
case "$STAGE" in
  inventory) run inventory "$@" ;;
  check) run check "$@" ;;
  ppsd) run ppsd --pressure-only --skip-uninventoried --workers max "$@" ;;
  plots) run plots --pressure-only --skip-uninventoried "$@" ;;
  all)
    run inventory
    check "$@"
    run ppsd --pressure-only --skip-uninventoried --workers max "$@"
    run plots --pressure-only --skip-uninventoried "$@"
    ;;
  *)
    printf 'Unknown stage %s; use inventory, check, ppsd, plots or all\n' "$STAGE" >&2
    exit 2 ;;
esac
printf 'Outputs in %s\n' "$OUT"
