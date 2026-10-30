#!/usr/bin/env bash
# Full-data coverage: read-only plan, dry run, explicit execution, report.
#   plan     inventory + plan; writes only a plan folder and header cache
#   dry-run  rebuild the plan; list tasks that would run; nothing executed
#   execute  run only missing tasks (campaign lock, logs, resume, protection)
#   report   coverage tables and pressure probability heat maps
# Never applies gains, assigns scores/thresholds or writes saved products.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run with bash; do not source this script.\n' >&2
  return 2
fi
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SOURCE_ROOT="${BRUITFM_SOURCE_ROOT:-$HOME/Projects/bruit-fm}"
CAMPAIGN="${BRUITFM_COVERAGE_CAMPAIGN:-full_data_coverage_v1}"
PYTHON="${BRUITFM_PYTHON:-/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python}"
MODE="${1:-}"
[[ $# -gt 0 ]] && shift
case "$MODE" in
  plan) ARGS=(plan) ;;
  dry-run) ARGS=(run) ;;
  execute) ARGS=(run --execute) ;;
  report) ARGS=(report) ;;
  *)
    printf 'Usage: bash scripts/run_full_coverage.sh {plan|dry-run|execute|report} [options]\n' >&2
    printf 'Options are passed to python -m compy_pressure.coverage (e.g. --workers 4,\n' >&2
    printf '%s\n' '--ppsd-scope all-channels, --only-kind ppsd, --only-seed NET.STA.LOC.CHA).' >&2
    exit 2 ;;
esac
cd "$ROOT"
# One numerical thread per PPSD worker process.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$ROOT/.."
exec "$PYTHON" -m compy_pressure.coverage "${ARGS[@]}" \
  --source-root "$SOURCE_ROOT" --campaign "$CAMPAIGN" "$@"
