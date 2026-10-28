#!/bin/bash
# Stage B (cleaning) for all stations in parallel, then Stage D (calibration).
cd "$(dirname "$0")/.."
PY=/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
printf "%s\n" RR52 RR28 RR29 RR34 RR36 RR38 RR40 RR50 | xargs -P 8 -I{} sh -c "$PY scripts/b_clean.py {} > logs/b_clean_{}.log 2>&1; echo B {} exit \$?"
printf "%s\n" RR28 RR29 RR34 RR36 RR38 RR40 RR50 | xargs -P 4 -I{} sh -c "$PY scripts/d_calibrate.py {} > logs/d_calibrate_{}.log 2>&1; echo D {} exit \$?"
echo ALL_DONE
