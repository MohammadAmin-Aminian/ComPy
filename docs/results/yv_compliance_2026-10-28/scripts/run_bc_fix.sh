#!/bin/bash
# Rerun cleaning (B) and compliance (C) for stations given as arguments, in parallel.
cd "$(dirname "$0")/.."
PY=/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for s in "$@"; do
  (rm -rf work/$s/clean
   $PY scripts/b_clean.py $s > logs/b_clean_$s.log 2>&1 \
     && $PY scripts/c_compliance.py $s > logs/c_compliance_$s.log 2>&1 \
     && $PY scripts/g_tiskit_check.py $s > logs/g_tiskit_$s.log 2>&1
   echo "BCG $s $?") &
done
wait
echo ALL_DONE
