#!/bin/bash
# Stage A for all stations (best vertical noise first), 14 workers.
cd "$(dirname "$0")/.."
PY=/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cd scripts
$PY a_preprocess.py CC11 CC05 CC06 EE04 WC03 CC04 WC02 EE01 CC09 WC05 CC08 EC04 WC04 WC01 CC03 EE03 EE02 \
    CC07 EC05 WW03 WW04 WW01 CC01 EC01 EC02 EC03 WW02 --workers 14 2> ../logs/a_preprocess.err | grep -v "FIR normalized"
echo A_DONE
