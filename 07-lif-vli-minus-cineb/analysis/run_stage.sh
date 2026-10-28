#!/bin/bash
# Reproduce a completed stage only after checking its input and current job state.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VASP="${VASP:?Set VASP to your licensed VTST-enabled vasp_std executable}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$ROOT/neb"
/usr/bin/time -f 'WALL_SECONDS=%e' -o "$ROOT/logs/${STAGE_NAME:?Set STAGE_NAME}.time" \
  mpirun --oversubscribe -np 6 "$VASP" > vasp.out 2>&1
