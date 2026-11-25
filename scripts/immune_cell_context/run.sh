#!/usr/bin/env bash
# Replay the frozen analyses with deposited inputs or full upstream sources.
set -euo pipefail
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  printf 'Usage: %s EXTRACTED_DATA NEW_OUTPUT_DIR [EXTERNAL_INPUT_ROOT]\n' "$0" >&2
  exit 2
fi
scripts="$(cd "$(dirname "$0")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONHASHSEED=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME PYTHONSTARTUP
exec "${PYTHON:-python3}" "$scripts/replay.py" "$@"
