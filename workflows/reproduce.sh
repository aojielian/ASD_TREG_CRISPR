#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
  printf 'Usage: %s VERIFIED_DATA_ROOT NEW_OUTPUT_ROOT [EXTERNAL_INPUT_ROOT]\n' "$0" >&2
  exit 2
fi

repo="$(cd "$(dirname "$0")/.." && pwd)"
data="$1"
out="$2"
if [[ "$#" -eq 2 ]]; then
  if [[ "$data" == *ASD_TREG_CRISPR_processed_data_v1.1.0.zip ]] || \
     [[ -f "$data/VERSION.txt" && "$(cat "$data/VERSION.txt")" == ASD_TREG_CRISPR_processed_data_v1.1.0 ]]; then
    exec bash "$repo/scripts/smoke_test.sh" "$data" "$out"
  fi
fi

export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONSTARTUP
python="${PYTHON:-python3}"
data="$(cd "$data" && pwd)"
out="$("$python" -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' "$out")"
case "$out/" in
  "$data/"*|"$repo/"*) printf 'Keep output outside the source and data directories.\n' >&2; exit 2 ;;
esac
if [[ "$#" -eq 3 ]]; then
  external="$(cd "$3" && pwd)"
  case "$out/" in
    "$external/"*) printf 'Keep output outside the original-source directory.\n' >&2; exit 2 ;;
  esac
fi
mkdir -p "$out"
if [[ -n "$(ls -A "$out")" ]]; then
  printf 'Use an empty output directory: %s\n' "$out" >&2
  exit 2
fi

compatibility_dir="$(mktemp -d "${TMPDIR:-/tmp}/asd-published-core-inputs.XXXXXX")"
trap 'rm -rf -- "$compatibility_dir"' EXIT
"$python" "$repo/scripts/release/verify_release.py" "$data" \
  --report "$out/TSV_SCHEMA_QA.tsv" --legacy-view "$compatibility_dir/published_baseline"
bash "$repo/scripts/smoke_test.sh" "$compatibility_dir/published_baseline" "$out/historical_core"
bash "$repo/scripts/functional_comparisons/run.sh" "$data" "$out/functional_comparisons"
if [[ "$#" -eq 3 ]]; then
  bash "$repo/scripts/immune_cell_context/run.sh" "$data" "$out/immune_cell_context" "$external"
  bash "$repo/scripts/conditional_knockout/run.sh" "$data" "$out/conditional_knockout" "$external" \
    --rebuild-samples --rebuild-signatures --rebuild-literature
else
  bash "$repo/scripts/immune_cell_context/run.sh" "$data" "$out/immune_cell_context"
  bash "$repo/scripts/conditional_knockout/run.sh" "$data" "$out/conditional_knockout"
fi
