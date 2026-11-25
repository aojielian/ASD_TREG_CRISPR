#!/usr/bin/env bash
# Replay directional and phenotype comparisons.
set -euo pipefail

if [ "$#" -ne 2 ]; then
  printf 'Usage: %s EXTRACTED_DATA NEW_OUTPUT_DIR\n' "$0" >&2
  exit 2
fi

scripts="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$scripts/../.." && pwd)"
data="$(cd "$1" && pwd)"
python="${PYTHON:-python3}"
rscript="${RSCRIPT:-Rscript}"

export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME PYTHONSTARTUP
output="$("$python" -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' "$2")"
case "$output/" in
  "$data/"*|"$repo/"*) printf 'Keep output outside the source and data directories.\n' >&2; exit 2 ;;
esac
if [ -e "$output" ] && { [ ! -d "$output" ] || [ -n "$(ls -A "$output")" ]; }; then
  printf 'Choose a new or empty output directory.\n' >&2
  exit 2
fi
mkdir -p "$output"/{results,checks}
"$python" "$scripts/baseline_release.py" "$data" "$repo"
if [[ -f "$data/data/manifests/compatibility/REPLAY_PATHS.tsv" ]]; then
  compatibility_dir="$(mktemp -d "${TMPDIR:-/tmp}/asd-functional-inputs.XXXXXX")"
  trap 'rm -rf -- "$compatibility_dir"' EXIT
  "$python" "$repo/scripts/release/layout.py" replay "$data" "$compatibility_dir/data"
  data="$compatibility_dir/data"
fi

status=0
"$python" "$scripts/lof_directional.py" "$data" "$repo" "$output/results" \
  --checks "$output/checks/LOF_DIRECTIONAL_CHECKS.tsv" || status=1
"$python" "$scripts/phenotype_preference.py" "$data" "$repo" "$output/results" || status=1
[ "$status" -eq 0 ] || exit "$status"

"$python" "$scripts/verify_exact_comparisons.py" "$data" "$output/checks" --production "$output/results"
"$rscript" --vanilla "$scripts/verify_discordant_firth.R" "$data" "$output/checks" "$output/results"
if [ -d "$data/data/r9/results" ]; then
  "$python" "$scripts/compare_tables.py" "$data/data/r9/results" "$output/results" "$output/checks/R9_REPLAY_PARITY.tsv"
fi
printf 'Functional comparisons and independent verification finished in %s\n' "$output"
