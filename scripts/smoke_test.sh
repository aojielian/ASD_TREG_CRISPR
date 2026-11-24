#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 2 ]; then
  printf 'Usage: %s DATA_ZIP_OR_DIRECTORY EMPTY_OUTPUT_DIRECTORY\n' "$0" >&2
  exit 2
fi
repo="$(cd "$(dirname "$0")/.." && pwd)"
if [ "$(uname -s)" = Darwin ] && [ "$(uname -m)" = x86_64 ]; then
  exec /usr/bin/arch -x86_64 "${PYTHON:-python3}" "$repo/scripts/smoke_test.py" "$1" --out-dir "$2" --verify-r --figures
fi
exec "${PYTHON:-python3}" "$repo/scripts/smoke_test.py" "$1" --out-dir "$2" --verify-r --figures
