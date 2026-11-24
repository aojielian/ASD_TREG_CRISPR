#!/usr/bin/env bash
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
exec bash "$repo/scripts/smoke_test.sh" "$@"
