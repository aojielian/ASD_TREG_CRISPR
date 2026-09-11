#!/bin/sh
set -eu
SNAPSHOT=${1:?usage: run_all_processed.sh EXTRACTED_ZENODO_ROOT [OUTPUT_DIR]}
OUTPUT=${2:-outputs/minimal_reproduction}
mkdir -p "$OUTPUT"
PYTHON=${PYTHON:-python3}
"$PYTHON" scripts/reproduce_from_processed/reproduce_core_statistics.py "$SNAPSHOT" "$OUTPUT/recomputed_core.tsv"
"$PYTHON" scripts/reproduce_from_processed/reproduce_patient_statistics.py "$SNAPSHOT" "$OUTPUT/recomputed_patient.tsv"
"$PYTHON" scripts/reproduce_from_processed/verify_against_expected.py "$SNAPSHOT/expected_results/EXPECTED_CORE_RESULTS.tsv" "$OUTPUT/recomputed_core.tsv" "$OUTPUT/recomputed_patient.tsv" "$OUTPUT/expected_results_verification.tsv"
status=$(awk -F '\t' 'NR>1 && $6 != "PASS" {bad=1} END {print bad ? "FAIL" : "PASS"}' "$OUTPUT/expected_results_verification.tsv")
failures=$(awk -F '\t' 'NR>1 && $6 != "PASS" {bad++} END {print bad+0}' "$OUTPUT/expected_results_verification.tsv")
printf 'EXPECTED_RESULTS_VERIFICATION=%s\n' "$status" | tee "$OUTPUT/EXPECTED_RESULTS_VERIFICATION.txt"
printf 'EXPECTED_RESULTS_FAILURES=%s\n' "$failures" | tee "$OUTPUT/EXPECTED_RESULTS_FAILURES.txt"
printf 'MINIMAL_REPRODUCTION_STATUS=%s\n' "$status" | tee "$OUTPUT/MINIMAL_REPRODUCTION_STATUS.txt"
[ "$status" = PASS ]
[ "$failures" -eq 0 ]
