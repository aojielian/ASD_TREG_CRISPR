#!/usr/bin/env python3
"""Run the published v1.0.0 numerical workflow against extracted release data."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import subprocess
import sys


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def flag(value: str) -> bool:
    return value.upper() in {"TRUE", "YES", "1"}


def contingency(
    rows: list[dict[str, str]], membership: str, outcome: str
) -> tuple[int, int, int, int]:
    counts = [0, 0, 0, 0]
    for row in rows:
        i = 0 if flag(row[membership]) else 2
        j = 0 if flag(row[outcome]) else 1
        counts[i + j] += 1
    a, b, c, d = counts
    return a, b, c, d


def fisher_enrichment(a: int, b: int, c: int, d: int) -> float:
    """One-sided exact upper-tail probability from the fixed margins."""
    hits = a + c
    cases = a + b
    controls = c + d
    count = sum(
        math.comb(cases, k) * math.comb(controls, hits - k)
        for k in range(max(a, hits - controls), min(cases, hits) + 1)
    )
    return count / math.comb(cases + controls, hits)


def run_independent_checks(
    data: Path, expected_path: Path, output: Path, rscript: str | None
) -> int:
    """Recompute values that the historical workflow carries from frozen rows."""
    expected = {row["result_id"]: row for row in read_table(expected_path)}
    comparisons: list[tuple[str, int | float, str]] = []
    removed_ids = {"FOXP1", "ZMYND8"}
    for screen, filename in [
        ("CHEN", "chen_genomewide_analysis_ready.tsv"),
        ("TREG", "umhoefer_primary_treg_analysis_ready.tsv"),
    ]:
        rows = read_table(data / filename)
        omitted = [row for row in rows if row["gene"] in removed_ids]
        if len(omitted) != len(removed_ids) or any(
            not (flag(row["sfari_score1"]) and flag(row["screen_hit"])) for row in omitted
        ):
            raise ValueError(f"{screen}: invalid frozen two-hit removal set")
        retained = [row for row in rows if row["gene"] not in removed_ids]
        cells = contingency(retained, "sfari_score1", "screen_hit")
        prefix = f"FOXP1_ZMYND8_{screen}"
        comparisons.extend(
            [
                (f"{prefix}_REMAINING_HITS", cells[0], "Rows after removal of two frozen hits"),
                (f"{prefix}_P", fisher_enrichment(*cells), "One-sided exact Fisher enrichment"),
            ]
        )

    if rscript is not None:
        mouse = read_table(data / "loo_mouse_analysis_ready.tsv")
        mouse_cells = contingency(mouse, "sfari_score1", "LOO_FOXP3_HIT")
        # Use the reference R implementation and its original root-finder precision.
        expr = (
            "cells <- as.integer(commandArgs(trailingOnly=TRUE)); "
            "fit <- fisher.test(matrix(cells, 2, byrow=TRUE), alternative='greater'); "
            "cat(sprintf('%.17g', unname(fit$estimate)))"
        )
        result = subprocess.run(
            [rscript, "--vanilla", "-e", expr, *map(str, mouse_cells)],
            check=True,
            capture_output=True,
            text=True,
        )
        comparisons.append(
            (
                "LOO_REFERENCE_R_CONDITIONAL_OR",
                float(result.stdout),
                "Base R fisher.test conditional odds ratio from gene-table fixed margins",
            )
        )

    moreno = read_table(data / "moreno_gene_level_analysis_ready.tsv")
    moreno_cells = contingency(moreno, "ASD14_MEMBER", "deg_fdr05")
    comparisons.append(
        (
            "MORENO_ASD14_FISHER_P",
            fisher_enrichment(*moreno_cells),
            "One-sided exact Fisher enrichment from the 14-gene membership",
        )
    )

    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["result_id", "expected_value", "observed_value", "tolerance", "status", "method"])
        failing = []
        for key, observed, source in comparisons:
            reference = expected[key]
            target = reference["expected_value"]
            tolerance = reference["tolerance"]
            passed = (
                int(observed) == int(target)
                if reference["value_type"] == "integer"
                else abs(float(target) - observed) <= float(tolerance)
            )
            if not passed:
                failing.append(key)
            writer.writerow([key, target, observed, tolerance, "PASS" if passed else "FAIL", source])
    if failing:
        raise ValueError(f"Legacy reference checks failed: {', '.join(failing)}")
    return len(comparisons)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path, help="Root of the extracted processed-data archive")
    parser.add_argument("output_dir", type=Path, help="Directory for numerical verification outputs")
    parser.add_argument(
        "--rscript", metavar="COMMAND", help="Also refit the historical R conditional odds ratio"
    )
    args = parser.parse_args()
    snapshot = args.data_root.resolve()
    output = args.output_dir.resolve()
    checkout = Path(__file__).resolve().parents[2]
    if output == snapshot or snapshot in output.parents:
        parser.error("choose an output directory outside the extracted data archive")
    old = checkout / "scripts" / "reproduce_from_processed"
    data = snapshot / "data" / "processed"
    expected = snapshot / "expected_results" / "EXPECTED_CORE_RESULTS.tsv"
    checkout_expected = checkout / "expected_results" / expected.name
    if not data.is_dir() or not expected.is_file():
        parser.error("data_root must contain data/processed/ and expected_results/EXPECTED_CORE_RESULTS.tsv")
    if expected.read_bytes() != checkout_expected.read_bytes():
        raise ValueError("Frozen v1.0.0 expected results differ between code and data")
    output.mkdir(parents=True, exist_ok=True)

    core_results = output / "recomputed_core.tsv"
    patient_results = output / "recomputed_patient.tsv"
    verification = output / "expected_results_verification.tsv"
    core_script = old / "reproduce_core_statistics.py"
    commands = [
        [core_script, snapshot, core_results],
        [old / "reproduce_patient_statistics.py", snapshot, patient_results],
        [old / "verify_against_expected.py", expected, core_results, patient_results, verification],
    ]
    for command in commands:
        subprocess.run([sys.executable, *map(str, command)], check=True)

    n_legacy = len(read_table(verification))
    n_check = run_independent_checks(
        data, expected, output / "independent_reference_checks.tsv", args.rscript
    )
    print(f"LEGACY_EXPECTED_RESULTS={n_legacy}")
    print(f"LEGACY_RECOMPUTED_REFERENCE_ROWS={n_check}")
    print(
        "LEGACY_R_CONDITIONAL_OR_REFIT="
        + ("PASS" if args.rscript is not None else "OPTIONAL_BASE_R_REQUIRED")
    )
    print("LEGACY_PROCESSED_REPRODUCTION_STATUS=PASS")


if __name__ == "__main__":
    main()
