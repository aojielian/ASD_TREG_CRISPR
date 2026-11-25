#!/usr/bin/env python3
"""Reproduce the gene-property analysis from the verified processed-data deposit."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

# The frozen fits use one BLAS thread for deterministic reduction and bounded RAM.
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                  "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_variable] = "1"

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True

import background_and_matching
import firth_models
from reconstruct_background import rebuild_background, rebuild_strata

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "release"))
from layout import REPLAY_MAP, replay_view
DEFAULT_EXPECTED = SCRIPT_DIR.parent.parent / "config/gene_property_expected_results.tsv"


def run_expected_checks(out_dir, fixture):
    checked = []
    cache = {}
    with fixture.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            rel = row["artifact"]
            if rel not in cache:
                path = (out_dir / rel).resolve(strict=True)
                if not path.is_relative_to(out_dir):
                    raise ValueError(f"Expected-results artifact escapes the output directory: {rel}")
                cache[rel] = firth_models.rd(path)
            frame = cache[rel]
            select = np.ones(len(frame), dtype=bool)
            for name, value in json.loads(row["selector"]).items():
                vals = frame[name]
                select &= vals.isna() if value is None else vals.eq(value)
            which = frame.loc[select, row["field"]]
            if len(which) != 1:
                raise ValueError(f"Expected one output row in {rel}: {row['selector']}")
            actual = which.iloc[0]
            expected, rule = row["value"], row["rule"]
            if rule == "missing":
                accepted = pd.isna(actual)
            elif rule == "text":
                accepted = not pd.isna(actual) and str(actual) == expected
            elif rule == "numeric":
                accepted = np.isclose(
                    float(actual), float(expected), rtol=float(row["relative_tolerance"]),
                    atol=float(row["absolute_tolerance"]), equal_nan=False,
                )
            else:
                raise ValueError(f"Unknown comparison rule: {rule}")
            checked.append(row | {"observed": actual, "status": "PASS" if accepted else "FAIL"})
    if not checked:
        raise ValueError("Expected-results fixture is empty")
    firth_models.write(checked, out_dir / "checks/EXPECTED_RESULTS_QA.tsv")
    failures = [row for row in checked if row["status"] == "FAIL"]
    if failures:
        first = failures[0]
        raise ValueError(
            f"{len(failures)} numeric or status checks failed. First failure in "
            f"{first['artifact']} {first['selector']} {first['field']}: "
            f"expected {first['value']!r}, observed {first['observed']!r}. "
            "See checks/EXPECTED_RESULTS_QA.tsv."
        )
    print(f"Verified {len(checked)} released numerical values and statuses.")


def reproduce(data_dir, out_dir, expected, verify_r=False):
    data_dir = data_dir.resolve(strict=True)
    out_dir = out_dir.resolve()
    if out_dir.is_relative_to(data_dir) or data_dir.is_relative_to(out_dir):
        raise ValueError("Use an output directory separate from the processed-data deposit")
    if out_dir.exists() and any(out_dir.iterdir()):
        raise ValueError("Specify a new or empty output directory")
    (out_dir / "results").mkdir(parents=True)
    (out_dir / "checks").mkdir()

    genes, historical = firth_models.read_inputs(data_dir)
    genes, state_counts, eligible = rebuild_background(data_dir, genes)
    strata = rebuild_strata(data_dir, genes, eligible)

    fields = ["screen", "approved_symbol", "hgnc_id", "is_measured", "asd", "hit",
              "primary_mean_log_expression", "hpa_mean_log_expression",
              "primary_status_0.5", "primary_status_1", "primary_status_2", "hpa_status_1"]
    gene_out = genes[fields].copy()
    gene_out["hit"] = gene_out["hit"].astype("Int64").where(gene_out["is_measured"])
    gene_out.to_csv(
        out_dir / "results/GENE_DETECTABILITY_RECALCULATION.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "mtime": 0}, float_format="%.13g",
    )
    firth_models.write(state_counts, out_dir / "checks/DETECTION_STATE_COUNTS.tsv")
    firth_models.write(strata, out_dir / "checks/R1_CONDITIONING_STRATA_RECALCULATION.tsv")
    print(f"Rebuilt detectability for {len(genes):,} gene-screen rows and {len(strata):,} exact strata.")

    firth_models.run(data_dir, out_dir, genes, historical)
    background_and_matching.main(data_dir, out_dir, genes, strata)
    run_expected_checks(out_dir, expected)
    if verify_r:
        rscript = shutil.which("Rscript")
        if not rscript:
            raise FileNotFoundError("The independent verification requires Rscript on PATH")
        subprocess.run(
            [rscript, "--vanilla", str(SCRIPT_DIR / "verify_primary_r.R"),
             str(data_dir), str(out_dir)],
            check=True,
        )
    print("Gene-property processed-data reproduction complete.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path,
                        help="Root of the extracted processed-data archive")
    parser.add_argument("--out-dir", required=True, type=Path, help="New work directory")
    parser.add_argument("--expected", type=Path, default=DEFAULT_EXPECTED,
                        help="Expected-results TSV (default is the versioned repository fixture)")
    parser.add_argument("--verify-r", action="store_true",
                        help="Also run the independent base-R primary-model check")
    args = parser.parse_args()
    data = args.data_dir.resolve()
    out = args.out_dir.resolve()
    if out.is_relative_to(data) or data.is_relative_to(out):
        parser.error("Use an output directory separate from the processed-data deposit")
    if (data / REPLAY_MAP).is_file():
        with replay_view(data) as compatibility:
            reproduce(compatibility, out, args.expected, args.verify_r)
    else:
        reproduce(data, out, args.expected, args.verify_r)


if __name__ == "__main__":
    main()
