#!/usr/bin/env python3
"""Replay immune-cell context computations using released and optional upstream inputs."""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

from baseline_recheck import baseline
from integrate_context_results import assemble, table, validate, write
from replay_inputs import records, verify_baseline_file, verify_file

SCRIPTS = Path(__file__).resolve().parent
REPO = SCRIPTS.parents[1]
sys.path.insert(0, str(REPO / "scripts/functional_comparisons"))
from compare_tables import compare, write_report
sys.path.insert(0, str(REPO / "scripts/release"))
from layout import replay_view

ATLAS_TABLES = frozenset("""
R10_EXPRESSION_DONOR_CONTRASTS.tsv
R10_EXPRESSION_MATCHED_NULL.tsv
R10_EXPRESSION_MATCH_BALANCE.tsv
R10_EXPRESSION_SENSITIVITY.tsv
R10_HUMAN_EXPRESSION_PANEL.tsv
R10_HUMAN_EXPRESSION_SENSITIVITY.tsv
R10_IMMGEN_CELLTYPE_CPM.tsv.gz
R10_IMMGEN_PANEL_LIBRARY_CPM.tsv.gz
R10_MONACO_CELLTYPE_MEDIAN_TPM.tsv.gz
R10_MONACO_EXPRESSION_HEATMAP_SOURCE.tsv
R10_MONACO_LINEAGE_HEATMAP_SOURCE.tsv
R10_MONACO_MATCH_MEMBERS.tsv.gz
R10_MONACO_MATCH_NULL_DRAWS.tsv.gz
R10_MONACO_MATCH_STRATA.tsv
R10_MOUSE_EXPRESSION_PANEL.tsv
R10_MOUSE_EXPRESSION_PLOT.tsv
R10_MOUSE_EXPRESSION_SENSITIVITY.tsv
R10_R9_BASELINE_RECHECK.tsv
""".split())


def run_command(output, name, command):
    print(f"Immune-cell context: {name}", flush=True)
    with (output / "08_logs" / f"{name}.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)


def run_python(output, name, *args):
    run_command(output, Path(name).stem, [sys.executable, SCRIPTS / name, *args])


def run_r(output, name, *args):
    run_command(output, Path(name).stem, [os.environ.get("RSCRIPT", "Rscript"), "--vanilla", SCRIPTS / name, *args])


def stage_verified_inputs(data, external, inputs, audit):
    rows = []
    for row in records("inputs.tsv"):
        category = row["mode"]
        if category == "external":
            if external is None:
                continue
            source, local = external / row["path"], inputs / row["path"]
        elif category == "packaged_source":
            source = data / row["path"]
            local = inputs / source.parent.name / source.name
        else:
            source, local = data / row["path"], None
        verify_file(source, row)
        rows.append({**row, "status": "PASS"})
        if local:
            local.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, local)
    with (audit / "R10_VERIFIED_INPUTS.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parity_report(data, output, full):
    frozen = data / "data/r10/results"
    results = output / "03_results"
    names = sorted(p.name for p in frozen.iterdir() if p.name.endswith((".tsv", ".tsv.gz")))
    if not ATLAS_TABLES.issubset(names):
        raise ValueError("The immune-context result inventory is incomplete")
    comparisons = []
    for name in names:
        if not full and name not in ATLAS_TABLES:
            comparisons.append((name, "REQUIRES_EXTERNAL", 0, 0, "Retained frozen result. Use the full upstream-input replay."))
            continue
        ignore = ("source",) if name == "R10_R9_BASELINE_RECHECK.tsv" else ()
        comparisons.append((name, *compare(frozen / name, results / name, rtol=1e-7, atol=1e-10, ignore=ignore)))
    write_report(output / "08_logs/R10_REPLAY_PARITY.tsv", comparisons)
    if any(row[1] in ("MISSING", "FAIL") for row in comparisons):
        for name, status, _, _, detail in comparisons:
            if status in ("MISSING", "FAIL"):
                print(f"{status}\t{name}\t{detail}")
        raise ValueError("One or more immune-context tables differ from the release")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="Extracted v1.2.0 Zenodo processed-data root")
    parser.add_argument("output", type=Path)
    parser.add_argument("external", type=Path, nargs="?", help="Directory laid out as in config/inputs.tsv")
    args = parser.parse_args()
    data, output = args.data.resolve(), args.output.resolve()
    external = args.external.resolve() if args.external else None
    if output.is_relative_to(data) or output.is_relative_to(REPO) or (external and output.is_relative_to(external)):
        parser.error("Keep output outside the source and data directories")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("Choose a new or empty output directory")
    with replay_view(data) as compatibility:
        replay(compatibility, data, output, external)


def replay(data, release_data, output, external):
    if not (data / "data/r10/results").is_dir():
        raise ValueError("Missing immune-context result tables; use the extracted v1.2.0 deposit")
    for subdir in ("01_audit", "03_results", "07_config", "08_logs"):
        (output / subdir).mkdir(parents=True, exist_ok=True)
    for source in (SCRIPTS / "config").glob("*.tsv"):
        shutil.copyfile(source, output / "07_config" / source.name)
    ref = data / "data/r10/reference"
    audit = output / "01_audit"
    results = output / "03_results"
    shifrut = ref / "SHIFRUT_CANONICAL_SCREEN.tsv.gz"

    try:
        for row in records("frozen_baseline_inputs.tsv"):
            verify_baseline_file(release_data, row["historical_alias"])
        with TemporaryDirectory(prefix="asd-immune-context-inputs-") as work:
            inputs = Path(work)
            stage_verified_inputs(release_data, external, inputs, audit)
            baseline(data, shifrut, results)
            monaco = ["--data", data, "--stage", output, "--inputs", inputs]
            if external is None:
                monaco.extend(("--frozen-mappings", ref))
            run_python(output, "monaco_expression.py", *monaco)
            run_python(output, "verify_monaco.py", "--stage", output, "--source-tpm",
                       inputs / "GSE107011/GSE107011_Processed_data_TPM.txt.gz")
            if external:
                run_python(output, "immgen_metadata.py", "--download-dir", inputs, "--out", output)
            else:
                for name in ("SAMPLE", "GROUP", "GENE"):
                    filename = f"R10_IMMGEN_{name}_CROSSWALK.tsv"
                    shutil.copyfile(ref / filename, audit / filename)
            run_r(output, "immgen_tmm.R", inputs, output)
            run_python(output, "immgen_expression.py", "--out", output)
            if external:
                run_python(output, "umhoefer_comparison.py", "--data", data, "--out", output,
                           "--download-dir", inputs / "GSE286470", "--firth-dir", REPO / "scripts/gene_property_adjustment",
                           "--shifrut-table", shifrut)
                run_r(output, "verify_umhoefer_numeric.R", output)
                run_python(output, "chen_perturb.py", inputs, data, results)
                assemble(output)
                validate(output)
            else:
                combine_atlas_sensitivities(results)
        parity_report(data, output, external is not None)
    except subprocess.CalledProcessError as exc:
        print(f"Immune-cell context subprocess failed with exit status {exc.returncode}; see {output / '08_logs'}", file=sys.stderr)
        raise SystemExit(1) from None
    except (OSError, ValueError) as exc:
        print(f"Immune-cell context replay failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from None
    mode = "full upstream-input" if external else "deposited-input atlas and baseline"
    print(f"Immune-cell context {mode} statistical replay finished in {output}")


def combine_atlas_sensitivities(results):
    import pandas as pd
    panels = [table(results / f"R10_{species}_EXPRESSION_SENSITIVITY.tsv") for species in ("HUMAN", "MOUSE")]
    write(pd.concat(panels, ignore_index=True), results / "R10_EXPRESSION_SENSITIVITY.tsv")


if __name__ == "__main__":
    main()
