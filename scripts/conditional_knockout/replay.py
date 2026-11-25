"""Recompute the frozen mouse perturbation analysis from the release matrices."""

import argparse
import csv
from hashlib import file_digest
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

from integrate_candidate_evidence import DOIS, KNOWN, PANEL
from verify_signatures import run as verify_signatures


SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS.parent / "release"))
from layout import replay_view
AUDIT_INPUTS = (
    "R11_SAMPLE_MAP.tsv",
    "R11_ID_MAP.tsv.gz",
    "R11_RNASEQ_DESIGN.tsv",
    "R11_SIGNATURE_DEFINITIONS.tsv",
)
FROZEN_LITERATURE = ("R11_14_GENE_LITERATURE_AUDIT.tsv", "R11_LITERATURE_GENE_SUMMARY.tsv")
CHEN_ARCHIVE = "Chen2025/41586_2025_8795_MOESM5_ESM.zip"
FIELD_NAMES = ("dataset", "contrast", "component", "program", "gene_id", "gene_symbol",
               "mouse_gene", "GSM", "sample_id", "pmid", "doi", "state", "genotype", "sex")


def tsv(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str, low_memory=False)


def sha256(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def verify_file(path, expected):
    if not path.is_file():
        raise FileNotFoundError(path)
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"{path.name} SHA-256 mismatch: expected {expected}; got {actual}")


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def verify_inputs(data, external, enabled):
    provenance = tsv(SCRIPTS / "inputs.tsv")
    verified = []
    for row in provenance.itertuples(index=False):
        if row.use not in enabled:
            continue
        root = external if row.use.startswith("rebuild_") and row.origin == "external" else data
        if row.origin == "code":
            root = SCRIPTS
        path = root / row.path
        try:
            verify_file(path, row.sha256)
        except FileNotFoundError as error:
            raise FileNotFoundError(f"{row.path} missing. Official source: {row.official_source}") from error
        verified.append((row.use, row.path, row.sha256, "PASS"))
    return verified


def verify_literature(results):
    evidence = tsv(results / "R11_14_GENE_LITERATURE_AUDIT.tsv")
    summary = tsv(results / "R11_LITERATURE_GENE_SUMMARY.tsv")
    if summary.gene.tolist() != PANEL or set(evidence.gene) != set(PANEL):
        raise ValueError("The literature records must cover the frozen 14-candidate panel.")
    cited = evidence[evidence.pmid.ne("NA")]
    if not (
        cited.pmid.str.fullmatch(r"\d+").all()
        and cited.doi.str.match(r"10\.").all()
        and cited.pubmed_url.eq("https://pubmed.ncbi.nlm.nih.gov/" + cited.pmid + "/").all()
        and cited.doi_url.eq("https://doi.org/" + cited.doi).all()
    ):
        raise ValueError("Incomplete or inconsistent publication identifiers in the frozen curation.")
    for gene, _, _, _, _, _, _, reference, _ in KNOWN:
        pmid, _ = reference.split(",", 1)
        rows = cited.loc[cited.gene.eq(gene) & cited.pmid.eq(pmid)]
        if rows.empty or not rows.doi.str.lower().eq(DOIS[pmid].lower()).all():
            raise ValueError(f"Candidate integration source absent from frozen curation: {gene}, PMID {pmid}.")
    return len(summary), len(evidence), cited.pmid.nunique()


def table_match(expected, actual, enforce_exact=False):
    if expected == actual or sha256(expected) == sha256(actual):
        frame = tsv(actual)
        return "BYTE_IDENTICAL", len(frame), len(frame.columns), 0.0
    reference = tsv(expected)
    observed = tsv(actual)
    if not reference.columns.equals(observed.columns) or reference.shape != observed.shape:
        raise ValueError(f"{actual.name}: output does not match the released table shape/schema")
    max_difference = 0.0
    for column in reference:
        left = reference[column]
        right = observed[column]
        if left.equals(right):
            continue
        if enforce_exact or column in FIELD_NAMES or not left.eq("NA").equals(right.eq("NA")):
            raise ValueError(f"{actual.name}: column {column} differs from the released table")
        measured = left.ne("NA")
        left_numeric = pd.to_numeric(left[measured], errors="coerce").to_numpy(dtype=float)
        right_numeric = pd.to_numeric(right[measured], errors="coerce").to_numpy(dtype=float)
        if not (np.isfinite(left_numeric).all() and np.isfinite(right_numeric).all()):
            raise ValueError(f"{actual.name}: non-numeric column {column} differs from the released table")
        if column.lower() in {"p", "pvalue", "p_adjusted", "adjusted_p", "camera_p", "padj"}:
            close = np.isclose(left_numeric, right_numeric, rtol=1e-6, atol=0)
            close &= (left_numeric < 0.05) == (right_numeric < 0.05)
        else:
            close = np.isclose(left_numeric, right_numeric, rtol=1e-7, atol=5e-12)
        if not np.all(close):
            indices = np.flatnonzero(~close)[:4].tolist()
            raise ValueError(f"{actual.name}: column {column} differs at measured row positions {indices}")
        max_difference = max(max_difference, float(np.max(np.abs(left_numeric - right_numeric))))
    return "NUMERIC_PARITY", len(reference), len(reference.columns), max_difference


def call(script, *args):
    python = sys.executable
    rscript = os.environ.get("RSCRIPT", "Rscript")
    interpreter = [rscript, "--vanilla"] if script.endswith(".R") else [python]
    print(f"Running {script}", flush=True)
    subprocess.run([*interpreter, str(SCRIPTS / script), *(str(arg) for arg in args)], check=True)


def replay(data, release_data, out, external, rebuild_samples, rebuild_signatures, rebuild_literature):
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError("Choose a new or empty output directory.")
    enabled = {"default"}
    for status, mode in [(rebuild_samples, "rebuild_samples"),
                         (rebuild_signatures, "rebuild_signatures"),
                         (rebuild_literature, "rebuild_literature")]:
        if status:
            enabled.add(mode)
    verified = verify_inputs(release_data, external, enabled)
    frozen = data / "data/r11"
    reference = frozen / "reference"
    expected = frozen / "results"
    audit = out / "01_audit"
    result = out / "03_results"
    qa = out / "07_qa"
    gmt_stage = out / "00_download/files/MSigDB"
    for directory in (audit, result, qa, gmt_stage, out / "08_logs"):
        directory.mkdir(parents=True, exist_ok=True)
    write_tsv(qa / "R11_RELEASE_INPUT_PARITY.tsv", ("use", "path", "sha256", "status"), verified)

    for field in AUDIT_INPUTS:
        shutil.copyfile(reference / field, audit / field)
    members = "R11_SIGNATURE_MEMBERS.tsv.gz"
    selected_gmt = "R11_MSIGDB_MOUSE_MH_2025.1_SELECTED.gmt"
    shutil.copyfile(expected / members, result / members)
    shutil.copyfile(reference / selected_gmt, gmt_stage / selected_gmt)

    counts = frozen / "geo"
    if rebuild_samples:
        combined = out / "00_download/files"
        for counts_file in counts.glob("GSE*/*.gz"):
            dest = combined / counts_file.relative_to(counts)
            dest.parent.mkdir(exist_ok=True)
            shutil.copyfile(counts_file, dest)
        for accession in ("GSE121279", "GSE263830", "GSE121251"):
            soft = Path(accession) / f"{accession}_family.soft.gz"
            shutil.copyfile(external / soft, combined / soft)
        call("prepare_rnaseq_metadata.R", combined, out)
        sample_parity = [(field, *table_match(reference / field, audit / field, enforce_exact=True))
                         for field in AUDIT_INPUTS[:3]]
        write_tsv(qa / "R11_RELEASE_SAMPLE_RECONSTRUCTION.tsv",
                  ("table", "comparison", "rows", "columns", "max_abs_difference"), sample_parity)

    mgi_report = None
    if rebuild_signatures:
        mgi_report = external / "annotation/HOM_ProteinCoding.rpt"
        call("prepare_signature_membership.R", data, out, external / "MSigDB/msigdb.2025.1.Hs.rds",
             external / "MSigDB/mh.all.v2025.1.Mm.symbols.gmt", mgi_report)
        for field in ("R11_SIGNATURE_DEFINITIONS.tsv", "R11_SIGNATURE_IMMGEN_GROUPS.tsv"):
            verify_file(audit / field, sha256(reference / field))
        verify_file(result / members, sha256(expected / members))
        verify_file(gmt_stage / selected_gmt, sha256(reference / selected_gmt))

    call("fit_differential_expression.R", counts, out)
    call("verify_differential_expression.R", counts, out)
    call("signature_analysis.R", counts, out)
    verify_signatures(data, out, counts, mgi_report)
    literature_results = expected
    if rebuild_literature:
        call("curate_literature.py", data, out, external / CHEN_ARCHIVE)
        literature_pairs = [(expected / table, result / table) for table in FROZEN_LITERATURE]
        literature_pairs.append((reference / "R11_LITERATURE_SCREENING.tsv", audit / "R11_LITERATURE_SCREENING.tsv"))
        literature_parity = [(actual.name, *table_match(frozen, actual, enforce_exact=True))
                             for frozen, actual in literature_pairs]
        write_tsv(qa / "R11_RELEASE_LITERATURE_PARITY.tsv",
                  ("table", "comparison", "rows", "columns", "max_abs_difference"), literature_parity)
        literature_results = result
    literature = verify_literature(literature_results)
    write_tsv(qa / "R11_RELEASE_LITERATURE_SOURCE_QA.tsv", ("candidates", "evidence_rows", "cited_pmids", "status"),
              [(*literature, "FROZEN_CURATION_CHECKED")])
    call("integrate_candidate_evidence.py", data, out)

    releases = tsv(SCRIPTS / "expected_results.tsv")
    actual = sorted(path.name for path in result.iterdir() if path.name not in {members, *FROZEN_LITERATURE})
    if sorted(releases.table.tolist()) != actual:
        raise ValueError("Unexpected set of conditional-knockout analysis results")
    parity = []
    for row in releases.itertuples(index=False):
        reference_path = release_data / row.path
        verify_file(reference_path, row.sha256)
        stat = table_match(reference_path, result / row.table)
        parity.append((row.table, *stat, sha256(result / row.table)))
    write_tsv(qa / "R11_RELEASE_STATISTIC_PARITY.tsv",
              ("table", "comparison", "rows", "columns", "max_abs_difference", "replay_sha256"), parity)
    print(f"Conditional knockout: {len(parity)} newly computed result tables verified against the frozen release.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="extracted processed-data archive root")
    parser.add_argument("out", type=Path, help="new or empty output directory")
    parser.add_argument("external", nargs="?", type=Path, help="root of optional source snapshots specified in inputs.tsv")
    parser.add_argument("--rebuild-samples", action="store_true", help="reconstruct the released sample and ID maps")
    parser.add_argument("--rebuild-signatures", action="store_true", help="reconstruct the frozen reference membership")
    parser.add_argument("--rebuild-literature", action="store_true", help="reconstruct literature tables from frozen curation and verified publisher workbooks")
    args = parser.parse_args()
    if (args.rebuild_samples or args.rebuild_signatures or args.rebuild_literature) and args.external is None:
        parser.error("provide EXTERNAL_INPUT_ROOT for upstream reconstruction")
    data, out = args.data.resolve(), args.out.resolve()
    external = args.external.resolve() if args.external else None
    if (out.is_relative_to(data) or out.is_relative_to(SCRIPTS.parents[1])
            or (external and out.is_relative_to(external))):
        parser.error("Keep output outside the source and data directories")
    with replay_view(data) as compatibility:
        replay(compatibility, data, out, external, args.rebuild_samples,
               args.rebuild_signatures, args.rebuild_literature)


if __name__ == "__main__":
    main()
