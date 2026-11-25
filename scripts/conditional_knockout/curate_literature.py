"""Join R11 paper curation to NLM citation records and original Chen workbooks."""
from argparse import ArgumentParser
from collections import Counter
import csv
from math import isclose
from pathlib import Path
from shutil import copyfileobj
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from openpyxl import load_workbook
from openpyxl.utils.cell import column_index_from_string


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def published_gene_rows(profile, supplements):
    path = supplements / "2023-09-16093D-Supplementary Table 1.xlsx"
    workbook = load_workbook(path, read_only=True, data_only=True)
    workbook_rows = workbook["MAGECK gene"].iter_rows(values_only=True)
    fields = next(workbook_rows)
    indices = {int(row["chen_source_row"]): row for row in profile}
    matched = {}
    for source_row, values in enumerate(workbook_rows, 2):
        if source_row not in indices:
            continue
        frozen = indices[source_row]
        gene = frozen["gene"]
        record = dict(zip(fields, values))
        assert record["id"] == ("SUV420H1" if gene == "KMT5B" else gene), (gene, record)
        assert isclose(record["pos|fdr"], float(frozen["chen_pos_fdr"]), abs_tol=1e-12), gene
        assert isclose(record["neg|fdr"], float(frozen["chen_neg_fdr"]), abs_tol=1e-12), gene
        lfc = record["pos|lfc"] if record["pos|fdr"] < record["neg|fdr"] else record["neg|lfc"]
        assert isclose(lfc, float(frozen["chen_author_lfc"]), abs_tol=1e-12), (gene, lfc)
        assert isclose(-lfc, float(frozen["chen_signed_foxp3_bin_lfc"]), abs_tol=1e-12), gene
        matched[gene] = (record["id"], record["neg|fdr"], record["pos|fdr"], lfc)
    workbook.close()
    assert set(matched) == {row["gene"] for row in profile}
    return matched


def chen_record(template, gene, grade, context, perturbation, locator, protein, other, relation, limitation):
    return dict.fromkeys(template, "NA") | {
        "gene": gene, "pmid": "40140585", "evidence_grade": grade,
        "independence": "Same study as discovery CRISPR", "species": "Homo sapiens",
        "context": context, "perturbation": perturbation, "treg_specific_cko": "NO",
        "figure_or_source_row": locator, "foxp3_protein_or_frequency": protein,
        "direct_suppression": "NO_CANDIDATE_SPECIFIC_ASSAY", "other_readout": other,
        "directional_relation_to_chen": relation, "source_access": "PUBLISHED_ORIGINAL_SOURCE_FROZEN",
        "primary_url": "https://www.nature.com/articles/s41586-025-08795-5",
        "limitation": limitation,
    }


def original_human_records(profile, supplements, protein, prior, template):
    matched = published_gene_rows(profile, supplements)
    records = []
    for row in profile:
        gene, source_row = row["gene"], row["chen_source_row"]
        symbol, neg_fdr, pos_fdr, original_lfc = matched[gene]
        effect = -original_lfc
        records.append(chen_record(
            template, gene, "SAME_STUDY_SCREEN", "Cord-blood-derived human iTreg induction",
            "Pooled CRISPR-Cas9 loss of target gene",
            f"Fig. 1A–C; Extended Data Fig. 1C; Supplementary Table 1, MAGECK gene row {source_row}, {symbol}",
            f"{'Higher' if effect > 0 else 'Lower'} FOXP3 protein-bin association; standardized signed LFC {effect:.10g}",
            f"Author low-versus-high bin LFC {original_lfc:.10g}; negative FDR {neg_fdr:.10g}; positive FDR {pos_fdr:.10g}",
            "Discovery reference",
            "Guide selection by protein bins differs from per-cell expression change and direct suppression",
        ))
    workbook_path = supplements / "2023-09-16093D-Supplementary Table 2.xlsx"
    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    rows = list(workbook["ProteinScore_Targets"].iter_rows(min_row=1, max_row=7, values_only=True))
    genes = {row["gene"]: row for row in profile}
    for row in protein:
        gene, score = row["gene"], row["score"]
        if gene not in genes or row["protein"] != "FOXP3" or score == "NA":
            continue
        source_row = int(float(row["source_row"]))
        column = column_index_from_string(row["source_column"]) - 1
        assert rows[0][column] == gene, (gene, rows[0][column])
        assert rows[source_row - 1][0] == "intra-B0041-FOXP3"
        value = rows[source_row - 1][column]
        assert isclose(value, float(score), abs_tol=1e-12), (gene, value, score)
        reference_sign = float(genes[gene]["chen_signed_foxp3_bin_lfc"])
        relation = "Concordant score sign" if value * reference_sign > 0 else "Discordant score sign"
        records.append(chen_record(
            template, gene, "SAME_STUDY_PROTEIN_SCORE", "Human iTreg induction Perturb-icCITE-seq",
            "Targeted CRISPR-Cas9 loss of target gene",
            f"Fig. 2B–E; Supplementary Table 2, ProteinScore_Targets {row['source_column']}{source_row}",
            f"Author-published intracellular FOXP3 protein score {value:.12g}",
            "Formula-defined zero non-targeting baseline",
            relation,
            "Published descriptive score without gene-specific donor-level inference or suppressive-function assay",
        ))
    workbook.close()
    for row in prior:
        gene = row["gene"]
        if gene not in genes or row["study"] != "Chen 2025" or row["assay"] != "Arrayed iTreg FOXP3 flow cytometry":
            continue
        records.append(chen_record(
            template, gene, "SAME_STUDY_ARRAYED_FLOW", row["cell_state"], row["perturbation"],
            row["source_locator"], row["effect"], row["biological_unit_and_n"], row["direction"],
            row["limitation"],
        ))
    return records


def augmented_rows(curated, metadata, original):
    rows = []
    for row in curated + original:
        pmid = row["pmid"]
        record = metadata.get(pmid, {})
        assert pmid == "NA" or record, pmid
        rows.append({
            "gene": row["gene"], "pmid": pmid,
            "doi": record.get("doi", "NA") or "NA",
            "title": record.get("title", "NA"),
            "first_author": record.get("first_author", "NA"),
            "journal": record.get("journal", "NA"),
            "publication_date": record.get("pub_date", "NA"),
            "publication_types": record.get("publication_types", "NA"),
            **{key: value for key, value in row.items() if key not in ("gene", "pmid")},
            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if record else "NA",
            "doi_url": f"https://doi.org/{record['doi']}" if record.get("doi") else "NA",
        })
    return sorted(rows, key=lambda row: row["gene"])


CONTEXT_ONLY = {
    ("PHIP", "42682783"): "Alternative chemical abbreviation; no perturbation of PHIP/DCAF14 in Treg",
    ("PHIP", "33102969"): "PHIP refers to phage immunoprecipitation; no perturbation of PHIP/DCAF14 in Treg",
    ("DYRK1A", "32295708"): "Corrigendum linked to the reported haploinsufficiency cohort",
    ("TLE3", "41699408"): "Contextual association under drug exposure; no TLE3-specific perturbation of Treg",
    ("TLE3", "32275817"): "Broad drug and Foxp3-related context; no gene-specific mature Treg intervention",
    ("EP300", "31182547"): "Combined CBP/p300 inhibitor and human Treg differentiation context; no paralog-specific suppression test",
    ("EP300", "38094248"): "Acetyltransferase network context; isolated EP300 perturbation in mature Treg with direct suppression not established",
    ("EP300", "42763295"): "AHR/cofactor context; isolated EP300 intervention in mature Treg not established",
}
MATURE_DIRECT_GRADES = {"SINGLE_GENE_FUNCTION", "SINGLE_GENE_FUNCTION_MIXED", "SINGLE_GENE_RESTORATION_FUNCTION"}


def screen_records(retrieved, evidence):
    by_gene_pmid = {}
    for row in evidence:
        by_gene_pmid.setdefault((row["gene"], row["pmid"]), []).append(row)
    screened = []
    for row in retrieved:
        key = row["gene"], row["pmid"]
        papers = by_gene_pmid.get(key, [])
        publication_types = row["publication_types"].split("|")
        if papers:
            access = " | ".join(dict.fromkeys(paper["source_access"] for paper in papers))
            verified = any(paper["source_access"].startswith(("PRIMARY_", "PUBLISHED_ORIGINAL_")) for paper in papers)
            classification = "INCLUDED_PRIMARY_OR_CONTEXT" if verified else "RETAINED_LIMITED_OR_INDIRECT"
            reason = " | ".join(dict.fromkeys(paper["evidence_grade"] for paper in papers))
        elif key in CONTEXT_ONLY:
            classification, reason, access = "EXCLUDED_CONTEXT_ONLY", CONTEXT_ONLY[key], "TITLE_ABSTRACT"
        elif any(value in publication_types for value in ("Review", "Comment", "Editorial", "Published Erratum")):
            classification, reason, access = "EXCLUDED_SECONDARY_OR_ERRATUM", "No independent candidate-Treg perturbation experiment assessed from this record", "TITLE_INDEX"
        else:
            classification, reason, access = "EXCLUDED_AT_TITLE_ABSTRACT", "Screening identified no eligible candidate-specific Treg experiment; full paper was not individually adjudicated", "TITLE_ABSTRACT_SCREEN"
        screened.append(row | {"screening_status": classification, "screening_reason": reason,
                               "evidence_access_level": access})
    return screened


def summarize(profile, search, screened, evidence):
    rows = []
    for assay in profile:
        gene = assay["gene"]
        gene_search = next(row for row in search if row["gene"] == gene)
        classifications = Counter(row["screening_status"] for row in screened if row["gene"] == gene)
        studies = [row for row in evidence if row["gene"] == gene and row["independence"].startswith("Independent")]
        mature = [
            row for row in studies
            if row["evidence_grade"] in MATURE_DIRECT_GRADES
            and row["treg_specific_cko"].startswith("YES")
            and row["direct_suppression"].startswith("YES")
        ]
        other_single_target = [
            row for row in studies if row["evidence_grade"] in MATURE_DIRECT_GRADES
            and row["direct_suppression"].startswith("YES") and row not in mature
        ]
        if mature:
            status = "VERIFIED_TREG_CONDITIONAL_AND_DIRECT_FUNCTION_OR_RESTORATION"
        elif other_single_target:
            status = "VERIFIED_SINGLE_TARGET_ASSAY_WITH_BROADER_OR_INDUCTION_CONTEXT"
        elif any(row["source_access"].startswith("PRIMARY") for row in studies):
            status = "NO_VERIFIED_SINGLE_GENE_DIRECT_MATURE_TREG_SUPPRESSION; PRIMARY_CONTEXT_RETAINED"
        else:
            status = "NO_VERIFIED_INDEPENDENT_DIRECT_TREG_FUNCTION_IN_BOUNDED_SEARCH"
        rows.append({
            "gene": gene, "pubmed_hits": gene_search["count_reported"],
            "curated_or_context_in_pubmed_search": classifications["INCLUDED_PRIMARY_OR_CONTEXT"],
            "retained_limited_or_indirect_in_pubmed_search": classifications["RETAINED_LIMITED_OR_INDIRECT"],
            "excluded_in_pubmed_search": sum(value for key, value in classifications.items() if key.startswith("EXCLUDED")),
            "independent_paper_pmids": "|".join(dict.fromkeys(row["pmid"] for row in studies)),
            "treg_conditional_direct_assay_or_restoration_pmids": "|".join(dict.fromkeys(row["pmid"] for row in mature)),
            "other_single_target_assay_pmids": "|".join(dict.fromkeys(row["pmid"] for row in other_single_target)),
            "status": status,
            "chen_screen_direction": assay["chen_direction"],
            "umhoefer_treg": assay["treg_screen_status"],
            "umhoefer_tconv": assay["tconv_screen_status"],
        })
    return rows


def qa_checks(profile, search, evidence, retrieved, screened, summary):
    claims = []
    for gene in {row["gene"] for row in profile}:
        indexed = len([row for row in screened if row["gene"] == gene])
        expected = int(next(row for row in search if row["gene"] == gene)["count_reported"])
        assert indexed == expected, gene
        claims.append({"check": f"PubMed count reconciled for {gene}", "status": "PASS", "observed": indexed})
    assert len(profile) == len(search) == len(summary) == 14
    assert len(retrieved) == len(screened) == len({(row["gene"], row["pmid"]) for row in screened})
    direct_evidence = [
        row for row in evidence
        if row["evidence_grade"] in MATURE_DIRECT_GRADES and row["direct_suppression"].startswith("YES")
    ]
    assert all(row["doi"] != "NA" and "Fig." in row["figure_or_source_row"]
               and row["source_access"].startswith("PRIMARY") for row in direct_evidence)
    assert all(row["doi"] != "NA" for row in evidence if row["pmid"] != "NA")
    mature = [row["gene"] for row in summary if row["treg_conditional_direct_assay_or_restoration_pmids"]]
    assert set(mature) == {"CREBBP", "EP300", "FOXP1", "KDM6B", "SIN3A", "SMARCA4"}
    for check, observed in (
        ("Frozen candidate genes reconciled with literature table and search", len(summary)),
        ("Query title/abstract records individually logged", len(screened)),
        ("Distinct PMIDs from systematic query", len({row["pmid"] for row in screened})),
        ("Candidate single-gene function rows linked to exact figure, DOI and original source", len(direct_evidence)),
        ("Distinct candidates with mature Treg conditional intervention or restoration and direct function assay", ", ".join(mature)),
        ("Chen author signed bin comparisons verified from original workbook", sum(row["evidence_grade"] == "SAME_STUDY_SCREEN" for row in evidence)),
        ("Chen author Perturb-icCITE FOXP3 protein scores verified from original workbook", sum(row["evidence_grade"] == "SAME_STUDY_PROTEIN_SCORE" for row in evidence)),
    ):
        claims.append({"check": check, "status": "PASS", "observed": observed})
    return sorted(claims, key=lambda row: row["check"])


def rebuild(data, stage, chen_archive):
    audit, results, qa = [stage / directory for directory in ("01_audit", "03_results", "07_qa")]
    ref = data / "data/r11/reference"
    frozen_result = data / "data/r11/results"
    r10 = data / "data/r10/results"
    profile = read_tsv(r10 / "R10_14_GENE_FUNCTION_PROFILE.tsv")
    protein = read_tsv(r10 / "R10_PERTURB_PROTEIN.tsv")
    retrieved = read_tsv(ref / "R11_LITERATURE_RETRIEVED_CITATIONS.tsv")
    linked = read_tsv(ref / "R11_LITERATURE_LINKED_SOURCE_CITATIONS.tsv")
    search = read_tsv(frozen_result / "R11_LITERATURE_SEARCH_LOG.tsv")
    metadata = {row["pmid"]: row for row in retrieved + linked}
    curated = read_tsv(ref / "R11_LITERATURE_EVIDENCE_CURATION.tsv")
    prior = read_tsv(Path(__file__).resolve().with_name("prior_chen_arrayed_flow.tsv"))
    with TemporaryDirectory(prefix="asd_chen_source_") as temp, ZipFile(chen_archive) as archive:
        supplements = Path(temp)
        for number in (1, 2):
            name = f"2023-09-16093D-Supplementary Table {number}.xlsx"
            with archive.open(f"2023-09-16093D-s5/{name}") as original, (supplements / name).open("wb") as target:
                copyfileobj(original, target)
        original = original_human_records(profile, supplements, protein, prior, curated[0])
    evidence = augmented_rows(curated, metadata, original)
    screened = screen_records(retrieved, evidence)
    summary = summarize(profile, search, screened, evidence)
    checks = qa_checks(profile, search, evidence, retrieved, screened, summary)
    write_tsv(results / "R11_14_GENE_LITERATURE_AUDIT.tsv", evidence)
    write_tsv(audit / "R11_LITERATURE_SCREENING.tsv", screened)
    write_tsv(results / "R11_LITERATURE_GENE_SUMMARY.tsv", summary)
    write_tsv(qa / "R11_LITERATURE_QA.tsv", checks)
    print(f"{len(summary)} candidates, {len(evidence)} evidence rows, {len(screened)} title/abstract records, {len(checks)} QA checks")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("stage", type=Path)
    parser.add_argument("chen_archive", type=Path)
    args = parser.parse_args()
    rebuild(args.data, args.stage, args.chen_archive)


if __name__ == "__main__":
    main()
