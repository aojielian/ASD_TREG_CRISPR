#!/usr/bin/env python3
"""Extract or verify candidate-level facts from the published Zhu et al. study."""

import argparse
import csv
import hashlib
import math
import re
import tempfile
from collections import defaultdict
from itertools import product
from pathlib import Path


REPOSITORY = "https://github.com/emdann/GWT_perturbseq_analysis_2025"
VIRTUAL_CELLS = "https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq"
DOI = "10.1016/j.cell.2026.08.002"
PMID = "42664972"
PMCID = "PMC13614707"
CONDITIONS = ("Rest", "Stim8hr", "Stim48hr")
POLARIZATION = {"Non-polarized": "nonpolarized", "Th1-polarized": "Th1", "Th2-polarized": "Th2"}
CANDIDATES = (
    "CREBBP", "CUL3", "DYRK1A", "EP300", "FOXP1", "KDM3B", "KDM6B", "KMT5B", "PHIP",
    "SIN3A", "SMARCA4", "TBL1XR1", "TLE3", "ZMYND8", "CTCF", "PTEN", "SATB1",
)
FLOW_ENDPOINTS = {
    "IFNG": ("IGG", "IFNG", "log2FC of percent cytokine-positive cells"),
    "IL5": ("IGG", "IL5", "log2FC of percent cytokine-positive cells"),
    "T-bet": ("TFRest", "T-bet", "log2FC of transcription-factor median fluorescence intensity"),
    "GATA3": ("TFRest", "GATA-3", "log2FC of transcription-factor median fluorescence intensity"),
}
FLOW_FIELDS = {"IFNG": "IFNG", "IL5": "IL5", "T-bet": "Tbet", "GATA3": "GATA3"}
RNA_ENDPOINTS = ("Th1 signature", "Th2 signature")
GUIDE_FILE = "guide_kd_efficiency.suppl_table.csv"
DE_FILE = "DE_stats.suppl_table.csv"
ARRAY_FILE = "Th1Th2_validation_summary.suppl_table.csv"
FLOW_FILE = "Th1Th2_validation_flow.csv"
SCREEN_METADATA = "sample_metadata.suppl_table.csv"
VALIDATION_METADATA = "validation_donor_metadata.csv"
COMPARISON_FILE = "K562_comparison.suppl_table.csv"
DOCUMENTATION = "data_sharing_readme.md"
CONDITION_FILE = "ZHU_CANDIDATE_17_BY_CONDITION.tsv"
SUMMARY_FILE = "ZHU_CANDIDATE_17_SUMMARY.tsv"
ARRAY_DERIVED_FILE = "ZHU_ARRAYED_VALIDATION_CANDIDATES.tsv"
BIBLIOGRAPHY_FILE = "ZHU_BIBLIOGRAPHY.tsv"
SOURCE_MANIFEST = "ZHU_SOURCE_FILES.tsv"

# Immutable source bytes and release roles, in reconstruction order.
SOURCE_SPECS = {
    DE_FILE: (7276533, "7011918976889aaa5d802e3fff35635d6e4b6757bb1af17306059c0be07fbdc0",
              "suppl_tables/", "Perturbation-condition DE coverage, counts and quality"),
    GUIDE_FILE: (13319852, "75929bef741f87a944838fe3638fe3eaeef881c14a6653bab96f4b5df899aee3",
                 "suppl_tables/", "Counts of distinct guide-quality records and author flags"),
    SCREEN_METADATA: (2966, "766134d11dabb5d63388e00b4d687a809c0254eb3c19d77f439a6bf9dad7abd4",
                      "suppl_tables/", "Aggregate original Perturb-seq study donor coverage"),
    COMPARISON_FILE: (969086, "b88941f5016da23f83e65475d5868597f63d964b05fd87f5a5eea5b6aeb880ec",
                      "suppl_tables/", "Source inventory only. K562 is outside this candidate evidence"),
    DOCUMENTATION: (41585, "bc2c104ff98e5810d3344745ab1b8a1f6a51527449edaa083afcbfb63e37b8fb",
                    "", "Author definitions of source fields and assays"),
    ARRAY_FILE: (20872, "adc8dfacb69dfae0ce3947ad67dcb3ef20a700cccfef60547025cbe278d5d5d2",
                 "suppl_tables/", "Author-reported arrayed-validation contrast statistics"),
    FLOW_FILE: (81076, "a8eef31227984b33da11a64d22be07268020d85a118052c92733e10174f7ff27",
                "", "Count unique same-batch, condition and endpoint-matched donors"),
    VALIDATION_METADATA: (1185, "a9ab7f1d83690b93433aa2ac8228c2b6dd58ea5002c0dcb3bc15db1e2c88ad66",
                          "suppl_tables/", "Validation cohort metadata inventory only. Flow matching uses within-assay labels"),
}
DE_COLUMNS = {
    "de_n_cells_target": "n_cells_target",
    "de_n_guides_used": "n_guides",
    "de_n_guides_significant_on_target": "guide_n_signif_ontarget",
    "de_on_target_significant": "ontarget_significant",
    "de_low_target_expression": "low_target_gex",
    "de_neighbouring_gene_knockdown": "neighboring_gene_KD",
    "de_distal_offtarget_flag": "distal_offtarget_flag",
    "de_single_guide_estimate": "single_guide_estimate",
    "de_n_up_genes": "n_up_genes",
    "de_n_down_genes": "n_down_genes",
    "de_n_total_genes": "n_total_de_genes",
    "de_n_downstream_genes": "n_downstream",
    "de_guide_correlation_hit_genes": "guide_correlation_signif",
    "de_mean_donor_correlation_all_genes": "donor_correlation_all_mean",
    "de_mean_donor_correlation_hit_genes": "donor_correlation_hits_mean",
}
CONDITION_FIELDS = (
    "gene", "culture_condition", "guide_quality_n_rows", "guide_quality_n_significant",
    "guide_quality_n_high_confidence_no_effect", "de_summary_status", "source_de_row", "source_de_key",
    *DE_COLUMNS,
)
SUMMARY_FIELDS = (
    "gene", "guide_quality_conditions", "n_guide_quality_rows", "n_guide_quality_significant_rows",
    "n_guide_quality_high_confidence_no_effect_rows", "n_de_conditions_observed", "observed_de_conditions",
    "missing_de_conditions", "n_de_conditions_significant_on_target", "arrayed_validation_status",
)
ARRAY_FIELDS = (
    "gene", "condition", "modality", "endpoint", "effect_measure", "source_effect", "source_sem",
    "source_p_value", "source_adjusted_p_value", "source_adjustment", "source_adjusted_p_lt_0_05",
    "n_paired_flow_donors", "n_bulk_rna_donors", "source_file", "source_row", "source_effect_column",
    "source_sem_column", "source_p_column", "source_adjusted_p_column", "flow_count_source_file",
    "flow_count_condition", "flow_count_measurement", "flow_count_batch", "publication_panel", "doi", "pmid",
)
SOURCE_FIELDS = ("filename", "bytes", "sha256", "source_url", "role", "redistribution")
BIBLIOGRAPHY_FIELDS = (
    "authors", "title", "journal", "publication_date", "online_date", "volume", "issue", "pages",
    "doi", "pmid", "pmcid", "article_license", "article_license_url", "figure_map_url", "processed_dataset",
    "processed_dataset_url", "raw_geo", "raw_sra", "screen_donors", "cell_context",
)


def check(condition, error):
    if not condition:
        raise ValueError(error)


def source_rows(path, required, delimiter=","):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = reader.fieldnames
        check(fields is not None and len(fields) == len(set(fields)) and set(required) <= set(fields),
              f"{path.name} has an unexpected header")
        for row in reader:
            check(None not in row and all(value is not None for value in row.values()),
                  f"{path.name} has a malformed record at line {reader.line_num}")
            yield reader.line_num, row


def expect_sources(source_dir):
    for name, (expected_bytes, expected_sha, _, _) in SOURCE_SPECS.items():
        path = source_dir / name
        check(path.is_file(), f"Missing study source {name}")
        with path.open("rb") as handle:
            actual_sha = hashlib.file_digest(handle, "sha256").hexdigest()
        check(path.stat().st_size == expected_bytes and actual_sha == expected_sha,
              f"Source bytes differ from the frozen snapshot for {name}")


def read_candidates(candidate_master):
    candidates = tuple(row["gene"] for _, row in source_rows(candidate_master, ("gene",), "\t"))
    check(candidates == CANDIDATES, "The candidate master does not match the 17 frozen genes in order")
    return candidates


def bibliography():
    return [{
        "authors": "Zhu R, Dann E, Yan J, Reyes Retana J, Goto R, Guitche RC, Brixi L, Ota M, Hartman A, "
                   "Roth TL, Satpathy AT, Pritchard JK, Marson A",
        "title": "Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators "
                 "of T cell programs and human immune traits",
        "journal": "Cell", "publication_date": "2026-10-29", "online_date": "2026-08-28",
        "volume": "189", "issue": "22", "pages": "6988-7012.e17", "doi": DOI, "pmid": PMID,
        "pmcid": PMCID, "article_license": "CC BY-NC 4.0",
        "article_license_url": "https://creativecommons.org/licenses/by-nc/4.0/",
        "figure_map_url": f"{REPOSITORY}/blob/master/metadata/figure_map.md",
        "processed_dataset": "Primary Human CD4+ T Cell Perturb-seq, processed v1.0",
        "processed_dataset_url": VIRTUAL_CELLS, "raw_geo": "GSE314342", "raw_sra": "SRP643211",
        "screen_donors": "4", "cell_context": "Primary human conventional CD4+ T cells",
    }]


def source_manifest():
    return [{
        "filename": name, "bytes": str(size), "sha256": digest,
        "source_url": f"{REPOSITORY}/blob/master/metadata/{subdirectory}{name}",
        "role": role, "redistribution": "Author source file linked only",
    } for name, (size, digest, subdirectory, role) in SOURCE_SPECS.items()]


def index_de(source_dir, candidates):
    records = {}
    required = ("index", "target_contrast_gene_name", "target_contrast", "culture_condition", *DE_COLUMNS.values())
    for line, row in source_rows(source_dir / DE_FILE, required):
        target, condition = row["target_contrast_gene_name"], row["culture_condition"]
        if target in candidates:
            key = (target, condition)
            check(condition in CONDITIONS and key not in records, f"Non-unique Perturb-seq key {key}")
            records[key] = (line, row)
    return records


def index_guides(source_dir, candidates):
    records = defaultdict(dict)
    required = ("", "perturbed_gene_id", "culture_condition", "signif_knockdown", "high_confidence_no_effect_guides")
    for _, row in source_rows(source_dir / GUIDE_FILE, required):
        guide = row[""]
        target, separator, suffix = guide.rpartition("-")
        if target in candidates and separator and suffix.isdigit():
            key = (target, row["culture_condition"])
            check(key[1] in CONDITIONS and guide not in records[key], f"Non-unique guide-quality record for {key}")
            records[key][guide] = row
    return records


def build_condition_rows(source_dir, candidates):
    de = index_de(source_dir, candidates)
    guides = index_guides(source_dir, candidates)
    result = []
    for gene, condition in product(candidates, CONDITIONS):
        key = (gene, condition)
        line, de_row = de.get(key, ("", {}))
        guide_rows = tuple(guides[key].values())
        if de_row:
            mapped_ensembl_ids = {row["perturbed_gene_id"] for row in guide_rows if row["perturbed_gene_id"]}
            check(mapped_ensembl_ids <= {de_row["target_contrast"]}, f"Source target ID disagrees for {key}")
        result.append({
            "gene": gene, "culture_condition": condition,
            "guide_quality_n_rows": str(len(guide_rows)),
            "guide_quality_n_significant": str(sum(row["signif_knockdown"] == "True" for row in guide_rows)),
            "guide_quality_n_high_confidence_no_effect": str(sum(
                row["high_confidence_no_effect_guides"] == "True" for row in guide_rows
            )),
            "de_summary_status": "observed" if de_row else "not_assessed",
            "source_de_row": str(line), "source_de_key": de_row.get("index", ""),
            **{column: de_row.get(original, "") for column, original in DE_COLUMNS.items()},
        })
    return result


def study_donor_categories(source_dir):
    by_condition = defaultdict(set)
    for _, row in source_rows(source_dir / SCREEN_METADATA, ("donor_id", "culture_condition")):
        check(row["donor_id"] != "", "Missing biological replicate identifier in screen metadata")
        by_condition[row["culture_condition"]].add(row["donor_id"])
    check(set(by_condition) == set(CONDITIONS), "Source Perturb-seq donor metadata has unexpected conditions")
    common = set.intersection(*by_condition.values())
    check(len(common) == 4 and all(donors == common for donors in by_condition.values()),
          "Expected four Perturb-seq study donors represented in each condition")


def flow_counts(source_dir, candidates):
    # Count overlap inside the flow table, where both arms use the same identifier namespace.
    # The separate validation metadata does not supply a crosswalk for these flow labels.
    donors = defaultdict(set)
    fields = ("Donor", "perturbation", "condition", "measurement", "batch", "value")
    for _, row in source_rows(source_dir / FLOW_FILE, fields):
        donor = row["Donor"]
        check(donor != "", "A source flow record lacks an assay donor label")
        target = row["perturbation"]
        is_control = target in ("NTC", "NTC1", "NTC2")
        if (is_control or target in candidates) and row["value"] and math.isfinite(float(row["value"])):
            key = (row["batch"], row["condition"], row["measurement"], "NTC" if is_control else target)
            donors[key].add(donor)

    return donors


def source_columns(endpoint):
    if endpoint in FLOW_FIELDS:
        protein = FLOW_FIELDS[endpoint]
        return f"flow_{protein}_log2FC", "", f"flow_{protein}_pval", f"flow_{protein}_fdr"
    signature = endpoint.split()[0]
    return (
        f"bulkRNA_{signature}_mean_zscore", f"bulkRNA_{signature}_sem_zscore",
        f"bulkRNA_{signature}_pvalue", f"bulkRNA_{signature}_adj_pvalue",
    )


def build_arrayed_rows(source_dir, candidates):
    donors = flow_counts(source_dir, candidates)
    endpoint_names = (*FLOW_ENDPOINTS, *RNA_ENDPOINTS)
    measurement_fields = tuple(column for endpoint in endpoint_names for column in source_columns(endpoint) if column)
    required = ("target_name", "condition", "flow_batch", "bulkRNA_n_donors", *measurement_fields)
    selected = {}
    for line, row in source_rows(source_dir / ARRAY_FILE, required):
        if row["target_name"] in candidates:
            key = (row["target_name"], row["condition"])
            check(key not in selected, f"Duplicate arrayed-validation source key {key}")
            selected[key] = (line, row)

    check(set(selected) == {("ZMYND8", condition) for condition in POLARIZATION},
          "Unexpected candidate/condition coverage in arrayed validation")
    result = []
    for condition in POLARIZATION:
        line, row = selected[("ZMYND8", condition)]
        for endpoint in endpoint_names:
            effect, sem, nominal, adjusted = source_columns(endpoint)
            is_flow = endpoint in FLOW_ENDPOINTS
            source_flow_condition = source_flow_measurement = flow_batch = paired = rna = ""
            if is_flow:
                condition_prefix, source_flow_measurement, measure = FLOW_ENDPOINTS[endpoint]
                source_flow_condition = f"{condition_prefix}_{POLARIZATION[condition]}"
                flow_batch = row["flow_batch"]
                common_key = (flow_batch, source_flow_condition, source_flow_measurement)
                paired_donors = donors[(*common_key, "ZMYND8")] & donors[(*common_key, "NTC")]
                paired = str(len(paired_donors))
            else:
                measure = "mean per-gene DE z-score of the signature"
                rna = str(int(float(row["bulkRNA_n_donors"])))
            result.append({
                "gene": "ZMYND8", "condition": condition, "modality": "flow" if is_flow else "bulk_RNA_seq",
                "endpoint": endpoint, "effect_measure": measure, "source_effect": row[effect],
                "source_sem": row[sem] if sem else "", "source_p_value": row[nominal],
                "source_adjusted_p_value": row[adjusted], "source_adjustment": "BH",
                "source_adjusted_p_lt_0_05": str(float(row[adjusted]) < 0.05),
                "n_paired_flow_donors": paired, "n_bulk_rna_donors": rna,
                "source_file": ARRAY_FILE, "source_row": str(line), "source_effect_column": effect,
                "source_sem_column": sem, "source_p_column": nominal, "source_adjusted_p_column": adjusted,
                "flow_count_source_file": FLOW_FILE if is_flow else "", "flow_count_condition": source_flow_condition,
                "flow_count_measurement": source_flow_measurement, "flow_count_batch": flow_batch,
                "publication_panel": "Figure 4G", "doi": DOI, "pmid": PMID,
            })
    return result


def summarize(candidates, condition_rows, arrayed_rows):
    arrayed = {row["gene"] for row in arrayed_rows}
    output = []
    for gene in candidates:
        rows = [row for row in condition_rows if row["gene"] == gene]
        guide_conditions = [row["culture_condition"] for row in rows if int(row["guide_quality_n_rows"])]
        observed = [row for row in rows if row["de_summary_status"] == "observed"]
        observed_conditions = [row["culture_condition"] for row in observed]
        output.append({
            "gene": gene, "guide_quality_conditions": "|".join(guide_conditions),
            "n_guide_quality_rows": str(sum(int(row["guide_quality_n_rows"]) for row in rows)),
            "n_guide_quality_significant_rows": str(sum(int(row["guide_quality_n_significant"]) for row in rows)),
            "n_guide_quality_high_confidence_no_effect_rows": str(sum(
                int(row["guide_quality_n_high_confidence_no_effect"]) for row in rows
            )),
            "n_de_conditions_observed": str(len(observed)), "observed_de_conditions": "|".join(observed_conditions),
            "missing_de_conditions": "|".join(cond for cond in CONDITIONS if cond not in observed_conditions),
            "n_de_conditions_significant_on_target": str(sum(
                row["de_on_target_significant"] == "True" for row in observed
            )) if observed else "",
            "arrayed_validation_status": "observed" if gene in arrayed else "not_assessed",
        })
    return output


def read_derivative(directory, file, fields):
    path = directory / file
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        check(reader.fieldnames == list(fields), f"Unexpected fields in {file}")
        rows = list(reader)
    check(all(None not in row and all(value is not None for value in row.values()) for row in rows),
          f"Malformed record in {file}")
    return rows


def verify(directory, candidates):
    condition_rows = read_derivative(directory, CONDITION_FILE, CONDITION_FIELDS)
    summary = read_derivative(directory, SUMMARY_FILE, SUMMARY_FIELDS)
    arrayed = read_derivative(directory, ARRAY_DERIVED_FILE, ARRAY_FIELDS)
    check(read_derivative(directory, SOURCE_MANIFEST, SOURCE_FIELDS) == source_manifest(),
          "Study source provenance changed")
    check(read_derivative(directory, BIBLIOGRAPHY_FILE, BIBLIOGRAPHY_FIELDS) == bibliography(),
          "Study bibliography changed")
    keys = [(row["gene"], row["culture_condition"]) for row in condition_rows]
    check(keys == list(product(candidates, CONDITIONS)), "Expected complete ordered 17 by 3 candidate grid")
    missing = {("CTCF", "Rest"), ("CTCF", "Stim48hr"), *(("TLE3", condition) for condition in CONDITIONS)}
    for row, key in zip(condition_rows, keys):
        guides = int(row["guide_quality_n_rows"])
        significant = int(row["guide_quality_n_significant"])
        high_conf_no_effect = int(row["guide_quality_n_high_confidence_no_effect"])
        check(guides == 2 and min(significant, high_conf_no_effect) >= 0
              and significant + high_conf_no_effect <= guides, f"Guide-quality counts disagree for {key}")
        if key in missing:
            check(row["de_summary_status"] == "not_assessed" and not any(
                row[column] for column in ("source_de_row", "source_de_key", *DE_COLUMNS)
            ), f"Absent DE summary unexpectedly contains measurements for {key}")
        else:
            check(row["de_summary_status"] == "observed" and int(row["source_de_row"]) >= 2
                  and re.fullmatch(rf"ENSG[0-9]+_{key[1]}", row["source_de_key"]),
                  f"Expected observed DE source row for {key}")
            check(row["de_on_target_significant"] in ("True", "False")
                  and int(row["de_n_up_genes"]) + int(row["de_n_down_genes"]) == int(row["de_n_total_genes"])
                  and 0 <= int(row["de_n_downstream_genes"]) <= int(row["de_n_total_genes"]),
                  f"DE count or on-target flag disagrees for {key}")

    check(summary == summarize(candidates, condition_rows, arrayed), "Candidate summary does not match detail rows")
    check(sum(int(row["n_guide_quality_significant_rows"]) for row in summary) == 76
          and sum(int(row["n_guide_quality_high_confidence_no_effect_rows"]) for row in summary) == 25
          and sum(int(row["n_de_conditions_significant_on_target"] or 0) for row in summary) == 36,
          "Published study QC counts for the candidate subset changed")
    endpoints = (*FLOW_ENDPOINTS, *RNA_ENDPOINTS)
    check([(row["gene"], row["condition"], row["endpoint"]) for row in arrayed]
          == [("ZMYND8", condition, endpoint) for condition, endpoint in product(POLARIZATION, endpoints)],
          "Expected all 18 ZMYND8 arrayed measurements in source condition order")
    source_confirmed_flow_n = {
        "Non-polarized": {"IFNG": "6", "IL5": "6", "T-bet": "5", "GATA3": "5"},
        "Th1-polarized": {"IFNG": "5", "IL5": "5", "T-bet": "3", "GATA3": "3"},
        "Th2-polarized": {"IFNG": "6", "IL5": "6", "T-bet": "6", "GATA3": "6"},
    }
    significant_endpoints = []
    for row in arrayed:
        condition, endpoint = row["condition"], row["endpoint"]
        effect, sem, nominal, adjusted = source_columns(endpoint)
        check([row[column] for column in ("source_effect_column", "source_sem_column", "source_p_column",
                                         "source_adjusted_p_column")] == [effect, sem, nominal, adjusted],
              f"Source field names disagree for {condition} {endpoint}")
        check(row["source_file"] == ARRAY_FILE and row["source_row"] == str(58 + list(POLARIZATION).index(condition))
              and row["doi"] == DOI and row["pmid"] == PMID and row["publication_panel"] == "Figure 4G",
              f"Arrayed provenance disagrees for {condition} {endpoint}")
        p_value, q_value = (float(row[name]) for name in ("source_p_value", "source_adjusted_p_value"))
        check(math.isfinite(float(row["source_effect"])) and 0 <= p_value <= q_value <= 1
              and row["source_adjustment"] == "BH", f"Invalid source statistic for {condition} {endpoint}")
        check(row["source_adjusted_p_lt_0_05"] == str(q_value < 0.05),
              f"Incorrect adjusted threshold flag for {condition} {endpoint}")
        if q_value < 0.05:
            significant_endpoints.append((condition, endpoint))
        if endpoint in FLOW_ENDPOINTS:
            condition_prefix, measurement, measure = FLOW_ENDPOINTS[endpoint]
            check(row["modality"] == "flow" and row["effect_measure"] == measure and not row["source_sem"]
                  and not row["n_bulk_rna_donors"] and row["n_paired_flow_donors"] == source_confirmed_flow_n[condition][endpoint]
                  and row["flow_count_source_file"] == FLOW_FILE and row["flow_count_measurement"] == measurement
                  and row["flow_count_condition"] == f"{condition_prefix}_{POLARIZATION[condition]}"
                  and row["flow_count_batch"] == "Diff084", f"Flow coverage/units disagree for {condition} {endpoint}")
        else:
            check(row["modality"] == "bulk_RNA_seq" and row["n_bulk_rna_donors"] == "6"
                  and math.isfinite(float(row["source_sem"])) and float(row["source_sem"]) >= 0
                  and all(not row[column] for column in ("n_paired_flow_donors", "flow_count_source_file",
                                                        "flow_count_condition", "flow_count_measurement",
                                                        "flow_count_batch")),
                  f"Bulk RNA-seq denominator disagrees for {condition} {endpoint}")
    lock = {("Th1-polarized", "GATA3"): ("0.32603825853458884", "0.023282001894810764"),
            ("Non-polarized", "GATA3"): ("0.2309517635528857", "0.05529254658887842"),
            ("Non-polarized", "IL5"): ("0.6475743291586045", "0.0985982210971664"),
            ("Th1-polarized", "Th2 signature"): ("0.5438308735895281", "0.08731345704907861")}
    observed_lock = {
        (row["condition"], row["endpoint"]): (row["source_effect"], row["source_adjusted_p_value"])
        for row in arrayed if (row["condition"], row["endpoint"]) in lock
    }
    check(observed_lock == lock, "Source-reported arrayed validation effects or adjusted P values changed")
    check(significant_endpoints == [("Th1-polarized", "GATA3")], "Arrayed validation threshold profile changed")


def write_tsv(directory, file, fields, rows):
    with (directory / file).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def extract(source_dir, candidate_master, output):
    check(not output.exists() and not output.is_relative_to(source_dir),
          "Extraction requires a new output directory outside the study source directory")
    expect_sources(source_dir)
    candidates = read_candidates(candidate_master)
    study_donor_categories(source_dir)
    conditions = build_condition_rows(source_dir, candidates)
    arrayed = build_arrayed_rows(source_dir, candidates)
    summary = summarize(candidates, conditions, arrayed)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".zhu-", dir=output.parent) as workspace:
        directory = Path(workspace) / "processed"
        directory.mkdir()
        write_tsv(directory, CONDITION_FILE, CONDITION_FIELDS, conditions)
        write_tsv(directory, SUMMARY_FILE, SUMMARY_FIELDS, summary)
        write_tsv(directory, ARRAY_DERIVED_FILE, ARRAY_FIELDS, arrayed)
        write_tsv(directory, SOURCE_MANIFEST, SOURCE_FIELDS, source_manifest())
        write_tsv(directory, BIBLIOGRAPHY_FILE, BIBLIOGRAPHY_FIELDS, bibliography())
        verify(directory, candidates)
        directory.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="mode", required=True)
    validate = commands.add_parser("verify", help="Check released derivative structure, coverage and source anchors")
    validate.add_argument("--derived-dir", type=Path, required=True)
    validate.add_argument("--candidate-master", type=Path)
    rebuild = commands.add_parser("extract", help="Rebuild derivatives from eight explicitly supplied study sources")
    rebuild.add_argument("--source-dir", type=Path, required=True)
    rebuild.add_argument("--candidate-master", type=Path, required=True)
    rebuild.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.mode == "verify":
            candidates = read_candidates(args.candidate_master) if args.candidate_master else CANDIDATES
            verify(args.derived_dir, candidates)
        else:
            extract(args.source_dir.resolve(), args.candidate_master, args.out_dir.resolve())
    except (OSError, ValueError) as error:
        parser.exit(1, f"Zhu evidence check failed: {error}\n")
    print("Zhu evidence verified: 17 candidate summaries, 51 condition rows (46 with DE), "
          "18 reported ZMYND8 contrasts. Th1 GATA3 flow n=3.")


if __name__ == "__main__":
    main()
