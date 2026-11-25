#!/usr/bin/env python3
"""Independently verify Monaco numerical and sample arithmetic from the source TPM."""
import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median, pstdev


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def close(value, expected, label):
    if not math.isclose(float(value), expected, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError(f"{label}: source {expected} versus result {value}")


def independent_check(stage, source):
    r = stage / "03_results"
    a = stage / "01_audit"
    groups = {item["cell_type"]: item for item in rows(stage / "07_config/monaco_cell_types.tsv")}
    panel = {item["gene"]: item for item in rows(r / "R10_HUMAN_EXPRESSION_PANEL.tsv")}
    check = []
    sample = {item["matrix_label"]: item for item in rows(a / "R10_MONACO_SAMPLE_CROSSWALK.tsv")}
    axis_by_label = defaultdict(list)
    for label, item in sample.items():
        if item["included_primary"] == "True":
            axis_by_label[(item["donor"], item["cell_type"])].append(label)
    counts = Counter(kind for _, kind in axis_by_label)
    if len(sample) != 127 or len(axis_by_label) != 114 or counts["CD4_TE"] != 2:
        raise ValueError("Sample counts shifted")
    wanted = {item["ensembl_gene_id"]: gene for gene, item in panel.items()}
    if len(wanted) != 18:
        raise ValueError("Source panel must have 18 unique and mapped rows")
    aggregate = {}
    donor_raw = {}
    with gzip.open(source, "rt", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        sample_col = {label: index for index, label in enumerate(header) if index > 0}
        if set(sample_col) != set(sample):
            raise ValueError("Raw header labels differ from sample crosswalk")
        for line in reader:
            if line[0] not in wanted:
                continue
            gene = wanted[line[0]]
            on_unit = {(donor, kind): median(float(line[sample_col[label]]) for label in labels)
                       for (donor, kind), labels in axis_by_label.items()}
            donor_raw[gene] = on_unit
            levels = {kind: median(value for (d, t), value in on_unit.items() if t == kind)
                      for kind in groups}
            aggregate[gene] = levels
            row = panel[gene]
            close(row["n_cell_types"], len(levels), f"{gene} non-missing cell types")
            close(row["treg_tpm"], levels["Treg"], f"{gene} Treg TPM")
            close(row["median_type_tpm"], median(levels.values()), f"{gene} median TPM")
            max_tpm = max(levels.values())
            close(row["tau"], sum(1-value/max_tpm for value in levels.values())/(len(levels)-1), f"{gene} Tau")
            rank = 1 + sum(v > levels["Treg"] for v in levels.values())
            rank += 0.5 * sum(v == levels["Treg"] for k, v in levels.items() if k != "Treg")
            close(row["treg_rank_desc"], rank, f"{gene} rank")
            close(row["treg_rank_percentile"], (len(levels)-rank)/(len(levels)-1), f"{gene} percentile")
            for th, name in ((0.5, "ge0_5"), (1, "ge1"), (2, "ge2")):
                detected = sum(value >= th for value in levels.values())
                close(row[f"detected_{name}"], detected, f"{gene} detected at {th:g} TPM")
                close(row[f"breadth_{name}"], detected/len(levels), f"{gene} breadth at {th:g} TPM")
            comparison = mean(value for kind, value in levels.items()
                              if groups[kind]["conventional_primary"] == "1")
            close(row["tconv_equal_type_tpm"], comparison, f"{gene} Tconv TPM")
            close(row["treg_tconv_log2_ratio"], math.log2((1+levels["Treg"])/(1+comparison)),
                  f"{gene} log2 Treg/Tconv")
    if set(aggregate) != set(panel):
        raise ValueError("Panel not recovered from TPM file")
    check.append(("raw_to_panel", "PASS",
                  "18 gene-level source rows independently re-aggregated; TPM, threshold metrics, Tau, rank and CD4 contrast match"))

    lineages = defaultdict(list)
    for kind, item in groups.items():
        lineages[item["lineage"]].append(kind)
    lineage_tpm = {gene: {lineage: mean(levels[kind] for kind in members) for lineage, members in lineages.items()}
                   for gene, levels in aggregate.items()}
    for label, matrix, filename in (
        ("all_29", aggregate, "R10_MONACO_EXPRESSION_HEATMAP_SOURCE.tsv"),
        ("ten_lineages", lineage_tpm, "R10_MONACO_LINEAGE_HEATMAP_SOURCE.tsv"),
    ):
        expected_keys = {(gene, kind) for gene, levels in matrix.items() for kind in levels}
        actual_keys = set()
        for item in rows(r / filename):
            gene, kind = item["gene"], item["cell_type"]
            if (gene, kind) in actual_keys:
                raise ValueError("Duplicate heatmap key")
            actual_keys.add((gene, kind))
            close(item["median_tpm"], matrix[gene][kind], f"{label} TPM {gene} {kind}")
            case_log = math.log2(1+matrix[gene][kind])
            close(item["log2_tpm_plus_1"], case_log, f"{label} log TPM {gene} {kind}")
            vector = [math.log2(1+v) for v in matrix[gene].values()]
            close(item["row_z"], (case_log-mean(vector))/pstdev(vector), f"{label} z-score {gene} {kind}")
            if label == "all_29":
                close(item["donor_n"], counts[kind], f"{gene} {kind} donor denominator")
        if actual_keys != expected_keys:
            raise ValueError(f"Incomplete heatmap {label}")
        check.append((f"heatmap_{label}", "PASS",
                      f"All {len(actual_keys)} displayed values match raw-derived TPM, log values and row z-scores"))

    n_by_contrast = Counter()
    for item in rows(r / "R10_EXPRESSION_DONOR_CONTRASTS.tsv"):
        gene, donor, contrast = item["gene"], item["donor"], item["contrast_set"]
        n_by_contrast[contrast] += 1
        n_types = 6 if contrast == "primary_balanced_six" else 7
        type_flag = "conventional_primary" if n_types == 6 else "conventional_sensitivity"
        tconv = mean(donor_raw[gene][(donor, kind)] for kind in groups if groups[kind][type_flag] == "1")
        treg = donor_raw[gene][(donor, "Treg")]
        close(item["treg_tpm"], treg, f"{gene} {donor} within-donor Treg TPM")
        close(item["tconv_equal_type_tpm"], tconv, f"{gene} {donor} within-donor conventional CD4 TPM")
        close(item["log2_ratio"], math.log2((1+treg)/(1+tconv)), f"{gene} {donor} within-donor log2 ratio")
    if n_by_contrast != {"primary_balanced_six": 72, "all_seven_complete_donors": 36}:
        raise ValueError("Paired donor-row count mismatch")
    check.append(("donor_ratios", "PASS", "All 108 primary and sensitivity paired rows match raw-derived TPM"))

    supported = defaultdict(list)
    for item in rows(r / "R10_MONACO_MATCH_MEMBERS.tsv.gz"):
        if item["common_support"] == "True":
            supported[(item["contrast"], item["stratum"])].append(item)
    margins = defaultdict(Counter)
    by_comparison = defaultdict(list)
    for (contrast, stratum), items in supported.items():
        counts = Counter(item["asd"] == "True" for item in items)
        if not counts[True] or not counts[False]:
            raise ValueError("Weight assigned to a one-sided stratum")
        for item in items:
            expected_weight = 1 if item["asd"] == "True" else counts[True]/counts[False]
            close(item["weight"], expected_weight, f"{contrast} {stratum} control reweight")
        margins[contrast].update(counts)
        by_comparison[contrast].extend(item for item in items if item["asd"] == "True")
    draw_values = defaultdict(list)
    for item in rows(r / "R10_MONACO_MATCH_NULL_DRAWS.tsv.gz"):
        draw_values[(item["contrast"], item["endpoint"])].append(float(item["null_median"]))
    for item in rows(r / "R10_EXPRESSION_MATCHED_NULL.tsv"):
        contrast, endpoint = item["contrast"], item["endpoint"]
        close(item["supported_asd_n"], margins[contrast][True], f"{contrast} ASD margin")
        close(item["supported_non_asd_n"], margins[contrast][False], f"{contrast} comparator margin")
        observed = median(float(value[endpoint]) for value in by_comparison[contrast])
        close(item["observed_asd_median"], observed, f"{contrast} observed median")
        draws = draw_values[(contrast, endpoint)]
        center = median(draws)
        close(item["null_median"], center, f"{contrast} simulated median")
        if len(draws) != 20_000:
            raise ValueError("Conditional draw count mismatch")
        if item["p"] == "NA":
            if item["status"] == "DEGENERATE_NULL" and any(value != center for value in draws):
                raise ValueError("Nonconstant simulated endpoint marked degenerate")
        else:
            extreme = sum(abs(value-center) >= abs(observed-center) - 1e-14 for value in draws)
            p = (extreme+1)/(len(draws)+1)
            close(item["p"], p, f"{contrast} two-sided plus-one permutation P")
            close(item["q_bh"], min(2*p, 1), f"{contrast} BH across both declared endpoints")
    check.append(("matched_margins_and_tail", "PASS",
                  "Stratum group counts and weights, 80,000 stored null values, status, P and declared BH verified"))
    with (a / "R10_MONACO_INDEPENDENT_QA.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("check", "status", "detail"))
        writer.writerows(check)
    for name, status, detail in check:
        print(f"{status} {name} {detail}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--source-tpm", type=Path, required=True)
    args = parser.parse_args()
    independent_check(args.stage, args.source_tpm)


if __name__ == "__main__":
    main()
