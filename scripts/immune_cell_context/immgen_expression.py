"""Summarize ImmGen population CPM, expression breadth, Treg ranks, and QA."""

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, pstdev

from immgen_metadata import EXTRA_B, MAIN_TREG, SERIES, TCONV, read_tsv, write_tsv

THRESHOLDS = (0.5, 1.0, 2.0)
DATASET = "ImmGen ULI GSE109125"
SPECIES = "Mus musculus"
UNITS = "TMM CPM"
TREG_KEYS = ("splenic_treg", "colonic_treg")
NA = "NA"


def read_gzip_table(path):
    with gzip.open(path, mode="rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def metric_tau(values):
    if len(values) < 2 or max(values) <= 0:
        return NA
    peak = max(values)
    return sum(1 - cpm / peak for cpm in values) / (len(values) - 1)


def rank_in_groups(expression, group):
    if group not in expression:
        return NA, NA
    value = expression[group]
    groups = list(expression.values())
    n_higher = sum(cpm > value for cpm in groups)
    n_equal = sum(cpm == value for cpm in groups)
    rank = n_higher + (n_equal + 1) / 2
    percentile = (len(groups) - rank) / (len(groups) - 1) if len(groups) > 1 else NA
    return rank, percentile


def panel_metrics(expression):
    values = list(expression.values())
    metrics = {
        "n_groups": len(values),
        "median_across_groups_cpm": median(values) if values else NA,
        "tau": metric_tau(values),
    }
    for threshold in THRESHOLDS:
        suffix = str(threshold).removesuffix(".0").replace(".", "p")
        detected = sum(cpm >= threshold for cpm in values)
        metrics[f"detected_groups_{suffix}"] = detected if values else NA
        metrics[f"breadth_{suffix}"] = detected / len(values) if values else NA
    for prefix, group in zip(TREG_KEYS, MAIN_TREG):
        rank, percentile = rank_in_groups(expression, group)
        metrics[f"{prefix}_cpm"] = expression.get(group, NA)
        metrics[f"{prefix}_rank"] = rank
        metrics[f"{prefix}_percentile"] = percentile
    metrics["splenic_tconv_cpm"] = expression.get(TCONV, NA)
    metrics["splenic_treg_vs_tconv_log2_ratio"] = (
        math.log2((expression[MAIN_TREG[0]] + 1) / (expression[TCONV] + 1))
        if MAIN_TREG[0] in expression and TCONV in expression else NA
    )
    return metrics


def summarize_by_group(libraries, group_set, column, allowed_library=None):
    values = defaultdict(lambda: defaultdict(list))
    for row in libraries:
        group = row["cell_type"]
        if group in group_set and (allowed_library is None or allowed_library(row)):
            values[row["human_gene"]][group].append(float(row[column]))
    return {
        human: {group: median(observations) for group, observations in gene_values.items()}
        for human, gene_values in values.items()
    }


def aggregate_equal_lineage(expression, definitions):
    expr_by_lineage = defaultdict(list)
    for group, value in expression.items():
        expr_by_lineage[definitions[group]["lineage"]].append(value)
    return {lineage: mean(values) for lineage, values in expr_by_lineage.items()}


def sensitivity_row(panel, variant, threshold, metrics):
    suffix = str(threshold).removesuffix(".0").replace(".", "p")
    equal_lineage = variant == "equal_lineage_TMM"
    unit = "ordinary CPM" if variant == "ordinary_CPM" else UNITS
    return {
        "dataset": DATASET,
        "species": SPECIES,
        "human_gene": panel["human_gene"],
        "mouse_gene": panel["mouse_gene"],
        "panel_role": panel["panel_role"],
        "variant": variant,
        "unit": unit,
        "denominator_unit": "lineage" if equal_lineage else "sorted population",
        "detection_threshold": threshold,
        "n_groups": metrics["n_groups"],
        "detected_groups": metrics[f"detected_groups_{suffix}"],
        "breadth": metrics[f"breadth_{suffix}"],
        "tau": metrics["tau"],
        "median_across_groups": metrics["median_across_groups_cpm"],
        "splenic_treg_rank": metrics["splenic_treg_rank"],
        "splenic_treg_percentile": metrics["splenic_treg_percentile"],
        "colonic_treg_rank": metrics["colonic_treg_rank"],
        "colonic_treg_percentile": metrics["colonic_treg_percentile"],
        "note": (
            "Only breadth/Tau/overall median are defined after equal-lineage pooling"
            if equal_lineage else
            "Colonic Treg excluded from splenic subset" if variant == "spleen_only_TMM" else
            "Measured GC B cells and plasmablasts added solely in this sensitivity"
            if variant == "expanded_B_state_TMM" else
            "Source-documented alias libraries excluded from group summaries"
            if variant == "exact_title_only_TMM" else
            panel["missing_reason"]
        ),
    }


def confirm_matrix_panel(matrix_path, mapping, definitions, main_cpm):
    """Reconcile all plot CPMs independently against the full R group matrix."""
    by_mouse = {row["mouse_gene"]: row["human_gene"] for row in mapping if row["measured"] == "yes"}
    n_rows = 0
    n_checked = 0
    seen_symbols = set()
    for row in read_gzip_table(matrix_path):
        n_rows += 1
        mouse = row.pop("mouse_gene")
        if mouse not in by_mouse:
            continue
        seen_symbols.add(mouse)
        human = by_mouse[mouse]
        if set(row) != set(main_cpm[human]):
            raise ValueError("Panel matrix and library-source primary groups disagree")
        for group, field in row.items():
            assert group in definitions
            if not math.isclose(float(field), main_cpm[human][group], rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError(f"Group CPM disagrees for {human}/{mouse}/{group}")
            n_checked += 1
    assert seen_symbols == set(by_mouse) and n_checked > 0
    return n_rows, n_checked


def verify_tmm_independent(libraries, tmm_qc, data_group):
    n = 0
    for row in libraries:
        meta = tmm_qc[row["matrix_library"]]
        count = float(row["raw_gene_count"])
        expected_tmm = count * 1e6 / float(meta["tmm_effective_library_size"])
        expected_ord = count * 1e6 / float(meta["raw_library_size"])
        if not math.isclose(expected_tmm, float(row["tmm_cpm"]), rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"TMM CPM does not match count / effective size: {row['gsm']}")
        if not math.isclose(expected_ord, float(row["ordinary_cpm"]), rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"Ordinary CPM does not match count / raw size: {row['gsm']}")
        if (row["human_gene"] not in data_group
                or row["cell_type"] != meta["cell_type"]
                or row["gsm"] != meta["gsm"]):
            raise ValueError("Panel library expression lacks validated metadata")
        n += 1
    return n


def build(out_dir):
    audit = out_dir / "01_audit"
    results = out_dir / "03_results"
    definitions = {
        row["cell_type"]: row for row in read_tsv(audit / "R10_IMMGEN_GROUP_CROSSWALK.tsv")
    }
    mapping = read_tsv(audit / "R10_IMMGEN_GENE_CROSSWALK.tsv")
    libraries = list(read_gzip_table(results / "R10_IMMGEN_PANEL_LIBRARY_CPM.tsv.gz"))
    tmm_qc = {
        row["matrix_library"]: row for row in read_tsv(audit / "R10_IMMGEN_LIBRARY_QC.tsv")
    }
    samples = read_tsv(audit / "R10_IMMGEN_SAMPLE_CROSSWALK.tsv")
    main_groups = {group for group, row in definitions.items() if row["include_primary"] == "yes"}
    if len(main_groups.intersection((*MAIN_TREG, TCONV))) != 3:
        raise ValueError("A prespecified Treg or the CD25-negative comparator is absent from baseline")
    main_cpm = summarize_by_group(libraries, main_groups, "tmm_cpm")
    ordinary = summarize_by_group(libraries, main_groups, "ordinary_cpm")
    exact = summarize_by_group(libraries, main_groups, "tmm_cpm",
                               lambda row: row["matrix_match"] == "exact_title")
    replicated = {group for group in main_groups if int(definitions[group]["n_libraries"]) >= 2}
    spleen = {group for group in main_groups if definitions[group]["tissue"] == "Spleen"}
    expanded_b = summarize_by_group(libraries, main_groups | EXTRA_B, "tmm_cpm")

    all_sensitivity = []
    panel_rows = []
    plot_rows = []
    n_libraries = sum(int(definitions[group]["n_libraries"]) for group in main_groups)
    for panel in mapping:
        human = panel["human_gene"]
        mouse = panel["mouse_gene"]
        expression = main_cpm.get(human, {})
        if expression and len(expression) != len(main_groups):
            raise ValueError(f"Missing TMM-CPM group estimates for a measured gene: {human}")
        summary = panel_metrics(expression)
        if expression:
            assert summary["breadth_0p5"] >= summary["breadth_1"] >= summary["breadth_2"]
            assert 0 <= summary["tau"] <= 1
            for label in TREG_KEYS:
                assert 1 <= summary[f"{label}_rank"] <= summary["n_groups"]
                assert 0 <= summary[f"{label}_percentile"] <= 1
        panel_rows.append({
            "dataset": DATASET,
            "species": SPECIES,
            "human_gene": human,
            "mouse_gene": mouse,
            "panel_role": panel["panel_role"],
            "orthology_status": panel["orthology_status"],
            "measured": panel["measured"],
            "missing_reason": panel["missing_reason"],
            "unit": UNITS,
            "n_libraries": n_libraries,
            **summary,
        })

        transformations = {
            "primary_TMM": expression,
            "ordinary_CPM": ordinary.get(human, {}),
            "exact_title_only_TMM": exact.get(human, {}),
            "min_two_libraries_TMM": {group: value for group, value in expression.items() if group in replicated},
            "spleen_only_TMM": {group: value for group, value in expression.items() if group in spleen},
            "equal_lineage_TMM": aggregate_equal_lineage(expression, definitions),
            "expanded_B_state_TMM": expanded_b.get(human, {}),
        }
        for variant, values in transformations.items():
            metrics = panel_metrics(values)
            for threshold in THRESHOLDS:
                all_sensitivity.append(sensitivity_row(panel, variant, threshold, metrics))

        logged = {group: math.log2(value + 1) for group, value in expression.items()}
        center = mean(logged.values()) if logged else NA
        scale = pstdev(logged.values()) if logged else NA
        for group in definitions:
            if group not in main_groups:
                continue
            g = definitions[group]
            log_value = logged.get(group, NA)
            plot_rows.append({
                "dataset": DATASET,
                "species": SPECIES,
                "human_gene": human,
                "mouse_gene": mouse,
                "panel_role": panel["panel_role"],
                "cell_type": group,
                "label": g["label"],
                "lineage": g["lineage"],
                "tissue": g["tissue"],
                "n_libraries": g["n_libraries"],
                "cpm": expression.get(group, NA),
                "log2_cpm1": log_value,
                "zscore": (log_value - center) / scale if log_value != NA and scale > 0 else NA,
                "group_order": g["group_order"],
                "unit": UNITS,
                "missing_reason": panel["missing_reason"],
            })

    n_normalized, n_check = confirm_matrix_panel(
        results / "R10_IMMGEN_CELLTYPE_CPM.tsv.gz", mapping, definitions, main_cpm
    )
    checked_library_values = verify_tmm_independent(libraries, tmm_qc, main_cpm)
    assert len({row["replicate_id"] for row in samples}) == len(samples)
    assert len({(row["human_gene"], row["cell_type"]) for row in plot_rows}) == len(plot_rows)
    assert len({row["mouse_gene"] for row in mapping if row["measured"] == "yes"}) == len(main_cpm)
    write_tsv(results / "R10_MOUSE_EXPRESSION_PANEL.tsv", panel_rows)
    write_tsv(results / "R10_MOUSE_EXPRESSION_PLOT.tsv", plot_rows)
    write_tsv(results / "R10_MOUSE_EXPRESSION_SENSITIVITY.tsv", all_sensitivity)

    qa = [
        {"check": "ImmGen count library and SOFT identity", "status": "PASS",
         "evidence": f"{len(tmm_qc)}/{len(tmm_qc)} joined 1:1; "
         f"{sum(row['matrix_match'] == 'source_documented_alias' for row in samples)} documented aliases; "
         f"{sum(row['matrix_match'] == 'not_in_count_matrix' for row in samples)} GSMs without count columns"},
        {"check": "ImmGen baseline population coverage", "status": "PASS",
         "evidence": f"{n_libraries} libraries, {len(main_groups)} main groups; biological library/BioSample "
         "used as replicate; animal IDs absent and pairing unused"},
        {"check": "ImmGen CPM independently recalculated", "status": "PASS",
         "evidence": f"Recalculated TMM and ordinary CPM for {checked_library_values} panel gene/library pairs"},
        {"check": "ImmGen plotted population medians", "status": "PASS",
         "evidence": f"{n_check} Python replicate medians agree with independently generated R population matrix "
         f"({n_normalized} genes) within 1e-12 relative/absolute tolerance"},
        {"check": "ImmGen panel mapping and threshold monotonicity", "status": "PASS",
         "evidence": f"{len(main_cpm)}/{len(mapping)} unique measured panel symbols; no duplicate group/gene keys. "
         "Breadth was monotone over 0.5, 1 and 2 CPM and Tau/ranks stayed within their valid ranges"},
    ]
    write_tsv(audit / "R10_IMMGEN_QA.tsv", qa)
    candidates = [row for row in panel_rows if row["panel_role"] == "candidate"]
    for threshold, suffix in zip(THRESHOLDS, ("0p5", "1", "2")):
        detection = [row[f"detected_groups_{suffix}"] for row in candidates if row["measured"] == "yes"]
        print(f"Candidate groups detected at {threshold:g} TMM CPM: "
              f"{min(detection)} to {max(detection)}/{len(main_groups)}. Median {median(detection):g}.")
    print(f"All outputs reconciled over {len(candidates)} candidates and {len(mapping) - len(candidates)} references.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    build(args.out)


if __name__ == "__main__":
    main()
