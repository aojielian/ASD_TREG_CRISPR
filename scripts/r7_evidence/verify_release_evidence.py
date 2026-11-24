#!/usr/bin/env python3
"""Reconcile the R7 figure sources, LoF sensitivity, and candidate evidence.

This verifier uses only the versioned processed-data directory from the Zenodo
archive. The exact LoF calculation uses the original conditional hypergeometric
null. Figure and evidence checks reuse the selections and aggregation rules of
the R7 source verifier, the R1 mechanism figure, and the R2 evidence builder.
"""

import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

CORE = {"CHEN": "CHEN_ITREG_NTC", "UMHOEFER_TREG": "CHEN_HTREG"}
SCREENS = (*CORE, "UMHOEFER_TCONV")
OLD_REFERENCES = {
    "unadjusted": "FROZEN_MEASURED",
    "class_conditioned": "FROZEN_MEASURED",
    "strict_joint": "FROZEN_GSE90600",
}
SCOPES = ("ASD_named", "ASD_feature_in_disorder", "NDD_named", "NDD_feature_in_disorder")
ASSAYED = {"DETECTABLE_HIT", "DETECTABLE_NONHIT"}
DIRECTIONS = {
    "SUPPORTS_FOXP3": "KO_REDUCES_FOXP3",
    "RESTRAINS_FOXP3": "KO_INCREASES_FOXP3",
    "NOT_SIGNIFICANT": "NO_SIGNIFICANT_KO_DIRECTION",
}
REFERENCES = {
    "Chen 2025": ("40140585", "10.1038/s41586-025-08795-5", "C"),
    "Umhoefer 2026": ("41237776", "10.1016/j.immuni.2025.10.020", "U"),
    "Schumann 2020": ("32989329", "10.1038/s41590-020-0784-4", "S"),
    "Kolb 2026": ("42286212", "10.1038/s44319-026-00812-6", "K"),
}
FOLLOW_UPS = {"SAME_STUDY_FOLLOW_UP", "INDEPENDENT_STUDY_FOLLOW_UP"}
AFFIRMATIVE = {"true", "yes", "1"}
NEGATIVE = {"false", "no", "0"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def truth(value):
    return value.lower() in AFFIRMATIVE


def equal_field(value, expected):
    if value == expected:
        return True
    if value.lower() in (AFFIRMATIVE | NEGATIVE) and expected.lower() in (AFFIRMATIVE | NEGATIVE):
        return truth(value) == truth(expected)
    try:
        return math.isclose(float(value), float(expected), rel_tol=1e-11, abs_tol=1e-13)
    except ValueError:
        return False


def compare(actual, expected, fields, context):
    for field in fields:
        require(field in actual and field in expected, f"{context}: missing field {field}")
        require(equal_field(actual[field], expected[field]),
                f"{context}: {field} differs ({actual[field]!r}, expected {expected[field]!r})")


def index(rows, fields):
    result = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        require(key not in result, f"Duplicate {fields}: {key}")
        result[key] = row
    return result


def select(rows, **criteria):
    found = [row for row in rows if all(equal_field(row[k], v) for k, v in criteria.items())]
    require(len(found) == 1, f"Expected one row for {criteria}, found {len(found)}")
    return found[0]


def joined(values):
    return " | ".join(dict.fromkeys(value for value in values if value))


class EvidenceAudit:
    def __init__(self, root):
        self.root = root.resolve(strict=True)
        self.checks = []

    def read(self, relative_path):
        path = (self.root / relative_path).resolve(strict=True)
        require(path.is_relative_to(self.root), f"Input escapes the released data root: {relative_path}")
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            header = reader.fieldnames or []
            require(header and "" not in header and len(header) == len(set(header)),
                    f"Empty or duplicate header in {relative_path}")
            records = []
            for row in reader:
                require(None not in row and None not in row.values(),
                        f"Unexpected column count in {relative_path}:{reader.line_num}")
                records.append(row)
        require(records, f"No data in {relative_path}")
        return records

    def record(self, group, check, evidence, status="PASS"):
        self.checks.append({"group": group, "check": check, "status": status, "evidence": evidence})


def figure_three(audit, inference):
    source = audit.read("figure_sources/Figure_3d_source.tsv")
    models = audit.read("results/LOEUF_ADJUSTED_MODELS.tsv")
    bridge = audit.read("results/R6_TO_R2_MODEL_BRIDGE.tsv")
    require(Counter(r["block"] for r in source) == {
        "enrichment": 12, "historical_regression": 2, "regression": 8, "regression_bridge": 4,
    }, "Figure 3d has unexpected blocks or row counts")
    inferential_fields = (
        "screen", "reference", "provenance", "threshold", "units", "analysis", "set",
        "total_n", "asd_hit", "asd_nonhit", "other_hit", "other_nonhit", "asd_n", "other_n",
        "original_asd_n", "informative_asd_n", "expected_hits", "effect_type", "effect",
        "ci_low", "ci_high", "p", "adjusted_p", "status", "original_strict_coverage_pass",
    )
    for screen, reference in CORE.items():
        for analysis, old_reference in OLD_REFERENCES.items():
            for context, ref, threshold in (
                ("original", old_reference, ""), ("CPM1", reference, "1"),
            ):
                expected = select(inference, screen=screen, reference=ref, threshold=threshold,
                                  set="SFARI_SCORE1", analysis=analysis)
                actual = select(source, screen=screen, block="enrichment", analysis=analysis,
                                display_context=context)
                label = f"{screen} {context} {analysis}"
                compare(actual, expected, inferential_fields, label)
                audit.record("Figure 3d", label,
                             f"{actual['effect_type']} {actual['effect']}, ASD {actual['asd_hit']}/{actual['asd_n']}")
        cohort_margins = set()
        for model in ("M0", "M1", "M2", "M3"):
            expected = select(models, screen=screen, reference=reference, threshold="1", units="CPM", model=model)
            actual = select(source, block="regression", screen=screen, model=model)
            label = f"{screen} CPM1 {model}"
            fields = (
                "screen", "reference", "threshold", "units", "model", "formula", "n_genes", "asd_n", "asd_hits",
                "other_hits", "other_n", "asd_nonhits", "other_nonhits", "OR", "CI_low", "CI_high",
                "P", "holm_p_m2", "converged", "interval_status", "estimator",
            )
            compare(actual, expected, fields, label)
            require(truth(actual["converged"]) and actual["interval_status"] == "ESTIMABLE",
                    f"{label}: an estimable interval was expected")
            n, an, ah, ac, cn, ch, cc = (int(actual[k]) for k in
                ("n_genes", "asd_n", "asd_hits", "asd_nonhits", "other_n", "other_hits", "other_nonhits"))
            require(n == an + cn and an == ah + ac and cn == ch + cc,
                    f"{label}: four-cell margins do not add up")
            point, low, high = (float(actual[k]) for k in ("OR", "CI_low", "CI_high"))
            require(0 < low < point < high, f"{label}: invalid interval")
            require(bool(actual["holm_p_m2"]) == (model == "M2"), f"{label}: correction family changed")
            cohort_margins.add((n, an, ah, cn, ch))
            audit.record("Figure 3d", label,
                         f"OR {actual['OR']}, 95% CI [{actual['CI_low']}, {actual['CI_high']}], {ah}/{an} ASD hits")
        require(len(cohort_margins) == 1, f"{screen}: CPM1 M0–M3 do not share the same margins")
        for era in ("OLD", "NEW"):
            analysis = f"CPM1_SHARED_COHORT_{era}_EXPRESSION"
            label = f"{screen} no-library expression bridge ({era})"
            expected = select(bridge, screen=screen, analysis=analysis)
            actual = select(source, block="regression_bridge", screen=screen, analysis=analysis)
            compare(actual, expected,
                    ("screen", "analysis", "n_genes", "asd_n", "OR", "CI_low", "CI_high",
                     "P", "interval_status", "wald_p_diagnostic"), label)
            audit.record("Figure 3d", label, f"OR {actual['OR']}, {actual['n_genes']} genes")
        expected = select(inference, screen=screen, reference="FROZEN_PREFLIGHT",
                          set="SFARI_SCORE1", analysis="preflight_M1")
        actual = select(source, block="historical_regression", screen=screen)
        label = f"{screen} historical adjusted estimate"
        compare(actual, expected, inferential_fields, label)
        if screen == "UMHOEFER_TREG":
            require(actual["status"] == "UNSTABLE_NOT_INFERENTIAL" and not (actual["ci_low"] or actual["ci_high"]),
                    "Historical Treg interval must remain unstable and blank")
        audit.record("Figure 3d", label, f"OR {actual['effect']}, {actual['status']}")


def mechanism(audit, background):
    master = audit.read("mechanism/ALLELIC_MECHANISM_MASTER.tsv")
    indexed = index(master, ("approved_symbol",))
    members = audit.read("figure_sources/Supplementary_Figure_S6_membership_source.tsv")
    plotted = index(members, ("approved_symbol",))
    require(len(master) == 245 and set(indexed) == set(plotted), "Incomplete or duplicated S6 membership")
    require(len(index(master, ("hgnc_id",))) == 245 and all(r["hgnc_id"] for r in master), "HGNC mapping is incomplete")
    fields = ("sfari_raw_symbol", "approved_symbol", "hgnc_id", "sfari_score",
              "mechanism_class", "lof_eligible", "lof_candidate_supports", "lof_candidate_scopes")
    for key, r in plotted.items():
        compare(r, indexed[key], fields, f"S6 membership {key[0]}")
        require(r["sfari_score"] == "1", f"{key[0]} is not a frozen SFARI score-1 gene")
    for screen in SCREENS:
        screen_genes = {gene: row for (sc, gene), row in background.items()
                        if sc == screen and truth(row["sfari_score1"])}
        require(set(screen_genes) == {k[0] for k in indexed}, f"{screen}: frozen SFARI IDs differ")
        for key, mr in indexed.items():
            row = screen_genes[key[0]]
            if "sfari_raw_symbol" in row:
                require(row["sfari_raw_symbol"] == mr["sfari_raw_symbol"],
                        f"{screen}/{key[0]}: SFARI symbol differs")
            if truth(row["measured"]):
                require(row["hgnc_id"] == mr["hgnc_id"], f"{screen}/{key[0]}: measured HGNC ID differs")
    classes = Counter(r["mechanism_class"] for r in master)
    require(classes == {"LoF-compatible": 130, "unknown": 99, "GoF": 8, "dominant negative": 4, "mixed": 4},
            f"S6 mechanism-class counts differ: {classes}")
    changes = index(audit.read("mechanism/LOF_ELIGIBILITY_CHANGE_LOG.tsv"), ("sfari_raw_symbol",))
    expect_dropped = {"GRIA2", "NLGN3", "PTEN", "WDFY3", "ZNF292"}
    require({k[0] for k in changes} == expect_dropped, "Primary-source exclusions differ")
    versions = Counter()
    scope_by_gene = {}
    support = Counter()
    disease_records = audit.read("mechanism/ALLELIC_MECHANISM_DISEASE_RECORDS.tsv")
    by_gene = defaultdict(list)
    for row in disease_records:
        by_gene[row["approved_symbol"]].append(row)
    require(set(by_gene) == {k[0] for k in indexed}, "Disease-record coverage differs")
    linked = [r for r in disease_records if r["g2p_id"]]
    require(len(disease_records) == 294 and len(linked) == 246, "Disease-record count differs")
    require(len(index(linked, ("g2p_id",))) == 246 and len({r["approved_symbol"] for r in linked}) == 197,
            "Disease-record identifiers or linked gene count differ")
    require(Counter(r["molecular_mechanism_support"] for r in linked) == {"inferred": 237, "evidence": 9},
            "G2P record-level support counts differ")
    valid_support = {"G2P_inferred_absent_product", "G2P_evidence"}
    for row in master:
        gene = row["approved_symbol"]
        records = by_gene[gene]
        mapped = [r for r in records if r["g2p_id"]]
        require(len(mapped) == int(row["g2p_record_n"]), f"{gene}: G2P count differs")
        require((bool(mapped) and len(mapped) == len(records)) or (not mapped and len(records) == 1),
                f"{gene}: absent-record handling differs")
        for record in records:
            compare(record, row, ("approved_symbol", "hgnc_id", "mechanism_class", "lof_eligible"),
                    f"{gene}/{record['g2p_id'] or 'no record'}")
        change = changes.get((row["sfari_raw_symbol"],))
        versions["v1"] += truth(row["lof_eligible_snapshot_v1"])
        versions["v2"] += truth(row["lof_eligible_v2_intermediate"])
        eligible = truth(row["lof_eligible"])
        if change:
            require(truth(change["lof_eligible_v1"]) and not truth(change["lof_eligible_v3_final"]) and not eligible,
                    f"{gene}: documented exclusion differs")
            require(truth(change["lof_eligible_v2_intermediate"]) == truth(row["lof_eligible_v2_intermediate"]),
                    f"{gene}: intermediate exclusion differs")
        if eligible:
            qualifying = {r["disease_scope"] for r in mapped if truth(r["lof_eligible_row"])}
            strongest = next((scope for scope in SCOPES if scope in qualifying), "")
            require(strongest and strongest == row["strongest_qualifying_lof_scope"],
                    f"{gene}: qualifying disease scope differs")
            require(row["lof_candidate_supports"] in valid_support,
                    f"{gene}: unknown eligible-gene evidence status")
            scope_by_gene[gene] = strongest
            support[row["lof_candidate_supports"]] += 1
        else:
            require(not row["strongest_qualifying_lof_scope"], f"{gene}: noneligible gene has an eligible scope")
    require(versions == {"v1": 109, "v2": 105} and len(scope_by_gene) == 104,
            "Mechanism eligibility history differs")
    require(support == {"G2P_inferred_absent_product": 100, "G2P_evidence": 4},
            f"Eligible-gene support differs: {support}")
    require(Counter(scope_by_gene.values()) == dict(zip(SCOPES, (7, 16, 63, 18))),
            "Disease-scope allocations differ")
    summary = index(audit.read("mechanism/ALLELIC_SCOPE_BREAKDOWN.tsv"), ("disease_scope",))
    scopes_in_records = set(r["disease_scope"] for r in disease_records)
    require(set(summary) == {(scope,) for scope in scopes_in_records}, "Disease record scopes differ")
    for scope in sorted(scopes_in_records):
        records = [r for r in disease_records if r["disease_scope"] == scope]
        expected = dict(
            record_or_missing_placeholder_n=len(records),
            unique_gene_n_at_this_scope=len({r["approved_symbol"] for r in records}),
            credible_gene_n_at_this_scope=len({r["approved_symbol"] for r in records
                                               if truth(r["credible_ASD_NDD_record"])}),
            final_LoF_n_allocated_to_this_scope=sum(v == scope for v in scope_by_gene.values()),
        )
        compare(summary[(scope,)], {k: str(v) for k, v in expected.items()}, expected.keys(),
                f"Mechanism association scope {scope}")
    audit.record("mechanism", "disease record scope breakdown", f"{len(summary)} reviewed or explicit missing scope groups")
    audit.record("Supplementary Figure S6", "frozen score-1 mapping and five mechanism classes",
                 "245 unique symbols and HGNC IDs; five class counts agree with the plot source")
    audit.record("mechanism", "disease records and eligible LoF scope",
                 "246 linked records; 48 explicit missing records; v1/v2/v3 counts 109/105/104; "
                 "100 inferred absent-product and 4 G2P-evidence genes; exclusive scope counts 7/16/63/18")
    return indexed, scope_by_gene


def conditional_lof_tail(strata):
    """Convolve the original per-stratum hypergeometric distributions."""
    n = hits = cases = observed = informative = 0
    expected = variance = 0.0
    pmf = [1.0]
    for stratum in strata:
        total, h, a, k = [int(stratum[key]) for key in ("N", "H", "A", "K")]
        lo = max(0, a + h - total)
        hi = min(a, h)
        require(total > 0 and 0 <= a <= total and 0 <= h <= total and lo <= k <= hi,
                f"Invalid stratum margins: {stratum['screen']}/{stratum['stratum_key']}")
        n, hits, cases, observed = n + total, hits + h, cases + a, observed + k
        expected += a * h / total
        if total > 1:
            variance += a * (h / total) * (1 - h / total) * ((total - a) / (total - 1))
        informative += a if 0 < a < total and 0 < h < total else 0
        denom = math.comb(total, a)
        updated = [0.0] * (len(pmf) + hi)
        for value in range(lo, hi + 1):
            probability = math.comb(h, value) * math.comb(total - h, a - value) / denom
            for i, earlier in enumerate(pmf):
                updated[i + value] += earlier * probability
        pmf = updated
    require(n > 0 and expected > 0 and math.isclose(sum(pmf), 1.0, abs_tol=1e-12),
            "Invalid combined LoF null distribution")
    return dict(total_n=n, asd_n=cases, asd_hit=observed, other_hit=hits-observed, other_n=n-cases,
                informative_asd_n=informative, expected_hits=expected, variance_hits=variance,
                p=sum(pmf[observed:]), effect=observed / expected)


def verify_lof_scope(audit, inference, background, scopes):
    confounder_strata = audit.read("inputs/R1_CONDITIONING_STRATA.tsv")
    matching = audit.read("inputs/R1_MATCHING_QC.tsv")
    final = audit.read("mechanism/LOF_GROUP_SCOPE_SENSITIVITIES.tsv")
    require(len(index(final, ("screen", "qualifying_disease_scope"))) == 10, "Expected 10 LoF scope rows")

    def is_lof(r):
        return (r["set"], r["analysis"], r["threshold"], r["units"]) == (
            "LOF_COMPATIBLE", "strict_joint", "1", "CPM")

    support = index((r for r in matching if is_lof(r) and r["record_type"] == "CASE_SUPPORT"),
                    ("screen", "approved_symbol"))
    exact = {}
    for screen, reference in CORE.items():
        strata = [r for r in confounder_strata if is_lof(r) and r["screen"] == screen]
        index(strata, ("stratum_key",))
        compare_to = select(inference, screen=screen, reference=reference, threshold="1", units="CPM",
                            analysis="strict_joint", set="LOF_COMPATIBLE")
        computed = {k: str(v) for k, v in conditional_lof_tail(strata).items()}
        compare(computed, compare_to, computed.keys(), f"{screen} exact LoF test")
        exact[screen] = float(computed["p"])
    holm = {}
    running = 0.0
    for i, screen in enumerate(sorted(exact, key=exact.get)):
        running = max(running, min(1.0, exact[screen] * (len(exact) - i)))
        holm[screen] = running

    for screen, reference in CORE.items():
        baseline = select(inference, screen=screen, reference=reference, threshold="1", units="CPM",
                          analysis="strict_joint", set="LOF_COMPATIBLE")
        require(equal_field(str(holm[screen]), baseline["adjusted_p"]), f"{screen}: LoF Holm P differs")
        for scope in ("ALL_104_LOF", *SCOPES):
            observed = select(final, screen=screen, qualifying_disease_scope=scope)
            genes = sorted(g for g, tag in scopes.items() if scope in ("ALL_104_LOF", tag))
            in_screen = [g for g in genes if truth(background[(screen, g)]["measured"])]
            expressed = [g for g in in_screen if background[(screen, g)]["primary_status_1"] in ASSAYED]
            strict = [g for g in expressed if support[(screen, g)]["stratum_key"]]
            informative = [g for g in strict if truth(support[(screen, g)]["informative_stratum"])]

            def hits(group):
                return [g for g in group if truth(background[(screen, g)]["author_hit"])]

            by_direction = Counter(background[(screen, g)]["FOXP3_REGULATOR_DIRECTION"] for g in hits(strict))
            computed = dict(frozen_case_genes=len(genes), screen_measured_case_genes=len(in_screen),
                            measured_case_hits=len(hits(in_screen)), detectable_case_genes=len(expressed),
                            detectable_case_hits=len(hits(expressed)), strict_complete_case_genes=len(strict),
                            strict_complete_case_hits=len(hits(strict)), informative_case_genes=len(informative),
                            informative_case_hits=len(hits(informative)),
                            strict_KO_reduces_FOXP3_hits=by_direction["SUPPORTS_FOXP3"],
                            strict_KO_increases_FOXP3_hits=by_direction["RESTRAINS_FOXP3"],
                            strict_hit_gene_symbols=joined(hits(strict)))
            compare(observed, {k: str(v) for k, v in computed.items()}, computed.keys(),
                    f"{screen}/{scope} LoF scope margins")
            if scope == "ALL_104_LOF":
                fields = {
                    "strict_other_gene_n": baseline["other_n"], "strict_other_gene_hits": baseline["other_hit"],
                    "R1_expected_hits": baseline["expected_hits"], "R1_effect_O_over_E": baseline["effect"],
                    "R1_exact_p": str(exact[screen]), "R1_two_screen_Holm_p": str(holm[screen]),
                    "R1_joint_gate_pass": baseline["original_strict_coverage_pass"],
                }
                compare(observed, fields, fields.keys(), f"{screen} LoF104 inference")
            else:
                require(observed["inference_role"] == "DESCRIPTIVE_SCOPE_NO_NEW_TEST" and
                        not (observed["R1_exact_p"] or observed["R1_two_screen_Holm_p"]),
                        f"{screen}/{scope}: descriptive subgroup contains a new inference")
            audit.record("mechanism", f"{screen} {scope}",
                         f"{len(hits(strict))}/{len(strict)} strict case hits"
                         + (f"; exact Holm P {holm[screen]:.10g}" if scope == "ALL_104_LOF" else "; descriptive"))


def candidate_directions(audit, master, background):
    figures = audit.read("figure_sources/Supplementary_Figure_S6_candidate_source.tsv")
    figure = index(figures, ("gene",))
    candidates = index(audit.read("mechanism/CANDIDATE_DIRECTION_MECHANISM.tsv"), ("gene",))
    evidence = index(audit.read("mechanism/CANDIDATE_SCREEN_DIRECTIONS.tsv"), ("gene", "screen"))
    union = {g for (screen, g), row in background.items() if screen in CORE
             and truth(row["sfari_score1"]) and truth(row["author_hit"])}
    require(len(union) == 17 and set(figure) == set(candidates) == {(g,) for g in union},
            "Candidate membership differs from the frozen hit union")
    require(set(evidence) == {(g, sc) for g in union for sc in SCREENS},
            "Expected 51 unique gene and screen contexts")
    measured_count = 0
    for gene in sorted(union):
        source, candidate, gene_master = figure[(gene,)], candidates[(gene,)], master[(gene,)]
        compare(source, gene_master, ("approved_symbol", "mechanism_class", "lof_eligible"),
                f"S6 candidate {gene} mechanism")
        require(
            candidate["frozen_G2P_gene_level_class"] == gene_master["mechanism_class"] and
            truth(candidate["frozen_LoF104_eligible"]) == truth(gene_master["lof_eligible"]),
            f"{gene}: candidate mechanism classification differs",
        )
        for screen in SCREENS:
            raw, row = background[(screen, gene)], evidence[(gene, screen)]
            measured = truth(raw["measured"])
            hit = truth(raw["author_hit"]) if measured else False
            measured_count += measured
            prefix = screen.lower()
            require(truth(row["screened"]) == measured and truth(candidate[f"{prefix}_screened"]) == measured,
                    f"{gene}/{screen}: measurement status differs")
            require(row["source_native_contrast"] == "FOXP3-low/high", f"{gene}/{screen}: native contrast differs")
            if measured:
                raw_to_output = {"raw_primary_id": "source_native_gene_id", "pos|lfc": "source_pos_LFC",
                                 "neg|lfc": "source_neg_LFC", "pos|fdr": "source_pos_FDR", "neg|fdr": "source_neg_FDR",
                                 "FOXP3_REGULATOR_DIRECTION": "frozen_regulator_label"}
                reconstructed = {dst: raw[src] for src, dst in raw_to_output.items()}
                compare(row, reconstructed, reconstructed.keys(), f"{gene}/{screen} native values")
                require(hit == truth(row["author_hit"]) == (min(float(row["source_pos_FDR"]),
                                                                float(row["source_neg_FDR"])) <= .05),
                        f"{gene}/{screen}: source FDR hit differs")
                standard = DIRECTIONS[row["frozen_regulator_label"]]
                if hit:
                    sign = 1 if standard == "KO_REDUCES_FOXP3" else -1
                    require(sign * float(row["source_pos_LFC"]) > 0 and sign * float(row["source_neg_LFC"]) > 0,
                            f"{gene}/{screen}: signed LFC and standardized KO direction disagree")
            else:
                require(not (row["author_hit"] or row["source_pos_LFC"] or row["source_neg_LFC"]),
                        f"{gene}/{screen}: unassayed result has screen statistics")
                standard = "NOT_ASSAYED"
            require(row["KO_effect_on_FOXP3_protein"] == standard ==
                    candidate[f"{prefix}_KO_effect_on_FOXP3_protein"],
                    f"{gene}/{screen}: standardized KO direction differs")
            require(row["expression_CPM1_status"] == raw["primary_status_1"],
                    f"{gene}/{screen}: detectable-universe label differs")
            fields = {
                "CHEN": ("Chen_hit", "Chen_author_LFC", "Chen_FOXP3_regulator_direction"),
                "UMHOEFER_TREG": ("Umhoefer_Treg_hit", "Umhoefer_Treg_LFC", "Umhoefer_Treg_direction"),
                "UMHOEFER_TCONV": ("Umhoefer_Tconv_hit", "", ""),
            }[screen]
            shown_hit = "YES" if hit else ("NO" if measured else "NOT_MEASURED")
            require(source[fields[0]] == shown_hit, f"S6 {gene}/{screen}: displayed hit differs")
            if fields[1]:
                require(equal_field(source[fields[1]], row["source_pos_LFC"]) and
                        source[fields[2]] == row["frozen_regulator_label"],
                        f"S6 {gene}/{screen}: displayed direction or native LFC differs")
    require(measured_count == 27, "Measured candidate-screen context count differs")
    audit.record("Supplementary Figure S6", "candidate measurement, source signs and standardized KO labels",
                 "17 candidates, 51 screen contexts, 27 measured. Native signs and KO direction agree.")
    return figure, candidates


def functional_evidence(audit, figures, candidates):
    details = audit.read("function/CANDIDATE_EXPERIMENT_DETAIL.tsv")
    designs = index(audit.read("function/FUNCTION_EXPERIMENT_DESIGNS.tsv"), ("experiment_id",))
    experiments = index(audit.read("function/EXPERIMENT_LEVEL_EVIDENCE.tsv"), ("experiment_id",))
    studies = index(audit.read("function/VERIFIED_PUBLICATION_IDENTIFIERS.tsv"), ("study",))
    cover = index(audit.read("function/CANDIDATE_EVIDENCE_COVERAGE.tsv"), ("gene", "study"))
    require(len(details) == 208 and len(designs) == 14 and set(designs) == set(experiments),
            "Evidence endpoint or experiment-design coverage differs")
    require({(g,) for g in {r["gene"] for r in details}} == set(figures), "Evidence candidates differ")
    require(set(studies) == {(name,) for name in REFERENCES}, "Publication identity table differs")
    for study, (pmid, doi, _) in REFERENCES.items():
        require((studies[(study,)]["pmid"], studies[(study,)]["doi"]) == (pmid, doi),
                f"{study}: DOI/PMID differs")
    endpoint_statuses = {"PUBLISHED_DESCRIPTIVE", "SCREEN_HIT", "SCREEN_NON_HIT",
                         "PUBLISHED_EFFECT", "PUBLISHED_NO_DECREASE", "DESIGN_ONLY"}
    generated_endpoint_ids = {}
    for key, design in designs.items():
        group = design["r1_experiment_group"]
        selection = set(design["selection"].split("|")) if design["selection"] else set()
        rows = [r for r in details if r["experiment_group"] == group and r["status"] != "NOT_TESTED"
                and (not selection or r["gene"] in selection)]
        actual = experiments[key]
        doi, pmid = REFERENCES[design["study"]][1], REFERENCES[design["study"]][0]
        expected = dict(study=design["study"], doi=doi, pmid=pmid, assay=design["assay"],
                        source_experiment_independence=design["independence"],
                        cell_state=design["cell_state"], biological_n=design["biological_n"],
                        technical_replication=design["technical_replication"],
                        level_of_evidence=design["evidence_level"], source_figure_or_table=design["source_figure"],
                        endpoint_rows_n=str(len(rows)), perturbed_candidate_genes=joined(r["gene"] for r in rows))
        compare(actual, expected, expected.keys(), f"Functional experiment {key[0]}")
        for r in rows:
            require(r["R1_evidence_file_row"] not in generated_endpoint_ids,
                    f"{key[0]}: endpoint is assigned to more than one assay")
            generated_endpoint_ids[r["R1_evidence_file_row"]] = key[0]
    for row in details:
        require(generated_endpoint_ids.get(row["R1_evidence_file_row"], "") == row["R2_experiment_id"]
                and bool(row["R2_experiment_id"]) == (row["status"] in endpoint_statuses),
                f"Endpoint-to-assay membership differs at row {row['R1_evidence_file_row']}")
        if row["study"] in REFERENCES:
            require((row["pmid"], row["doi"]) == REFERENCES[row["study"]][:2],
                    f"Endpoint bibliography differs at row {row['R1_evidence_file_row']}")
    require(len(index(details, ("R1_evidence_file_row",))) == 208, "Evidence row IDs are not unique")
    direct = [k[0] for k, r in experiments.items() if r["level_of_evidence"] == "ARRAYED_KO_DIRECT_SUPPRESSION"]
    require(direct == ["KOLB_SUPPRESSION"], "Unexpected direct suppression assay")
    kolb = experiments[("KOLB_SUPPRESSION",)]
    require(kolb["perturbed_candidate_genes"] == "SATB1" and kolb["endpoint_rows_n"] == "6"
            and "4 biological replicates" in kolb["biological_n"] and "9 donor-labelled columns" in kolb["technical_replication"],
            "SATB1 assay, sample-size boundary, or workbook limitation differs")
    require(all(r["status"] in ("NOT_FOUND", "NOT_TESTED") for r in details
                if r["evidence_tier"] in ("VARIANT_SPECIFIC_FUNCTION", "PATIENT_EVIDENCE")),
            "Variant-specific or patient evidence tier differs from the frozen review")
    require(set(cover) == {(g[0], study) for g in figures for study in REFERENCES},
            "Candidate and study coverage is incomplete")
    for (gene,), fig in figures.items():
        gene_records = [r for r in details if r["gene"] == gene]
        for study in REFERENCES:
            r = cover[(gene, study)]
            rows = [d for d in gene_records if d["study"] == study and
                    d["status"] in endpoint_statuses and d["status"] != "DESIGN_ONLY"]
            require(r["experiment_ids"] == joined(d["R2_experiment_id"] for d in rows),
                    f"{gene}/{study}: assay coverage differs")
            has_suppression = any(d["evidence_tier"] == "SUPPRESSION_FUNCTION" and d["status"] == "PUBLISHED_EFFECT"
                                  for d in rows)
            require((r["direct_Treg_suppression"] == "YES") == has_suppression,
                    f"{gene}/{study}: assay is not verified as direct suppression")
            require(r["ASD_variant_human_Treg_assay"] == r["ASD_patient_human_Treg_function"] ==
                    "Not established in the reviewed sources", f"{gene}/{study}: evidence boundary differs")
        for short, tier in (("expr", "FOXP3_EXPRESSION"), ("phen", "TREG_PHENOTYPE"),
                            ("sup", "SUPPRESSION_FUNCTION")):
            rows = [d for d in gene_records if d["evidence_tier"] == tier and
                    d["status"].startswith("PUBLISHED_") and d["analysis_role"] in FOLLOW_UPS]
            labels = ", ".join(sorted({REFERENCES[d["study"]][2] for d in rows})) or "–"
            require(fig[f"display_{short}"] == labels and
                    fig[f"sources_{short}"] == " | ".join(d["source_locator"] for d in rows) and
                    fig[f"status_{short}"] == " | ".join(d["status"] for d in rows),
                    f"S6 {gene}/{tier}: displayed readout availability differs")
        followups = [r for r in gene_records if r["R2_experiment_id"] and r["analysis_role"] in FOLLOW_UPS
                     and r["status"] != "DESIGN_ONLY"]
        same = joined(r["R2_experiment_id"] for r in followups if r["study"] in ("Chen 2025", "Umhoefer 2026"))
        external = joined(r["R2_experiment_id"] for r in followups if r["study"] in ("Schumann 2020", "Kolb 2026"))
        candidate = candidates[(gene,)]
        require(candidate["same_paper_followup_experiment_IDs"] == same and
                candidate["independent_followup_experiment_IDs"] == external,
                f"{gene}: same-study versus external evidence differs")
        audit.record("Supplementary Figure S6", f"{gene} published assay cells",
                     f"FOXP3 {fig['display_expr']}, phenotype {fig['display_phen']}, suppression {fig['display_sup']}")
    audit.record("published evidence", "publication and experiment provenance",
                 "4 verified DOI/PMID pairs, 14 experiment contexts, 208 endpoint/coverage rows, 68 candidate-study rows")
    audit.record("published evidence", "SATB1 assay scope",
                 "Kolb PMID 42286212 is the sole direct suppression context; 5C/5D workbook labels and "
                 "published biological replication remain distinct; source reports no ASD-variant or patient assay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path, help="unpacked archive's data/r7 directory")
    parser.add_argument("--report", type=Path, help="optional output for audit TSV")
    args = parser.parse_args()
    audit = EvidenceAudit(args.data_root)
    inference = audit.read("results/R1_FROZEN_SENSITIVITY_GRID.tsv")
    figure_three(audit, inference)
    background = index(audit.read("inputs/R7_GENE_BACKGROUND.tsv.gz"), ("screen", "approved_symbol"))
    master, scopes = mechanism(audit, background)
    verify_lof_scope(audit, inference, background, scopes)
    figures, candidates = candidate_directions(audit, master, background)
    functional_evidence(audit, figures, candidates)
    audit.record("Figure 3d/S6", "original vector/raster and full upstream rebuild",
                 "The release verifies the figure-source values and evidence from processed inputs. "
                 "Original print figures, earlier submission artwork, and separately hosted annotations "
                 "and experiment source materials are separate upstream inputs.", "UPSTREAM_DEPENDENT")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("group", "check", "status", "evidence"),
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(audit.checks)
    print(f"R7 evidence: {sum(r['status'] == 'PASS' for r in audit.checks)} processed-data checks passed.")
    print("Original figure geometry and a rebuild from author-hosted sources are upstream-dependent.")


if __name__ == "__main__":
    main()
