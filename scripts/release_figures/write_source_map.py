"""Build the public panel-to-source map for the final figures."""
import argparse
import csv
import sys
from hashlib import file_digest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "release"))
from layout import canonical_path


CONTEXT_SOURCES = "figure_sources/R10/"
KNOCKOUT_RESULTS = "data/r11/results/"
REFERENCE = "data/r11/reference/"
LEGACY = "figure_sources/supp/"
MAIN_LEGEND = "figure_sources/R11/Figure_Legends.md"
SUPP_LEGEND = "figure_sources/R11/Supplementary_Figure_Legends.md"

CHEN = "Chen et al., PMID 40140585"
UMHOEFER = "Umhoefer et al., PMID 41237776, GSE286470"
SHIFRUT = "Shifrut et al., PMID 30449619"
MONACO = "Monaco et al., PMID 30726743, GSE107011"
IMMGEN = "Yoshida et al., PMID 30686579, ImmGen GSE109125"
HILL = "Hill et al., PMID 18024188, GSE7460"
MOUSE = "Konopacki et al., PMID 30643266, GSE121279 | Christensen et al., PMID 39156895, GSE263830 | Ren et al., PMID 31125332, GSE121251"
MOUSE_MH = "MSigDB 2025.1.Mm Hallmark | Liberzon et al., PMID 26771021"
MORENO = "Moreno et al., PMID 41703554"
YAN = "Yan et al., PMID 41383128"
FROZEN_PAGES = "figure_sources/frozen_r10/R10_Supplementary_Figures.pdf"
HUMAN_CANDIDATE_PROFILE = CONTEXT_SOURCES + "R10_14_GENE_FUNCTION_PROFILE.tsv"
F3_NULL = (CONTEXT_SOURCES + "F3_foxp3_tcr.tsv", CONTEXT_SOURCES + "F3_foxp3_tcr_null.tsv.gz")
TREG_TCONV_HITS = (CONTEXT_SOURCES + "F3C_hit_classes.tsv", CONTEXT_SOURCES + "F3C_discordant_fisher.tsv", CONTEXT_SOURCES + "R10_TREG_TCONV_ADJUSTED.tsv")
TREG_TCONV_CONTINUOUS = (CONTEXT_SOURCES + "F3D_rank_difference.tsv.gz", CONTEXT_SOURCES + "F3D_continuous_contrast.tsv")
BREADTH = (CONTEXT_SOURCES + "F4B_breadth.tsv", CONTEXT_SOURCES + "F4A_human_expression.tsv", CONTEXT_SOURCES + "R10_MOUSE_EXPRESSION_PLOT.tsv")
LEGACY_S2 = LEGACY + "S2_source.tsv"
LEGACY_S3 = LEGACY + "S3_source.tsv"
LEGACY_S4 = LEGACY + "S4_source.tsv"
LEGACY_S5 = LEGACY + "S5_source.tsv"
RNA_ENDPOINT = KNOCKOUT_RESULTS + "R11_FOXP3_ENDPOINT.tsv"
SCORES = KNOCKOUT_RESULTS + "R11_SIGNATURE_SAMPLE_SCORES.tsv"
CONTRASTS = KNOCKOUT_RESULTS + "R11_SIGNATURE_CONTRASTS.tsv"


def panels():
    rows = []

    def add(panel, readout, files, selection, *citations):
        if isinstance(files, str):
            files = (files,)
        rows.extend((panel, readout, source, selection, " | ".join(citations)) for source in files)

    add("1A", "Induced-Treg screen", CONTEXT_SOURCES + "F1_chen.tsv.gz",
        "screen_gene rows; native MAGeCK low/high log fold change and -log10(minimum raw tail P); SFARI score 1 and author hit flags", CHEN)
    add("1B", "Induced-Treg screen enrichment", CONTEXT_SOURCES + "F1_chen.tsv.gz",
        "effect row; case/control hit margins, conditional OR, exact CI and one-sided P", CHEN)
    add("1C", "14 discovery effects", CONTEXT_SOURCES + "F1_chen.tsv.gz",
        "screen_gene rows with author hit and SFARI score 1; native effect and standardized direction", CHEN)
    add("2A", "Primary-Treg screen enrichment", CONTEXT_SOURCES + "F2_cross_screen.tsv",
        "proportion and effect rows; hit counts and denominators, conditional OR, exact CI and one-sided P", UMHOEFER)
    add("2B", "Shared score-1 measurements", CONTEXT_SOURCES + "F2_cross_screen.tsv",
        "common_asd_gene rows; show all five genes with an author hit in either Treg screen", CHEN, UMHOEFER)
    add("2C", "Gene-property adjustment", CONTEXT_SOURCES + "F2_models.tsv",
        "matched complete-case M0–M3 Firth model rows per screen; OR, CI_low, CI_high and Holm P for M2", CHEN, UMHOEFER)
    add("2D", "LoF directional observed-to-expected ratio", CONTEXT_SOURCES + "F2_lof.tsv",
        "STRICT_CONDITIONAL_EXACT; retain four direction-by-screen tests; effect, case_hits, expected_case_hits and p_holm", CHEN, UMHOEFER)
    add("3A", "FOXP3 vs proliferation exact preference", F3_NULL,
        "FULL_MEASURED_UNADJUSTED and COMPLETE_CASE_CONDITIONED; observed S, full null probabilities, one-sided raw_P", CHEN, SHIFRUT)
    add("3B", "Exclusive primary Treg vs Tconv hits", TREG_TCONV_HITS,
        "common measured background and discordant author hits; exact Fisher OR/CI/P and property-conditioned and Firth sensitivities", UMHOEFER)
    add("3C", "Continuous primary Treg vs Tconv comparison", TREG_TCONV_CONTINUOUS + (CONTEXT_SOURCES + "R10_TREG_TCONV_ADJUSTED.tsv",),
        "all comparable genes, absolute-effect rank differences; label permutation P, bootstrap CI and property-adjusted sensitivity", UMHOEFER)
    add("3D", "Immune atlas expression breadth", BREADTH + (CONTEXT_SOURCES + "R10_MOUSE_EXPRESSION_PANEL.tsv",),
        "14 candidates, human detection ≥1 TPM and mouse ≥1 TMM CPM; detected/evaluated counts; mouse ≥2 CPM reported as sensitivity", MONACO, IMMGEN)
    add("4A", "Foxp3 RNA abundance", (KNOCKOUT_RESULTS + "R11_SAMPLE_EXPRESSION.tsv", KNOCKOUT_RESULTS + "R11_QC_SAMPLES.tsv"),
        "Foxp3 log2_norm_count_plus1 for every experimental library; per-group medians", MOUSE)
    add("4B", "Foxp3 RNA effects", RNA_ENDPOINT,
        "all five contrasts; unshrunken log2FC, Wald CI and adjusted_p; separate Holm families", MOUSE)
    add("4C", "Hill identity score per library", SCORES,
        "program=Hill_LN, component=score; all biological libraries and per-group medians", MOUSE, HILL)
    add("4D", "Hill identity score effects", CONTRASTS,
        "program=Hill_LN, component=score; all five estimates, t intervals and adjusted_P; separate Holm families", MOUSE, HILL)

    for p, selection in (
        ("a", "33 shared score-1 genes; source-native log fold changes and author hit calls"),
        ("b", "cross-screen categories, tested populations and common directional effects"),
        ("c", "Tconv enrichment; measured hit counts, conditional OR, exact CI and one-sided P"),
    ):
        add("S1" + p, "Common FOXP3 assay measurements", LEGACY + "S1_source.tsv", selection, CHEN, UMHOEFER)
    add("S2a", "Calibration and shared-hit exclusions", LEGACY_S2, "summary rows; calibrated hit calls and exclusion of FOXP1 plus ZMYND8", CHEN, UMHOEFER)
    add("S2b", "Class-conditioned recalibration", LEGACY_S2, "conditioned summary rows; observed_expected_ratio and exact P", CHEN, UMHOEFER)
    add("S2c", "Leave-one-gene-out", LEGACY_S2, "one row per removed gene and screen; one_sided_fisher_p", CHEN, UMHOEFER)
    add("S3a", "Regulator class margins", LEGACY_S3, "strata; class_label, N, H, A, K and screen", CHEN, UMHOEFER)
    add("S3b", "Expression/attribute-conditioned summaries", (LEGACY_S3, LEGACY + "S3_conditioning.tsv", LEGACY + "S3_regression.tsv", LEGACY + "S3_model_sensitivities.tsv"),
        "reference and model sensitivities; O/E, within-stratum exact P and two-screen Holm P", CHEN, UMHOEFER)
    add("S3c", "Complete-case and informative coverage", (LEGACY_S3, LEGACY + "S3_conditioning.tsv"),
        "measured candidate denominators, complete-case and informative counts, reference background and coverage gates", CHEN, UMHOEFER)
    add("S3d", "Auxiliary mouse screen", LEGACY_S3, "Loo screen rows; measured and fitness-exclusion odds ratios, CI and P", "Loo et al., mouse Foxp3 screen, figure legend")
    add("S4a", "Paired whole-blood program scores", LEGACY_S4,
        "20 discordant sibling pairs; mean ASD-minus-sibling differences, CIs, exact P and BH q", "Tomaiuolo et al., PMID 36646776, GSE212645", CHEN)
    add("S4b", "Descriptive T-cell program overlap", LEGACY_S4,
        "Yan T_CD4_c03_ANXA2_FOXP3 cluster, annotated Chen program proportions and reported overlaps", YAN, CHEN)
    add("S5a", "Purified-Treg program coverage and effects", LEGACY_S5,
        "10 prespecified programs, gene list coverage, eligible distributional tests, BH family", MORENO, CHEN)
    add("S5b", "Exact-gene overlap in purified Tregs", LEGACY_S5,
        "Moreno exact-gene row; discovery genes measured, overlap within FDR-defined set", MORENO, CHEN)
    add("S5c", "Exact-gene overlap in reported cluster list", LEGACY_S5,
        "Yan row; overlap in the published list, measured background unavailable", YAN, CHEN)
    add("S6a", "Gene2Phenotype mechanisms", LEGACY + "S6_membership.tsv",
        "245 SFARI score-1 genes; mechanism classification, eligible LoF and competing/disease scope fields", "Gene2Phenotype and SFARI snapshots")
    add("S6b", "Candidate-level perturbation and assays", LEGACY + "S6_candidates.tsv",
        "17 candidate rows; standardized direction and published expression, phenotype and suppression evidence", CHEN, UMHOEFER, "Primary citations in source table and legend")
    for p, selection in (
        ("a", "all ten programs; measured coverage and test eligibility"),
        ("b", "five eligible mean signed-effect contrasts; 95% CIs and supported program flag"),
        ("c", "paired empirical cumulative curves for the four BH-supported programs and other features"),
    ):
        add("S7" + p, "Moreno Treg transcriptomic comparison", CONTEXT_SOURCES + "S7_moreno.tsv.gz", selection, MORENO, CHEN)
    add("S8A", "Human immune expression abundance", CONTEXT_SOURCES + "R10_MONACO_EXPRESSION_HEATMAP_SOURCE.tsv",
        "18 genes and 29 sorted cell types; log2_tpm_plus_1", MONACO)
    add("S8B", "Human within-gene relative expression", CONTEXT_SOURCES + "R10_MONACO_EXPRESSION_HEATMAP_SOURCE.tsv",
        "same 18 genes and 29 cell types; row_z", MONACO)
    add("S9", "Mouse baseline abundance", CONTEXT_SOURCES + "R10_MOUSE_EXPRESSION_PLOT.tsv",
        "18 one-to-one orthologs and 56 baseline groups; log2_cpm1", IMMGEN)
    add("S10", "Mouse within-gene relative expression", CONTEXT_SOURCES + "R10_MOUSE_EXPRESSION_PLOT.tsv",
        "same 18 orthologs and 56 groups; zscore", IMMGEN)
    add("S11A", "Human paired expression values", CONTEXT_SOURCES + "R10_EXPRESSION_DONOR_CONTRASTS.tsv",
        "primary_balanced_six; Treg and equal-type-weighted Tconv TPM per donor; candidate plus reference panel", MONACO)
    add("S11B", "Human within-donor expression contrasts", CONTEXT_SOURCES + "R10_EXPRESSION_DONOR_CONTRASTS.tsv",
        "primary_balanced_six; log2_ratio per donor", MONACO)
    add("S12", "Intracellular protein scores", CONTEXT_SOURCES + "R10_PERTURB_PROTEIN.tsv",
        "compartment=intra; 7 measured candidates × 27 marker scores", CHEN)
    add("S13", "Surface protein scores", CONTEXT_SOURCES + "R10_PERTURB_PROTEIN.tsv",
        "compartment=surface; 7 measured candidates × 61 marker scores", CHEN)
    add("S14A", "Human Treg expression ranks", CONTEXT_SOURCES + "R10_HUMAN_EXPRESSION_PANEL.tsv",
        "candidate and reference genes; treg_rank_desc in 29 cell types", MONACO)
    add("S14B", "Mouse Treg expression ranks", CONTEXT_SOURCES + "R10_MOUSE_EXPRESSION_PANEL.tsv",
        "candidate and reference orthologs; splenic_treg_rank and colonic_treg_rank among 56 groups", IMMGEN)
    add("S15A", "FOXP3 vs proliferation preference", F3_NULL,
        "same shared-background and conditioned null probability display as 3A", CHEN, SHIFRUT)
    add("S15B", "Discordant Firth phenotype contrast", CONTEXT_SOURCES + "F3_foxp3_tcr.tsv",
        "DISCORDANT_UNADJUSTED_FIRTH and DISCORDANT_ADJUSTED_FIRTH; profile OR, CI and P", CHEN, SHIFRUT)
    add("S15C", "Primary Treg vs Tconv exclusive hits", TREG_TCONV_HITS,
        "same common-gene, discordant Fisher and property sensitivity evidence as 3B", UMHOEFER)
    add("S15D", "Primary Treg minus Tconv rank contrast", TREG_TCONV_CONTINUOUS,
        "all 1,349 measured genes; ECDFs, mean difference, bootstrap CI and label-permutation P", UMHOEFER)
    add("S16A", "Human atlas expression heatmap", CONTEXT_SOURCES + "R10_MONACO_EXPRESSION_HEATMAP_SOURCE.tsv",
        "candidate plus reference genes, all 29 types, donor-median log2_tpm_plus_1", MONACO)
    add("S16B", "Atlas expression breadth", BREADTH,
        "detection/evaluated ratios for the 14 candidates; human ≥1 TPM and mouse ≥1 TMM CPM", MONACO, IMMGEN)
    add("S16C", "Four selected protein markers", CONTEXT_SOURCES + "F4C_protein_scores.tsv",
        "7 measured candidates by the four selected reference protein scores; relative to non-targeting control", CHEN)
    add("S17A", "Within-cohort principal components", KNOCKOUT_RESULTS + "R11_QC_SAMPLE_PCA.tsv",
        "32 libraries; PC1, PC2 and variance fractions from centered log2(TMM CPM + 0.5)", MOUSE)
    add("S17B", "Raw library sizes", KNOCKOUT_RESULTS + "R11_QC_SAMPLES.tsv",
        "32 samples; raw_library_size by cohort and genotype", MOUSE)
    add("S18A", "Sample correlation", (KNOCKOUT_RESULTS + "R11_QC_SAMPLE_DISTANCES.tsv.gz", REFERENCE + "R11_SAMPLE_MAP.tsv"),
        "within-cohort pairwise Pearson correlations and GEO sample order", MOUSE)
    add("S18B", "Sample distance", (KNOCKOUT_RESULTS + "R11_QC_SAMPLE_DISTANCES.tsv.gz", REFERENCE + "R11_SAMPLE_MAP.tsv"),
        "within-cohort pairwise Euclidean distances and GEO sample order", MOUSE)
    add("S19", "Count distribution by library", KNOCKOUT_RESULTS + "R11_QC_SAMPLES.tsv",
        "32 libraries; raw q25, median, q75, q95 and q99 displayed on log2(count + 1) scale", MOUSE)
    add("S20A", "Prespecified marker effects", KNOCKOUT_RESULTS + "R11_MARKER_EFFECTS.tsv",
        "13 marker genes × five contrasts; unshrunken log2FC, BH status and unestimable state", MOUSE)
    add("S20B", "Target gene RNA abundance effects", KNOCKOUT_RESULTS + "R11_TARGET_QC.tsv",
        "five target transcript effects; unshrunken log2FC and 95% Wald intervals", MOUSE)
    add("S21", "Treg identity reference effects", CONTRASTS,
        "component=score; both Hill_LN and ImmGen_spleen, all five contrasts, estimate/t CI/adjusted_P", MOUSE, HILL, IMMGEN)
    add("S22", "IL-2/STAT5 and broad cell program effects", CONTRASTS,
        "component=score; six fixed Hallmark programs, all five contrasts, estimate/t CI/BH q", MOUSE, MOUSE_MH)
    add("S23A", "Human candidate perturbation evidence", HUMAN_CANDIDATE_PROFILE,
        "FOXP1 and SIN3A; original iTreg, Treg, Tconv and proliferation calls; distinguish unmeasured and nonhit", CHEN, UMHOEFER, SHIFRUT)
    add("S23B", "Mouse RNA and signed Hill estimates", (RNA_ENDPOINT, CONTRASTS),
        "five contrasts; unshrunken Foxp3 RNA log2FC, Hill score estimate and Holm support for each endpoint", MOUSE, HILL)
    add("S23C", "Published direct suppression assays", REFERENCE + "R11_FIGURE_FUNCTION_ASSAYS.tsv",
        "source-author outcome, experimental context, PMID and exact figure locator",
        "Konopacki et al., PMID 30643266 | Ren et al., PMID 31125332 | Christensen et al., PMID 39156895 | Dong et al., PMID 40794436")
    return rows


def sha(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def write_map(data, output):
    rows = []
    for panel, readout, path, fields, citation in panels():
        if panel.startswith("S"):
            number = int("".join(c for c in panel if c.isdigit()))
            legend = SUPP_LEGEND
            art = (FROZEN_PAGES if number <= 14 else
                   f"figure_sources/frozen_r10/Figure_{number - 12}.pdf" if number <= 16 else "")
            mode = ("reuse vector page" if number <= 14 else
                    "reuse vector artwork and append caption" if number <= 16 else "plot from source statistics")
        else:
            number = int(panel[0])
            legend = MAIN_LEGEND
            art = f"figure_sources/frozen_r10/Figure_{number}.pdf" if number <= 2 else ""
            mode = "reuse vector artwork" if number <= 2 else "plot from source statistics"
        source = canonical_path(data, path)
        artwork = canonical_path(data, art) if art else None
        legend_path = canonical_path(data, legend)
        rows.append((panel, readout, source.relative_to(data).as_posix(), sha(source),
                     fields, citation, artwork.relative_to(data).as_posix() if artwork else "",
                     sha(artwork) if artwork else "", mode, "scripts/release_figures/render.sh",
                     legend_path.relative_to(data).as_posix()))
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("panel", "readout", "public_source", "source_sha256", "selection_and_fields", "source_citation",
                         "immutable_vector_source", "immutable_vector_sha256", "render_mode", "replay_entrypoint", "legend_path"))
        writer.writerows(rows)
    print(f"Mapped {len(set(row[0] for row in rows))} panels to {len(rows)} public source records.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    write_map(args.data.resolve(), args.output)


if __name__ == "__main__":
    main()
