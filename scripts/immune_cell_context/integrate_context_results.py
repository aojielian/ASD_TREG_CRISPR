"""Build shared metadata and run independent checks on R10 summary outputs."""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact

from baseline_recheck import PANEL


def table(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def write(frame, path):
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def assemble(stage):
    audit, res = stage / "01_audit", stage / "03_results"
    human = table(audit / "R10_MONACO_SAMPLE_CROSSWALK.tsv").rename(columns={
        "accession": "GSM", "original_title": "original_label", "matrix_label": "matrix_sample_label",
        "original_biological_replicate": "replicate", "included_primary": "included",
    })
    mouse = table(audit / "R10_IMMGEN_SAMPLE_CROSSWALK.tsv").rename(columns={
        "gsm": "GSM", "source_title": "original_label", "matrix_library": "matrix_sample_label",
        "donor_id": "donor", "replicate_id": "replicate", "include_primary": "included",
    })
    um = table(audit / "R10_UMHOEFER_SAMPLE_METADATA.tsv")
    columns = ["dataset", "GSM", "original_label", "matrix_sample_label", "donor", "replicate", "species",
               "cell_type", "lineage", "tissue", "state", "included", "reason"]
    samples = pd.concat([d[columns] for d in (human, mouse, um)], ignore_index=True)
    samples["included"] = samples.included.astype(str).str.lower().isin(["true", "yes", "1"])
    assert samples.GSM.is_unique and samples.GSM.notna().all()
    write(samples, audit / "R10_SAMPLE_METADATA.tsv")
    geno = table(res / "R10_HUMAN_EXPRESSION_PANEL.tsv")[["gene", "ensembl_gene_id", "gencode_symbol", "hgnc_id",
                                                          "mapping_route", "missing_reason"]].rename(columns={
        "gene": "human_gene", "missing_reason": "monaco_missing_reason",
    })
    ortho = table(audit / "R10_IMMGEN_GENE_CROSSWALK.tsv")
    gene = geno.merge(ortho.drop(columns="dataset"), on="human_gene", how="outer", validate="one_to_one")
    assert set(PANEL) <= set(gene.human_gene)
    assert (gene.loc[gene.human_gene.isin(PANEL), "hgnc_id"]
            .eq(gene.loc[gene.human_gene.isin(PANEL), "human_hgnc"])).all()
    write(gene, audit / "R10_GENE_ID_CROSSWALK.tsv")

    schema = pd.concat([table(audit / filename) for filename in [
        "R10_MONACO_SCHEMA.tsv", "R10_IMMGEN_SCHEMA.tsv", "R10_UMHOEFER_DATASET_SCHEMA.tsv",
    ]] + [table(res / "R10_PERTURB_SOURCE_SHEET_AUDIT.tsv").assign(dataset="Chen2025")], ignore_index=True)
    write(schema, audit / "R10_DATASET_SCHEMA.tsv")
    sensitivity = pd.concat([table(res / f"R10_{sp}_EXPRESSION_SENSITIVITY.tsv") for sp in ["HUMAN", "MOUSE"]],
                            ignore_index=True)
    write(sensitivity, res / "R10_EXPRESSION_SENSITIVITY.tsv")


def validate(stage):
    res, audit, logs = stage / "03_results", stage / "01_audit", stage / "08_logs"
    checks = []

    def check(name, ok, evidence):
        checks.append({"check": name, "status": "PASS" if ok else "FAIL", "evidence": evidence})

    sample = table(audit / "R10_SAMPLE_METADATA.tsv")
    check("GSM uniqueness and metadata", len(sample) == 411 and sample.GSM.is_unique and
          sample.groupby("dataset").included.sum().to_dict() == {"GSE107011": 114, "GSE109125": 115, "GSE286470": 10},
          "127 Monaco + 274 ImmGen + 10 Umhoefer metadata; one row per GSM")
    ids = table(audit / "R10_GENE_ID_CROSSWALK.tsv")
    check("Candidate mappings", set(PANEL) <= set(ids.human_gene) and ids.human_gene.is_unique and
          ids[ids.human_gene.isin(PANEL)].measured.eq("yes").all(), "Panel GENCODE v26/HGNC–MGI keys")
    source = table(res / "R10_14_GENE_FUNCTION_PROFILE.tsv")
    check("Candidate freeze", list(source.gene) == PANEL and source.chen_hit.eq(1).all(),
          "Exact ordered 14-gene panel, all original Chen hits")
    for prefix, nm in [("HUMAN", "gene"), ("MOUSE", "human_gene")]:
        expression = table(res / f"R10_{prefix}_EXPRESSION_PANEL.tsv").set_index(nm).loc[PANEL]
        if prefix == "HUMAN":
            ok = expression.n_cell_types.eq(29).all() and expression.detected_ge1.value_counts().to_dict() == {29: 13, 28: 1}
            evidence = "At 1 TPM, 13 genes meet 29/29 and one gene meets 28/29"
        else:
            ok = expression.n_groups.eq(56).all() and expression.detected_groups_2.eq(56).all()
            evidence = "All 14 orthologs meet 56/56 at 2 TMM CPM"
        check(prefix + " breadth", ok, evidence)
    u = table(res / "R10_TREG_TCONV_UNIVERSE.tsv.gz")
    counts = table(res / "R10_TREG_TCONV_COUNTS.tsv")
    categories = counts[counts.section.eq("shared_hit_categories")].set_index("category")
    primary = table(res / "R10_TREG_TCONV_PRIMARY.tsv").iloc[0]
    observed = categories.loc[["Treg-only", "Tconv-only"], ["asd_n", "non_asd_n"]].to_numpy().T
    check("Umhoefer partition and denominator", len(u) == 1349 and categories.total_n.sum() == len(u)
          and categories.asd_n.sum() == 33 and observed.sum() == primary.n,
          "Shared categories sum to 1,349 genes/33 ASD and 47 discordant/3 ASD")
    check("Umhoefer Fisher p", np.isclose(fisher_exact(observed).pvalue, primary.P, atol=1e-14) and
          primary.CI_low <= primary.OR <= primary.CI_high, "Independent two-sided fixed-margin recomputation")
    sens = table(res / "R10_TREG_TCONV_SENSITIVITY.tsv")
    check("Umhoefer BH family", np.allclose(false_discovery_control(sens.P), sens.P_BH, atol=1e-14) and
          np.allclose(sens.target_fraction, [.05, .10, .20]), "One family of three preregistered ranks")
    continuous = table(res / "R10_TREG_TCONV_CONTINUOUS.tsv").set_index("analysis").loc["abs_rank_mean"]
    check("Umhoefer randomized p", continuous.B >= 100_000 and
          np.isclose((continuous.mc_extreme_n + 1) / (continuous.B + 1), continuous.P, atol=1e-14),
          "Finite Monte Carlo correction and 100,000 fixed-seed allocations")
    cov = table(res / "R10_PERTURB_COVERAGE.tsv")
    protein = table(res / "R10_PERTURB_PROTEIN.tsv")
    check("Chen coverage and NTC scoring", list(cov.gene) == PANEL and cov.targeted_in_input_library.sum() == 8
          and cov.protein_score_available.sum() == 7 and protein.score.notna().sum() == 7 * 88 and
          protein.loc[protein.score.notna(), "reference_score"].eq(0).all(),
          "7 × 88 author NTC-relative values, one additional input-library-only gene, six absent targets")
    cx = table(res / "R10_PERTURB_EXPLORATORY.tsv").query("analysis == 'protein'")
    check("Chen protein BH family", len(cx) == 88 and np.allclose(false_discovery_control(cx.p), cx.q_BH_88, atol=1e-14),
          "One selected-target exploratory family of 88 epitopes with severe remaining LOEUF imbalance")
    combined = pd.DataFrame(checks)
    write(combined, logs / "R10_INTEGRATED_NUMERIC_QA.tsv")
    assert combined.status.eq("PASS").all(), combined[combined.status.eq("FAIL")]

