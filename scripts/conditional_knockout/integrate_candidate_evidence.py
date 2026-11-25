"""Integrate typed candidate endpoints and calculate descriptive gene-effect concordance."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PANEL = "CREBBP CUL3 DYRK1A EP300 FOXP1 KDM3B KDM6B KMT5B PHIP SIN3A SMARCA4 TBL1XR1 TLE3 ZMYND8".split()


def tsv(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def write(data, path):
    pd.DataFrame(data).to_csv(path, sep="\t", index=False, na_rep="NA")

CONTEXT = {
    "F1N": ("Konopacki et al.", "Naïve Treg", "GSE121279"),
    "F1A": ("Konopacki et al.", "Activated Treg", "GSE121279"),
    "S3": ("Christensen et al.", "Sorted Treg; sex-adjusted", "GSE263830"),
    "F2R": ("Ren et al.", "Resting Treg", "GSE121251"),
    "F2A": ("Ren et al.", "Activated Treg", "GSE121251"),
}
KNOWN = [
    ("FOXP1", "Konopacki et al.", "Naïve Treg", "F1N", "Foxp3 protein", "Intracellular flow cytometry", "down",
     "30643266, Fig. 2e–f", "Distinct published protein readout"),
    ("FOXP1", "Konopacki et al.", "Activated Treg", "F1A", "Foxp3 protein", "Intracellular flow cytometry", "no clear change",
     "30643266, Fig. 2e–f", "Protein abundance had no detected decrease in this state despite the R11 RNA result"),
    ("FOXP1", "Konopacki et al.", "Adoptive transfer", "F1", "Immune suppression", "In vivo transfer assay",
     "impaired", "30643266, Fig. 5e", "Published weaker protection against inflammation; independent endpoint in this study"),
    ("FOXP1", "Konopacki et al.", "Naïve Treg, low vs high IL-2", "F1", "STAT5 phosphorylation", "Flow cytometry",
     "lower at low IL-2", "30643266, Fig. 7d", "Low-dose response fell; high-dose response did not show the same reduction"),
    ("FOXP1", "Konopacki et al.", "Naïve and activated Treg", "F1", "Foxp3 chromatin occupancy", "ChIP-seq", "subset lost",
     "30643266, Fig. 4a–e", "Binding differences do not measure suppression"),
    ("FOXP1", "Ren et al.", "Resting Treg", "F2R", "Foxp3 RNA and protein", "qPCR, immunoblot, flow cytometry", "down",
     "31125332, Fig. 6a–c", "Different assay replicates and endpoints from the R11 RNA-seq analysis"),
    ("FOXP1", "Ren et al.", "Activated Treg", "F2A", "Foxp3 RNA and protein", "qPCR, immunoblot, flow cytometry",
     "no clear change", "31125332, Fig. 6a–c", "Published baseline, separate from challenge conditions"),
    ("FOXP1", "Ren et al.", "Challenge and induced Treg", "F2", "Foxp3 stability", "Transfer and stimulation follow-up",
     "impaired", "31125332, Fig. 6e–g", "Challenge data are a separate cell-state context"),
    ("FOXP1", "Ren et al.", "Induced Treg in colitis transfer", "F2", "Immune suppression", "In vivo cotransfer protection",
     "impaired", "31125332, Fig. 5b–d", "Published cotransfer result; phenotype is context dependent"),
    ("SIN3A", "Christensen et al.", "Conditional Treg deletion", "S3", "Foxp3 protein", "Intracellular flow cytometry",
     "down", "39156895, Fig. 2h–j", "Distinct published readout with developmental conditional deletion"),
    ("SIN3A", "Christensen et al.", "Conditional Treg deletion", "S3", "Foxp3 RNA", "Published RNA-seq and ancillary measurements",
     "down", "39156895, Fig. 4c/S16B", "Author analysis of the study's public RNA data; separate classification from new R11 fit"),
    ("SIN3A", "Christensen et al.", "Conditional Treg deletion", "S3", "Immune suppression",
     "Responder-proliferation suppression assay", "impaired", "39156895, Fig. 5a", "Published direct suppression endpoint"),
    ("SIN3A", "Christensen et al.", "Conditional Treg deletion", "S3", "Lineage stability", "Fate tracing",
     "more ex-Treg", "39156895, Fig. 6a", "Published fate tracing; RNA population analysis is distinct"),
]
DOIS = {"30643266": "10.1038/s41590-018-0291-z", "31125332": "10.1371/journal.pbio.3000270",
        "39156895": "10.3389/fimmu.2024.1444937"}


def direction(v):
    if pd.isna(v):
        return "unavailable"
    return "up" if v > 0 else "down" if v < 0 else "zero"


def join_table(data, stage):
    frozen = data / "data/r10/results"
    res = stage / "03_results"
    human = tsv(frozen / "R10_14_GENE_FUNCTION_PROFILE.tsv").set_index("gene")
    protein = tsv(frozen / "R10_PERTURB_PROTEIN.tsv").query("assay_label == 'intra-B0041-FOXP3'").set_index("gene")
    foxp3 = tsv(res / "R11_FOXP3_ENDPOINT.tsv")
    sig = tsv(res / "R11_SIGNATURE_CONTRASTS.tsv")
    assert list(human.index) == PANEL and len(foxp3) == 5
    rows = []

    def record(gene, origin, species, study, context, contrast, endpoint, assay, direction_label,
               status, source, explanation, n=None, estimate=None, low=None, high=None,
               p=None, adjusted=None, comparison=None, unit=None, pmid=None):
        rows.append(dict(
            gene=gene, evidence_origin=origin, species=species, study=study, context=context, contrast=contrast,
            endpoint=endpoint, assay=assay, n=n, estimate=estimate, ci_low=low, ci_high=high, p=p,
            p_adjusted=adjusted, direction=direction_label, comparison_to_Chen=comparison, status=status,
            source=source, interpretation=explanation, unit=unit, pmid=pmid, doi=DOIS.get(pmid),
        ))

    for gene, a in human.iterrows():
        sign = direction(a.chen_signed_foxp3_bin_lfc)
        record(gene, "Frozen human screen", "Human", "Chen et al.", "Induced Treg", "CHEN",
               "FOXP3 protein sorting phenotype", "Pooled CRISPR KO / intracellular FOXP3 bins", sign,
               a.chen_screen_status, "R10_14_GENE_FUNCTION_PROFILE.tsv; Chen Supplementary Table 1",
               "Signed bin enrichment on the operational KO scale. RNA/protein concentration and suppression are separate.",
               estimate=a.chen_signed_foxp3_bin_lfc, adjusted=min(a.chen_pos_fdr, a.chen_neg_fdr),
               unit="Signed MAGeCK bin log2 fold change", pmid="40140585")
        for prefix, name in [("treg", "Primary Treg"), ("tconv", "Stimulated Tconv")]:
            measured = a[prefix + "_measured"] == 1
            effect = a[prefix + "_signed_foxp3_bin_lfc"]
            record(gene, "Frozen human screen", "Human", "Umhoefer et al.", name, prefix.upper(),
                   "FOXP3 protein sorting phenotype", "Targeted CRISPR KO / intracellular FOXP3 bins",
                   direction(effect) if measured else "unavailable", a[prefix + "_screen_status"],
                   "R10_14_GENE_FUNCTION_PROFILE.tsv; GSE286470 Combined_gene_summary",
                   "The other cell context has a distinct donor set. Unmeasured targeted-library rows remain missing.",
                   estimate=effect if measured else None, adjusted=min(a[prefix + "_pos_fdr"], a[prefix + "_neg_fdr"])
                   if measured else None, unit="Signed MAGeCK bin log2 fold change", pmid="41237776")
        record(gene, "Frozen human screen", "Human", "Shifrut et al.", "Primary CD8 T cell", "SHIFRUT",
               "TCR-stimulated proliferation", "Pooled CRISPR KO / proliferation", a.shifrut_direction,
               a.shifrut_screen_status, "R10_14_GENE_FUNCTION_PROFILE.tsv; Shifrut Table S4",
               "Separate phenotype. Direction corresponds to the published TCR screen, with no cross-endpoint comparison.",
               adjusted=min(a.shifrut_pos_fdr, a.shifrut_neg_fdr) if a.shifrut_measured else None, pmid="30449619")
        b = protein.loc[gene]
        if pd.notna(b.score):
            record(gene, "Frozen within-Chen follow-up", "Human", "Chen et al.", "Induced Treg", "IC_CITE",
                   "FOXP3 protein", "Perturb-icCITE-seq", direction(b.score), "Author score available",
                   f"Chen Table 2, ProteinScore_Targets, cell {b.source_column}{int(b.source_row)}",
                   "Selected-target author score, descriptive; no available donor-level CI.",
                   estimate=b.score, unit="Mean protein difference divided by NTC mean", pmid="40140585")

    for row in foxp3.itertuples(index=False):
        study, context, accession = CONTEXT[row.contrast]
        gene = row.target.upper()
        s = direction(row.log2FC)
        comparison = ("concordant decrease across distinct RNA and protein-bin readouts" if s == "down"
                      else "estimate points opposite; uncertainty spans zero" if s == "up"
                      else "no directional estimate")
        record(gene, "R11 public-RNA reanalysis", "Mouse", study, context, row.contrast,
               "Foxp3 RNA", "DESeq2, KO minus control", s, row.status, "R11_FOXP3_ENDPOINT.tsv; " + accession,
               "Independent from human screening. Gene knockout is context-specific and differs from patient variants.",
               n=row.n_KO+row.n_control, estimate=row.log2FC, low=row.CI_low, high=row.CI_high,
               p=row.pvalue, adjusted=row.adjusted_p, comparison=comparison, unit="Unshrunken log2 fold change")
        for v in sig[(sig.contrast.eq(row.contrast)) & sig.component.eq("score")].itertuples(index=False):
            record(gene, "R11 public-RNA reanalysis", "Mouse", study, context, row.contrast,
                   v.program, "External fixed transcriptional program score", direction(v.estimate), v.status,
                   "R11_SIGNATURE_CONTRASTS.tsv; " + accession,
                   "Sample-level RNA score. Treg identities, IL-2 response and broader programs are distinct RNA features.",
                   n=v.n_KO+v.n_control, estimate=v.estimate, low=v.CI_low, high=v.CI_high,
                   p=v.P, adjusted=v.adjusted_P, unit=v.unit)
    for gene, study, context, contrast, endpoint, assay, sign, source, explanation in KNOWN:
        pmid = source.split(",")[0]
        record(gene, "Published original experiment", "Mouse", study, context, contrast, endpoint, assay, sign,
               "Verified primary source", source + "; https://pubmed.ncbi.nlm.nih.gov/" + pmid + "/", explanation, pmid=pmid)
    write(rows, res / "R11_CROSS_STUDY_EVIDENCE.tsv")


def effects(stage):
    res = stage / "03_results"
    tables = {}
    for key in CONTEXT:
        d = tsv(res / f"R11_DE_{key}.tsv.gz")
        d = d[d.symbol_unique.eq(True) & d.gene_symbol.notna() & np.isfinite(d.log2FC) & np.isfinite(d.pvalue)]
        assert d.gene_symbol.is_unique
        tables[key] = d.set_index("gene_symbol")[["log2FC", "padj"]]
    comparisons = [(x, "S3") for x in ("F1N", "F1A", "F2R", "F2A")] + [("F1N", "F2R"), ("F1A", "F2A")]
    rows = []
    for x, y in comparisons:
        d = tables[x].join(tables[y], lsuffix="_a", rsuffix="_b", how="inner", validate="one_to_one")
        a, b = d.log2FC_a, d.log2FC_b
        nonzero = a.ne(0) & b.ne(0)
        eligible = d.padj_a.notna() & d.padj_b.notna()
        shared = (d.padj_a.lt(.05) & d.padj_b.lt(.05))
        same = np.sign(a) == np.sign(b)
        rows.append(dict(
            contrast_a=x, contrast_b=y, tested_uniquely_named_a=len(tables[x]), tested_uniquely_named_b=len(tables[y]),
            common_tested_unique_mouse_symbols=len(d), finite_genome_BH_on_both=int(eligible.sum()),
            spearman_effect_rho=spearmanr(a, b).statistic, nonzero_effect_symbols=int(nonzero.sum()),
            consistent_direction_n=int((same & nonzero).sum()),
            consistent_direction_fraction=float(same[nonzero].mean()),
            both_BH_below_0p05=int(shared.sum()), both_BH_below_0p05_same_direction=int((shared & same & nonzero).sum()),
            both_BH_below_0p05_opposite_direction=int((shared & ~same & nonzero).sum()),
            method="Exact unique mouse gene symbols, separate DESeq2 Wald models; descriptive gene-level comparison",
        ))
    write(rows, res / "R11_CROSS_PERTURBATION_CONCORDANCE.tsv")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data", type=Path)
    p.add_argument("stage", type=Path)
    args = p.parse_args()
    join_table(args.data, args.stage)
    effects(args.stage)
