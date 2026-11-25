# Methods and statistics

## Screen ascertainment and primary statistics

The frozen SFARI 2026 Q2 score-1 set contains 245 genes. Human-screen hits use the authors' union of positive and negative FDR ≤0.05. Chen measured 18,966 induced-Treg genes with 14 hits among 242 measured score-1 genes. Umhoefer measured 1,349 primary-Treg genes with 5 hits among 33 measured score-1 genes. The original enrichments compare each study's complete measured universe with one-sided Fisher exact tests. Tconv, Shifrut proliferation and Loo mouse Treg each have their own universe and readout. Source-native LFCs and interpreted knockout directions remain separate.

## Detection and gene-property sensitivity

The primary detectable background is CPM ≥1. Chen iTreg NTC requires two of three libraries at threshold in each guide group. The primary Umhoefer Treg reference uses Chen peripheral hTreg with at least three of four libraries at threshold. HPA nTPM sensitivity uses its own reference. Counts are normalized before canonical-symbol mapping. The model covariate is the mean per-library log2(CPM + 1), using log2(nTPM + 1) in the HPA scenario. Measurement, reference missingness and threshold eligibility are separate fields.

The primary detectable case counts are 174 in Chen and 18 in primary Treg. Their expressed-background unadjusted ORs are approximately 3.95 and 10.86. Class conditioning uses fixed TF and EpiFactors classes. Stricter scenarios and matching retain distinct covariate representation, comparator support, and coverage and balance diagnostics. The two primary detectable-background class-conditioned tests each yield a Holm P of 0.0160. The frozen NDD list and overlap-stripped sensitivity comparisons are included.

Within each screen, reference and threshold, the same genes enter M0–M3. The primary complete-case cohorts are 10,463 total and 173 score-1 genes in Chen and 698 total and 18 score-1 genes in primary Treg. M0 uses the ASD indicator. M1 adds catalogue flags and Chen targeted-library membership. M2 adds expression and gnomAD v2.1.1 LOEUF. M3 adds log10 genomic span and guide count. Firth regression uses the Jeffreys-penalized binomial logit likelihood. Main ORs exponentiate the ASD coefficients, with penalized profile CIs and P values. The two primary M2 tests form one Holm family and give 0.3599574937658 and 0.1006914496552. Exploratory M3 and historical Wald results remain distinguishable. [Frozen model methods](../config/gene_property_methods.md)

## Direction and phenotype comparisons

- The directional LoF set contains 104 scoped LoF-compatible genes (100 inferred and four evidence-tagged). The two screens and two knock-out directions form four exact tests with fixed matched covariate strata and Holm adjustment. Opposite-direction hits and ambiguous author directions remain distinguishable. Full coverage, separate patient disease scopes, OR and O/E, and all sensitivity records are supplied.
- The FOXP3–TCR contrast uses 18,615 common measured Chen–Shifrut genes, including 242 score-1 genes. The exclusive-hit statistic is +1 for FOXP3-only, −1 for TCR-stimulated proliferation-only, and 0 for both or neither. The unconditioned one-sided exact P is 0.000318886582681224. The property-conditioned complete-case null uses class, Chen library membership, frozen healthy-CD4 expression and LOEUF strata. Its exact P is 0.01189315432641426. The result is flagged for limited case coverage and residual LOEUF imbalance. The adjusted discordant-case Firth sensitivity estimates exclusive-hit type in the reduced cohort, giving P = 0.063994632231608506. [Exact and Firth specifications](../scripts/functional_comparisons/METHODS.md)
- The primary Treg–Tconv comparison uses a two-sided Fisher exact test on 47 discordant genes within the common screened universe. Its P is 0.56706753006475485. Property-conditioned, continuous-effect and additional threshold sensitivities use their own defined estimands. [Umhoefer analysis methods](../scripts/immune_cell_context/config/umhoefer_methods.md)

## Atlas expression and perturbation

Monaco median human TPM and ImmGen TMM-normalized mouse CPM are separate abundance measures. Figure 3D uses ≥1 TPM and ≥1 TMM CPM. The mouse ≥2 CPM sensitivity likewise shows broad detectability. Human donor contrasts and expression-matched comparisons, orthologue mappings, and atlas cell-type and lineage breadth are distributed. Chen Perturb-icCITE-seq provides defined coverage for protein and program endpoints with explicit unmeasured states. [Monaco methods](../scripts/immune_cell_context/config/monaco_methods.md) · [ImmGen methods](../scripts/immune_cell_context/config/immgen_methods.md)

## Mouse differential expression and transcriptional programs

The primary DESeq2 models use genotype by state for GSE121279, a sex-adjusted Sin3a fit for GSE263830 and the separately matched auxiliary Foxp1 design for GSE121251. The Sin3a design has six pooled biological libraries per genotype, three female and three male. There are 32 included RNA-seq libraries. Three core and two auxiliary contrasts form separate Holm families for Foxp3 and the Hill identity score. Prespecified signatures use within-library ranks independently of the knockout estimates. CAMERA competitive sensitivity and gene-wise DESeq2 adjusted P values have their own families. The supplied evidence tables retain auxiliary and discordant effects.

## Patient-program and published evidence context

Eligibility in Moreno requires at least 20 measured symbols and at least 70% membership coverage across 5,851 finite effect feature rows. The two-sided rank comparison uses a ten-program BH family, with P = 1 for ineligible programs in the correction calculation. Whole-blood program differences use exact two-sided sign flips across 20 paired ASD-minus-sibling contrasts and BH across ten programs. Moreno has zero of 14 discovery genes among its DE hits. Yan remains descriptive due to its source universe.

Published experiment-level evidence retains expression, protein, transcriptional identity and directly assayed suppression or in vivo phenotypes alongside experiment design. The Zhu conventional CD4 analysis follows source-reported gene and assay coverage, including the distinct three matched flow donor records and six-donor bulk RNA series. The ZMYND8 arrayed endpoint is Th1 GATA3 median fluorescence intensity. [Zhu](https://doi.org/10.1016/j.cell.2026.08.002) · [Source and citation guide](SOURCES.md)

Source paths, reproduction entries and output identifiers for the primary analyses are indexed in the [analysis map](analysis_map.tsv).
