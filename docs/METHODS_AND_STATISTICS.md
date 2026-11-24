# Methods and statistics

## Original processed analysis

The frozen SFARI 2026 Q2 score-1 set contains 245 genes. Human screen hits use the author union of positive and negative FDR ≤0.05. The original comparisons use each study's full frozen measured gene universe, and one-sided Fisher exact tests compare score-1 genes to other measured genes. The Chen human induced-Treg screen measured 18,966 genes and had 14 hits among 242 measured score-1 genes. The Umhoefer targeted primary human Treg screen measured 1,349 genes and had 5 hits among 33 measured score-1 genes. The two sorts measure intracellular FOXP3 expression. The measured Tconv, Shifrut proliferation and mouse Loo contexts have their own screened universes and biological endpoints. The original conditional odds ratios use the conditional exact implementation, which differs from the sample cross-product odds ratio in some sensitivity tables.

For direct phenotype preference, Chen and Shifrut share 18,615 measured genes, including 242 score-1 genes. Each gene enters one of four hit combinations. The signed statistic among score-1 genes counts FOXP3-only hits minus proliferation-only hits. The exact upper tail is calculated under sampling 242 labels without replacement, with shared and neither hits pooled as neutral for the signed statistic. The observed difference is +3. The exact one-sided P is 0.000318886582681224.

Frozen Chen M1–M10 perturbation-program memberships are tested against 5,851 finite Moreno patient Treg RNA effect **feature rows**. Formal program eligibility requires at least 20 measured symbols and at least 70% coverage. The two-sided Mann–Whitney analysis adjusts across all ten programs, with P = 1 for the ineligible programs in the correction family. The associated ECDFs use the actual signed feature effects. This is a feature-level program analysis. Whole-blood program scores instead use the exact two-sided sign flips of 20 paired ASD-minus-sibling contrasts and BH across ten programs. The paired-t interval is used only for display. The exact discovery-gene list returned zero of 14 differential-expression hits in the Moreno analysis; the cohort, endpoint and tested universe remain recorded with each overlap result.

## Post-review sensitivity analysis dated 23 November 2026

The original screen comparisons and all original thresholds remain in the legacy tables. The added sensitivity analyses use expression thresholds fixed for the 23 November 2026 post-review analysis. The primary expression restriction is CPM ≥1. Chen iTreg uses the study's baseline non-targeting-guide expression libraries, with at least two of three libraries meeting the CPM cutoff in each guide group. The Umhoefer primary-Treg reference is Chen's peripheral hTreg expression, with at least three of four libraries meeting the cutoff. Other cutoffs and the HPA nTPM external reference are distinct scenarios. Library normalization uses counts before canonical-symbol mapping. The model expression covariate within each reference is the mean per-library log2(CPM + 1), or log2(nTPM + 1) for the HPA scenario. The scaled terms and expression summaries are released.

Under primary detectability filtering, the case-gene denominators shift to 174 for Chen and 18 for primary Treg. The expressed-background unadjusted OR is approximately 3.95 and 10.86, respectively. Exact conditional analysis stratifies by fixed TF and EpiFactors classes (and by the specified variables in strict variants). It evaluates observed hits against the sum of stratum-wise hypergeometric counts with coverage and identifiability gates. In the primary detectable-background class-conditioned comparison, the two-screen Holm P is 0.0160 in each screen. Strict exact conditioning, matched comparison and Firth regression use different coverage and comparisons. NDD disease-control sets also intersect FOXP3 regulation; the source comparisons and overlap-stripped tests remain available.

## Same-complete-case model grid

Within each screen/reference/threshold scenario, all four models use the identical complete-case cohort for the declared M0–M3 covariate set. The primary CPM ≥1 model cohorts are 10,463 total genes / 173 score-1 genes in Chen and 698 / 18 in Umhoefer primary Treg.

| Model | Terms beyond intercept |
| --- | --- |
| M0 | ASD score-1 indicator |
| M1 | M0 plus TF, EpiFactors and targeted-library membership for Chen |
| M2 | M1 plus continuous expression and gnomAD v2.1.1 LOEUF |
| M3 | M2 plus log10 genomic span and guide count |

Firth mean-bias reduction uses the Jeffreys-penalized binomial logit likelihood. Main odds ratios are exponentiated ASD-indicator coefficients. Intervals and P values use the penalized profile likelihood. The two primary M2 tests defined for the post-review analysis share a Holm family. M0, M1, M3, external references and matching remain exploratory or diagnostic as documented. Original model Wald diagnostics are retained separately in the historical bridge. Estimation or coverage failures remain explicitly flagged and numeric cells can be missing.

| Model and screen | OR | Lower 95% CI | Upper 95% CI | Nominal P | Two-screen M2 Holm P |
| --- | ---: | ---: | ---: | ---: | ---: |
| M2 Chen iTreg | 1.331614172154 | 0.7061893091742 | 2.345216980248 | 0.3599574937658 | 0.3599574937658 |
| M2 Umhoefer Treg | 3.837761029373 | 0.9978377796148 | 13.09664415239 | 0.05034572482759 | 0.1006914496552 |
| Exploratory M3 Umhoefer Treg | 5.757611149375 | 1.40996735197 | 21.61374911821 | 0.01610907384996 | — |

Both M2 tests fall above the family-adjusted 5% threshold. The nominal M3 result uses the same small case cohort and has exploratory status. Coverage, balance and stage-specific margins accompany the numerical results.

## Allelic mechanism and published functional evidence

Source-native positive and negative LFC values remain separate from standardized biological knockout direction. An interpreted KO effect refers to the protein phenotype in that screen. Across the two human FOXP3 screens, the direction for FOXP1 is concordant and the direction for ZMYND8 is context-dependent. Allelic-mechanism records preserve G2P confidence and relevant ASD or broader neurodevelopmental disease scope. Of 245 genes, 130 have a broad LoF-compatible classification. The final scoped LoF group contains 104 genes, comprising 100 database-inferred and four evidence-tagged classifications. Gene-wide perturbations and individual patient alleles have distinct meanings.

The functional evidence tables separate Chen same-study follow-up, targeted Umhoefer screening, Schumann Treg measurements and Kolb SATB1 suppression. Zhu et al. adds selected-candidate context from conventional CD4 T-cell CRISPRi/Perturb-seq and arrayed testing. Among 17 candidates, all have guide-quality records and 16 have DE summary rows across 46 observed gene-condition entries. TLE3 has no reported DE summary. In Zhu Fig. 4G / Supplementary Table 18, Th1 ZMYND8 CRISPR interference has a larger GATA3 median fluorescence intensity, log2FC +0.32603825853458884, source Welch P 0.0033260002706872517 and BH-adjusted P 0.023282001894810764. Three matched flow donor records describe coverage, and a separate bulk RNA series reports six donors. Other available reported adjusted arrayed endpoints remain in the source-derived table. Zhu evidence concerns this conventional CD4 endpoint and selected genes. The FOXP3 model inputs and M2 family are unchanged. [Zhu et al., Cell](https://doi.org/10.1016/j.cell.2026.08.002)

See the [frozen configuration](../config/r7_method.md), [data dictionary](DATA_DICTIONARY.md) and [source guide](SOURCES.md).
