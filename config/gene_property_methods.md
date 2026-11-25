# Gene-property analysis configuration

The frozen score-1 list is SFARI Gene Q2 2026. It contains 245 members. A hit uses the source author union of positive or negative FDR ≤0.05. The measured universe, five direction/status labels and technical definitions are stock release inputs. The source screen readout is intracellular FOXP3 expression by FACS.

Expression-restricted comparisons include CPM ≥0.5, ≥1 and ≥2 in study/state references, and an external HPA ≥1 nTPM sensitivity. CPM libraries use the full valid source gene rows as denominators. Chen iTreg requires two of three NTC libraries for each guide. The Chen peripheral references each require three of four libraries. HPA Treg uses the T-reg aggregate, and its Tconv proxy uses the maximum of the naïve and memory CD4 aggregates when both are available. Strict exact conditioning and its coverage gates are retained. The old strict strata do not include LOEUF.

M0–M3 fit the same complete cases within each screen and expression scenario. M0 is unadjusted; M1 adds frozen TF and EpiFactors flags and the targeted-library indicator for Chen; M2 adds continuous mean log expression and LOEUF; M3 adds log span and guide count. Firth mean bias reduction is achieved by the Jeffreys penalty for the binomial logit model. The two primary CPM1 M2 penalized profile P values are adjusted with Holm. M0–M3, the Chen library sensitivity, matching, historical bridge, and other reference/threshold combinations retain their individual raw or diagnostic outcomes.

LoF classification distinguishes the phenotype scope and G2P evidence tag per record. The final 104-member sensitivity includes multiple NDD scopes. Mechanism labels and source-native directions are linked but gene-wide knockout and specific patient variants are not equivalent. Refer to `data/gene_property_adjustment/allelic_mechanisms` for source-correction ancestry and the scope breakdown.

These sensitivities are secondary to the original primary analysis. Coverage or estimation failures are carried through without inventing intervals or imputing results.
