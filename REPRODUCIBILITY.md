# Reproducibility

Processed-data minimal reproduction is supported from the companion Zenodo analysis-ready snapshot, DOI [10.5281/zenodo.22704659](https://doi.org/10.5281/zenodo.22704659). The row-level tables contain the complete inputs needed by the two standard-library scripts, including the Chen and Umhoefer universes, common comparator universe, regulator-class table, strict joint complete-case tables, M1–M10 memberships, Loo rows, Moreno rows, and whole-blood pair/program differences.

The workflow recomputes contingency counts, Fisher/rank-sum P values, exact conditioning observed K, expected K, O/E, and direct two-test Holm-adjusted P from the recomputed exact P values. It does not import Holm values from a summary table. It also checks class-conditioning and joint-conditioning coverage, including 189/242 and 22/33, and explicitly verifies FOXP1/ZMYND8 as BOTH_HIT and CTCF/PTEN/SATB1 as TREG_ONLY among the 33 common ASD score-1 genes.

For the Loo boundary, row-level counts and Fisher P are recomputed. The conditional MLE odds ratio is retained as `LOO_REFERENCE_R_CONDITIONAL_OR` from the released result table; it is a reference estimate rather than a newly recomputed OR.

The package does not claim a full raw-source rebuild because the exact historical upstream ingest and canonicalization scripts are not included. The minimal scripts use only Python 3 standard-library modules.
