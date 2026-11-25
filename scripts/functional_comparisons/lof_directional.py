"""Directional LoF-compatible FOXP3 screen enrichment on frozen v1.1.0 inputs.

Usage: python lof_directional.py DATA_ROOT CODE_REPOSITORY OUTPUT_DIRECTORY
The optional --checks argument writes a separate numerical verification TSV.
"""
from argparse import ArgumentParser
import gzip
import io
from math import fsum
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

VERSION = "v1.1.0"
DOI = "10.5281/zenodo.24487407"
SCOPES = {
    "ASD_named": 7,
    "ASD_feature_in_disorder": 16,
    "NDD_named": 63,
    "NDD_feature_in_disorder": 18,
}
DIRECTIONS = ("KO_REDUCES_FOXP3", "KO_INCREASES_FOXP3")
AMBIGUOUS = "BIDIRECTIONAL_AMBIGUOUS"
REGULATOR_TO_DIRECTION = {
    "SUPPORTS_FOXP3": DIRECTIONS[0],
    "RESTRAINS_FOXP3": DIRECTIONS[1],
    "NOT_SIGNIFICANT": "MEASURED_NOT_HIT",
    AMBIGUOUS: AMBIGUOUS,
}
ALL = "ALL_104_LOF"
COUNTS = ["N", "H", "A", "K"]
DETECTED = ("DETECTABLE_HIT", "DETECTABLE_NONHIT")
CLASSES = ("BOTH_LISTED", "TF_ONLY", "EPIFACTOR_ONLY", "NEITHER_LISTED")


def require(ok, condition):
    if not ok:
        raise ValueError(condition)


def tsv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    opts = dict(sep="\t", index=False, encoding="utf-8", na_rep="NA",
                lineterminator="\n", float_format="%.17g")
    if path.suffix != ".gz":
        frame.to_csv(path, **opts)
        return
    with path.open("wb") as out, gzip.GzipFile(
        filename="", fileobj=out, mode="wb", mtime=0
    ) as archive, io.TextIOWrapper(archive, encoding="utf-8", newline="") as stream:
        frame.to_csv(stream, **opts)


def same(actual, expected, key):
    try:
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False,
                                      rtol=1e-12, atol=1e-14)
    except AssertionError as error:
        raise ValueError(f"{key} disagrees with v1.1.0 inputs: {error}") from error


def import_frozen(repo):
    location = repo / "scripts/gene_property_adjustment"
    for name in ("firth_models", "reconstruct_background", "background_and_matching",
                 "release_inputs"):
        require((location / f"{name}.py").is_file(), f"Missing release module: scripts/gene_property_adjustment/{name}.py")
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(location))
    import firth_models as firth
    import reconstruct_background as background
    import background_and_matching as matching
    return firth, background, matching


def validate_directions(genes, data, firth, matching):
    fdr = genes[["pos|fdr", "neg|fdr"]]
    measured = genes["is_measured"]
    require(np.isfinite(fdr.loc[measured]).all(axis=None), "A measured gene lacks a finite author FDR")
    require(((fdr.loc[measured] >= 0) & (fdr.loc[measured] <= 1)).all(axis=None),
            "An author FDR is outside [0, 1]")
    pos, neg = fdr["pos|fdr"].le(.05), fdr["neg|fdr"].le(.05)
    expected = pd.Series(np.select(
        [~measured, pos & neg, pos, neg],
        ["NOT_MEASURED", AMBIGUOUS, "SUPPORTS_FOXP3", "RESTRAINS_FOXP3"],
        default="NOT_SIGNIFICANT",
    ), index=genes.index)
    # Both frozen human contrasts use FOXP3-low/high. The directional FDR choice
    # identifies the enriched guide bin and reproduces the operational labels.
    native = genes["FOXP3_REGULATOR_DIRECTION"].fillna("NOT_MEASURED")
    bad = expected.ne(native)
    require(not bad.any(), f"Frozen direction disagrees with native FDR for "
                          f"{genes.loc[bad, ['screen', 'approved_symbol']].to_dict('records')}")
    require(genes.loc[measured, "hit"].eq((pos | neg).loc[measured]).all(),
            "The author hit flag differs from the union of directional FDRs")
    require(genes.loc[~measured, ["author_hit", "FOXP3_REGULATOR_DIRECTION",
                                  "pos|fdr", "neg|fdr"]].isna().all(axis=None),
            "An unmeasured gene has a frozen hit or direction observation")

    # Reuse the released R2 comparison, then independently cross-check the other
    # released candidate history including its FDRs and detectability field.
    standard = matching.direction_table(data, genes)
    curated = firth.rd(data / "data/r7/mechanism/CANDIDATE_SCREEN_DIRECTIONS.tsv")
    right = genes.loc[genes["approved_symbol"].isin(standard["approved_symbol"]), [
        "approved_symbol", "screen", "raw_primary_id", "is_measured", "author_hit",
        "pos|lfc", "neg|lfc", "pos|fdr", "neg|fdr", "FOXP3_REGULATOR_DIRECTION",
        "primary_status_1",
    ]].copy()
    right["KO_effect_on_FOXP3_protein"] = right["FOXP3_REGULATOR_DIRECTION"].map(
        REGULATOR_TO_DIRECTION | {"NOT_SIGNIFICANT": "NO_SIGNIFICANT_KO_DIRECTION"}
    ).fillna("NOT_ASSAYED")
    alias = {
        "approved_symbol": "gene", "raw_primary_id": "source_native_gene_id",
        "is_measured": "screened", "pos|lfc": "source_pos_LFC",
        "neg|lfc": "source_neg_LFC", "pos|fdr": "source_pos_FDR",
        "neg|fdr": "source_neg_FDR", "FOXP3_REGULATOR_DIRECTION": "frozen_regulator_label",
        "primary_status_1": "expression_CPM1_status",
    }
    right = right.rename(columns=alias)
    key = ["gene", "screen"]
    require(len(curated) == len(right) == len(standard) == 51
            and not curated.duplicated(key).any() and not right.duplicated(key).any(),
            "Direction audit requires 51 unique candidate × screen records in each table")
    require(curated["source_native_contrast"].eq("FOXP3-low/high").all(),
            "Curated candidate file has an unexpected native contrast")
    same(right.sort_values(key, ignore_index=True),
         curated[right.columns].sort_values(key, ignore_index=True),
         "Candidate mechanism directions, measurement, FDR and detectability")
    return native.map(REGULATOR_TO_DIRECTION).fillna("NOT_MEASURED")


def validate_identity(genes, master, firth):
    require(len(genes) == 22091, "Expected 22,091 frozen screen × gene rows")
    eligible = master.loc[firth.truth(master["lof_eligible"])]
    require(len(master) == 245 and master["approved_symbol"].is_unique,
            "Expected 245 uniquely identified frozen score-1 genes")
    require(master["sfari_score"].eq(1).all(), "Mechanism master has a non-score-1 record")
    require(master["hgnc_id"].is_unique
            and master["hgnc_id"].str.fullmatch(r"HGNC:[1-9][0-9]*", na=False).all(),
            "Frozen score-1 mechanism rows require unique HGNC accessions")
    require(eligible.groupby("strongest_qualifying_lof_scope").size().to_dict() == SCOPES,
            "Frozen eligible disease scopes must sum to 104 with seven ASD-named genes")
    for screen, g in genes.groupby("screen", sort=True):
        # Unmeasured score-1 placeholders have a frozen approved symbol but no
        # screen identifier. Require their symbol and separately audit IDs only
        # on observed records or wherever a frozen screen ID exists.
        observed_ids = g.loc[g["is_measured"], "hgnc_id"]
        require(g["approved_symbol"].is_unique and observed_ids.is_unique
                and observed_ids.str.fullmatch(r"HGNC:[1-9][0-9]*", na=False).all(),
                f"Released gene names and HGNC accessions are not uniquely reconciled for {screen}")
        audit = g.loc[g["asd"].eq(1), ["approved_symbol", "hgnc_id", "is_measured"]].merge(
            master[["approved_symbol", "hgnc_id"]], on="approved_symbol", how="outer",
            validate="one_to_one", suffixes=("_background", "_mechanism"), indicator=True,
        )
        compare = audit["is_measured"].fillna(False) | audit["hgnc_id_background"].notna()
        require(audit["_merge"].eq("both").all() and
                audit.loc[compare, "hgnc_id_background"].eq(audit.loc[compare, "hgnc_id_mechanism"]).all(),
                f"Score-1 HGNC identities differ between frozen sources for {screen}")
    return eligible.set_index("approved_symbol")["strongest_qualifying_lof_scope"]


def build_population(genes, screen, scope):
    members = genes.loc[genes["screen"].eq(screen)].copy()
    members["frozen_disease_scope"] = members["approved_symbol"].map(scope)
    members["lof_case"] = members["frozen_disease_scope"].notna()
    members["analysis_group"] = np.select(
        [members["lof_case"], members["asd"].eq(1)],
        ["LOF_CASE", "OTHER_SCORE1_EXCLUDED"], default="NON_SCORE1_CONTROL",
    )
    keep = members["asd"].eq(0) | members["lof_case"]
    members["in_measured_comparison"] = members["is_measured"] & keep
    members["in_detectable_comparison"] = (
        members["in_measured_comparison"] & members["primary_status_1"].isin(DETECTED)
    )
    joint_fields = (
        members["class_label"].isin(CLASSES)
        & np.isfinite(members["primary_mean_log_expression"])
        & np.isfinite(members["log10_span"])
        & np.isfinite(members["guide_count"])
        & members["guide_count"].gt(0)
    )
    if screen == "CHEN":
        joint_fields &= members["library_member"].isin([False, True])
    members["in_strict_comparison"] = members["in_detectable_comparison"] & joint_fields
    members["comparison_status"] = np.select(
        [~members["is_measured"], ~keep, ~members["in_detectable_comparison"],
         ~members["in_strict_comparison"]],
        ["UNMEASURED", "OTHER_SCORE1_OMITTED", "NOT_CPM1_DETECTABLE", "MISSING_STRICT_PROPERTY"],
        default="IN_STRICT_COMPARISON",
    )
    strict = members.loc[members["in_strict_comparison"]].copy()
    require(len(strict) > 0, f"No strict complete cases for {screen}")
    strict["A_flag"] = strict["lof_case"].astype(int)
    for value, label in (("primary_mean_log_expression", "expression_tercile"),
                         ("log10_span", "span_tercile")):
        # R1 ranks the filtered full case/control population. Ties share a rank.
        strict[label] = np.minimum(3, np.ceil(
            3 * strict[value].rank(method="average") / len(strict)
        )).astype(int)
    strict["guide_class"] = pd.cut(
        strict["guide_count"], [0, 3, 4, np.inf],
        labels=["1_TO_3", "EXACTLY_4", "5_OR_MORE"],
    ).astype(str)
    within = strict["class_label"].astype(str)
    if screen == "CHEN":
        within += "|" + strict["targeted_library"].map({1: "INSIDE", 0: "OUTSIDE"})
    strict["stratum_key"] = (within + "|" + strict["expression_tercile"].astype(str)
                             + "|" + strict["span_tercile"].astype(str)
                             + "|" + strict["guide_class"])
    assignment = ["stratum_key", "expression_tercile", "span_tercile", "guide_class"]
    members = members.join(strict[assignment], validate="one_to_one")
    for direction in DIRECTIONS:
        members[direction + "_Y"] = pd.Series(pd.NA, index=members.index, dtype="Int8")
        members.loc[strict.index, direction + "_Y"] = strict["ko_direction"].eq(direction).astype("Int8")
    return members, strict


def margins(d, y="hit"):
    case = d["A_flag"]
    return d.assign(case_hit=case * d[y]).groupby("stratum_key", sort=True).agg(
        N=(y, "size"), H=(y, "sum"), A=("A_flag", "sum"), K=("case_hit", "sum"),
    )


def local_null(N, H, A):
    lo, hi = max(0, H + A - N), min(A, H)
    k = np.arange(lo, hi + 1)
    p = hypergeom.pmf(k, N, H, A)
    mean = A * H / N
    var = A * (H / N) * (1 - H / N) * (N - A) / (N - 1) if N > 1 else 0
    verify_null(k, p, mean, var, f"N={N}, H={H}, A={A}")
    return k, p, mean, var


def verify_null(k, p, mean, var, name):
    mass = fsum(p)
    m = float(np.dot(k, p))
    v = float(np.dot((k - m) ** 2, p))
    require(np.isfinite(p).all() and np.all(p >= 0) and
            np.allclose([mass, m, v], [1, mean, var], rtol=2e-11, atol=1e-12),
            f"Exact PMF normalization or analytic moments disagree: {name}")


def add_null(records, metadata, level, key, k, p, mean, var, observed):
    total = fsum(p)
    tail = np.minimum(1.0, np.cumsum(p[::-1])[::-1])
    cdf = np.minimum(1.0, np.cumsum(p))
    tail[0] = cdf[-1] = 1.0
    for hit_count, prob, low_p, high_p in zip(k, p, cdf, tail):
        records.append(metadata | dict(
            distribution_level=level, stratum_key=key, case_hit_count=int(hit_count),
            probability=float(prob), cdf=float(low_p), upper_tail=float(high_p),
            pmf_mass=total, null_mean=mean, null_variance=var, observed_k=observed,
        ))


def convolve(table, metadata=None, records=None):
    offset = 0
    pmf = np.array([1.0])
    expected = variance = 0.0
    support = []
    for row in table.itertuples():
        k, p, m, v = local_null(row.N, row.H, row.A)
        require(k[0] <= row.K <= k[-1], f"Observed K is outside the null support: {row.Index}")
        offset += int(k[0])
        pmf = np.convolve(pmf, p)
        expected += m
        variance += v
        support.append((int(k[0]), int(k[-1]), fsum(p), m, v))
        if records is not None:
            add_null(records, metadata, "STRATUM", row.Index, k, p, m, v, int(row.K))
    k = np.arange(offset, offset + len(pmf))
    verify_null(k, pmf, expected, variance, "across strict strata")
    observed = int(table["K"].sum())
    require(k[0] <= observed <= k[-1], "Observed total is outside the combined null support")
    lo, hi = k[np.searchsorted(np.cumsum(pmf), [.025, .975])]
    if records is not None:
        add_null(records, metadata, "COMBINED", "ALL_STRATA", k, pmf, expected, variance, observed)
    p_raw = np.nan
    if variance > 0:
        p_raw = 1.0 if observed == k[0] else min(1.0, fsum(pmf[k >= observed]))
    return dict(expected_case_hits=expected, null_variance=variance, null_count_lo=int(lo),
                null_count_hi=int(hi), null_pmf_mass=fsum(pmf), p_raw=p_raw), support


def baseline_verification(data, genes, firth, checks):
    # Reconcile the two original 2 x 2 hit tables before choosing a LoF subset.
    margin = {"CHEN": (14, 228, 243, 18481), "UMHOEFER_TREG": (5, 28, 26, 1290)}
    grid = firth.rd(data / "data/r7/results/R1_FROZEN_SENSITIVITY_GRID.tsv")
    for screen in firth.CORE:
        d = genes.loc[genes["screen"].eq(screen) & genes["is_measured"]]
        expected = grid.loc[grid["screen"].eq(screen) & grid["analysis"].eq("unadjusted") &
                            grid["set"].eq("SFARI_SCORE1") & grid["reference"].eq("FROZEN_MEASURED")]
        cols = ["asd_hit", "asd_nonhit", "other_hit", "other_nonhit"]
        require(len(expected) == 1 and tuple(expected.iloc[0][cols].astype(int)) ==
                firth.margins(d) == margin[screen], f"Frozen baseline hit margins disagree: {screen}")
        a, b, c, e = margin[screen]
        checks.append((f"{screen}_frozen_full_measured", f"{a}/{a+b} vs {c}/{c+e}", "PASS"))
    return grid


def verify_strict(screen, strict, frozen, prior, scope_table, firth):
    table = margins(strict)
    target = frozen.loc[frozen["screen"].eq(screen)].set_index("stratum_key").sort_index()
    same(table, target[COUNTS], f"{screen} complete strict stratum N/H/A/K")
    target = target.loc[table.index]
    classes = strict.drop_duplicates("stratum_key").set_index("stratum_key").loc[table.index]
    cols = ["class_label", "expression_tercile", "span_tercile", "guide_class"]
    same(classes[cols], target[cols], f"{screen} strict stratum labels")
    if screen == "CHEN":
        require(classes["targeted_library"].eq(firth.truth(target["library_member"])).all(),
                "Frozen Chen strict stratum library assignments disagree")
    shared = table["A"].gt(0) & table["A"].lt(table["N"])
    informative = shared & table["H"].gt(0) & table["H"].lt(table["N"])
    require(shared.eq(firth.truth(target["randomizable_labels"])).all() and
            (table["H"].gt(0) & table["H"].lt(table["N"])).eq(
                firth.truth(target["nondegenerate_hit_distribution"])
            ).all(), f"Frozen strict support flags disagree: {screen}")

    expected = {"CHEN": (74, 10357, 10, 240, 8, 2), "UMHOEFER_TREG": (7, 687, 2, 24, 2, 0)}
    n, hits, cases, case_hits = table[COUNTS].sum()
    split = tuple(int((strict["lof_case"] & strict["ko_direction"].eq(d)).sum()) for d in DIRECTIONS)
    require((cases, n-cases, case_hits, hits-case_hits, *split) == expected[screen],
            f"Frozen ALL_104_LOF strict totals or directional case split disagree for {screen}")

    historical = prior.loc[prior["screen"].eq(screen)]
    require(len(historical) == 1, f"Expected one frozen exact LoF primary row for {screen}")
    old = historical.iloc[0]
    actual = [n, cases, n-cases, case_hits, hits-case_hits,
              table.loc[informative, "A"].sum(), table.loc[shared, "A"].sum(), informative.sum()]
    eq = ["total_n", "asd_n", "other_n", "asd_hit", "other_hit", "informative_asd_n",
          "comparator_supported_asd_n", "informative_strata"]
    require(actual == list(old[eq].astype(int)), f"Strict total, hit, and support margins differ from R1 for {screen}")
    stats, _ = convolve(table)
    require(np.allclose(
        [stats["expected_case_hits"], stats["null_variance"], stats["p_raw"],
         stats["null_count_lo"], stats["null_count_hi"], case_hits/stats["expected_case_hits"]],
        old[["expected_hits", "variance_hits", "p", "null_lo", "null_hi", "effect"]].astype(float),
        rtol=1e-9, atol=1e-13,
    ), f"Frozen all-hit exact tail, moments, or quantiles disagree for {screen}")
    scope_row = scope_table.loc[scope_table["screen"].eq(screen) &
                                scope_table["qualifying_disease_scope"].eq(ALL)].iloc[0]
    require(np.allclose(
        old[["expected_hits", "effect", "p", "adjusted_p"]].astype(float),
        scope_row[["R1_expected_hits", "R1_effect_O_over_E",
                   "R1_exact_p", "R1_two_screen_Holm_p"]].astype(float),
        rtol=1e-10,
    ), f"Two frozen all-hit reference tables disagree for {screen}")
    return table, shared, informative, stats["p_raw"]


def summarize_population(d, case):
    counts = dict(case_n=int(case.sum()), case_all_hits=int(d.loc[case, "hit"].sum()))
    for prefix, group in (("case", d.loc[case]), ("control", d.loc[~case])):
        counts |= {
            f"{prefix}_reduce_hits": int(group["ko_direction"].eq(DIRECTIONS[0]).sum()),
            f"{prefix}_increase_hits": int(group["ko_direction"].eq(DIRECTIONS[1]).sum()),
            f"{prefix}_ambiguous_hits": int(group["ko_direction"].eq(AMBIGUOUS).sum()),
        }
    return counts | dict(control_n=int((~case).sum()), control_all_hits=int(d.loc[~case, "hit"].sum()))


def coverage_rows(members, strict, shared, informative, scope_table, metadata, firth):
    frozen_cases = members.loc[members["lof_case"]]
    measured = members.loc[members["in_measured_comparison"]]
    detectable = members.loc[members["in_detectable_comparison"]]
    stage_sets = {
        "FROZEN_CASE_SET": frozen_cases,
        "MEASURED": measured,
        "CPM1_DETECTABLE": detectable,
        "STRICT_COMPLETE_CASE": strict,
        "STRICT_COMPARATOR_SUPPORTED": strict.loc[strict["stratum_key"].isin(shared[shared].index)],
        "STRICT_ALL_HIT_INFORMATIVE": strict.loc[strict["stratum_key"].isin(informative[informative].index)],
    }
    rows = []
    for scope in (ALL, *SCOPES):
        all_scopes = scope == ALL
        x = scope_table.loc[scope_table["qualifying_disease_scope"].eq(scope)
                            & scope_table["screen"].eq(metadata["screen"])]
        require(len(x) == 1, f"Missing or duplicate scope row: {metadata['screen']} {scope}")
        frozen = x.iloc[0]
        current = {}
        for name, pop in stage_sets.items():
            choose = pop["lof_case"] if all_scopes else pop["frozen_disease_scope"].eq(scope)
            if not all_scopes:
                pop = pop.loc[choose]
                choose = pop["lof_case"]
            stats = summarize_population(pop, choose)
            if name == "FROZEN_CASE_SET":
                stats = dict(case_n=stats["case_n"])
            if not all_scopes or name == "FROZEN_CASE_SET":
                stats = {k: v for k, v in stats.items() if not k.startswith("control_")}
            current[name] = stats
            record = metadata | dict(
                disease_scope=scope, stage=name, inference_role="DESCRIPTIVE_COVERAGE",
                case_fraction_of_frozen=stats["case_n"] / int(frozen["frozen_case_genes"]),
                case_fraction_of_measured=np.nan if name == "FROZEN_CASE_SET" else
                (stats["case_n"] / int(frozen["screen_measured_case_genes"])
                 if int(frozen["screen_measured_case_genes"]) else np.nan),
            ) | stats
            if all_scopes and name == "STRICT_COMPLETE_CASE":
                complete_fraction = record["case_fraction_of_measured"]
                informative_fraction = int(frozen["informative_case_genes"]) / stats["case_n"]
                gate = complete_fraction >= .8 and informative_fraction >= .7
                require(gate == bool(firth.truth(pd.Series([frozen["R1_joint_gate_pass"]])).iloc[0]),
                        f"Frozen prior broad gate setting differs for {metadata['screen']}")
                record.update(prior_all_hit_gate_pass=gate, prior_all_hit_informative_fraction=informative_fraction,
                              prior_all_hit_expected=frozen["R1_expected_hits"],
                              prior_all_hit_effect=frozen["R1_effect_O_over_E"],
                              prior_all_hit_p=frozen["R1_exact_p"], prior_all_hit_holm=frozen["R1_two_screen_Holm_p"])
            rows.append(record)

        count_keys = {
            "FROZEN_CASE_SET": ["frozen_case_genes"],
            "MEASURED": ["screen_measured_case_genes", "measured_case_hits"],
            "CPM1_DETECTABLE": ["detectable_case_genes", "detectable_case_hits"],
            "STRICT_COMPLETE_CASE": ["strict_complete_case_genes", "strict_complete_case_hits"],
            "STRICT_ALL_HIT_INFORMATIVE": ["informative_case_genes", "informative_case_hits"],
        }
        for name, names in count_keys.items():
            values = current[name]
            got = [values["case_n"]] if len(names) == 1 else [values["case_n"], values["case_all_hits"]]
            require(got == list(frozen[names].astype(int)),
                    f"Frozen {scope} case counts disagree for {metadata['screen']} {name}")
        stats = current["STRICT_COMPLETE_CASE"]
        require([stats["case_reduce_hits"], stats["case_increase_hits"]] ==
                list(frozen[["strict_KO_reduces_FOXP3_hits", "strict_KO_increases_FOXP3_hits"]].astype(int)),
                f"Frozen {scope} strict directional case counts disagree for {metadata['screen']}")
        if all_scopes:
            require([stats["control_n"], stats["control_all_hits"]] ==
                    list(frozen[["strict_other_gene_n", "strict_other_gene_hits"]].astype(int)),
                    f"Frozen strict control counts disagree for {metadata['screen']}")
        mask = strict["lof_case"] if all_scopes else strict["frozen_disease_scope"].eq(scope)
        hit_symbols = sorted(strict.loc[mask & strict["hit"].eq(1), "approved_symbol"])
        frozen_symbols = sorted(str(frozen["strict_hit_gene_symbols"]).split(" | ")) if (
            pd.notna(frozen["strict_hit_gene_symbols"])
        ) else []
        require(hit_symbols == frozen_symbols,
                f"Frozen {scope} strict case-hit identities differ for {metadata['screen']}")
    return rows


def directional_tests(strict, all_hits, metadata, measured_cases, matching):
    results, strata, nulls = [], [], []
    ambig = margins(strict.assign(y=strict["ko_direction"].eq(AMBIGUOUS).astype(int)), "y")
    endpoints = {
        label: margins(strict.assign(y=strict["ko_direction"].eq(label).astype(int)), "y")
        for label in DIRECTIONS
    }
    require((endpoints[DIRECTIONS[0]][["H", "K"]] + endpoints[DIRECTIONS[1]][["H", "K"]]
             + ambig[["H", "K"]]).eq(all_hits[["H", "K"]]).all(axis=None),
            f"Directions do not sum to the released all-hit margins for {metadata['screen']}")
    class_tags = strict.drop_duplicates("stratum_key").set_index("stratum_key")
    for label, table in endpoints.items():
        other_label = next(d for d in DIRECTIONS if d != label)
        spec = metadata | dict(direction=label)
        exact, support = convolve(table, spec, nulls)
        n, h, a, k = (int(value) for value in table[COUNTS].sum())
        shared = table["A"].gt(0) & table["A"].lt(table["N"])
        direction_varies = table["H"].gt(0) & table["H"].lt(table["N"])
        informative = shared & direction_varies
        comparator_cases = int(table.loc[shared, "A"].sum())
        informative_cases = int(table.loc[informative, "A"].sum())
        require(bool(informative.any()) == (exact["null_variance"] > 0),
                f"Informative-stratum flags disagree with null variance for {metadata['screen']} {label}")
        broad_gate = a / measured_cases >= .8 and informative_cases / a >= .7
        coverage = dict(
            frozen_case_n=104, measured_case_n=measured_cases,
            strict_fraction_of_measured=a/measured_cases, randomizable_strata=int(shared.sum()),
            randomizable_case_n=comparator_cases, randomizable_case_fraction=comparator_cases/a,
            informative_strata=int(informative.sum()), informative_case_n=informative_cases,
            informative_case_fraction=informative_cases/a, directional_coverage_gate_pass=broad_gate,
            coverage_status="WITHIN_R1_COVERAGE_BENCHMARKS" if broad_gate else
            "LIMITED_COVERAGE_OR_INFORMATION",
        )
        common = spec | coverage | dict(
            total_n=n, case_n=a, control_n=n-a, case_hits=k, control_hits=h-k,
            case_nonhits=a-k, control_nonhits=n-a-h+k,
            case_other_direction_hits=int(endpoints[other_label]["K"].sum()),
            control_other_direction_hits=int((endpoints[other_label]["H"] -
                                              endpoints[other_label]["K"]).sum()),
            case_ambiguous_hits=int(ambig["K"].sum()),
            control_ambiguous_hits=int((ambig["H"] - ambig["K"]).sum()),
            population="COMPARABLE_MEASURED_CPM1_STRICT_COMPLETE_CASES",
            denominator_rule="LoF-compatible score-1 cases and non-score-1 controls",
            outcome_rule="Specified direction=1; all other measured genes=0",
            holm_family_size=4,
        )
        oe = k / exact["expected_case_hits"] if exact["expected_case_hits"] > 0 else np.nan
        results.append(common | exact | dict(
            analysis="STRICT_CONDITIONAL_EXACT", inference_role="PRIMARY_DIRECTIONAL_CONDITIONAL",
            effect_type="OBSERVED_OVER_EXPECTED", effect=oe, ci_low=np.nan, ci_high=np.nan,
            status="ESTIMABLE" if informative.any() else "NON_INFORMATIVE_ZERO_NULL_VARIANCE",
            holm_family="FOUR_DIRECTIONAL_CONDITIONAL_TESTS",
        ))
        f = matching.fisher_stats(strict.assign(
            asd=strict["A_flag"], hit=strict["ko_direction"].eq(label).astype(int)
        ))
        require([f["ASD_hits"], f["ASD_nonhits"], f["other_hits"], f["other_nonhits"]] ==
                [k, a-k, h-k, n-a-h+k], "Fisher table must use the same strict population")
        expected_fisher = hypergeom.sf(k-1, n, h, a)
        if np.isfinite(f["P"]):
            require(np.isclose(f["P"], expected_fisher, rtol=1e-12, atol=1e-15),
                    "One-sided Fisher P differs from an independent hypergeometric tail")
        results.append(common | dict(
            analysis="SAME_POPULATION_FISHER", inference_role="UNCONDITIONED_SENSITIVITY",
            effect_type="CONDITIONAL_MAXIMUM_LIKELIHOOD_ODDS_RATIO", effect=f["OR"],
            ci_low=f["CI_low"], ci_high=f["CI_high"], p_raw=f["P"], status=f["status"],
            holm_family="FOUR_SAME_POPULATION_FISHER_TESTS",
        ))
        alternate = endpoints[other_label]
        for row, (lo, hi, mass, mean, variance) in zip(table.itertuples(), support):
            key = row.Index
            tag = class_tags.loc[key]
            strata.append(spec | dict(
                stratum_key=key, class_label=tag["class_label"],
                library_member=bool(tag["targeted_library"]) if spec["screen"] == "CHEN" else np.nan,
                expression_tercile=int(tag["expression_tercile"]), span_tercile=int(tag["span_tercile"]),
                guide_class=tag["guide_class"], N=row.N, H=row.H, A=row.A, K=row.K,
                all_hit_H=int(all_hits.loc[key, "H"]), all_hit_K=int(all_hits.loc[key, "K"]),
                other_direction_H=int(alternate.loc[key, "H"]), other_direction_K=int(alternate.loc[key, "K"]),
                ambiguous_H=int(ambig.loc[key, "H"]), ambiguous_K=int(ambig.loc[key, "K"]),
                randomizable_labels=bool(shared.loc[key]), direction_outcome_variation=bool(direction_varies.loc[key]),
                informative=bool(informative.loc[key]), pmf_min=lo, pmf_max=hi,
                pmf_mass=mass, expected_case_hits=mean, null_variance=variance,
            ))
    return results, strata, nulls


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("data_root", type=Path)
    p.add_argument("repo_root", type=Path)
    p.add_argument("output_dir", type=Path)
    p.add_argument("--checks", type=Path)
    args = p.parse_args()
    data, repo, out = (path.resolve() for path in (args.data_root, args.repo_root, args.output_dir))
    for path in [out, *([args.checks.resolve()] if args.checks else [])]:
        require(not path.is_relative_to(data) and not path.is_relative_to(repo),
                "Output cannot be written inside frozen data or code")
    from baseline_release import verify_release
    verify_release(data)
    firth, background, matching = import_frozen(repo)
    checks = []
    genes, _ = firth.read_inputs(data)
    genes, _, eligible = background.rebuild_background(data, genes)
    master = firth.rd(data / "data/r7/mechanism/ALLELIC_MECHANISM_MASTER.tsv")
    scope = validate_identity(genes, master, firth)
    require(set(scope.index) == eligible, "The independent LoF eligible-member sets differ")
    checks.append(("frozen_mapping_and_LoF_membership",
                   "22,091 unique screen × approved symbol records; measured HGNC accessions match; 245 score-1; 104 LoF-compatible incl. 7 ASD-named", "PASS"))
    genes["ko_direction"] = validate_directions(genes, data, firth, matching)
    ambiguous_n = int(genes["ko_direction"].eq(AMBIGUOUS).sum())
    checks.append(("author_hit_and_direction",
                   f"Full measured hit/FDR and direction agreement; both 51-row candidate comparisons; ambiguous hits={ambiguous_n}",
                   "PASS"))
    grid = baseline_verification(data, genes, firth, checks)
    original = background.rebuild_strata(data, genes, eligible)
    mask = lambda d: (d["set"].eq("LOF_COMPATIBLE") & d["analysis"].eq("strict_joint")
                      & d["threshold"].eq(1) & d["units"].eq("CPM") & d["screen"].isin(firth.CORE))
    frozen = original.loc[mask(original)]
    prior = grid.loc[mask(grid)]
    require(len(frozen) == 103, "Expected 103 archived LoF strict stratum keys")
    scope_table = firth.rd(data / "data/r7/mechanism/LOF_GROUP_SCOPE_SENSITIVITIES.tsv")
    tables = {"ENRICHMENT": [], "STRATA": [], "MEMBERS": [], "COVERAGE": [], "NULL": []}
    baseline = []
    for screen in firth.CORE:
        meta = dict(source_version=VERSION, source_doi=DOI, screen=screen,
                    reference=firth.REF[screen], threshold=1, units="CPM", case_set="LOF_COMPATIBLE")
        members, strict = build_population(genes, screen, scope)
        all_hits, shared, informative, p_value = verify_strict(screen, strict, frozen, prior, scope_table, firth)
        baseline.append(p_value)
        coverage = coverage_rows(members, strict, shared, informative, scope_table, meta, firth)
        measured = next(r for r in coverage if r["disease_scope"] == ALL and r["stage"] == "MEASURED")
        tests, strat, pmfs = directional_tests(strict, all_hits, meta, measured["case_n"], matching)
        for name, records in (("ENRICHMENT", tests), ("STRATA", strat), ("COVERAGE", coverage), ("NULL", pmfs)):
            tables[name].extend(records)
        member_cols = [
            "approved_symbol", "hgnc_id", "sfari_score1", "lof_case", "frozen_disease_scope",
            "analysis_group", "is_measured", "author_hit", "ko_direction", "FOXP3_REGULATOR_DIRECTION",
            "primary_status_1", "in_measured_comparison", "in_detectable_comparison", "in_strict_comparison",
            "comparison_status", "class_label", "library_member", "primary_mean_log_expression",
            "GENE_SPAN_BP", "log10_span", "guide_count", "expression_tercile", "span_tercile",
            "guide_class", "stratum_key", *[d + "_Y" for d in DIRECTIONS],
        ]
        tables["MEMBERS"].append(
            members[member_cols].sort_values("approved_symbol").assign(**meta)
                   .loc[:, [*meta, *member_cols]]
        )
        n, h, a, k = all_hits[COUNTS].sum()
        checks.append((f"{screen}_all_hit_R1_strata",
                       f"{len(all_hits)} full N/H/A/K keys; {k}/{a} vs {h-k}/{n-a}; previous exact tail {p_value:.16g}",
                       "PASS"))
        checks.append((f"{screen}_scope_and_direction_partition",
                       "All 5 frozen scope rows; published strict membership and support; both directions plus ambiguity sum to all hit per stratum",
                       "PASS"))

    previous_adjusted = firth.p_adjust_two(baseline)
    old_holm = prior.set_index("screen").loc[list(firth.CORE), "adjusted_p"]
    require(np.allclose(previous_adjusted, old_holm, rtol=1e-10, atol=1e-14),
            "Archived all-hit Holm P values did not reproduce")
    require(np.allclose(firth.p_adjust_two([.01, .02, .5, np.nan]), [.04, .06, 1, np.nan], equal_nan=True),
            "Holm implementation must preserve a missing fourth family member")
    all_out = pd.DataFrame(tables["ENRICHMENT"])
    require(len(all_out) == 8, "Expected four exact and four Fisher tests")
    for name, family in all_out.groupby("holm_family", sort=False):
        require(len(family) == 4, f"Each directional Holm family must have four rows: {name}")
        all_out.loc[family.index, "p_holm"] = firth.p_adjust_two(family["p_raw"])
    sample = pd.DataFrame({"asd": [1, 1, 0, 0], "hit": [0, 0, 1, 1]})
    lower, upper = (matching.fisher_stats(x) for x in (sample, sample.assign(hit=1-sample["hit"])))
    require(lower["OR"] == lower["CI_low"] == 0 and lower["CI_high"] > 0 and lower["P"] == 1
            and np.isinf(upper["OR"]) and np.isinf(upper["CI_high"]) and upper["CI_low"] > 0,
            "The reused Fisher implementation must preserve zero and infinite MLE OR boundaries")
    checks.extend([
        ("directional_null_normalization_and_moments",
         f"206 local directional PMFs and 4 combined PMFs; mass and first two moments agree with analytic hypergeometric moments",
         "PASS"),
        ("four_member_separate_Holm_families",
         "Separate exact and Fisher groups of four, including verification of the unestimable-row convention", "PASS"),
        ("fisher_boundaries", "One-sided tails cross-checked against hypergeometric survival; tested zero and infinite conditional odds ratios", "PASS"),
    ])
    tables["ENRICHMENT"] = all_out
    tables["MEMBERS"] = pd.concat(tables["MEMBERS"], ignore_index=True)
    for label, rows in tables.items():
        suffix = ".tsv.gz" if label == "MEMBERS" else ".tsv"
        tsv(pd.DataFrame(rows), out / f"LOF_DIRECTIONAL_{label}{suffix}")
    if args.checks:
        tsv(pd.DataFrame(checks, columns=["check", "observed", "status"]), args.checks)


if __name__ == "__main__":
    main()
