#!/usr/bin/env python3
"""Reproduce the signed FOXP3–proliferation comparison and its R9 adjustments.

Usage: python phenotype_preference.py DATA_ROOT CODE_ROOT OUTPUT_DIR
DATA_ROOT is the extracted v1.1.0 deposit and CODE_ROOT the corresponding code.
"""
from argparse import ArgumentParser
from collections import Counter
from importlib import import_module
from itertools import combinations
import gzip
import io
import math
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

RELEASE = "v1.1.0"
FULL = "FULL_MEASURED_UNADJUSTED"
COMPLETE = "COMPLETE_CASE_UNADJUSTED"
CONDITIONED = "COMPLETE_CASE_CONDITIONED"
UNADJUSTED_FIRTH = "DISCORDANT_UNADJUSTED_FIRTH"
ADJUSTED_FIRTH = "DISCORDANT_ADJUSTED_FIRTH"
CASE_COVERAGE_MIN = .8
INFORMATIVE_FRACTION_MIN = .7
RESIDUAL_SMD_LIMIT = .1
CLASS_BY_FLAGS = {
    (0, 0): "NEITHER_LISTED",
    (1, 0): "TF_ONLY",
    (0, 1): "EPIFACTOR_ONLY",
    (1, 1): "BOTH_LISTED",
}
PREDICTORS = ("TF", "EpiFactors", "expression", "LOEUF", "targeted_library")
ATTRIBUTES = {
    "healthy_cd4_mean_logcpm": "mean log2(CPM+1), frozen healthy CD4 GSE90600",
    "LOEUF_v2_1_1": "gnomAD v2.1.1 upper observed/expected pLoF bound",
}
TERTILE_RULE = "min(3,ceil(3*average_rank/n)); rank among all complete cases; ties retained"
SMD_RULE = (
    "Mean difference (score1 minus other) divided by sqrt(mean of group unweighted "
    "sample variances); weighted row uses case weight 1 and control weight A/(N-A) "
    "within strata containing cases and controls"
)
POSITIVE_FLAGS = {"true", "1", "1.0", "yes"}
NEGATIVE_FLAGS = {"false", "0", "0.0", "no"}


def check(test, condition, observation, checks):
    if not bool(condition):
        raise ValueError(f"{test} disagrees with the analysis definition: {observation}")
    checks.append(dict(branch="B", check=test, status="PASS", observation=observation))


def rd(root, name, **kwargs):
    return pd.read_csv(root / name, sep="\t", low_memory=False, **kwargs)


def tsv(path, rows):
    d = pd.DataFrame(rows)
    options = dict(sep="\t", index=False, lineterminator="\n", float_format="%.17g")
    if path.suffix == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                    d.to_csv(stream, **options)
    else:
        d.to_csv(path, encoding="utf-8", **options)


def flags(raw, label, allow_unknown=False):
    encoded = raw.astype("string").str.strip().str.lower()
    yes = encoded.isin(POSITIVE_FLAGS)
    no = encoded.isin(NEGATIVE_FLAGS)
    known = yes | no
    if not allow_unknown and not known.all():
        raise ValueError(f"Missing or invalid Boolean values in {label}: {(~known).sum()}")
    result = yes.astype("Int8")
    return result.where(known)


def keys(d, column, label, checks):
    check(f"{label}_keys",
          d[column].notna().all() and d[column].ne("").all() and d[column].is_unique,
          f"{len(d):,} nonempty unique {column} keys", checks)


def load_universe(data, release_firth, checks):
    background, _ = release_firth.read_inputs(data)
    check("r7_background_and_r6_covariates", len(background) == 22_091,
          "22,091 screen/gene records and released read_inputs() covariate reconciliation", checks)
    for screen, expected in (("CHEN", (14, 242, 243, 18724)),
                             ("UMHOEFER_TREG", (5, 33, 26, 1316))):
        measured = background[background.screen.eq(screen) & background.is_measured]
        case = measured.asd.eq(1)
        fdr_hit = measured["neg|fdr"].le(.05) | measured["pos|fdr"].le(.05)
        margins = (int(measured.loc[case, "hit"].sum()), int(case.sum()),
                   int(measured.loc[~case, "hit"].sum()), int((~case).sum()))
        check(f"{screen}_measured_baseline", margins == expected and
              measured["hit"].eq(fdr_hit.astype(int)).all(),
              f"score-1 hits/N and other hits/N = {margins}; author union FDR <=.05", checks)

    common = rd(data, "data/processed/chen_shifrut_common_universe.tsv",
                dtype={"approved_symbol": "string"})
    keys(common, "approved_symbol", "shared_measured_universe", checks)
    common["score1"] = flags(common.SFARI_SCORE1, "shared SFARI_SCORE1").astype(int)
    common["chen_hit"] = flags(common.CHEN_HIT_BOOL, "shared CHEN_HIT_BOOL").astype(int)
    common["shifrut_hit"] = flags(common.SHIFRUT_HIT, "shared SHIFRUT_HIT").astype(int)
    D = common.chen_hit - common.shifrut_hit
    expected_category = D.map({1: "FOXP3_ONLY", 0: "CONCORDANT", -1: "TCR_ONLY"})
    check("signed_phenotype_and_shifrut_hit_directions",
          D.eq(common.D).all() and common.derived_category.eq(expected_category).all() and
          common.SHIFRUT_DIRECTION.eq("NO_HIT").eq(common.shifrut_hit.eq(0)).all() and
          common.SHIFRUT_DIRECTION.isin(
              ["NO_HIT", "NEG_TCR_PHENOTYPE", "POS_TCR_PHENOTYPE"]).all(),
          "D recomputed as Chen hit minus Shifrut hit; sign and source categories agree", checks)
    common["D"] = D

    cols = ["approved_symbol", "hgnc_id", "asd", "hit", "FOXP3_REGULATOR_DIRECTION",
            "TF", "EpiFactors", "class_label", "CLASS_ANNOTATION_UNCERTAIN", "annotation_status",
            "library_member", "HEALTHY_CD4_MEAN_LOGCPM", "LOEUF_v2_1_1"]
    chen = background.loc[background.screen.eq("CHEN") & background.is_measured, cols]
    keys(chen, "approved_symbol", "chen_measured", checks)
    keys(chen, "hgnc_id", "chen_frozen_hgnc", checks)
    check("chen_hgnc_id_format", chen.hgnc_id.str.fullmatch(r"HGNC:\d+").all(),
          "Frozen HGNC identifiers match the canonical ID format", checks)
    # The historical shared table serializes symbols in upper case. Match the
    # representation only; preserve the original values and frozen HGNC spelling.
    exact_matches = common.approved_symbol.isin(chen.approved_symbol)
    common = common.rename(columns={"approved_symbol": "legacy_common_symbol"})
    common["symbol_join_key"] = common.legacy_common_symbol.str.upper()
    chen = chen.assign(symbol_join_key=chen.approved_symbol.str.upper())
    keys(common, "symbol_join_key", "shared_uppercase_symbol", checks)
    keys(chen, "symbol_join_key", "chen_uppercase_symbol", checks)
    common = common.merge(chen, how="left", on="symbol_join_key",
                          validate="one_to_one", indicator=True, sort=False)
    common["case_only_symbol_match"] = common.legacy_common_symbol.ne(common.approved_symbol)
    case_only = common.loc[common.case_only_symbol_match.fillna(False)]
    check("legacy_uppercase_serialization",
          common["_merge"].eq("both").all() and common.legacy_common_symbol.eq(common.symbol_join_key).all() and
          (~exact_matches).sum() == len(case_only) == 192 and
          case_only.score1.eq(0).all() and case_only.D.value_counts().to_dict() == {0: 186, -1: 6},
          "192 capitalization-only matches (186 neutral, 6 TCR-only controls); all 18,615 match; zero aliases",
          checks)
    check("shared_join_and_hit_mapping",
          common["_merge"].eq("both").all() and common.score1.eq(common.asd).all() and
          common.chen_hit.eq(common.hit).all() and
          common.CHEN_DIRECTION.eq(common.FOXP3_REGULATOR_DIRECTION).all() and
          common.CHEN_DIRECTION.isin(
              ["NOT_SIGNIFICANT", "SUPPORTS_FOXP3", "RESTRAINS_FOXP3"]).all() and
          common.CHEN_DIRECTION.eq("NOT_SIGNIFICANT").eq(common.chen_hit.eq(0)).all(),
          f"One-to-one shared-to-measured Chen join; {len(common):,} symbols, cases and hit directions agree",
          checks)

    frozen_class = rd(data, "data/processed/regulator_class_gene_level.tsv",
                      dtype={"approved_symbol": "string"})
    keys(frozen_class, "approved_symbol", "class_source", checks)
    columns = {
        "class_label": "class_label", "annotation_status": "annotation_status",
        "TF": "tf_listed", "EpiFactors": "epifactors_genes_listed",
        "CLASS_ANNOTATION_UNCERTAIN": "class_annotation_uncertain",
        "library_member": "chen_library_member_umhoefer_treg",
        "HEALTHY_CD4_MEAN_LOGCPM": "healthy_cd4_mean_logcpm",
    }
    paired = common.merge(frozen_class, how="left", on="approved_symbol",
                          suffixes=("", "_class_source"), validate="one_to_one", indicator="class_merge")
    equal = paired.class_merge.eq("both").all()
    for left, right in columns.items():
        cr = right + "_class_source" if left == right else right
        if left == "HEALTHY_CD4_MEAN_LOGCPM":
            equal &= np.allclose(paired[left], paired[cr], rtol=1e-12, atol=0, equal_nan=True)
        elif left == "library_member":
            equal &= flags(paired[left], left, True).fillna(-1).eq(
                flags(paired[cr], cr, True).fillna(-1)).all()
        else:
            equal &= paired[left].eq(paired[cr]).fillna(False).where(
                paired[left].notna() | paired[cr].notna(), True).all()
    equal &= common.score1.eq(flags(paired.chen_sfari_score1, "class source score-1")).all()
    equal &= common.chen_hit.eq(flags(paired.chen_screen_hit, "class source Chen hit")).all()
    check("frozen_class_expression_projection", equal,
          "Regulator flags, targeted library, healthy CD4 expression, Chen hits and case membership agree", checks)

    uncertain = flags(common.CLASS_ANNOTATION_UNCERTAIN, "class uncertainty", True)
    class_label_usable = (common.class_label.isin(CLASS_BY_FLAGS.values()) & uncertain.eq(0) &
                          common.annotation_status.eq("CLASS_ANNOTATION_USABLE")).fillna(False)
    specified = pd.Series([CLASS_BY_FLAGS.get(key) for key in zip(common.TF, common.EpiFactors)],
                          dtype="string", index=common.index)
    check("class_usable_flags_agree",
          common.loc[class_label_usable, "class_label"].eq(specified.loc[class_label_usable]).all(),
          "Usable four-level regulator class agrees with both binary catalogue flags", checks)
    common["missing_class"] = ~class_label_usable
    common["library_member"] = flags(common.library_member, "targeted library", True)
    common["missing_library"] = common.library_member.isna()
    common["healthy_cd4_mean_logcpm"] = common.HEALTHY_CD4_MEAN_LOGCPM
    common["missing_expression"] = ~np.isfinite(common.healthy_cd4_mean_logcpm)
    common["missing_loeuf"] = ~np.isfinite(common.LOEUF_v2_1_1)
    reasons = {
        "missing_class": "REGULATOR_CLASS", "missing_library": "TARGETED_LIBRARY",
        "missing_expression": "HEALTHY_CD4_EXPRESSION", "missing_loeuf": "LOEUF",
    }
    complete = ~common[list(reasons)].any(axis=1)
    common["complete_case"] = complete
    common["exclusion_reasons"] = [
        "|".join(label for label, bad in zip(reasons.values(), flags_) if bad)
        for flags_ in common[list(reasons)].itertuples(index=False, name=None)
    ]
    common["discordant_complete_case"] = complete & common.D.ne(0)
    return common, reasons


def margins(d):
    n_counts = d.D.value_counts()
    a_counts = d.loc[d.score1.eq(1), "D"].value_counts()
    return dict(N=len(d), A=int(d.score1.sum()),
                n_plus=int(n_counts.get(1, 0)), n_zero=int(n_counts.get(0, 0)), n_minus=int(n_counts.get(-1, 0)),
                a_plus=int(a_counts.get(1, 0)), a_zero=int(a_counts.get(0, 0)), a_minus=int(a_counts.get(-1, 0)))


def validate_margins(m):
    if not (0 <= m["A"] <= m["N"] and
            sum(m[f"n_{d}"] for d in ("plus", "zero", "minus")) == m["N"] and
            sum(m[f"a_{d}"] for d in ("plus", "zero", "minus")) == m["A"] and
            all(0 <= m[f"a_{d}"] <= m[f"n_{d}"] for d in ("plus", "zero", "minus"))):
        raise ValueError(f"Invalid signed-score margins: {m}")


def signed_pmf(N, A, n_plus, n_zero, n_minus):
    """Shifted, dense PMF. First draw discordant count then draw its FOXP3-only count."""
    if A == 0:
        return 0, np.array([1.0])
    if A == N:
        return n_plus - n_minus, np.array([1.0])
    discordant = n_plus + n_minus
    if not discordant:
        return 0, np.array([1.0])
    least_minus = min(A, n_minus)
    least_S = max(0, A - least_minus - n_zero) - least_minus
    most_plus = min(A, n_plus)
    most_S = most_plus - max(0, A - most_plus - n_zero)
    mass = np.zeros(most_S - least_S + 1)
    counts = np.arange(max(0, A - n_zero), min(A, discordant) + 1)
    weights = hypergeom.pmf(counts, N, discordant, A)
    for k, weight in zip(counts, weights):
        positive = np.arange(max(0, k - n_minus), min(k, n_plus) + 1)
        conditional = hypergeom.pmf(positive, discordant, n_plus, k)
        mass[2 * positive - k - least_S] += weight * conditional
    return least_S, mass


def null_moments(m):
    N, A = m["N"], m["A"]
    ep, em = A * m["n_plus"] / N, A * m["n_minus"] / N
    mu = (m["n_plus"] - m["n_minus"]) / N
    s2 = (m["n_plus"] + m["n_minus"]) / N - mu * mu
    var = A * (N - A) / (N - 1) * s2 if N > 1 else 0.
    return ep, em, var


def pmf_stats(offset, mass, expected, variance, label):
    score = np.arange(offset, offset + len(mass))
    total = math.fsum(mass)
    mean_error = abs(float(np.dot(score, mass)) - expected)
    variance_error = abs(float(np.dot((score - expected) ** 2, mass)) - variance)
    if not (np.all(np.isfinite(mass)) and np.all(mass >= 0) and
            abs(total - 1) <= 1e-12 and mean_error <= 1e-10 and variance_error <= 1e-10):
        raise ArithmeticError(f"Non-normalized or inconsistent exact null: {label}, mass={total}, "
                              f"mean_error={mean_error}, variance_error={variance_error}")
    q_lo, q_hi = offset + np.searchsorted(np.cumsum(mass), (.025, .975), side="left")
    return dict(pmf_mass=total, pmf_mean_error=mean_error, pmf_variance_error=variance_error,
                null_low95=int(q_lo), null_high95=int(q_hi))


def exact_tests(data, analysis, checks):
    work = data if analysis == CONDITIONED else data.assign(stratum_key="ALL")
    strata = []
    offset = 0
    mass = np.array([1.0])
    for key, genes in work.groupby("stratum_key", sort=True):
        m = margins(genes)
        validate_margins(m)
        n, a = m["N"], m["A"]
        observed = m["a_plus"] - m["a_minus"]
        ep, em, variance = null_moments(m)
        k_offset, k_mass = signed_pmf(n, a, m["n_plus"], m["n_zero"], m["n_minus"])
        stats = pmf_stats(k_offset, k_mass, ep - em, variance, key)
        first = genes.iloc[0]
        strata.append(dict(analysis=analysis, stratum_key=key, **m,
                           class_label=first.class_label if analysis == CONDITIONED else None,
                           library_member=first.library_member if analysis == CONDITIONED else None,
                           expression_tertile=first.expression_tertile if analysis == CONDITIONED else None,
                           LOEUF_tertile=first.LOEUF_tertile if analysis == CONDITIONED else None,
                           expected_ASD_plus=ep, expected_ASD_minus=em,
                           observed_S=observed, expected_S=ep-em, var_S=variance,
                           stratum_has_controls=n > a, comparator_supported=0 < a < n,
                           informative=variance > 0, **stats))
        offset += k_offset
        mass = np.convolve(mass, k_mass)
    m = margins(data)
    validate_margins(m)
    check(f"{analysis}_strata_margins",
          all(sum(x[col] for x in strata) == value for col, value in m.items()),
          f"{len(strata)} occupied strata; pooled and summed margins agree", checks)
    ep = math.fsum(s["expected_ASD_plus"] for s in strata)
    em = math.fsum(s["expected_ASD_minus"] for s in strata)
    var = math.fsum(s["var_S"] for s in strata)
    observed = m["a_plus"] - m["a_minus"]
    stats = pmf_stats(offset, mass, ep - em, var, analysis)
    check(f"{analysis}_null_mass_and_moments", True,
          f"Every stratum and combined null passed mass and analytic moment checks; combined {stats}", checks)
    if not offset <= observed < offset + len(mass):
        raise ArithmeticError(f"Observed score outside the exact support for {analysis}: {observed}")
    cond_cov = dict(
        n_strata=len(strata), overlapping_strata=sum(s["comparator_supported"] for s in strata),
        informative_strata=sum(s["informative"] for s in strata),
        comparator_supported_ASD=sum(s["A"] for s in strata if s["comparator_supported"]),
        informative_ASD=sum(s["A"] for s in strata if s["informative"]),
    )
    svals = range(offset, offset + len(mass))
    pmf_table = [dict(analysis=analysis, S=S, probability=p) for S, p in zip(svals, mass)]
    p = math.fsum(mass[observed - offset:]) if var > 0 else np.nan
    result = dict(analysis=analysis, source_version=RELEASE,
                  estimand="Signed exclusive-hit preference in the jointly measured gene universe",
                  test="One-sided exact multivariate hypergeometric convolution" if analysis == CONDITIONED
                  else "One-sided exact multivariate hypergeometric tail", **m,
                  observed_S=observed, expected_S=ep-em, observed_minus_expected_S=observed-(ep-em),
                  var_S=var, expected_ASD_plus=ep, expected_ASD_minus=em, raw_P=p,
                  exact_status="ESTIMABLE" if var > 0 else "NON_INFORMATIVE_NULL",
                  multiplicity="Single fixed conditioned test; no adjustment" if analysis == CONDITIONED
                  else "Context only; no adjustment", **stats, **cond_cov)
    return result, strata, pmf_table


def verify_baseline(data, full, checks):
    summary = rd(data, "data/processed/direct_phenotype_preference_summary.tsv")
    check("legacy_direct_summary_unique", len(summary) == 1,
          "Unique frozen direct-phenotype summary", checks)
    prior = summary.iloc[0]
    for result, source, expected in (
        ("N", "common_universe", 18615), ("A", "asd_score1_common", 242),
        ("n_plus", "d_plus", 179), ("n_zero", "d_zero", 17562), ("n_minus", "d_minus", 874),
        ("observed_S", "S_observed", 3),
    ):
        check(f"historical_{result}", full[result] == prior[source] == expected,
              f"recomputed/released/expected {full[result]}/{prior[source]}/{expected}", checks)
    check("historical_asd_signed_margin",
          (full["a_plus"], full["a_zero"], full["a_minus"]) == (11, 223, 8),
          "Score-1 plus/zero/minus=11/223/8", checks)
    check("historical_total_hit_flags",
          (prior.chen_hits_common, prior.shifrut_hits_common) == (257, 952),
          "257 Chen and 952 Shifrut hits within shared genes", checks)
    check("historical_exact_tail_and_mass",
          math.isclose(full["raw_P"], prior.exact_null_p, abs_tol=1e-15, rel_tol=1e-12) and
          math.isclose(full["pmf_mass"], prior.exact_null_mass, abs_tol=1e-12),
          f"Exact P {full['raw_P']:.17g}; frozen {prior.exact_null_p:.17g}", checks)


def toy_checks(checks):
    # Enumerate sampled labelled genes to verify signed offsets and multiple limiting cases.
    outcomes = [((-1, -1, 0, 1, 1), 3), ((-1, -1, 0, 1, 1), 5),
                ((-1, -1, 0), 2), ((-1, -1), 1), ((0, 0, 0), 2), ((-1, 0, 1), 0)]
    for genes, draw in outcomes:
        counts = Counter(map(sum, combinations(genes, draw)))
        offset, mass = signed_pmf(len(genes), draw, genes.count(1), genes.count(0), genes.count(-1))
        expected = np.array([counts[s] / sum(counts.values()) for s in range(offset, offset + len(mass))])
        if not (min(counts) == offset and max(counts) == offset + len(mass) - 1 and
                np.allclose(mass, expected, rtol=0, atol=1e-14)):
            raise ArithmeticError(f"Signed PMF disagrees with exhaustive toy test: {genes}, A={draw}")
    check("pmf_exhaustive_toys", True, "Six labelled-toy subset enumerations including signed and deterministic supports",
          checks)


def set_tertiles(common, checks):
    d = common.loc[common.complete_case].copy()
    if d.empty:
        raise ValueError("No common genes retain the four prespecified covariates")
    rows = []
    for source, target in (("healthy_cd4_mean_logcpm", "expression"), ("LOEUF_v2_1_1", "LOEUF")):
        rank = d[source].rank(method="average")
        name = target + "_tertile"
        d[name] = np.minimum(3, np.ceil(3 * rank / len(d))).astype(int)
        d[target + "_average_rank"] = rank
        check(f"{target}_tertile_ties", d.groupby(source)[name].nunique().max() == 1,
              "Equal observed covariate values never split across tertiles", checks)
        for tertile, group in d.groupby(name, sort=True):
            rows.append(dict(section="TERTILE_BOUNDARY", field=source, subset=f"TERTILE_{tertile}", N=len(group),
                             A=int(group.score1.sum()), other_n=int((1 - group.score1).sum()),
                             low_value=group[source].min(), high_value=group[source].max(),
                             low_average_rank=group[target + "_average_rank"].min(),
                             high_average_rank=group[target + "_average_rank"].max(),
                             all_complete_case_N=len(d), definition=TERTILE_RULE))
    d["stratum_key"] = (d.class_label + "|" + d.library_member.map({0: "OUTSIDE", 1: "INSIDE"}) +
                        "|" + d.expression_tertile.astype(str) + "|" + d.LOEUF_tertile.astype(str))
    check("fixed_class_library_tertile_strata",
          d.stratum_key.notna().all() and d.stratum_key.nunique() <= 4 * 2 * 3 * 3,
          f"{d.stratum_key.nunique()} occupied fixed strata, set on all {len(d):,} complete cases", checks)
    for name in ["expression_tertile", "expression_average_rank", "LOEUF_tertile",
                 "LOEUF_average_rank", "stratum_key"]:
        common[name] = d[name].reindex(common.index)
    return d, rows


def coverage_row(section, field, subset, genes, basis, definition):
    N = len(genes)
    A = int(genes.score1.sum())
    To = len(basis) - int(basis.score1.sum())
    Ta = int(basis.score1.sum())
    return dict(section=section, field=field, subset=subset,
                N=N, A=A, other_n=N-A, denominator_N=len(basis), denominator_A=Ta,
                denominator_other_n=To, score1_fraction=A/Ta if Ta else np.nan,
                other_fraction=(N-A)/To if To else np.nan, definition=definition)


def balance_row(phase, d, field, weights):
    case = d.score1.eq(1)
    a, o = d.loc[case, field], d.loc[~case, field]
    wa, wo = weights.loc[case], weights.loc[~case]
    mu_a = float(np.average(a, weights=wa)) if len(a) else np.nan
    mu_o = float(np.average(o, weights=wo)) if len(o) else np.nan
    denominator = math.sqrt((a.var(ddof=1) + o.var(ddof=1)) / 2)
    delta = mu_a - mu_o
    smd = delta / denominator if denominator > 0 else np.nan
    if denominator == 0 and delta == 0:
        smd = 0.
    return dict(section="BALANCE", field=field, subset=phase, N=len(d), A=len(a), other_n=len(o),
                score1_mean=mu_a, other_mean=mu_o, mean_difference=delta,
                pooled_unweighted_SD=denominator, signed_SMD=smd, abs_SMD=abs(smd),
                score1_effective_n=float(wa.sum() ** 2 / (wa ** 2).sum()) if len(wa) else np.nan,
                other_effective_n=float(wo.sum() ** 2 / (wo ** 2).sum()) if len(wo) else np.nan,
                score1_sum_weights=float(wa.sum()), other_sum_weights=float(wo.sum()),
                residual_imbalance_flag=(abs(smd) > RESIDUAL_SMD_LIMIT) if np.isfinite(smd) else None,
                definition=f"{ATTRIBUTES[field]}. {SMD_RULE}")


def support_balance(common, d, strata, checks, reasons):
    by_key = pd.DataFrame(strata).set_index("stratum_key")
    is_supported = d.stratum_key.map(by_key.comparator_supported)
    is_informative = d.stratum_key.map(by_key.informative)
    common["comparator_supported"] = is_supported.astype("boolean").reindex(common.index)
    common["informative_stratum"] = is_informative.astype("boolean").reindex(common.index)
    supported = d.loc[is_supported].copy()
    info = d.loc[is_informative]
    unsupported = d.loc[~is_supported]
    comparator_weight = by_key.A / (by_key.N - by_key.A).where(by_key.N > by_key.A)
    control_weights = supported.stratum_key.map(comparator_weight)
    weights = control_weights.where(supported.score1.eq(0), 1.)
    check("common_support_weights", math.isclose(weights.loc[supported.score1.eq(0)].sum(),
                                                 supported.score1.sum(), abs_tol=1e-11),
          "Controls weighted A/(N−A) reproduce the case stratum frequencies", checks)
    rows = [
        coverage_row("COVERAGE", "common_universe", "FULL_MEASURED", common, common, "Measured in both screens"),
        coverage_row("COVERAGE", "common_universe", "COMPLETE_CASE", d, common,
                     "All four fixed gene attributes usable; no imputation"),
        coverage_row("COVERAGE", "common_universe", "EXCLUDED", common.loc[~common.complete_case], common,
                     "At least one missing required attribute"),
        coverage_row("COVERAGE", "common_universe", "DISCORDANT_COMPLETE_CASE", d.loc[d.D.ne(0)], d,
                     "Same complete-case population restricted to exclusive hits"),
        coverage_row("POSITIVITY", "strata", "MIXED_CASE_CONTROL", supported, d,
                     "Strata contain score-1 cases and other genes"),
        coverage_row("POSITIVITY", "strata", "NO_MIXED_CASE_CONTROL", unsupported, d,
                     "Complement of positivity subset; all strata retained in exact test"),
        coverage_row("POSITIVITY", "strata", "INFORMATIVE_SIGNED_NULL", info, d,
                     "Within-stratum signed null has positive variance"),
        coverage_row("POSITIVITY", "strata", "NON_INFORMATIVE_SIGNED_NULL", d.loc[~is_informative], d,
                     "Zero within-stratum variance; strata retained in exact test"),
    ]
    for column, label in reasons.items():
        rows.append(coverage_row("MISSINGNESS", label, "OVERLAPPING_REASON", common.loc[common[column]], common,
                                 "Missing or unusable attribute; row may have multiple reasons"))
    for labels, subset in common.groupby("exclusion_reasons", sort=True):
        rows.append(coverage_row("MISSINGNESS_PATTERN", labels or "NONE", "EXCLUSIVE_PATTERN", subset, common,
                                 "Exclusive combination of missing attributes"))
    smds = {}
    for field in ATTRIBUTES:
        for phase, part, w in (("COMPLETE_CASE_UNWEIGHTED", d, pd.Series(1., index=d.index)),
                               ("COMMON_SUPPORT_UNWEIGHTED", supported, pd.Series(1., index=supported.index)),
                               ("WITHIN_STRATUM_WEIGHTED", supported, weights)):
            b = balance_row(phase, part, field, w)
            rows.append(b)
            if phase == "WITHIN_STRATUM_WEIGHTED":
                smds[field] = b["signed_SMD"]
    return rows, smds


def firth_tests(d, release_firth, checks):
    d = d.loc[d.D.ne(0)].copy()
    d["asd"] = d.score1
    d["hit"] = d.D.eq(1).astype(int)
    d["targeted_library"] = d.library_member.astype(int)
    d["expression"] = d.healthy_cd4_mean_logcpm
    d["LOEUF"] = d.LOEUF_v2_1_1
    fits, diagnostics = [], []
    for analysis, predictors in ((UNADJUSTED_FIRTH, ()), (ADJUSTED_FIRTH, PREDICTORS)):
        retained = [name for name in predictors if d[name].nunique() > 1]
        dropped = [name for name in predictors if name not in retained]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                fit, coefs = release_firth.fit(d, retained)
            except (ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
                fit, coefs = (dict(OR=np.nan, CI_low=np.nan, CI_high=np.nan, P=np.nan, converged=False,
                                   interval_status="NUMERIC_FAILURE", warnings=str(exc)), [])
        actual_warnings = " | ".join(dict.fromkeys(str(w.message) for w in caught))
        context = dict(analysis=analysis, source_version=RELEASE,
                       estimand="Odds of FOXP3-only rather than TCR/proliferation-only, conditional on discordance",
                       test="Firth/Jeffreys penalized logistic regression; two-sided penalized profile LR",
                       formula="FOXP3_only ~ score1" + "".join(f" + {x}" for x in retained),
                       dropped_constant_terms="|".join(dropped),
                       multiplicity="Supporting sensitivity; no adjustment with the exact estimand",
                       **margins(d))
        fit = dict(fit)
        fit["raw_P"] = fit.pop("P")
        diagnostics.append(context | fit | dict(record_type="MODEL", emitted_runtime_warnings=actual_warnings))
        for term in coefs:
            diagnostics.append(context | term | dict(
                record_type="COEFFICIENT", emitted_runtime_warnings="",
                scaling_mean=float(d[term["term"]].mean()) if term["term"] in ("expression", "LOEUF") else None,
                scaling_population_SD=float(d[term["term"]].std(ddof=0)) if term["term"] in ("expression", "LOEUF") else None))
        status = "SENSITIVITY_ESTIMABLE" if fit["interval_status"] == "ESTIMABLE" else "SENSITIVITY_NOT_FULLY_ESTIMABLE"
        fits.append(context | dict(OR=fit["OR"], CI_low=fit["CI_low"], CI_high=fit["CI_high"],
                                   raw_P=fit["raw_P"], model_status=fit["interval_status"],
                                   support_status=status, warnings=fit["warnings"]))
    check("firth_discordant_same_subset",
          len(fits) == 2 and all(fits[0][k] == fits[1][k] for k in margins(d)),
          f"Both fits use {len(d):,} discordant complete cases and {d.score1.sum():,} score-1 genes", checks)
    for model in fits:
        checks.append(dict(branch="B", check=model["analysis"] + "_model_status", status=(
            "PASS" if model["model_status"] == "ESTIMABLE" else "NON_ESTIMABLE"),
            observation=model["model_status"] + (": " + model["warnings"] if model["warnings"] else "")))
    return fits, diagnostics


def run(data, code, out):
    from baseline_release import verify_release
    verify_release(data)
    sys.path.insert(0, str((code / "scripts/gene_property_adjustment").resolve(strict=True)))
    release_firth = import_module("firth_models")
    checks = []
    toy_checks(checks)
    universe, reasons = load_universe(data, release_firth, checks)
    full, full_strata, full_pmf = exact_tests(universe, FULL, checks)
    verify_baseline(data, full, checks)
    d, bins = set_tertiles(universe, checks)
    cc, cc_strata, cc_pmf = exact_tests(d, COMPLETE, checks)
    adjusted, adj_strata, adj_pmf = exact_tests(d, CONDITIONED, checks)
    cov, smds = support_balance(universe, d, adj_strata, checks, reasons)

    for row in (full, cc, adjusted):
        row["score1_complete_over_full"] = row["A"] / full["A"]
        row["score1_informative_over_analysis"] = row["informative_ASD"] / row["A"] if row["A"] else np.nan
        row["score1_comparator_over_analysis"] = row["comparator_supported_ASD"] / row["A"] if row["A"] else np.nan
    full["support_status"] = "HISTORICAL_UNADJUSTED"
    cc["support_status"] = "COMPLETE_CASE_RESTRICTION_CONTEXT"
    coverage_pass = (adjusted["score1_complete_over_full"] >= CASE_COVERAGE_MIN and
                     adjusted["score1_informative_over_analysis"] >= INFORMATIVE_FRACTION_MIN)
    balanced = all(np.isfinite(smd) and abs(smd) <= RESIDUAL_SMD_LIMIT for smd in smds.values())
    adjusted.update(coverage_gate_pass=coverage_pass, residual_balance_pass=balanced,
                    weighted_expression_SMD=smds["healthy_cd4_mean_logcpm"],
                    weighted_LOEUF_SMD=smds["LOEUF_v2_1_1"])
    if adjusted["exact_status"] != "ESTIMABLE":
        support = "NON_INFORMATIVE_NULL"
    elif not (coverage_pass and balanced):
        support = "LIMITED_SUPPORT"
    elif adjusted["raw_P"] < .05:
        support = "CONDITIONAL_PREFERENCE_SUPPORTED"
    else:
        support = "CONDITIONAL_PREFERENCE_NOT_SUPPORTED"
    adjusted["support_status"] = support
    models, diagnostic_rows = firth_tests(d, release_firth, checks)

    columns = ["approved_symbol", "hgnc_id", "legacy_common_symbol", "case_only_symbol_match",
               "score1", "chen_hit", "CHEN_DIRECTION",
               "shifrut_hit", "SHIFRUT_DIRECTION", "D", "derived_category", "TF", "EpiFactors",
               "class_label", "annotation_status", "library_member", *ATTRIBUTES,
               *reasons, "complete_case", "exclusion_reasons", "expression_average_rank",
               "LOEUF_average_rank", "expression_tertile", "LOEUF_tertile", "stratum_key",
               "comparator_supported", "informative_stratum", "discordant_complete_case"]
    out.mkdir(parents=True, exist_ok=True)
    tsv(out / "PHENOTYPE_PREFERENCE_R9.tsv", [full, cc, adjusted, *models])
    tsv(out / "PHENOTYPE_PREFERENCE_STRATA.tsv", [*full_strata, *cc_strata, *adj_strata])
    tsv(out / "PHENOTYPE_PREFERENCE_MEMBERS.tsv.gz", universe[columns].sort_values("approved_symbol"))
    tsv(out / "PHENOTYPE_PREFERENCE_NULL.tsv", [*full_pmf, *cc_pmf, *adj_pmf])
    tsv(out / "PHENOTYPE_PREFERENCE_COVERAGE_BALANCE.tsv", [*cov, *bins])
    tsv(out / "PHENOTYPE_PREFERENCE_MODEL_DIAGNOSTICS.tsv", diagnostic_rows)
    tsv(out / "PHENOTYPE_PREFERENCE_CHECKS.tsv", checks)
    print(f"Shared preference historical P={full['raw_P']:.16g}, complete-case P={cc['raw_P']:.16g}, "
          f"conditioned P={adjusted['raw_P']:.16g} ({support}).")
    for result in models:
        print(f"{result['analysis']}: {result['model_status']}; "
              f"OR={result['OR']:.10g}, 95% profile CI [{result['CI_low']:.10g}, {result['CI_high']:.10g}], "
              f"P={result['raw_P']:.10g} (N={result['N']}, ASD={result['A']}).")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path)
    parser.add_argument("code_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    run(args.data_root, args.code_root, args.output_dir)


if __name__ == "__main__":
    main()
