#!/usr/bin/env python3
"""Run primary, continuous and property-conditioned GSE286470 comparisons."""

import argparse
from hashlib import sha256
import math
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, hypergeom
from scipy.stats.contingency import odds_ratio

from umhoefer_sources import (
    GENE_SHEET, GUIDE_SHEET, PROPERTIES, SOURCE_LINK,
    SOURCE_STUDY, candidate_guides, functional_panel, release_tables, shifrut_panel,
    soft_samples, verify, workbook_genes, write_tsv,
)

FIRTH_HASH = "e679b82442e6c711c835e8a830a338b8458eec88dc82199428b331c599610692"
CATEGORIES = ("Treg-only", "Tconv-only", "Both", "Neither")
CLASS_BY_FLAGS = {
    (0, 0): "NEITHER_LISTED", (1, 0): "TF_ONLY",
    (0, 1): "EPIFACTOR_ONLY", (1, 1): "BOTH_LISTED",
}
PREDICTORS = ("TF", "EpiFactors", "expression", "LOEUF", "targeted_library")
N_PERM = 100_000
N_BOOT = 20_000
RANKS = ("delta_abs_rank", "delta_signed_rank")


def ranked(effect):
    return (effect.rank(method="average") - 1) / (len(effect) - 1)


def build_universe(screens, checks):
    measure_map = {
        "source_row": "source_row", "id": "raw_id", "pos|lfc": "author_pos_lfc",
        "neg|lfc": "author_neg_lfc", "pos|fdr": "pos_fdr", "neg|fdr": "neg_fdr",
        "pos|p-value": "pos_p", "neg|p-value": "neg_p", "pos|rank": "pos_rank", "neg|rank": "neg_rank",
        "num": "guide_n", "hit": "hit", "direction": "effect_sign", "mapping_route": "mapping_route",
        "signed_foxp3_bin_lfc": "signed_foxp3_bin_lfc", "sfari_score1": "raw_asd",
        "hgnc_id": "hgnc_id", "FOXP3_REGULATOR_DIRECTION": "frozen_direction",
    }
    attrs = list(PROPERTIES.values())
    sides = []
    for context, measured in screens.items():
        frame = measured[list(measure_map) + attrs].rename(columns=measure_map)
        sides.append(frame.add_prefix(context.lower() + "_"))
    common = sides[0].join(sides[1], how="inner", validate="one_to_one").sort_index()
    properties_agree = all(
        (common[f"treg_{v}"].eq(common[f"tconv_{v}"])
         | (common[f"treg_{v}"].isna() & common[f"tconv_{v}"].isna())).all() for v in attrs
    )
    verify(common["treg_raw_asd"].eq(common["tconv_raw_asd"]).all() and properties_agree
           and common["treg_hgnc_id"].eq(common["tconv_hgnc_id"]).all(),
           "common_asd_id_and_attributes",
           "Both measured screen rows use the identical frozen ASD label, HGNC key and gene-property values",
           checks)
    common["asd"] = common.pop("treg_raw_asd").astype(int)
    common = common.drop(columns="tconv_raw_asd")
    common = common.rename(columns={"treg_hgnc_id": "hgnc_id"}).drop(columns="tconv_hgnc_id")
    common = common.rename(columns={f"treg_{v}": v for v in attrs}).drop(columns=[f"tconv_{v}" for v in attrs])
    paired = 2 * common["treg_hit"].astype(int) + common["tconv_hit"].astype(int)
    common["category"] = paired.map({0: "Neither", 1: "Tconv-only", 2: "Treg-only", 3: "Both"})
    common["discordant"] = common["category"].isin(CATEGORIES[:2])
    for key in ("treg", "tconv"):
        effect = common[f"{key}_signed_foxp3_bin_lfc"]
        common[f"{key}_abs_rank"] = ranked(effect.abs())
        common[f"{key}_signed_rank"] = ranked(effect)
    common["delta_abs_rank"] = common["treg_abs_rank"] - common["tconv_abs_rank"]
    common["delta_signed_rank"] = common["treg_signed_rank"] - common["tconv_signed_rank"]
    expected_class = pd.Series(
        [CLASS_BY_FLAGS.get((a, b)) for a, b in zip(common["TF"], common["EpiFactors"])], index=common.index)
    class_valid = common["class_label"].isin(CLASS_BY_FLAGS.values()) & common["class_uncertain"].eq(False) & (
        common["annotation_status"].eq("CLASS_ANNOTATION_USABLE"))
    verify(expected_class[class_valid].eq(common.loc[class_valid, "class_label"]).all(),
           "frozen_regulator_class_parity", "Usable four-level class matches TF and EpiFactors catalogue flags", checks)
    missing_fields = {
        "regulator_class": ~class_valid, "healthy_CD4_expression": ~np.isfinite(common["expression"]),
        "LOEUF": ~np.isfinite(common["LOEUF"]), "targeted_library": common["targeted_library"].isna(),
    }
    for name, flag in missing_fields.items():
        common[f"missing_{name}"] = flag
    common["adjustment_exclusions"] = [
        "|".join(name for name, bad in zip(missing_fields, bits) if bad)
        for bits in zip(*missing_fields.values())
    ]
    common["complete_case"] = common["adjustment_exclusions"].eq("")
    cc = common["complete_case"]
    for col in ("expression", "LOEUF"):
        ranks = common.loc[cc, col].rank(method="average")
        common.loc[cc, f"{col}_tertile"] = np.ceil(3 * ranks / cc.sum()).clip(upper=3).astype(int)
        checks.append(dict(check=f"{col}_tertile_ties", status="PASS",
                           observation="Identical values tied by their average rank"))
        assert common[cc].groupby(col)[f"{col}_tertile"].nunique().max() == 1
    k = common.loc[cc]
    common.loc[cc, "stratum"] = (
        k["class_label"] + "|" + k["targeted_library"].astype(int).astype(str) + "|" +
        k["expression_tertile"].astype(int).astype(str) + "|" + k["LOEUF_tertile"].astype(int).astype(str))
    both = common["category"].eq("Both")
    directions_available = ~common[["treg_signed_foxp3_bin_lfc", "tconv_signed_foxp3_bin_lfc"]].eq(0).any(axis=1)
    same_dir = np.sign(common["treg_signed_foxp3_bin_lfc"]).eq(np.sign(common["tconv_signed_foxp3_bin_lfc"]))
    common.loc[both, "both_hit_direction"] = np.select(
        [same_dir[both] & directions_available[both], directions_available[both]],
        ["same", "opposite"], default="zero estimate")
    for k in ("treg", "tconv"):
        common[f"{k}_hit"] = common[f"{k}_hit"].astype(int)
    return common


def grouped_count(label, mask, asd, **fields):
    cases = int(asd[mask].sum())
    count = int(mask.sum())
    return dict(category=label, asd_n=cases, non_asd_n=count - cases, total_n=count, **fields)


def table_counts(common, screens, all_cases, checks):
    out = [grouped_count(label, common["category"].eq(label), common["asd"], section="shared_hit_categories")
           for label in CATEGORIES]
    verify(sum(r["total_n"] for r in out) == len(common)
           and sum(r["asd_n"] for r in out) == common["asd"].sum(),
           "four_categories_partition_U", "Four hit categories sum to shared U and preserve the ASD denominator",
           checks)
    for ctx, measured in screens.items():
        asd = measured["sfari_score1"].astype(int)
        out.append(grouped_count(ctx + " measured", asd.notna(), asd, section="measurement"))
        out.append(grouped_count(ctx + " measured only", ~measured.index.isin(common.index),
                                 asd, section="measurement"))
        out.append(dict(category=ctx + " frozen ASD missing from experiment",
                        asd_n=len(all_cases) - int(asd.sum()), non_asd_n=None, total_n=None,
                        section="measurement"))
    out.append(grouped_count("Measured in both", common["asd"].notna(), common["asd"], section="measurement"))
    out.append(dict(category="Frozen ASD absent from shared U", asd_n=len(all_cases) - int(common["asd"].sum()),
                    non_asd_n=None, total_n=None, section="measurement"))
    for direction in ("same", "opposite", "zero estimate"):
        out.append(grouped_count("Both hits, " + direction, common["both_hit_direction"].eq(direction),
                                 common["asd"], section="both_hit_directions"))
    return pd.DataFrame(out)


def two_by_two(asd, outcome):
    row = np.asarray(asd, dtype=bool)
    col = np.asarray(outcome, dtype=bool)
    return np.array([[(row & col).sum(), (row & ~col).sum()],
                     [(~row & col).sum(), (~row & ~col).sum()]], dtype=int)


def describe_fisher(mat):
    case_out, case_other, ctrl_out, ctrl_other = mat.ravel()
    details = dict(n=mat.sum(), asd_n=mat[0].sum(), non_asd_n=mat[1].sum(),
                   asd_hit=case_out, asd_other=case_other, non_asd_hit=ctrl_out, non_asd_other=ctrl_other,
                   OR=np.nan, CI_low=np.nan, CI_high=np.nan, P=np.nan, sample_OR=np.nan, upper_tail_P=np.nan,
                   status="NOT_ESTIMABLE")
    if min(mat.sum(axis=1).min(), mat.sum(axis=0).min()) == 0:
        return details
    result = fisher_exact(mat, alternative="two-sided")
    exact = odds_ratio(mat, kind="conditional")
    ci = exact.confidence_interval(confidence_level=.95, alternative="two-sided")
    return details | dict(OR=exact.statistic, CI_low=ci.low, CI_high=ci.high, P=result.pvalue,
                          sample_OR=result.statistic, upper_tail_P=fisher_exact(mat, alternative="greater").pvalue,
                          status="ESTIMABLE")


def main_tests(common, screens, checks):
    d = common[common["discordant"]]
    box = two_by_two(d["asd"], d["treg_hit"])
    a, b, c, e = box.ravel()
    method = "Two-sided Fisher exact; conditional-MLE OR with equal-tail exact 95% CI"
    rows = [dict(analysis="discordant_primary", role="primary", background="Shared U, author-hit discordant D",
                 outcome="Treg-only vs Tconv-only", method=method, **describe_fisher(box),
                 asd_treg_only=a, asd_tconv_only=b, non_asd_treg_only=c, non_asd_tconv_only=e,
                 interpretation="Relative exclusive FOXP3 hit type among common discordant genes")]
    expected = {"Treg": (1349, 33, 5, 26), "Tconv": (1349, 33, 4, 34)}
    for ctx, tab in screens.items():
        full = describe_fisher(two_by_two(tab["sfari_score1"], tab["hit"]))
        observed = tuple(int(full[x]) for x in ("n", "asd_n", "asd_hit", "non_asd_hit"))
        verify(observed == expected[ctx], "original_full_margins_" + ctx,
               f"Observed N, ASD N, ASD hits, other hits = {observed}; exact frozen parity", checks)
        for label, data, info in (
            ("individual", full, "Original actual measured " + ctx + " screen universe"),
            ("shared_U", describe_fisher(two_by_two(common["asd"], common[f"{ctx.lower()}_hit"])),
             "Shared U without discordance restriction"),
        ):
            rows.append(dict(analysis=f"{ctx.lower()}_enrichment_{label}", role="descriptive",
                             background=info, outcome=f"{ctx} author hit vs non-hit",
                             method=method, interpretation="Within one screen; does not compare the two cell contexts",
                             **data))
    return pd.DataFrame(rows)


def bh(p_values):
    v = np.asarray(p_values)
    order = np.argsort(v)
    q_sorted = np.minimum.accumulate((len(v) * v[order] / np.arange(1, 1 + len(v)))[::-1])[::-1]
    result = np.empty(len(v))
    result[order] = np.minimum(q_sorted, 1)
    return result


def threshold_tests(common):
    rows = []
    for frac in (.05, .10, .20):
        k = math.ceil(len(common) * frac)
        scores = [common[f"{ctx}_signed_foxp3_bin_lfc"].abs() for ctx in ("treg", "tconv")]
        cuts = [s.nlargest(k).min() for s in scores]
        a, b = [s.ge(cut) for s, cut in zip(scores, cuts)]
        d = a != b
        counts = two_by_two(common.loc[d, "asd"], a[d])
        rows.append(dict(analysis=f"equal_fraction_top_{int(frac * 100)}pct", role="sensitivity",
                         background="Shared U, top-fraction discordant genes",
                         target_fraction=frac, common_U=len(common), target_k=k, tie_rule="include all cutoff ties",
                         treg_cut=cuts[0], tconv_cut=cuts[1], treg_n=a.sum(), tconv_n=b.sum(),
                         treg_actual_fraction=a.mean(), tconv_actual_fraction=b.mean(),
                         both_n=(a & b).sum(), neither_n=(~a & ~b).sum(),
                         method="Two-sided Fisher exact; |signed selection LFC| ranking, not author FDR hits",
                         interpretation="Relative exclusive absolute-effect ranking; a threshold sensitivity",
                         **describe_fisher(counts)))
    frame = pd.DataFrame(rows)
    frame["P_BH"] = np.where(frame["P"].notna(), bh(frame["P"].fillna(1)), np.nan)
    return frame


def sampled_null(values, a, rng):
    """Gene-label permutations for both rank endpoints, preserving sample counts."""
    n = len(values)
    arr = np.empty((N_PERM, values.shape[1]))
    raw_sum = values.sum(axis=0)
    for start in range(0, N_PERM, 1024):
        stop = min(N_PERM, start + 1024)
        if a == 1:
            index = rng.integers(0, n, size=(stop - start, 1))
        else:
            permutations = np.broadcast_to(np.arange(n), (stop - start, n))
            index = rng.permuted(permutations, axis=1)[:, :a]
        subtotal = values[index].sum(axis=1)
        arr[start:stop] = subtotal / a - (raw_sum - subtotal) / (n - a)
    return arr


def matched_bootstrap(values, flags, rng, both_statistics):
    """Percentile bootstrap of the two observed gene-label groups at fixed group sizes."""
    first, second = values[flags], values[~flags]
    length = 2 if both_statistics else 1
    out = np.empty((N_BOOT, values.shape[1], length))
    for start in range(0, N_BOOT, 1024):
        stop = min(N_BOOT, start + 1024)
        a = first[rng.integers(0, len(first), size=(stop - start, len(first)))]
        b = second[rng.integers(0, len(second), size=(stop - start, len(second)))]
        out[start:stop, :, 0] = a.mean(axis=1) - b.mean(axis=1)
        if both_statistics:
            out[start:stop, :, 1] = np.median(a, axis=1) - np.median(b, axis=1)
    return out


def mc_summary(null, observed):
    b = int(np.sum(np.abs(null) >= abs(observed) - 1e-12))
    p = (b + 1) / (len(null) + 1)
    return dict(P=p, mc_extreme_n=b, B=len(null), mc_se=np.sqrt(p * (1 - p) / len(null)),
                null_q025=np.quantile(null, .025), null_q975=np.quantile(null, .975))


def continuous_tests(common, seed_perm, seed_boot):
    data = common[list(RANKS)].to_numpy()
    flags = common["asd"].to_numpy(dtype=bool)
    n_case = int(flags.sum())
    null = sampled_null(data, n_case, seed_perm)
    boot = matched_bootstrap(data, flags, seed_boot, True)
    rows = []
    for i, endpoint in enumerate(RANKS):
        for j, (label, estimator) in enumerate((("mean", np.mean), ("median", np.median))):
            group_stats = {}
            for name, sample in (("asd", data[flags, i]), ("non_asd", data[~flags, i])):
                group_stats.update({name + "_value": estimator(sample),
                                    name + "_q25": np.quantile(sample, .25),
                                    name + "_q75": np.quantile(sample, .75)})
            observed = group_stats["asd_value"] - group_stats["non_asd_value"]
            low, high = np.quantile(boot[:, i, j], (.025, .975))
            method = (("Two-sided ASD-label permutation; " if j == 0 else "") +
                      f"{N_BOOT:,} gene-level percentile bootstrap resamples within labels")
            rows.append(dict(analysis=f"{endpoint.removeprefix('delta_')}_{label}", n=len(common), asd_n=n_case,
                             non_asd_n=len(common) - n_case, effect=observed, CI_low=low, CI_high=high,
                             effect_unit="Treg minus Tconv rank on 0–1 scale, ASD minus non-ASD",
                             method=method, background="Shared U; both effects available", role="supportive",
                             interpretation=("Within-screen absolute gene effect ranks" if i == 0
                                             else "Within-screen signed FOXP3 selection ranks"),
                             **group_stats, **(mc_summary(null[:, i], observed) if j == 0 else {})))
    return pd.DataFrame(rows), pd.DataFrame(null, columns=RANKS).rename(columns=lambda c: "unadjusted_" + c)


def property_histogram(d, checks):
    """Conditional exact fixed-margin hypergeometric convolution on discordant genes."""
    mass = np.array([1.0])
    offset = 0
    expected, variance = 0.0, 0.0
    observed = 0
    informative_asd = 0
    for _, g in d.groupby("stratum", sort=True):
        n = len(g)
        a = int(g["asd"].sum())
        successes = int(g["treg_hit"].sum())
        k = int((g["treg_hit"] * g["asd"]).sum())
        lo = max(0, a - (n - successes))
        hi = min(a, successes)
        local = hypergeom.pmf(np.arange(lo, hi + 1), n, successes, a)
        mass = np.convolve(mass, local)
        offset += lo
        observed += k
        expected += a * successes / n
        if lo < hi:
            informative_asd += a
            variance += a * (successes / n) * (1 - successes / n) * (n - a) / (n - 1)
    axis = np.arange(offset, offset + len(mass))
    verify(np.isclose(mass.sum(), 1) and np.isclose(axis @ mass, expected)
           and np.isclose((axis - expected) ** 2 @ mass, variance),
           "conditioned_discordant_exact_pmf",
           "Convolved exact PMF has total mass 1 and agrees with analytical null mean and variance",
           checks)
    if variance > 0:
        p = mass[np.abs(axis - expected) >= abs(observed - expected) - 1e-12].sum()
        upper = mass[axis >= observed].sum()
    else:
        p, upper = np.nan, np.nan
    lo, hi = axis[np.minimum(np.searchsorted(np.cumsum(mass), [.025, .975]), len(axis) - 1)]
    summary = dict(observed_count=observed, expected_count=expected, effect=observed - expected,
                   effect_unit="observed minus null-expected Treg-only ASD genes",
                   observed_over_expected=observed / expected if expected > 0 else np.nan,
                   null_variance=variance, null_q025=lo, null_q975=hi, P=p, upper_tail_P=upper,
                   informative_asd_n=informative_asd, status="ESTIMABLE" if variance > 0 else "NOT_ESTIMABLE")
    return summary, pd.DataFrame({"ASD_Treg_only_count": axis, "null_mass": mass})


def coverage_and_balance(common):
    cc = common[common["complete_case"]]
    sets = {"continuous": (common, cc),
            "binary_discordant": (common[common["discordant"]], cc[cc["discordant"]])}
    rows, strata = [], []
    for analysis, (intended, complete) in sets.items():
        groups = complete.groupby("stratum", sort=True)
        supported = groups["asd"].transform("nunique").eq(2)
        continuous_support = complete[supported]
        weights = pd.Series(np.nan, index=continuous_support.index)
        outcome_informative = pd.Series(False, index=complete.index)
        for key, g in groups:
            count = len(g)
            a = int(g["asd"].sum())
            hits = int(g["treg_hit"].sum()) if analysis == "binary_discordant" else None
            eligible = 0 < a < count
            informative = eligible and (analysis == "continuous" or 0 < hits < count)
            outcome_informative.loc[g.index] = informative
            if eligible:
                weights.loc[g.index] = np.where(g["asd"], 1., a / (count - a))
            strata.append(dict(analysis=analysis, stratum=key, n=count, asd_n=a, non_asd_n=count - a,
                               tregrank_only_n=hits,
                               asd_treg_only_n=(g["asd"] * g["treg_hit"]).sum() if hits is not None else None,
                               comparator_supported=eligible, outcome_informative=informative,
                               expression_min=g["expression"].min(), expression_max=g["expression"].max(),
                               LOEUF_min=g["LOEUF"].min(), LOEUF_max=g["LOEUF"].max()))
        for name, g in (("full", intended), ("complete_case", complete), ("supported_strata", continuous_support),
                        ("informative_strata", complete[outcome_informative])):
            a = int(g["asd"].sum())
            rows.append(dict(section="coverage", analysis=analysis, subset=name, n=len(g), asd_n=a,
                             non_asd_n=len(g) - a, retained_fraction=len(g) / len(intended),
                             ASD_retained_fraction=a / intended["asd"].sum(),
                             total_strata=g["stratum"].nunique()))
        for col in [c for c in intended if c.startswith("missing_")]:
            for label in (0, 1):
                n_label = intended["asd"].eq(label).sum()
                count = int((intended[col] & intended["asd"].eq(label)).sum())
                rows.append(dict(section="overlapping_missingness", analysis=analysis,
                                 subset="ASD" if label else "non-ASD", field=col,
                                 group_n=n_label, missing_n=count, missing_fraction=count / n_label))
        for reason, g in intended.groupby("adjustment_exclusions", dropna=False):
            rows.append(dict(section="exclusive_missingness", analysis=analysis,
                             subset=reason or "complete_case", n=len(g), asd_n=g["asd"].sum(),
                             non_asd_n=len(g) - g["asd"].sum()))
        for phase, data, w in (("complete_case_unweighted", complete, None),
                               ("supported_unweighted", continuous_support, None),
                               ("supported_stratum_weighted", continuous_support, weights)):
            for field in PREDICTORS:
                a = data[data["asd"].eq(1)][field].astype(float)
                b = data[data["asd"].eq(0)][field].astype(float)
                if w is None:
                    means = a.mean(), b.mean()
                elif len(a) and len(b):
                    means = np.average(a, weights=w.loc[a.index]), np.average(b, weights=w.loc[b.index])
                else:
                    means = np.nan, np.nan
                den = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
                smd = (means[0] - means[1]) / den if den > 0 else np.nan
                rows.append(dict(section="covariate_balance", analysis=analysis, subset=phase, field=field,
                                 asd_n=len(a), non_asd_n=len(b), ASD_mean=means[0], non_ASD_mean=means[1],
                                 ASD_sample_SD=a.std(ddof=1), non_ASD_sample_SD=b.std(ddof=1), SMD=smd,
                                 flagged_abs_SMD_gt_0_1=abs(smd) > .1 if np.isfinite(smd) else None,
                                 note="SMD undefined for zero variance or too few cases" if np.isnan(smd) else ""))
    return pd.DataFrame(rows), pd.DataFrame(strata)


def firth_comparisons(d, fitter):
    rows, coefs = [], []
    d = d.copy().assign(hit=d["treg_hit"])
    for model, candidates in (("discordant_complete_case_Firth_ASD_only", ()),
                              ("discordant_complete_case_Firth_adjusted", PREDICTORS)):
        used = [name for name in candidates if d[name].nunique() > 1]
        dropped = [name for name in candidates if name not in used]
        fit, coeff = fitter(d, used)
        for param in coeff:
            coefs.append(dict(analysis=model, **param))
        rows.append(dict(analysis=model, n=len(d), asd_n=d["asd"].sum(), non_asd_n=len(d) - d["asd"].sum(),
                         background="Discordant common U with complete frozen property fields",
                         role="supportive binary sensitivity",
                         method="R7/R9 Jeffreys-penalized logistic model; two-sided profile likelihood ratio P",
                         formula="Treg-only ~ ASD" + "".join(" + " + p for p in used),
                         dropped_constant_terms="|".join(dropped), status=fit["interval_status"],
                         interpretation="Conditional odds of Treg-only vs Tconv-only among discordant genes",
                         **fit))
    return rows, coefs


def adjusted_continuous(complete, rng_perm, rng_boot):
    groups = []
    for _, g in complete.groupby("stratum", sort=True):
        if g["asd"].nunique() > 1:
            groups.append(g)
    if not groups:
        return [dict(analysis="continuous_delta_abs_rank_conditioned", status="NOT_ESTIMABLE",
                     note="No property stratum supports both label groups")], pd.DataFrame()
    supported = pd.concat(groups)
    a_total = int(supported["asd"].sum())
    mean = np.zeros(len(RANKS))
    boot = np.zeros((N_BOOT, len(RANKS)))
    null = np.zeros((N_PERM, len(RANKS)))
    for g in groups:
        values = g[list(RANKS)].to_numpy()
        flags = g["asd"].to_numpy(dtype=bool)
        a = int(flags.sum())
        w = a / a_total
        mean += w * (values[flags].mean(axis=0) - values[~flags].mean(axis=0))
        boot += w * matched_bootstrap(values, flags, rng_boot, False)[:, :, 0]
        null += w * sampled_null(values, a, rng_perm)
    rows = []
    for i, endpoint in enumerate(RANKS):
        low, high = np.quantile(boot[:, i], (.025, .975))
        rows.append(dict(analysis=f"continuous_{endpoint}_conditioned", status="ESTIMABLE",
                         n=len(supported), asd_n=a_total, non_asd_n=len(supported) - a_total,
                         effect=mean[i], CI_low=low, CI_high=high,
                         effect_unit="ASD minus non-ASD within-stratum average Δrank",
                         method="Case-count-weighted within-stratum mean contrasts; two-sided fixed-margin label permutations; gene bootstrap",
                         background="Covariate-complete U in strata containing both label groups",
                         informative_strata=len(groups), role="supportive continuous sensitivity",
                         interpretation="Conditional cross-screen within-gene rank difference",
                         **mc_summary(null[:, i], mean[i])))
    return rows, pd.DataFrame(null, columns=RANKS).rename(columns=lambda c: "adjusted_" + c)


def numeric_validation(primary, sensitivity, continuous, adjusted, checks):
    for name, table in (("primary", primary), ("sensitivity", sensitivity),
                        ("continuous", continuous), ("adjusted", adjusted)):
        p = table["P"].dropna()
        verify(p.between(0, 1).all(), name + "_p_valid", "All estimable P values lie in [0,1]", checks)
    for name, frame in (("primary", primary), ("equal_fraction", sensitivity)):
        bounds = frame[frame["status"].eq("ESTIMABLE")]
        verify((bounds["CI_low"].le(bounds["OR"]) & bounds["OR"].le(bounds["CI_high"])).all(),
               name + "_exact_ci_valid", "Reported conditional ORs lie within exact 95% confidence intervals", checks)
    verify(sensitivity["P_BH"].between(0, 1).all()
           and sensitivity["P_BH"].ge(sensitivity["P"] - 1e-12).all() and len(sensitivity) == 3,
           "equal_fraction_BH_family", "Three fixed sensitivity fractions only enter the BH family", checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--download-dir", type=Path, required=True)
    parser.add_argument("--firth-dir", type=Path, required=True)
    parser.add_argument("--shifrut-table", type=Path, required=True)
    args = parser.parse_args()
    seed = 20261125
    release, firth_path = args.data, args.firth_dir
    shifrut_table, download = args.shifrut_table, args.download_dir
    output = args.out.resolve()
    audit, results, logs = [output / d for d in ("01_audit", "03_results", "08_logs")]
    for dirpath in (audit, results, logs):
        dirpath.mkdir(parents=True, exist_ok=True)
    streams = [np.random.default_rng(k) for k in np.random.SeedSequence(seed).spawn(4)]
    checks = []
    sys.path.insert(0, str(firth_path))
    from firth_models import fit
    try:
        if sha256((firth_path / "firth_models.py").read_bytes()).hexdigest() != FIRTH_HASH:
            raise ValueError("Firth-model implementation differs from the R9-locked public script")
        inputs = release_tables(release, checks)
        gene_path = download / "GSE286470_Combined_gene_summary.xlsx"
        guide_path = download / "GSE286470_Combined_grna_summary.xlsx"
        soft = download / "GSE286470_family.soft.gz"
        samples = soft_samples(soft, checks)
        write_tsv(samples, audit / "R10_UMHOEFER_SAMPLE_METADATA.tsv")
        screens, controls, genes, all_cases = workbook_genes(gene_path, inputs, checks)
        write_tsv(controls, audit / "R10_UMHOEFER_CONTROL_EXCLUSIONS.tsv")
        guides, guide_rows, guide_frame, guide_disagreements = candidate_guides(guide_path, screens, checks)
        write_tsv(guides, results / "R10_UMHOEFER_GUIDE_AUDIT.tsv")
        write_tsv(guide_rows, results / "R10_UMHOEFER_CANDIDATE_GUIDE_STATS.tsv")
        write_tsv(guide_disagreements, audit / "R10_UMHOEFER_GUIDE_DISCREPANCIES.tsv")
        schema = [
            dict(dataset="GSE286470", file=gene_path.name, sheet=GENE_SHEET, row_unit="gene/screen pair",
                 data_rows=len(genes) + len(controls), n_cols=len(genes.columns) - 3, key="celltype, id",
                 columns="|".join(col for col in genes if col not in ("source_row", "hit", "signed_foxp3_bin_lfc")),
                 value_unit="MAGeCK gene summary LFC, P/FDR, guide count and ranks", source=SOURCE_STUDY),
            dict(dataset="GSE286470", file=guide_path.name, sheet=GUIDE_SHEET, row_unit="guide/screen pair",
                 data_rows=len(guide_frame), n_cols=len(guide_frame.columns) - 1, key="celltype, sgrna",
                 columns="|".join(col for col in guide_frame if col != "source_row"),
                 value_unit="MAGeCK normalized counts (unlabelled positions), pooled guide LFC, P/FDR", source=SOURCE_STUDY),
            dict(dataset="GSE286470", file=soft.name, row_unit="GSM (FOXP3 sort bin within cell-context donor)",
                 data_rows=len(samples), key="GSM", value_unit="metadata", source=SOURCE_STUDY),
        ]
        write_tsv(schema, audit / "R10_UMHOEFER_DATASET_SCHEMA.tsv")
        shifrut = shifrut_panel(shifrut_table, inputs, checks)
        common = build_universe(screens, checks)
        write_tsv(common.reset_index(), results / "R10_TREG_TCONV_UNIVERSE.tsv.gz")
        panel = functional_panel(screens, inputs, shifrut, common, checks)
        write_tsv(panel, results / "R10_14_GENE_FUNCTION_PROFILE.tsv")
        coverage = table_counts(common, screens, all_cases, checks)
        write_tsv(coverage, results / "R10_TREG_TCONV_COUNTS.tsv")
        primary = main_tests(common, screens, checks)
        write_tsv(primary, results / "R10_TREG_TCONV_PRIMARY.tsv")
        sensitivity = threshold_tests(common)
        write_tsv(sensitivity, results / "R10_TREG_TCONV_SENSITIVITY.tsv")
        raw_continuous, null = continuous_tests(common, *streams[:2])
        write_tsv(raw_continuous, results / "R10_TREG_TCONV_CONTINUOUS.tsv")
        balance, strata = coverage_and_balance(common)
        write_tsv(balance, results / "R10_TREG_TCONV_MATCH_BALANCE.tsv")
        write_tsv(strata, results / "R10_TREG_TCONV_CONDITIONING_STRATA.tsv")
        cc = common[common["complete_case"]]
        d = cc[cc["discordant"]]
        exact, mass = property_histogram(d, checks)
        write_tsv(mass, results / "R10_TREG_TCONV_ADJUSTED_EXACT_NULL.tsv")
        adjusted_rows = [dict(analysis="discordant_complete_case_conditioned_exact", n=len(d),
                              asd_n=d["asd"].sum(), non_asd_n=len(d) - d["asd"].sum(),
                              method="Exact within-stratum hypergeometric convolution; two-sided distance from null mean",
                              role="supportive conditional exact sensitivity",
                              background="Discordant common U with complete frozen property fields",
                              interpretation="Conditional ASD Treg-only count when stratum margins are fixed", **exact)]
        firth, coefficients = firth_comparisons(d, fit)
        adjusted_rows.extend(firth)
        write_tsv(coefficients, results / "R10_TREG_TCONV_FIRTH_COEFFICIENTS.tsv")
        adjusted, adj_null = adjusted_continuous(cc, *streams[2:])
        adjusted_rows.extend(adjusted)
        adjusted_tab = pd.DataFrame(adjusted_rows)
        write_tsv(adjusted_tab, results / "R10_TREG_TCONV_ADJUSTED.tsv")
        write_tsv(pd.concat([null, adj_null], axis=1), results / "R10_TREG_TCONV_CONTINUOUS_NULL.tsv.gz")
        numeric_validation(primary, sensitivity, raw_continuous, adjusted_tab, checks)
        paths = [gene_path, guide_path, soft, firth_path / "firth_models.py", shifrut_table,
                 Path(__file__).with_name("strict_xlsx.py")]
        provenance = [dict(path=p.name, bytes=p.stat().st_size, sha256=sha256(p.read_bytes()).hexdigest())
                      for p in paths]
        write_tsv(provenance, audit / "R10_UMHOEFER_INPUT_HASHES.tsv")
        print(f"Python {platform.python_version()} ({platform.machine()}), pandas {pd.__version__}, "
              f"NumPy {np.__version__}, SciPy {scipy.__version__}; seed={seed}")
        print(f"Source verification: {SOURCE_LINK}; GSE286470 donor and bin metadata from local SOFT")
        print(primary.iloc[0].to_json())
        print(raw_continuous[["analysis", "effect", "P"]].to_string(index=False))
        print(adjusted_tab[["analysis", "n", "asd_n", "effect", "OR", "P", "status"]].to_string(index=False))
    finally:
        write_tsv(checks, logs / "R10_UMHOEFER_CHECKS.tsv")


if __name__ == "__main__":
    main()
