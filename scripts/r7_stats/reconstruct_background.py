"""Rebuild state calls and exact strata from the released author-derived snapshot.

Expression references in this module begin at frozen per-gene reference
summaries. They follow the original R1 grouping and threshold definitions.
"""
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from firth_models import CORE, DESIGNS, DETECTED, REF, p_adjust_two, rd, truth
from release_inputs import released_file

EXACT_KEYS = ["screen", "reference", "threshold", "units", "set", "analysis"]
STRATA_COUNTS = ["N", "H", "A", "K"]
STATE_GROUPS = {
    "CHEN": ("itreg", (2, 2)),
    "UMHOEFER_TREG": ("htreg", (3,)),
    "UMHOEFER_TCONV": ("htconv", (3,)),
}
DETECTION_STATES = [
    "NOT_MEASURED", "EXPRESSION_REFERENCE_MISSING", "BELOW_THRESHOLD",
    "DETECTABLE_HIT", "DETECTABLE_NONHIT",
]


def detection_status(measured, known, detected, hit):
    return np.select(
        [~measured, ~known, ~detected, hit],
        DETECTION_STATES[:-1],
        default=DETECTION_STATES[-1],
    )


def _assert_values(actual, frozen, label, rtol=0, atol=0):
    if not np.allclose(actual, frozen, rtol=rtol, atol=atol, equal_nan=True):
        raise ValueError(f"Recalculation disagrees with the released value: {label}")


def rebuild_background(data_dir, genes):
    rebuilt = genes.copy()
    report = []
    for screen, (prefix, minimum_votes) in STATE_GROUPS.items():
        d = genes.loc[genes["screen"].eq(screen)]
        measured = d["is_measured"]
        expression = d[f"{prefix}_mean_log_expression"]
        _assert_values(expression, d["primary_mean_log_expression"], f"{screen} CPM expression")
        hpa_col = "hpa_tconv_max_nTPM" if screen == "UMHOEFER_TCONV" else "hpa_T-reg"
        hpa = d[hpa_col]
        hpa_log = np.log2(hpa + 1)
        _assert_values(hpa_log, d["hpa_mean_log_expression"], f"{screen} HPA expression",
                       rtol=1e-12, atol=1e-14)
        fdr_hit = d["pos|fdr"].le(.05) | d["neg|fdr"].le(.05)
        if not d.loc[measured, "hit"].eq(fdr_hit.loc[measured].astype(int)).all():
            raise ValueError(f"Frozen measured hit flags disagree with directional FDR: {screen}")
        for threshold, units in DESIGNS:
            if units == "CPM":
                ref = REF[screen]
                status_name = f"primary_status_{threshold:g}"
                group_pass = [
                    d[f"{prefix}_group{i}_n_ge_{threshold:g}"].ge(votes)
                    for i, votes in enumerate(minimum_votes, start=1)
                ]
                detected = pd.concat(group_pass, axis=1).all(axis=1)
                known = expression.notna()
                frozen_vote_result = truth(d[f"{prefix}_detect_{threshold:g}"])
                if not detected.eq(frozen_vote_result).all():
                    raise ValueError(f"Frozen votes disagree with detection indicator: {screen} CPM {threshold:g}")
            else:
                ref = "HPA_IMMUNE"
                status_name = "hpa_status_1"
                detected = hpa.ge(1)
                known = hpa.notna()
            state = detection_status(measured, known, detected, fdr_hit)
            if not d[status_name].eq(state).all():
                raise ValueError(f"Frozen detection states disagree with vote reconstruction: {screen} {ref} {threshold:g}")
            rebuilt.loc[d.index, status_name] = state
            for label in DETECTION_STATES:
                mask = state == label
                report.append({
                    "screen": screen, "reference": ref, "threshold": threshold,
                    "units": units, "detection_status": label, "n_genes": int(mask.sum()),
                    "asd_n": int(d.loc[mask, "asd"].sum()),
                    "asd_hits": int((d.loc[mask, "asd"] * d.loc[mask, "hit"]).sum()),
                    "other_hits": int(((1 - d.loc[mask, "asd"]) * d.loc[mask, "hit"]).sum()),
                })
    if set(rebuilt["screen"]) != set(STATE_GROUPS):
        raise ValueError("The background must contain only the three frozen screen universes")
    membership = rd(released_file(data_dir, "mechanism/ALLELIC_MECHANISM_MASTER.tsv"))
    if len(membership) != 245 or not membership["approved_symbol"].is_unique:
        raise ValueError("Expected 245 unique frozen score-1 members in the mechanism table")
    symbols = set(rebuilt.loc[rebuilt["asd"].eq(1), "approved_symbol"])
    if symbols != set(membership["approved_symbol"]):
        raise ValueError("Score-1 membership differs between the background and mechanism table")
    eligible = set(membership.loc[truth(membership["lof_eligible"]), "approved_symbol"])
    if len(eligible) != 104:
        raise ValueError("LoF sensitivity requires the 104 frozen qualifying members")
    return rebuilt, pd.DataFrame(report), eligible


def rebuild_strata(data_dir, genes, lof_eligible):
    frozen = rd(released_file(data_dir, "inputs/R1_CONDITIONING_STRATA.tsv"))
    output = []
    for key, table in frozen.groupby(EXACT_KEYS, dropna=False, sort=False):
        screen, ref, threshold, units, gene_set, method = key
        x = genes.loc[genes["screen"].eq(screen) & genes["is_measured"]].copy()
        if ref == "FROZEN_GSE90600":
            expression = "HEALTHY_CD4_MEAN_LOGCPM"
        elif ref == "FROZEN_MEASURED":
            expression = "HEALTHY_CD4_MEAN_LOGCPM"
        else:
            expression = "hpa_mean_log_expression" if units == "nTPM" else "primary_mean_log_expression"
            status = "hpa_status_1" if units == "nTPM" else f"primary_status_{threshold:g}"
            x = x.loc[x[status].isin(DETECTED)].copy()
        if gene_set == "LOF_COMPATIBLE":
            eligible = x["approved_symbol"].isin(lof_eligible)
            x = x.loc[x["asd"].eq(0) | eligible].copy()
            x["asd"] = eligible.loc[x.index].astype(int)
        elif gene_set != "SFARI_SCORE1":
            raise ValueError(f"Unexpected analysis set: {gene_set}")

        x = x.loc[x["class_label"].notna()].copy()
        joint = method == "strict_joint"
        if joint:
            x = x.loc[np.isfinite(x[expression]) & np.isfinite(x["log10_span"]) &
                      np.isfinite(x["guide_count"]) & x["guide_count"].gt(0)].copy()
        elif method != "class_conditioned":
            raise ValueError(f"Unexpected conditional analysis: {method}")
        stratum = x["class_label"].astype(str)
        if screen == "CHEN":
            stratum += "|" + x["targeted_library"].map({1: "INSIDE", 0: "OUTSIDE"})

        if joint:
            # R1 uses average ranks and a right-closed 1–3/4/≥5 guide grouping.
            for col in (expression, "log10_span"):
                tercile = np.minimum(3, np.ceil(3 * x[col].rank(method="average") / len(x)))
                stratum += "|" + tercile.astype(int).astype(str)
            guide_class = pd.cut(
                x["guide_count"], [0, 3, 4, np.inf],
                labels=["1_TO_3", "EXACTLY_4", "5_OR_MORE"],
            )
            stratum += "|" + guide_class.astype(str)
        x["stratum_key"] = stratum
        x["asd_hit"] = x["asd"] * x["hit"]
        counts = x.groupby("stratum_key", sort=True).agg(
            N=("hit", "size"), H=("hit", "sum"), A=("asd", "sum"), K=("asd_hit", "sum"),
        )
        target = table.set_index("stratum_key")[STRATA_COUNTS].sort_index()
        try:
            pd.testing.assert_frame_equal(counts, target, check_dtype=False)
        except AssertionError as error:
            raise ValueError(f"Stratum counts disagree with frozen R1 assignment: {key}") from error
        merged = table.drop(columns=STRATA_COUNTS).join(counts, on="stratum_key", validate="one_to_one")
        row = merged.loc[:, table.columns]
        if not row["randomizable_labels"].pipe(truth).eq(row["A"].gt(0) & row["A"].lt(row["N"])).all():
            raise ValueError(f"Randomizable-label indicators disagree: {key}")
        if not row["nondegenerate_hit_distribution"].pipe(truth).eq(
            row["H"].gt(0) & row["H"].lt(row["N"])
        ).all():
            raise ValueError(f"Nondegenerate-hit indicators disagree: {key}")
        output.append(row)
    return pd.concat(output, ignore_index=True)


def recalculate_exact_tests(strata, frozen):
    rows = []
    for key, d in strata.groupby(EXACT_KEYS, dropna=False):
        pmf = np.array([1.])
        for row in d.itertuples():
            prob = hypergeom.pmf(np.arange(min(row.A, row.H) + 1), row.N, row.H, row.A)
            _assert_values(prob.sum(), 1, f"{key} stratum PMF", atol=1e-12)
            pmf = np.convolve(pmf, prob)
        _assert_values(pmf.sum(), 1, f"{key} combined PMF", atol=1e-12)
        index = np.arange(len(pmf))
        observed = int(d["K"].sum())
        expectation = np.dot(index, pmf)
        variance = np.dot((index - expectation) ** 2, pmf)
        informative = truth(d["randomizable_labels"]) & truth(d["nondegenerate_hit_distribution"])
        shared = d["A"].gt(0) & d["A"].lt(d["N"])
        lower, upper = np.searchsorted(np.cumsum(pmf), [.025, .975])
        rows.append(dict(
            zip(EXACT_KEYS, key),
            total_n=int(d["N"].sum()), asd_n=int(d["A"].sum()), asd_hit=observed,
            other_n=int((d["N"] - d["A"]).sum()), other_hit=int((d["H"] - d["K"]).sum()),
            informative_asd_n=int(d.loc[informative, "A"].sum()),
            informative_strata=int(informative.sum()),
            comparator_supported_asd_n=int(d.loc[shared, "A"].sum()),
            expected_hits=expectation, variance_hits=variance,
            effect=observed / expectation if expectation > 0 else np.nan,
            null_lo=lower, null_hi=upper,
            p=pmf[observed:].sum() if informative.any() else np.nan,
            adjusted_p=np.nan,
        ))
    out = pd.DataFrame(rows)
    if len(out) != 30:
        raise ValueError("The released exact-test grid must contain all 30 frozen comparisons")
    family = out["reference"].where(out["units"].ne("CPM"), "CHEN_CPM")
    for _, group in out.loc[out["screen"].isin(CORE)].groupby(
        [family, "threshold", "units", "set", "analysis"], dropna=False
    ):
        if set(group["screen"]) != set(CORE) or len(group) != 2:
            raise ValueError("Each Holm family must contain exactly the two original core screens")
        out.loc[group.index, "adjusted_p"] = p_adjust_two(group["p"])
    checked = out.merge(frozen, on=EXACT_KEYS, suffixes=("", "_frozen"), how="left",
                        validate="one_to_one", indicator=True)
    if not checked["_merge"].eq("both").all():
        raise ValueError("The exact-test grid does not match the released reference rows")
    for col in out.columns.difference(EXACT_KEYS):
        _assert_values(checked[col], checked[f"{col}_frozen"], col, rtol=1e-9, atol=1e-12)
    return out
