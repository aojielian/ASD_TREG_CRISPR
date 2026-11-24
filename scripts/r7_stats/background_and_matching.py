"""Recalculate frozen R1 enrichment, R7 coverage and R7 matching."""
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import binomtest, fisher_exact
from scipy.stats.contingency import odds_ratio

from firth_models import ALL_FIELDS, CORE, DESIGNS, REF, margins, rd, reference_data, truth, write
from release_inputs import released_file
from reconstruct_background import recalculate_exact_tests

CLASSES = ["BOTH_LISTED", "TF_ONLY", "EPIFACTOR_ONLY", "NEITHER_LISTED"]


def fisher_stats(d):
    a, b, c, e = margins(d)
    row = dict(total_n=len(d), ASD_hits=a, ASD_nonhits=b, ASD_n=a+b,
               other_hits=c, other_nonhits=e, other_n=c+e, OR=np.nan, CI_low=np.nan, CI_high=np.nan,
               P=np.nan)
    if min(a+b, c+e, a+c, b+e) == 0:
        row["status"] = "NO_OUTCOME_OR_GROUP_VARIATION"
        return row
    table = [[a,b], [c,e]]
    ft = fisher_exact(table, alternative="greater")
    estimate = odds_ratio(table, kind="conditional")
    ci = estimate.confidence_interval(.95)
    row.update(P=ft.pvalue, OR=estimate.statistic, CI_low=ci.low, CI_high=ci.high,
               status="ESTIMABLE" if 0 < estimate.statistic < np.inf else "BOUNDARY_OR")
    return row


def rowsame(df, spec):
    select = np.ones(len(df), dtype=bool)
    for k, v in spec.items():
        select &= df[k].isna() if pd.isna(v) else df[k].eq(v)
    return df[select]


def stage_row(spec, stage, members, source):
    a, b, c, e = margins(members)
    return spec | dict(stage=stage, total_n=len(members), ASD_n=a+b, ASD_hits=a,
                       other_hits=c, other_n=c+e, source=source)


def subset_stratum_rows(s, spec, kind):
    return rowsame(s, spec | {"set":"SFARI_SCORE1", "analysis":kind})


def aggregate_strata_rows(spec, kind, s, informative=False):
    if informative:
        s = s[truth(s["randomizable_labels"]) & truth(s["nondegenerate_hit_distribution"])]
    n,h,a,k = s[["N","H","A","K"]].sum()
    return spec | dict(stage=kind, total_n=n, ASD_n=a, ASD_hits=k, other_hits=h-k, other_n=n-a,
                       source="R1 CONDITIONING_STRATA.tsv; verified exact N/H/A/K")


def balance(s, phase, spec):
    rows = []
    for col in ALL_FIELDS:
        x = s.loc[s.asd.eq(1), col].dropna()
        y = s.loc[s.asd.eq(0), col].dropna()
        delta = x.mean() - y.mean()
        denom = np.sqrt((x.var(ddof=0) + y.var(ddof=0)) / 2)
        smd = delta / denom if denom > 0 else (0 if abs(delta) < 1e-12 else np.nan)
        rows.append(spec | dict(phase=phase, covariate=col, asd_observed=len(x), other_observed=len(y),
                                asd_mean=x.mean(), other_mean=y.mean(), SMD=smd,
                                weight_rule="Unit weight per observed gene", source="R1 gene table and R2 eligibility"))
    return rows


def match(cases, measured, spec):
    features = ["expression", "LOEUF", "log10_span"]
    d = cases.sort_values("approved_symbol").copy()
    z = (d[features] - d[features].mean()) / d[features].std(ddof=0)
    a = d[d.asd.eq(1)]
    c = d[d.asd.eq(0)]
    exact_names = ["class_label"]
    if spec["screen"] == "CHEN":
        exact_names.append("targeted_library")
    feasible = np.ones((len(a), len(c)), dtype=bool)
    for key in exact_names:
        feasible &= a[key].to_numpy()[:, None] == c[key].to_numpy()[None, :]
    distance = np.zeros(feasible.shape)
    for feature in features:
        diff = z.loc[a.index, feature].to_numpy()[:, None] - z.loc[c.index, feature].to_numpy()[None,:]
        feasible &= np.abs(diff) <= .5
        distance += diff ** 2
    # Valid edge distance is <=sqrt(.75). Dummy penalty exceeds total valid cost.
    dummy = len(a) + 1
    cost = np.hstack([np.where(feasible, np.sqrt(distance), 2 * dummy),
                      np.full((len(a), len(a)), dummy)])
    indices, assigned = linear_sum_assignment(cost)
    pairs = {}
    for i,j in zip(indices,assigned):
        if j < len(c) and feasible[i,j]:
            pairs[i] = j
    positivity = []
    for i, case in enumerate(a.itertuples()):
        j = pairs.get(i)
        control = None if j is None else c.iloc[j]
        matches = int(feasible[i].sum())
        positivity.append(spec | dict(
            ASD_gene=case.approved_symbol, ASD_hit=case.hit, eligible_controls=matches,
            control_gene="" if j is None else control.approved_symbol,
            control_hit=np.nan if j is None else control.hit,
            weight=0 if j is None else 1,
            status="MATCHED" if j is not None else "NO_CALIPER_SUPPORT" if not matches else "NO_ASSIGNMENT_WITHOUT_REPLACEMENT",
            distance=np.nan if j is None else np.sqrt(distance[i,j]),
            expression_delta_z=np.nan if j is None else z.loc[a.index[i], "expression"] - z.loc[c.index[j], "expression"],
            loeuf_delta_z=np.nan if j is None else z.loc[a.index[i], "LOEUF"] - z.loc[c.index[j], "LOEUF"],
            span_delta_z=np.nan if j is None else z.loc[a.index[i], "log10_span"] - z.loc[c.index[j], "log10_span"]))
    x = pd.concat([a.iloc[list(pairs)], c.iloc[list(pairs.values())]])
    bal = balance(x, "R2_MATCHED", spec)
    smd = np.array([r["SMD"] for r in bal])
    kept = len(pairs)
    original_fraction = kept / measured.asd.sum()
    complete_fraction = kept / len(a)
    n10 = sum(a.iloc[i].hit == 1 and c.iloc[j].hit == 0 for i,j in pairs.items())
    n01 = sum(a.iloc[i].hit == 0 and c.iloc[j].hit == 1 for i,j in pairs.items())
    pass_ = original_fraction >= .8 and complete_fraction >= .7 and np.isfinite(smd).all() and np.max(np.abs(smd)) <= .1
    summary = spec | dict(common_case_asd_n=len(a), original_measured_asd_n=int(measured.asd.sum()),
        matched_asd=kept, matched_controls=kept, unmatched_ASD=len(a)-kept,
        no_caliper_support=int((feasible.sum(axis=1) == 0).sum()),
        max_abs_SMD=np.nanmax(np.abs(smd)), coverage_original=original_fraction,
        coverage_model_complete=complete_fraction, ASD_hit_control_nonhit=n10,
        ASD_nonhit_control_hit=n01, matched_p=binomtest(n10,n10+n01,.5,alternative="greater").pvalue
        if pass_ and n10+n01 else np.nan, status="USABLE" if pass_ and n10+n01 else "P_WITHHELD_INSUFFICIENT_OVERLAP_OR_BALANCE",
        weight_rule="Deterministic 1:1 without replacement; unit weights")
    return positivity, bal, summary


def direction_table(data_dir, genes):
    expected = rd(released_file(data_dir, "results/DIRECTION_ORIGINAL_TO_MANUSCRIPT.tsv"))
    fields = ["screen", "approved_symbol", "raw_primary_id", "pos|lfc", "neg|lfc",
              "author_hit", "is_measured", "FOXP3_REGULATOR_DIRECTION",
              "standardized_protein_sort_direction"]
    candidates = sorted(expected["approved_symbol"].unique())
    if len(candidates) != 17 or len(expected) != 51:
        raise ValueError("Expected 17 released candidate symbols across three screen contexts")
    native = genes.loc[genes["approved_symbol"].isin(candidates), fields[:-1]].copy()
    tags = {"SUPPORTS_FOXP3": "KO_REDUCES_FOXP3", "RESTRAINS_FOXP3": "KO_INCREASES_FOXP3",
            "NOT_SIGNIFICANT": "MEASURED_NOT_HIT", "BIDIRECTIONAL_AMBIGUOUS": "BIDIRECTIONAL_AMBIGUOUS"}
    native[fields[-1]] = native["FOXP3_REGULATOR_DIRECTION"].map(tags).fillna("NOT_MEASURED")
    key = ["approved_symbol", "screen"]
    native = native.sort_values(key, ignore_index=True)
    archived = expected[fields].sort_values(key, ignore_index=True)
    for col in ("author_hit", "is_measured"):
        native[col] = native[col].astype("boolean")
        archived[col] = archived[col].astype("boolean")
    pd.testing.assert_frame_equal(native, archived, check_dtype=False, rtol=1e-12)
    return native


def main(data_dir, out, genes, strata):
    old = rd(released_file(data_dir, "results/R1_FROZEN_SENSITIVITY_GRID.tsv"))
    mqc = rd(released_file(data_dir, "inputs/R1_MATCHING_QC.tsv"))
    model = rd(out / "results/LOEUF_ADJUSTED_MODELS.tsv")
    s = strata
    flow, recalc, breakdown, balances, coverage_qc, composition, positivity, support = ([] for _ in range(8))
    for screen in REF:
        measured0 = genes[genes.screen.eq(screen) & genes.is_measured].copy()
        smeta = dict(screen=screen,reference="FROZEN_MEASURED",threshold=np.nan,units=np.nan)
        original = fisher_stats(measured0)
        recalc.append(smeta | original)
        assert len(measured0) == (18966 if screen == "CHEN" else 1349)
        missing_score = genes.loc[genes.screen.eq(screen) & ~genes.is_measured & genes.asd.eq(1),"approved_symbol"]
        coverage_qc.append(dict(screen=screen, check="measured canonical", total=len(measured0),
                                ASD_n=int(measured0.asd.sum()), ASD_hits=int((measured0.asd*measured0.hit).sum()),
                                mapping_route="R1 accepted canonical, original raw design",
                                detail=("|".join(missing_score) if screen == "CHEN" else
                                        f"{len(missing_score)} absent from targeted library")))
        for threshold, unit in DESIGNS:
            ref, measured, detected, cases, status = reference_data(genes, screen, threshold, unit)
            spec = dict(screen=screen,reference=ref,threshold=threshold,units=unit)
            flow.append(stage_row(spec,"ASSAYED",measured,"R1 measured canonical QC"))
            reference_present = measured[~measured[status].eq("EXPRESSION_REFERENCE_MISSING")]
            flow.append(stage_row(spec, "REFERENCE_MATCHED", reference_present, "R1 reference ID mapping plus state value"))
            for subs, label in (("EXPRESSION_REFERENCE_MISSING","REFERENCE_MISSING"),("BELOW_THRESHOLD","BELOW_THRESHOLD")):
                flow.append(stage_row(spec,label,measured[measured[status].eq(subs)],"R1 frozen detection calls"))
            flow.append(stage_row(spec,"DETECTABLE",detected,"R1 frozen detection calls"))
            recalc.append(spec | fisher_stats(detected))
            strict = subset_stratum_rows(s,spec,"strict_joint")
            flow.append(aggregate_strata_rows(spec,"R1_STRICT_COMPLETE",strict))
            flow.append(aggregate_strata_rows(spec,"R1_STRICT_INFORMATIVE",strict,True))
            if screen in CORE:
                flow.append(stage_row(spec,"R2_COMMON_MODEL_COMPLETE",cases,"Frozen HGNC, TF/EpiFactors, expression, LOEUF, span and guides"))
            for label in CLASSES:
                all_d = measured[measured.class_label.eq(label)]
                exp_d = detected[detected.class_label.eq(label)]
                for name,d in [("MEASURED",all_d),("DETECTABLE",exp_d)]:
                    breakdown.append(spec | dict(regulator_class=label, background=name,
                                                  measured_class_n=len(all_d), detectable_class_n=len(exp_d)) | fisher_stats(d))
            map_field = "hpa_reference_mapping" if unit == "nTPM" else {
                "CHEN":"itreg_reference_mapping","UMHOEFER_TREG":"htreg_reference_mapping",
                "UMHOEFER_TCONV":"htconv_reference_mapping"}[screen]
            z = measured.assign(reference_map=measured[map_field].fillna("UNAVAILABLE"),
                                status=measured[status], cclass=measured.class_label.fillna("UNCLASSIFIED"))
            grouped = z.groupby(["mapping_route","reference_map","cclass","targeted_library","status"],dropna=False)
            for key, d in grouped:
                composition.append(spec | dict(zip(["screen_map","reference_map","class","targeted_library","status"],key)) |
                    dict(total_n=len(d), ASD_n=int(d.asd.sum()), total_hits=int(d.hit.sum()), ASD_hits=int((d.asd*d.hit).sum())))
            if threshold == 1 and unit == "CPM":
                balances += balance(measured,"MEASURED",spec) + balance(detected,"DETECTABLE",spec)
                if screen in CORE:
                    balances += balance(cases,"R2_MODEL_COMPLETE",spec)
                    got, bs, summ = match(cases,measured,spec)
                    positivity.extend(got)
                    balances.extend(bs)
                    support.append(summ)
                found = rowsame(mqc,spec | {"analysis":"strict_joint","set":"SFARI_SCORE1",
                                             "record_type":"BALANCE"})
                for q in found.itertuples():
                    cov = {"gene_span":"log10_span","LOEUF_v2_1_1":"LOEUF"}.get(q.covariate,q.covariate)
                    balances.append(spec | dict(phase="R1_STRICT_"+q.phase, covariate=cov,
                        asd_observed=q.asd_n, other_observed=q.control_n,
                        asd_mean=q.asd_mean, other_mean=q.control_mean,
                        SMD=q.standardized_mean_diff, weight_rule=q.weight_rule, source="R1 R3_MATCHING_QC.tsv"))
    for row in recalc:
        key = {k:row[k] for k in ("screen","reference","threshold","units")}
        a = rowsame(old,key | {"analysis":"unadjusted","set":"SFARI_SCORE1"}).iloc[0]
        assert (row["ASD_hits"],row["ASD_n"],row["other_hits"],row["other_n"]) == (a.asd_hit,a.asd_n,a.other_hit,a.other_n)
        assert np.isclose(row["P"],a.p,rtol=1e-9)
        assert np.allclose([row["OR"],row["CI_low"],row["CI_high"]],[a.effect,a.ci_low,a.ci_high],rtol=.002)
        row.update(R1_OR=a.effect,R1_CI_low=a.ci_low,R1_CI_high=a.ci_high,R1_P=a.p)
    exact = recalculate_exact_tests(s, old)
    write(exact, out / "results/R1_EXACT_RECALCULATION.tsv")
    delta = []
    pairs = {"CHEN":(14,242,11348,174,3.95,1.41,.0896,173,-.855,.174),
             "UMHOEFER_TREG":(5,33,726,18,10.86,2.33,.0389,18,-.836,.670)}
    for screen, e in pairs.items():
        ref = REF[screen]
        oldraw = rowsame(old,dict(screen=screen,reference="FROZEN_MEASURED",analysis="unadjusted")).iloc[0]
        raw, clas, strict = [rowsame(old,dict(screen=screen,reference=ref,threshold=1, units="CPM",
                 analysis=a,set="SFARI_SCORE1")).iloc[0] for a in ["unadjusted","class_conditioned","strict_joint"]]
        b = pd.DataFrame(balances)
        smd = b[b.screen.eq(screen) & b.phase.eq("R1_STRICT_AFTER_STRATUM_WEIGHTING_SHARED_SUPPORT")].set_index("covariate").SMD
        checks = {"measured ASD hits":(e[0],oldraw.asd_hit,0),
            "measured ASD":(e[1],oldraw.asd_n,0),"detectable genes":(e[2],raw.total_n,0),
            "detectable ASD":(e[3],raw.asd_n,0), "detectable ASD hits":(e[0],raw.asd_hit,0),
            "detectable OR":(e[4],raw.effect,.005), "class Holm":(.0160,clas.adjusted_p,.00005),
            "strict O/E":(e[5],strict.effect,.005), "strict Holm":(e[6],strict.adjusted_p,.00005),
            "strict ASD complete cases":(e[7],strict.asd_n,0),
            "strict LOEUF SMD":(e[8],smd["LOEUF"],.001),
            "strict log10 span SMD":(e[9],smd["log10_span"],.001)}
        for measure,(expected,actual,tol) in checks.items():
            assert abs(actual-expected) <= tol, (screen,measure,expected,actual)
            delta.append(dict(screen=screen,metric=measure, rounded_reference_value=expected, validated_frozen_value=actual,
                              status="AGREES_AT_REFERENCE_PRECISION", R2_scope="R1 reference result retained"))
        mm = rowsame(model,dict(screen=screen,units="CPM",threshold=1,model="M2")).iloc[0]
        for metric, value in [("M2 OR",mm.OR),("M2 profile P",mm.P),("M2 Holm",mm.holm_p_m2),
                              ("M2 genes",mm.n_genes),("M2 ASD complete cases",mm.asd_n)]:
            delta.append(dict(screen=screen,metric=metric,rounded_reference_value=np.nan,validated_frozen_value=np.nan,
                              R2_value=value,status="R2_NEW_LOEUF_INCLUSIVE_MODEL",R2_scope="Same M0–M3 complete cases"))
    direction = direction_table(data_dir, genes)
    write(direction, out / "results/CANDIDATE_SCREEN_DIRECTIONS_RECALCULATION.tsv")
    for rows,folder,file in [(flow,"results","MEASURED_TO_DETECTABLE_FLOW.tsv"),
        (recalc,"results","R1_FISHER_INDEPENDENT_RECALCULATION.tsv"),
        (breakdown,"results","REGULATOR_CLASS_BREAKDOWN.tsv"),
        (balances,"results","PROPERTY_BALANCE_BEFORE_AFTER.tsv"),
        (coverage_qc,"checks","GENE_UNIVERSE_QC.tsv"),
        (composition,"checks","REFERENCE_MAPPING_CLASS_COMPOSITION.tsv"),
        (positivity,"results","MATCH_CASE_SUPPORT.tsv"),
        (support,"results","MATCHING_FEASIBILITY.tsv"),
        (delta,"checks","R1_TO_R2_NUMERIC_DELTA.tsv")]:
        write(rows,out/folder/file)
    print(pd.DataFrame(support)[["screen","matched_asd","unmatched_ASD","no_caliper_support","max_abs_SMD",
                                  "ASD_hit_control_nonhit","ASD_nonhit_control_hit","status"]].to_string(index=False))
    print(f"Recomputed {len(recalc)} Fisher tests and {len(exact)} exact convolution tails.")

