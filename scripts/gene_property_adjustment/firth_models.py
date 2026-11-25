"""Frozen R7 bias-reduced binomial fits and penalized likelihood profiles.

This release adapts the R2 input and output paths. The estimator, support rules,
optimizer, model grid, and historical bridge formulas are unchanged.
"""
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq, minimize
from scipy.special import expit, logit
from scipy.stats import chi2, norm

from release_inputs import released_file

CORE = ("CHEN", "UMHOEFER_TREG")
REF = {"CHEN": "CHEN_ITREG_NTC", "UMHOEFER_TREG": "CHEN_HTREG", "UMHOEFER_TCONV": "CHEN_HTCONV"}
DESIGNS = [(0.5, "CPM"), (1, "CPM"), (2, "CPM"), (1, "nTPM")]
DETECTED = ("DETECTABLE_HIT", "DETECTABLE_NONHIT")
ALL_FIELDS = ["TF", "EpiFactors", "expression", "LOEUF", "log10_span", "guide_count", "targeted_library"]
MODEL_VARIABLES = {"M0": [], "M1": ["TF", "EpiFactors"],
                   "M2": ["TF", "EpiFactors", "expression", "LOEUF"],
                   "M3": ["TF", "EpiFactors", "expression", "LOEUF", "log10_span", "guide_count"],
                   "M2_NO_LIBRARY": ["TF", "EpiFactors", "expression", "LOEUF"]}
CONTINUOUS = ["expression", "LOEUF", "log10_span", "guide_count"]
SCORE_TOL = 1e-5


def rd(path, **kwargs):
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def write(df, path):
    pd.DataFrame(df).to_csv(path, sep="\t", index=False, float_format="%.13g")


def truth(s):
    return s.astype(str).str.lower().isin(["true", "1", "1.0", "yes"])


def read_inputs(data_dir):
    g = rd(released_file(data_dir, "inputs/R7_GENE_BACKGROUND.tsv.gz"))
    g["asd"] = truth(g["sfari_score1"]).astype(int)
    g["hit"] = truth(g["author_hit"]).astype(int)
    g["is_measured"] = truth(g["measured"])
    g["TF"] = g["TF_LISTED"]
    g["EpiFactors"] = g["EPIFACTORS_GENES_LISTED"]
    g["targeted_library"] = truth(g["library_member"]).astype(int)
    g["LOEUF"] = g["LOEUF_v2_1_1"]
    g["log10_span"] = np.log10(g["GENE_SPAN_BP"])
    g["guide_count"] = g["num"]
    if g[["screen", "approved_symbol"]].isna().any(axis=None) or g.duplicated(["screen", "approved_symbol"]).any():
        raise ValueError("Released gene records require unique screen and approved-symbol keys")
    sets = g.loc[g["asd"].eq(1)].groupby("screen")["approved_symbol"].agg(frozenset)
    if set(sets.index) != set(REF) or len(set(sets)) != 1 or len(sets.iloc[0]) != 245:
        raise ValueError("All three screen universes must preserve the same 245 frozen score-1 symbols")
    master = rd(released_file(data_dir, "inputs/R6_GENE_COVARIATES.tsv"))
    master["hit"] = truth(master["hit"]).astype(int)
    pairs = g[g["is_measured"] & g["screen"].isin(CORE)].merge(
        master, left_on=["screen", "approved_symbol"], right_on=["screen", "gene_symbol"],
        validate="one_to_one", suffixes=("", "_old"))
    assert len(pairs) == 18966 + 1349
    assert np.allclose(pairs["LOEUF"], pairs["LOEUF_v2_1_1_old"], equal_nan=True)
    assert pairs["hit"].eq(pairs["hit_old"]).all() and pairs["asd"].eq(pairs["ASD_score1"]).all()
    return g, master


def reference_data(g, screen, threshold, units):
    measured = g[g["screen"].eq(screen) & g["is_measured"]].copy()
    ref = REF[screen] if units == "CPM" else "HPA_IMMUNE"
    field = "primary_mean_log_expression" if units == "CPM" else "hpa_mean_log_expression"
    label = f"primary_status_{threshold:g}" if units == "CPM" else "hpa_status_1"
    measured["expression"] = measured[field]
    detectable = measured[measured[label].isin(DETECTED)].copy()
    eligible = detectable[np.isfinite(detectable[ALL_FIELDS]).all(axis=1) &
                          detectable["hgnc_id"].notna() & detectable["class_label"].notna()]
    return ref, measured, detectable, eligible.copy(), label


def margins(d):
    a = int((d.asd * d.hit).sum())
    b = int(d.asd.sum()) - a
    c = int(d.hit.sum()) - a
    return a, b, c, len(d) - a - b - c


def p_adjust_two(values):
    a = np.asarray(values, dtype=float)
    order = np.argsort(np.where(np.isfinite(a), a, 1.0))
    adjusted = np.zeros(len(a), dtype=float)
    sorted_p = np.where(np.isfinite(a[order]), a[order], 1.0)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(sorted_p * np.arange(len(a), 0, -1)))
    return np.where(np.isfinite(a), adjusted, np.nan)


def design(d, predictors):
    cols = {}
    for name in predictors:
        s = d[name].astype(float)
        if name in CONTINUOUS:
            sd = s.std(ddof=0)
            s = (s - s.mean()) / sd if sd > 0 else s - s.mean()
        cols[name] = s.to_numpy()
    cols = {"intercept": np.ones(len(d)), "asd": d.asd.to_numpy(), **cols}
    return np.column_stack(list(cols.values())), list(cols)


class PenalizedLogit:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def evaluate(self, beta):
        x, y = self.x, self.y
        eta = x @ beta
        prob = expit(eta)
        weights = expit(eta) * expit(-eta)
        information = x.T @ (weights[:, None] * x)
        try:
            factor = cho_factor(information, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            return -np.inf, np.full(x.shape[1], np.nan), information
        inv_info = cho_solve(factor, np.eye(x.shape[1]), check_finite=False)
        h = weights * np.einsum("ij,ij->i", x @ inv_info, x)
        loglik = float(np.sum(y * eta - np.logaddexp(0, eta)) + np.log(np.diag(factor[0])).sum())
        score = x.T @ (y - prob + h * (.5 - prob))
        return loglik, score, information

    def solve(self, init, fixed=None):
        beta = np.asarray(init).copy()
        free = np.arange(len(beta))
        if fixed is not None:
            beta[1] = fixed
            free = free[free != 1]

        def neg(values):
            beta[free] = values
            ll, sco, _ = self.evaluate(beta)
            return -ll, -sco[free]

        original_free = beta[free].copy()
        result = minimize(neg, original_free, jac=True, method="BFGS",
                          options={"gtol": 1e-7, "maxiter": 300})
        beta[free] = result.x
        ll, score, info = self.evaluate(beta)
        residual = np.linalg.norm(score[free], ord=np.inf)
        messages = [str(result.message)]
        if residual > SCORE_TOL or not np.isfinite(ll):
            before_alternate = beta.copy()
            alternate = minimize(neg, original_free, jac=True, method="L-BFGS-B",
                                 options={"gtol": 1e-9, "ftol": 1e-15, "maxiter": 1200, "maxls": 80})
            alternate_beta = beta.copy()
            alternate_beta[free] = alternate.x
            value, gradient, matrix = self.evaluate(alternate_beta)
            candidate_residual = np.linalg.norm(gradient[free], ord=np.inf)
            if np.isfinite(value) and (value > ll or not np.isfinite(ll) or
                                       (value >= ll - 1e-10 and candidate_residual < residual)):
                beta, ll, score, info, residual = alternate_beta, value, gradient, matrix, candidate_residual
            else:
                beta[:] = before_alternate
            messages.append("L-BFGS-B recovery: " + str(alternate.message))
        if residual > SCORE_TOL and np.isfinite(ll):
            # A fixed Fisher-scoring refinement addresses BFGS precision warnings.
            for _ in range(100):
                step = np.linalg.solve(info[np.ix_(free, free)], score[free])
                for i in range(20):
                    proposal = beta.copy()
                    proposal[free] += (0.5 ** i) * step
                    q, qs, qi = self.evaluate(proposal)
                    if np.isfinite(q) and q >= ll - 1e-11:
                        beta, ll, score, info = proposal, q, qs, qi
                        break
                residual = np.linalg.norm(score[free], ord=np.inf)
                if residual <= SCORE_TOL:
                    break
        return beta, ll, residual, info, " | ".join(messages)

    def profile(self, estimate, best):
        def loss(b):
            _, ll, residual, _, _ = self.solve(estimate, fixed=b)
            if not np.isfinite(ll) or residual > SCORE_TOL or ll > best + 1e-6:
                raise ArithmeticError(f"Profile verification failed at beta={b:.8g}, score={residual:.4g}")
            return 2 * (best - ll)

        target = chi2.ppf(.95, 1)
        endpoints = []
        middle = estimate[1]
        for sign in (-1, 1):
            steps = (.25, .5, 1, 2, 4, 8, 12, 20)
            candidates = (middle + sign * step for step in steps)
            bracket = next((b for b in candidates if loss(b) >= target), None)
            if bracket is None:
                raise ArithmeticError("Profile failed to bracket a 95% endpoint within 20 log-odds units")
            endpoints.append(brentq(lambda b: loss(b) - target,
                                    *sorted([middle, bracket]), xtol=1e-8))
        return np.exp(endpoints), chi2.sf(max(0, loss(0.0)), 1)


def fit(d, predictors):
    blank = dict(OR=np.nan, CI_low=np.nan, CI_high=np.nan, P=np.nan, wald_p_diagnostic=np.nan,
                 converged=False, interval_status="NOT_ESTIMABLE", warnings="", score_residual=np.nan,
                 second_start_residual=np.nan, second_start_delta=np.nan, design_condition=np.nan,
                 information_condition=np.nan, design_rank=np.nan, max_vif=np.nan)
    x, terms = design(d, predictors)
    n, k = x.shape
    if n <= k or d.asd.nunique() < 2 or d.hit.nunique() < 2:
        return blank | {"warnings": "Insufficient rows, outcome variability or score-1 variability"}, []
    rank = np.linalg.matrix_rank(x)
    cond = np.linalg.cond(x)
    blank.update(design_rank=rank, design_condition=cond)
    if rank < k or cond > 1e10:
        return blank | {"warnings": "Rank-deficient or ill-conditioned design"}, []
    vifs = [1.0] if k == 2 else np.diag(np.linalg.inv(np.corrcoef(x[:, 1:], rowvar=False)))
    engine = PenalizedLogit(x, d.hit.to_numpy())
    init = np.zeros(k)
    init[0] = logit((d.hit.sum() + .5) / (n + 1))
    shifted = init.copy()
    shifted[1] = 1.0
    est, ll, score, info, message = engine.solve(init)
    alt, all_, ascore, _, amessage = engine.solve(shifted)
    delta = float(np.max(np.abs(alt - est)))
    icond = float(np.linalg.cond(info))
    fits_agree = bool(np.isfinite(ll) and np.isfinite(all_) and score <= SCORE_TOL
                      and ascore <= SCORE_TOL and abs(ll - all_) <= 1e-7 and delta <= .002
                      and icond <= 1e10)
    blank.update(score_residual=score, second_start_residual=ascore, second_start_delta=delta,
                 information_condition=icond, max_vif=max(vifs), converged=fits_agree,
                 optimizer_message=message, second_start_optimizer_message=amessage)
    if not fits_agree:
        return blank | {"warnings": "Adjusted-score, conditioning, or second-start agreement check failed"}, []
    if k == 2:
        a, b, c, e = margins(d)
        expected = ((a + .5) * (e + .5)) / ((b + .5) * (c + .5))
        assert np.isclose(np.exp(est[1]), expected, rtol=1e-5)
    se = np.sqrt(np.diag(np.linalg.inv(info)))
    blank.update(OR=np.exp(est[1]), wald_p_diagnostic=2 * norm.sf(abs(est[1] / se[1])))
    try:
        ci, p = engine.profile(est, ll)
        if np.isfinite(ci).all() and ci[0] < np.exp(est[1]) < ci[1] and np.isfinite(p):
            blank.update(CI_low=ci[0], CI_high=ci[1], P=p, interval_status="ESTIMABLE")
        else:
            raise ArithmeticError("Nonfinite or invalid profile interval")
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as err:
        blank.update(interval_status="UNSTABLE_PROFILE", warnings=str(err))
    coef = [dict(term=t, beta=b, expected_information_SE=s,
                 vif=np.nan if i == 0 else vifs[i-1]) for i, (t, b, s) in enumerate(zip(terms, est, se))]
    return blank, coef


def count_columns(measured, detected, cases):
    a, b, c, d = margins(cases)
    return dict(n_genes=len(cases), asd_n=a+b, asd_hits=a, other_hits=c, other_n=c+d,
                asd_nonhits=b, other_nonhits=d, total_hits=a+c,
                measured_n=len(measured), measured_asd_n=int(measured.asd.sum()),
                detectable_n=len(detected), detectable_asd_n=int(detected.asd.sum()),
                excluded_covariate_n=len(detected)-len(cases),
                excluded_covariate_asd_n=int(detected.asd.sum())-a-b,
                loeuf_missing_detectable=int(detected.LOEUF.isna().sum()),
                loeuf_missing_detectable_asd=int(detected.loc[detected.LOEUF.isna(), "asd"].sum()),
                loeuf_missing_measured=int(measured.LOEUF.isna().sum()),
                loeuf_missing_measured_asd=int(measured.loc[measured.LOEUF.isna(), "asd"].sum()),
                complete_over_measured_asd=(a+b) / measured.asd.sum(),
                complete_over_detectable_asd=(a+b) / detected.asd.sum())


def run(data_dir, out, genes, historical):
    models, coefficients, exclusions, cohort_rows = [], [], [], []
    for screen in CORE:
        for threshold, units in DESIGNS:
            reference, measured, detected, cases, _ = reference_data(genes, screen, threshold, units)
            spec = {"screen": screen, "reference": reference, "threshold": threshold, "units": units}
            missing = detected[~detected.approved_symbol.isin(cases.approved_symbol)].copy()
            for row in missing.itertuples():
                reasons = [col for col in ALL_FIELDS if pd.isna(getattr(row, col))]
                if pd.isna(row.class_label):
                    reasons.append("class_label")
                if pd.isna(row.hgnc_id):
                    reasons.append("frozen_HGNC")
                exclusions.append(spec | dict(gene=row.approved_symbol, asd=row.asd,
                                               hit=row.hit, missing="|".join(reasons)))
            for col in ALL_FIELDS:
                v = cases[col]
                cohort_rows.append(spec | dict(covariate=col, common_mean=v.mean(),
                                                common_sd=v.std(ddof=0), min=v.min(), max=v.max()))
            for method, names in MODEL_VARIABLES.items():
                if method == "M2_NO_LIBRARY" and (screen, units, threshold) != ("CHEN", "CPM", 1):
                    continue
                predictors = list(names)
                if screen == "CHEN" and method in ("M1", "M2", "M3"):
                    predictors.insert(2, "targeted_library")
                model_fit, coef = fit(cases, predictors)
                model = spec | {"model": method, "formula": "hit ~ " + " + ".join(["asd", *predictors])}
                models.append(model | count_columns(measured, detected, cases) | model_fit)
                coefficients.extend(model | c for c in coef)
    model_results = pd.DataFrame(models)
    primary = model_results.model.eq("M2") & model_results.threshold.eq(1) & model_results.units.eq("CPM")
    assert primary.sum() == 2
    model_results["holm_p_m2"] = np.nan
    model_results.loc[primary, "holm_p_m2"] = p_adjust_two(model_results.loc[primary, "P"])
    model_results["estimator"] = "Firth/Jeffreys mean bias reduction; two-sided penalized profile LR"
    write(model_results, out / "results/LOEUF_ADJUSTED_MODELS.tsv")
    write(coefficients, out / "results/LOEUF_MODEL_COEFFICIENTS.tsv")
    write(exclusions, out / "results/LOEUF_COMPLETE_CASE_EXCLUSIONS.tsv")
    write(cohort_rows, out / "results/LOEUF_SCALING.tsv")
    rows = []
    old_result = rd(released_file(data_dir, "results/R6_HISTORICAL_FIRTH_ASD_SUMMARY.tsv"))
    for screen in CORE:
        _, _, _, cases, _ = reference_data(genes, screen, 1, "CPM")
        o = old_result[old_result.screen.eq(screen) & old_result.model.eq("M1")].iloc[0]
        rows.append(dict(screen=screen, analysis="R6_M1_PUBLISHED_HISTORICAL", n_genes=o.complete_case_N,
                         asd_n=o.ASD_score1_complete_case_N, OR=o.adjusted_OR, CI_low=o.profile_likelihood_CI_low,
                         CI_high=o.profile_likelihood_CI_high, wald_p_diagnostic=o.p_value,
                         converged=o.convergence,
                         interval_status="REPORTED_HISTORICAL" if o.profile_CI_success else "UNSTABLE_IN_R6"))
        oz = historical[historical.screen.eq(screen)].rename(columns={"ASD_score1":"asd",
            "TF_member":"TF", "EpiFactors_member":"EpiFactors", "LOEUF_z":"LOEUF",
            "expression_z":"expression"})
        for col in ("TF", "EpiFactors"):
            oz[col] = truth(oz[col]).astype(int)
        preds = ["TF", "EpiFactors", "expression", "LOEUF"]
        old = oz.dropna(subset=preds).copy()
        assert len(old) == o.complete_case_N and old.asd.sum() == o.ASD_score1_complete_case_N
        shared = cases.merge(oz[["gene_symbol", "expression"]].rename(columns={"expression":"old_expression"}),
                             left_on="approved_symbol", right_on="gene_symbol", validate="one_to_one")
        shared = shared.dropna(subset=["old_expression"])
        for name, d in (("R6_DESIGN_REFIT_OLD_COHORT", old),
                        ("CPM1_SHARED_COHORT_OLD_EXPRESSION", shared.assign(expression=shared.old_expression)),
                        ("CPM1_SHARED_COHORT_NEW_EXPRESSION", shared)):
            z, _ = fit(d, preds)
            a, b, c, e = margins(d)
            if name == "R6_DESIGN_REFIT_OLD_COHORT":
                assert np.isclose(z["OR"], o.adjusted_OR, rtol=1e-3), (screen, z["OR"], o.adjusted_OR)
            rows.append(dict(screen=screen, analysis=name, n_genes=len(d), asd_n=a+b,
                             asd_hits=a, other_hits=c) | z)
    write(rows, out / "results/R6_TO_R2_MODEL_BRIDGE.tsv")
    print(model_results.loc[model_results.threshold.eq(1) & model_results.units.eq("CPM"),
          ["screen","model","n_genes","asd_n","asd_hits","other_hits","OR","CI_low","CI_high","P","holm_p_m2","interval_status"]].to_string(index=False))
    print(f"Finished {len(model_results)} locked model rows and {len(rows)} historical bridge rows.")

