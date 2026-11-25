#!/usr/bin/env python3
"""Recompute R9 exact tests with integer combinatorics from the v1.1.0 release.

Uses only the Python standard library. Each probability starts as an exact
subset count divided by the exact product of choose(N, A). No production
hypergeometric, ranking, convolution, or Holm code is imported.
"""
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import gzip
from itertools import combinations, groupby
import math
from pathlib import Path
from statistics import mean, variance

SCREENS = ("CHEN", "UMHOEFER_TREG")
REFERENCES = {"CHEN": "CHEN_ITREG_NTC", "UMHOEFER_TREG": "CHEN_HTREG"}
CLASSES = {"NEITHER_LISTED", "TF_ONLY", "EPIFACTOR_ONLY", "BOTH_LISTED"}
DIRECTION = {
    "SUPPORTS_FOXP3": "KO_REDUCES_FOXP3",
    "RESTRAINS_FOXP3": "KO_INCREASES_FOXP3",
}
A_DIRECTIONS = tuple(DIRECTION.values())
EXPRESSION = "HEALTHY_CD4_MEAN_LOGCPM"
LOEUF = "LOEUF_v2_1_1"
B_LABELS = (
    "FULL_MEASURED_UNADJUSTED",
    "COMPLETE_CASE_UNADJUSTED",
    "COMPLETE_CASE_CONDITIONED",
)


def load(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def index(rows, keys):
    indexed = {tuple(row[key] for key in keys): row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"Duplicate key in fields {keys}")
    return indexed


def flag(s):
    flags = {"true": True, "false": False, "1": True, "0": False}
    if s.lower() not in flags:
        raise ValueError(f"Unrecognized flag: {s!r}")
    return flags[s.lower()]


def number(s):
    return float(s) if s not in ("", "NA", "NaN", "nan") else math.nan


def check_close(known, observed, descriptor, rtol=5e-9, atol=2e-12):
    if not math.isclose(float(known), number(str(observed)), rel_tol=rtol, abs_tol=atol):
        raise ValueError(f"{descriptor}: verified {known}, comparison {observed}")


def save(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def choose_list(n, max_k):
    return [math.comb(n, k) for k in range(min(n, max_k) + 1)]


def score_weights(plus, zero, minus, cases):
    """Sum distinct favorable subsets of size A for each signed score."""
    n = plus + zero + minus
    if not 0 <= cases <= n:
        raise ValueError((plus, zero, minus, cases))
    positive = choose_list(plus, cases)
    neutral = choose_list(zero, cases)
    negative = choose_list(minus, cases)
    weight = Counter()
    for ap in range(max(0, cases - zero - minus), len(positive)):
        m_start = max(0, cases - ap - zero)
        m_end = min(minus, cases - ap)
        for am in range(m_start, m_end + 1):
            az = cases - ap - am
            weight[ap - am] += positive[ap] * neutral[az] * negative[am]
    base = math.comb(n, cases)
    if sum(weight.values()) != base:
        raise ArithmeticError((plus, zero, minus, cases, "exact mass != 1"))
    return weight, base


def combine(local_pmfs):
    pmf, denominator = Counter({0: 1}), 1
    for local, base in local_pmfs:
        next_pmf = Counter()
        for k, count in local.items():
            for old_k, old_count in pmf.items():
                next_pmf[k + old_k] += count * old_count
        pmf, denominator = next_pmf, denominator * base
        if sum(pmf.values()) != denominator:
            raise ArithmeticError("exact combined mass != 1")
    return pmf, denominator


def exact_stats(weights, base, observed):
    expectation = Fraction(sum(score * count for score, count in weights.items()), base)
    second = Fraction(sum(score * score * count for score, count in weights.items()), base)
    cumulative = 0
    lower, upper = None, None
    for score, count in sorted(weights.items()):
        cumulative += count
        if lower is None and 40 * cumulative >= base:
            lower = score
        if 40 * cumulative >= 39 * base:
            upper = score
            break
    return {
        "expected": expectation,
        "variance": second - expectation * expectation,
        "lower": lower,
        "upper": upper,
        "p": Fraction(sum(count for score, count in weights.items() if score >= observed), base),
    }


def analytic_moments(plus, zero, minus, cases):
    n = plus + zero + minus
    population_mean = Fraction(plus - minus, n)
    population_var = Fraction(plus + minus, n) - population_mean * population_mean
    correction = Fraction(cases * (n - cases), n - 1) if n > 1 else 0
    return cases * population_mean, correction * population_var


def toy_check():
    small = [
        ((1, 1, 0, 0, -1, -1, -1), 3),
        ((1, -1, 1), 2),
        ((1, -1, 0), 3),
        ((0, -1, 0), 0),
    ]
    exhaustive = []
    exact = []
    for values, cases in small:
        observed = Counter(sum(z) for z in combinations(values, cases))
        result = score_weights(values.count(1), values.count(0), values.count(-1), cases)
        assert result == (observed, sum(observed.values()))
        exhaustive.append(observed)
        exact.append(result)
    straight = Counter()
    for i, a in exhaustive[0].items():
        for j, b in exhaustive[1].items():
            for k, c in exhaustive[2].items():
                for l, d in exhaustive[3].items():
                    straight[i + j + k + l] += a * b * c * d
    combined, total = combine(exact)
    assert combined == straight
    return len(small), total


def avg_rank_tertiles(rows, value_column):
    # For tied ranks r...s, use the exact integer identity
    # ceil(3*(r+s)/(2*n)) to avoid floating point at cut boundaries.
    result = {}
    first = 1
    n = len(rows)
    for _, tied in groupby(sorted(rows, key=lambda row: number(row[value_column])),
                           key=lambda row: number(row[value_column])):
        block = list(tied)
        last = first + len(block) - 1
        tercile = min(3, (3 * (first + last) + 2 * n - 1) // (2 * n))
        result.update({row["approved_symbol"]: tercile for row in block})
        first = last + 1
    assert first == n + 1
    return result


def validate_background(bg):
    keys = index(bg, ("screen", "approved_symbol"))
    assert len(keys) == 22091
    ids_by_symbol = {}
    symbols_by_id = {}
    all_sets = defaultdict(set)
    baseline = defaultdict(Counter)
    for (screen, symbol), row in keys.items():
        assert symbol
        case = flag(row["sfari_score1"])
        if case:
            all_sets[screen].add(symbol)
        if not flag(row["measured"]):
            assert row["author_hit"] == "" and row["primary_status_1"] == "NOT_MEASURED"
            continue
        hit = flag(row["author_hit"])
        id_ = row["hgnc_id"]
        assert id_.startswith("HGNC:")
        assert ids_by_symbol.setdefault(symbol, id_) == id_
        assert symbols_by_id.setdefault(id_, symbol) == symbol
        fdr_hit = min(number(row["pos|fdr"]), number(row["neg|fdr"])) <= .05
        assert hit == fdr_hit, (screen, symbol)
        baseline[screen][case, hit] += 1
    assert len(set(map(frozenset, all_sets.values()))) == 1
    score1, = set(map(frozenset, all_sets.values()))
    assert len(score1) == 245
    for screen, expected in zip(SCREENS, ((14, 228, 243, 18481), (5, 28, 26, 1290))):
        n = baseline[screen]
        actual = tuple(n[x] for x in [(True, True), (True, False), (False, True), (False, False)])
        assert actual == expected, (screen, actual, expected)
    return keys, score1


def category_counts(rows):
    total = Counter(r["D"] for r in rows)
    cases = Counter(r["D"] for r in rows if r["case"])
    return {
        "N": len(rows), "A": sum(cases.values()), "n_plus": total[1],
        "n_zero": total[0], "n_minus": total[-1], "a_plus": cases[1],
        "a_zero": cases[0], "a_minus": cases[-1], "observed_S": cases[1] - cases[-1],
    }


def phenotype_margins(data, keys, checks_dir):
    common = index(load(data / "data/processed/chen_shifrut_common_universe.tsv"), ("approved_symbol",))
    assert len(common) == 18615
    assert len({symbol.upper() for (symbol,) in common}) == len(common)
    measured = [row for row in keys.values() if row["screen"] == "CHEN" and flag(row["measured"])]
    normalized = {row["approved_symbol"].upper(): row for row in measured}
    assert len(normalized) == len(measured), "Nonunique uppercase measured Chen approved symbol"
    with_covariates = []
    casing_only = []
    for (symbol,), frozen in common.items():
        if symbol.upper() not in normalized:
            raise ValueError(f"Frozen phenotype symbol has no case-normalized Chen match: {symbol}")
        chen = normalized[symbol.upper()]
        if symbol != chen["approved_symbol"]:
            casing_only.append((symbol, chen["approved_symbol"], chen["hgnc_id"]))
        case = flag(frozen["SFARI_SCORE1"])
        fo = flag(frozen["CHEN_HIT_BOOL"])
        tcr = flag(frozen["SHIFRUT_HIT"])
        D = int(fo) - int(tcr)
        assert D == int(frozen["D"])
        assert case == flag(chen["sfari_score1"])
        assert fo == flag(chen["author_hit"])
        assert frozen["CHEN_DIRECTION"] == chen["FOXP3_REGULATOR_DIRECTION"]
        available = (
            chen["class_label"] in CLASSES and chen["library_member"].lower() in ("true", "false")
            and math.isfinite(number(chen[EXPRESSION])) and math.isfinite(number(chen[LOEUF]))
        )
        with_covariates.append(chen | {"case": case, "D": D, "complete": available,
                                       "frozen_common_symbol": symbol})
    assert len(casing_only) == 192
    save(checks_dir / "INDEPENDENT_PHENOTYPE_CASE_ONLY_JOIN.tsv",
         ("frozen_common_symbol", "frozen_chen_approved_symbol", "frozen_chen_hgnc_id"), casing_only)
    full_counts = category_counts(with_covariates)
    assert full_counts == dict(N=18615, A=242, n_plus=179, n_zero=17562, n_minus=874,
                               a_plus=11, a_zero=223, a_minus=8, observed_S=3)
    complete = [r for r in with_covariates if r["complete"]]
    strata = defaultdict(list)
    expression_bins = avg_rank_tertiles(complete, EXPRESSION)
    loeuf_bins = avg_rank_tertiles(complete, LOEUF)
    for row in complete:
        symbol = row["approved_symbol"]
        inside = "INSIDE" if flag(row["library_member"]) else "OUTSIDE"
        stratum = "|".join((row["class_label"], inside, str(expression_bins[symbol]), str(loeuf_bins[symbol])))
        row["stratum_key"] = stratum
        strata[stratum].append(row)
    return with_covariates, complete, dict(sorted(strata.items()))


def phenotype_balance(strata, complete):
    n_supported, n_informative = 0, 0
    case_weight, control_weight = {}, {}
    for genes in strata.values():
        counts = category_counts(genes)
        a, b = counts["A"], counts["N"] - counts["A"]
        if a and b:
            n_supported += a
            if len({r["D"] for r in genes}) > 1:
                n_informative += a
            for r in genes:
                symbol = r["approved_symbol"]
                if r["case"]:
                    case_weight[symbol] = 1.0
                else:
                    control_weight[symbol] = a / b
    rows = []
    for variable in (EXPRESSION, LOEUF):
        pre_a = [number(r[variable]) for r in complete if r["case"]]
        pre_b = [number(r[variable]) for r in complete if not r["case"]]
        supported_a = [number(r[variable]) for r in complete if r["approved_symbol"] in case_weight]
        supported_b = [number(r[variable]) for r in complete if r["approved_symbol"] in control_weight]
        denominator = math.sqrt((variance(supported_a) + variance(supported_b)) / 2)
        weighted_b = sum(number(r[variable]) * control_weight.get(r["approved_symbol"], 0)
                         for r in complete) / sum(control_weight.values())
        pre_smd = (mean(pre_a) - mean(pre_b)) / math.sqrt((variance(pre_a) + variance(pre_b)) / 2)
        after_smd = (mean(supported_a) - weighted_b) / denominator
        rows.append((variable, pre_smd, after_smd, len(supported_a), len(supported_b)))
    return n_supported, n_informative, rows


def compare_pmf(weights, denominator, rows, score_field, descriptor):
    supplied = {int(r[score_field]): number(r["probability"]) for r in rows}
    assert len(rows) == len(supplied), descriptor
    score_set = set(weights).union(supplied)
    maximum = 0.0
    for score in score_set:
        exact = weights.get(score, 0) / denominator
        actual = supplied.get(score, 0)
        check_close(exact, actual, f"{descriptor} score {score}", rtol=2e-8, atol=3e-14)
        if exact and score not in supplied:
            raise ValueError(f"{descriptor} missing representable score {score}")
        maximum = max(maximum, abs(actual - exact))
    check_close(1, math.fsum(supplied.values()), f"{descriptor} released mass")
    return maximum


def run_b(data, keys, production, checks_dir):
    full, cc, strata = phenotype_margins(data, keys, checks_dir)
    supported, informative, balance = phenotype_balance(strata, cc)
    runs = (
        (B_LABELS[0], {"ALL": full}, full),
        (B_LABELS[1], {"ALL": cc}, cc),
        (B_LABELS[2], strata, cc),
    )
    prod_stats, prod_strata, prod_null = {}, {}, []
    if production:
        prod_stats = index(load(production / "PHENOTYPE_PREFERENCE_R9.tsv"), ("analysis",))
        prod_strata = index(load(production / "PHENOTYPE_PREFERENCE_STRATA.tsv"), ("analysis", "stratum_key"))
        prod_null = load(production / "PHENOTYPE_PREFERENCE_NULL.tsv")
        members = index(load(production / "PHENOTYPE_PREFERENCE_MEMBERS.tsv.gz"), ("approved_symbol",))
        balance_prod = index([r for r in load(production / "PHENOTYPE_PREFERENCE_COVERAGE_BALANCE.tsv")
                              if r["section"] == "BALANCE"], ("field", "subset"))
        for variable, raw, weighted, ac, bc in balance:
            field = "healthy_cd4_mean_logcpm" if variable == EXPRESSION else variable
            check_close(raw, balance_prod[field, "COMPLETE_CASE_UNWEIGHTED"]["signed_SMD"], field)
            answer = balance_prod[field, "WITHIN_STRATUM_WEIGHTED"]
            check_close(weighted, answer["signed_SMD"], field)
            assert (int(answer["A"]), int(answer["other_n"])) == (ac, bc)
        assert set(members) == {(r["approved_symbol"],) for r in full}
        for row in full:
            a = members[row["approved_symbol"],]
            assert bool(int(a["score1"])) == row["case"] and int(a["D"]) == row["D"]
            assert flag(a["complete_case"]) == row["complete"], row["approved_symbol"]
            assert a["stratum_key"] == row.get("stratum_key", ""), row["approved_symbol"]
            assert a["legacy_common_symbol"] == row["frozen_common_symbol"], row["approved_symbol"]
            assert flag(a["case_only_symbol_match"]) == (row["approved_symbol"] != row["frozen_common_symbol"])
    stats_rows = []
    strata_rows = []
    pmf_rows = []
    for label, groups, records in runs:
        comparison_strata = {(a, k) for (a, k) in prod_strata if a == label}
        if production:
            assert comparison_strata == {(label, k) for k in groups}, label
        overall = category_counts(records)
        local_pmfs = []
        ex, var = Fraction(), Fraction()
        ex_plus, ex_minus = Fraction(), Fraction()
        for name, group in groups.items():
            counts = category_counts(group)
            margs = tuple(counts[k] for k in ("n_plus", "n_zero", "n_minus", "A"))
            local_pmfs.append(score_weights(*margs))
            ex_i, var_i = analytic_moments(*margs)
            ex += ex_i
            var += var_i
            ex_plus += Fraction(counts["n_plus"] * counts["A"], counts["N"])
            ex_minus += Fraction(counts["n_minus"] * counts["A"], counts["N"])
            strata_rows.append((label, name, *counts.values(), float(ex_i), float(var_i)))
            if production:
                details = prod_strata[label, name]
                for col, val in counts.items():
                    assert int(details[col]) == val, (label, name, col)
        weights, denom = combine(local_pmfs)
        stat = exact_stats(weights, denom, overall["observed_S"])
        assert stat["expected"] == ex and stat["variance"] == var, label
        p = float(stat["p"])
        if label == B_LABELS[0]:
            historic = load(data / "data/processed/direct_phenotype_preference_summary.tsv")[0]
            check_close(p, historic["exact_null_p"], "frozen exact full-universe tail", atol=1e-15)
        maximum = ""
        if production:
            actual = prod_stats[label,]
            for col, val in overall.items():
                assert int(actual[col]) == val, (label, col)
            for independent, production_field in (("expected", "expected_S"), ("variance", "var_S"),
                                                   ("lower", "null_low95"), ("upper", "null_high95"),
                                                   ("p", "raw_P")):
                check_close(float(stat[independent]), actual[production_field], (label, production_field))
            for col, val in (("expected_ASD_plus", ex_plus), ("expected_ASD_minus", ex_minus),
                             ("observed_minus_expected_S", overall["observed_S"] - ex)):
                check_close(float(val), actual[col], (label, col))
            if label == B_LABELS[2]:
                assert int(actual["comparator_supported_ASD"]) == supported
                assert int(actual["informative_ASD"]) == informative
            maximum = compare_pmf(weights, denom, [r for r in prod_null if r["analysis"] == label], "S", label)
        method = "integer subset counts; exact mass/mean/variance"
        parity = "PASS" if production else "FROZEN_INDEPENDENT"
        stats_rows.append((label, len(groups), *overall.values(), float(ex), float(var), stat["lower"],
                           stat["upper"], p, maximum, parity, method))
        for score in sorted(weights):
            pmf_rows.append((label, score, weights[score] / denom))
    count_header = tuple(category_counts(full))
    save(checks_dir / "INDEPENDENT_PHENOTYPE_EXACT.tsv",
         ("analysis", "strata", *count_header, "expected_S", "variance", "null_low95", "null_high95",
          "exact_P", "largest_probability_difference", "status", "method"), stats_rows)
    save(checks_dir / "INDEPENDENT_PHENOTYPE_STRATA.tsv",
         ("analysis", "stratum_key", *count_header, "expected_S", "variance"), strata_rows)
    save(checks_dir / "INDEPENDENT_PHENOTYPE_NULL.tsv", ("analysis", "S", "probability"), pmf_rows)
    save(checks_dir / "INDEPENDENT_PHENOTYPE_BALANCE.tsv",
         ("covariate", "cc_unweighted_smd", "after_case_weighted_smd", "supported_cases", "supported_controls"),
         balance)
    for row in stats_rows:
        print(row[0], "N", row[2], "A", row[3], "observed", row[10], "expected", row[11],
              "variance", row[12], "p", row[15])
    print("Case support and informative cases:", supported, informative)


def author_direction(row):
    if flag(row["author_hit"]):
        return DIRECTION.get(row["FOXP3_REGULATOR_DIRECTION"], "AMBIGUOUS_HIT")
    return "MEASURED_NONHIT"


def validate_directions(data, keys):
    candidate = index(load(data / "data/r7/mechanism/CANDIDATE_SCREEN_DIRECTIONS.tsv"), ("screen", "gene"))
    audit = index(load(data / "data/r7/results/DIRECTION_ORIGINAL_TO_MANUSCRIPT.tsv"),
                  ("screen", "approved_symbol"))
    assert set(candidate) == set(audit) and len(candidate) == 51
    for key, source in candidate.items():
        row = keys[key]
        is_measured = flag(row["measured"])
        assert is_measured == flag(source["screened"]) == flag(audit[key]["is_measured"]), key
        if is_measured:
            label = row["FOXP3_REGULATOR_DIRECTION"]
            assert label == source["frozen_regulator_label"] == audit[key]["FOXP3_REGULATOR_DIRECTION"], key
            if flag(row["author_hit"]):
                direction = author_direction(row)
                assert direction == source["KO_effect_on_FOXP3_protein"]
                assert direction == audit[key]["standardized_protein_sort_direction"]
    return len(candidate)


def all_hit_counts(rows):
    return (len(rows), sum(r["hit"] for r in rows), sum(r["case"] for r in rows),
            sum(r["case"] and r["hit"] for r in rows))


def lof_assign(data, keys, score1):
    master = index(load(data / "data/r7/mechanism/ALLELIC_MECHANISM_MASTER.tsv"), ("approved_symbol",))
    assert set(master) == {(s,) for s in score1}
    eligible = {symbol for (symbol,), row in master.items() if flag(row["lof_eligible"])}
    assert len(eligible) == 104
    assert sum(master[s,]["strongest_qualifying_lof_scope"] == "ASD_named" for s in eligible) == 7
    frozen = index([r for r in load(data / "data/r7/inputs/R1_CONDITIONING_STRATA.tsv")
                    if r["set"] == "LOF_COMPATIBLE" and r["analysis"] == "strict_joint"
                    and r["units"] == "CPM" and number(r["threshold"]) == 1],
                   ("screen", "stratum_key"))
    scope = index([r for r in load(data / "data/r7/mechanism/LOF_GROUP_SCOPE_SENSITIVITIES.tsv")
                   if r["qualifying_disease_scope"] == "ALL_104_LOF"], ("screen",))
    assigned = {}
    historical = {}
    for screen in SCREENS:
        strict = []
        measured_cases = 0
        detected_cases = 0
        for (state, symbol), row in keys.items():
            if state != screen or not flag(row["measured"]) or (symbol in score1 and symbol not in eligible):
                continue
            case = symbol in eligible
            measured_cases += case
            if row["primary_status_1"] not in ("DETECTABLE_HIT", "DETECTABLE_NONHIT"):
                continue
            detected_cases += case
            span = number(row["GENE_SPAN_BP"])
            guides = number(row["num"])
            expr = number(row["primary_mean_log_expression"])
            if row["class_label"] not in CLASSES or not all(map(math.isfinite, (span, guides, expr))):
                continue
            if span <= 0 or guides <= 0:
                continue
            strict.append(row | {"case": case, "hit": flag(row["author_hit"]), "direction": author_direction(row)})
        expr_bins = avg_rank_tertiles(strict, "primary_mean_log_expression")
        span_bins = avg_rank_tertiles(strict, "GENE_SPAN_BP")
        groups = defaultdict(list)
        for row in strict:
            symbol = row["approved_symbol"]
            factors = [row["class_label"]]
            if screen == "CHEN":
                factors.append("INSIDE" if flag(row["library_member"]) else "OUTSIDE")
            factors.extend((str(expr_bins[symbol]), str(span_bins[symbol])))
            guides = number(row["num"])
            factors.append("1_TO_3" if guides <= 3 else "EXACTLY_4" if guides <= 4 else "5_OR_MORE")
            groups["|".join(factors)].append(row)
        assert set(groups) == {s for h, s in frozen if h == screen}, screen
        for stratum, rows in groups.items():
            r = frozen[screen, stratum]
            assert r["reference"] == REFERENCES[screen]
            assert all_hit_counts(rows) == tuple(int(r[c]) for c in ("N", "H", "A", "K")), (screen, stratum)
        n, h, a, k = all_hit_counts(strict)
        old = scope[screen,]
        assert (n - a, h - k, a, k, measured_cases, detected_cases) == tuple(int(old[c]) for c in (
            "strict_other_gene_n", "strict_other_gene_hits", "strict_complete_case_genes",
            "strict_complete_case_hits", "screen_measured_case_genes", "detectable_case_genes"))
        dirs = Counter(r["direction"] for r in strict if r["case"] and r["hit"])
        assert (dirs[A_DIRECTIONS[0]], dirs[A_DIRECTIONS[1]]) == tuple(int(old[c]) for c in (
            "strict_KO_reduces_FOXP3_hits", "strict_KO_increases_FOXP3_hits"))
        weights, denominator = combine([
            score_weights(h_j, n_j - h_j, 0, a_j)
            for n_j, h_j, a_j, k_j in map(all_hit_counts, groups.values())
        ])
        baseline = exact_stats(weights, denominator, k)
        check_close(float(baseline["expected"]), old["R1_expected_hits"], (screen, "frozen all-hit expectation"))
        check_close(float(baseline["p"]), old["R1_exact_p"], (screen, "frozen all-hit P"))
        historical[screen] = baseline
        assigned[screen] = dict(sorted(groups.items()))
    historic_holm = holm({key: row["p"] for key, row in historical.items()})
    for key, adjusted in historic_holm.items():
        check_close(float(adjusted), scope[key,]["R1_two_screen_Holm_p"], (key, "historic Holm P"))
    return assigned, len(frozen)


def holm(items):
    # Sort unestimable hypotheses as 1 within the fixed family, then mask.
    order = sorted(items, key=lambda key: items[key] if items[key] is not None else Fraction(1))
    cur = Fraction()
    results = {}
    for rank, key in enumerate(order):
        value = items[key]
        cur = max(cur, min(Fraction(1), (len(items) - rank) * (Fraction(1) if value is None else value)))
        results[key] = cur if value is not None else None
    return results


def exact_odds_ratio(weights, observed):
    scores = sorted(weights)
    if len(scores) < 2:
        return math.nan, 0.0, math.inf
    log_weights = [math.log(weights[k]) for k in scores]

    def summaries(theta):
        logs = [a + k * theta for k, a in zip(scores, log_weights)]
        m = max(logs)
        relative = [math.exp(x - m) for x in logs]
        mass = math.fsum(relative)
        average = math.fsum(k * v for k, v in zip(scores, relative)) / mass
        up = math.fsum(v for k, v in zip(scores, relative) if k >= observed) / mass
        down = math.fsum(v for k, v in zip(scores, relative) if k <= observed) / mass
        return average, up, down

    def solve(function):
        lo, hi = -55.0, 55.0
        assert function(lo) <= 0 <= function(hi)
        for _ in range(80):
            mid = (lo + hi) / 2
            if function(mid) > 0:
                hi = mid
            else:
                lo = mid
        return math.exp((lo + hi) / 2)

    low = observed == min(scores)
    high = observed == max(scores)
    estimate = 0.0 if low else math.inf if high else solve(lambda z: summaries(z)[0] - observed)
    lower = 0.0 if low else solve(lambda z: summaries(z)[1] - .025)
    upper = math.inf if high else solve(lambda z: .025 - summaries(z)[2])
    return estimate, lower, upper


def run_a(data, keys, score1, production, checks_dir):
    n_candidates = validate_directions(data, keys)
    assigned, n_frozen = lof_assign(data, keys, score1)
    prod_stats, prod_strata, prod_null = {}, {}, defaultdict(list)
    if production:
        prod_stats = index(load(production / "LOF_DIRECTIONAL_ENRICHMENT.tsv"),
                           ("analysis", "screen", "direction"))
        prod_strata = index(load(production / "LOF_DIRECTIONAL_STRATA.tsv"), ("screen", "direction", "stratum_key"))
        for row in load(production / "LOF_DIRECTIONAL_NULL.tsv"):
            key = tuple(row[c] for c in ("screen", "direction", "distribution_level", "stratum_key"))
            prod_null[key].append(row)
        member_prod = index(load(production / "LOF_DIRECTIONAL_MEMBERS.tsv.gz"), ("screen", "approved_symbol"))
        strict = {key for key, row in member_prod.items() if flag(row["in_strict_comparison"])}
        assert strict == {(screen, r["approved_symbol"]) for screen, groups in assigned.items()
                          for members in groups.values() for r in members}
        for screen, groups in assigned.items():
            for name, members in groups.items():
                for r in members:
                    row = member_prod[screen, r["approved_symbol"]]
                    assert row["stratum_key"] == name and flag(row["lof_case"]) == r["case"]
                    assert (int(row[A_DIRECTIONS[0] + "_Y"]), int(row[A_DIRECTIONS[1] + "_Y"])) == (
                        r["direction"] == A_DIRECTIONS[0], r["direction"] == A_DIRECTIONS[1])
    results = {}
    for screen, groups in assigned.items():
        for direction in A_DIRECTIONS:
            strata = []
            local_pmfs = []
            for name, rows in groups.items():
                n, all_hit_h, a, all_hit_k = all_hit_counts(rows)
                h = sum(r["direction"] == direction for r in rows)
                k = sum(r["case"] and r["direction"] == direction for r in rows)
                local = score_weights(h, n - h, 0, a)
                local_pmfs.append(local)
                s_local = exact_stats(*local, k)
                assert (s_local["expected"], s_local["variance"]) == analytic_moments(h, n - h, 0, a)
                if production:
                    details = prod_strata[screen, direction, name]
                    for field, value in zip(("N", "H", "A", "K", "all_hit_H", "all_hit_K"),
                                            (n, h, a, k, all_hit_h, all_hit_k)):
                        assert int(details[field]) == value, (screen, direction, name, field)
                    assert h + int(details["other_direction_H"]) + int(details["ambiguous_H"]) == all_hit_h
                    assert k + int(details["other_direction_K"]) + int(details["ambiguous_K"]) == all_hit_k
                    compare_pmf(*local, prod_null[screen, direction, "STRATUM", name], "case_hit_count",
                                (screen, direction, name))
                strata.append((n, h, a, k))
            total_n, total_h, total_a, total_k = map(sum, zip(*strata))
            weights, denominator = combine(local_pmfs)
            stat = exact_stats(weights, denominator, total_k)
            moments = [analytic_moments(h, n - h, 0, a) for n, h, a, k in strata]
            assert stat["expected"] == sum(m[0] for m in moments)
            assert stat["variance"] == sum(m[1] for m in moments)
            fisher_weights, fisher_base = score_weights(total_h, total_n - total_h, 0, total_a)
            fisher = exact_stats(fisher_weights, fisher_base, total_k)
            odds = exact_odds_ratio(fisher_weights, total_k)
            check_max = ""
            if production:
                local_rows = prod_null[screen, direction, "COMBINED", "ALL_STRATA"]
                check_max = compare_pmf(weights, denominator, local_rows, "case_hit_count", (screen, direction))
            results[screen, direction] = dict(
                n=total_n, h=total_h, a=total_a, k=total_k, stat=stat, fisher=fisher, odds=odds, maxdiff=check_max)
    holm_cond = holm({key: x["stat"]["p"] if x["stat"]["variance"] else None for key, x in results.items()})
    holm_fisher = holm({key: x["fisher"]["p"] if x["fisher"]["variance"] else None for key, x in results.items()})
    rows = []
    for (screen, direction), r in results.items():
        s = r["stat"]
        p1 = holm_cond[screen, direction]
        p2 = holm_fisher[screen, direction]
        effect = r["k"] / s["expected"] if s["expected"] else math.nan
        if production:
            for kind, p, adj in (("STRICT_CONDITIONAL_EXACT", s["p"], p1),
                                 ("SAME_POPULATION_FISHER", r["fisher"]["p"], p2)):
                prod = prod_stats[kind, screen, direction]
                for field, val in (("case_n", r["a"]), ("control_n", r["n"] - r["a"]),
                                   ("case_hits", r["k"]), ("control_hits", r["h"] - r["k"])):
                    assert int(prod[field]) == val, (kind, screen, direction, field)
                if adj is None:
                    assert math.isnan(number(prod["p_raw"])) and math.isnan(number(prod["p_holm"]))
                else:
                    check_close(float(p), prod["p_raw"], kind)
                    check_close(float(adj), prod["p_holm"], kind)
                if kind == "SAME_POPULATION_FISHER":
                    for col, val in zip(("effect", "ci_low", "ci_high"), r["odds"]):
                        if math.isnan(val):
                            assert math.isnan(number(prod[col]))
                        else:
                            check_close(val, prod[col], (screen, direction, col), rtol=4e-7)
                else:
                    for col, val in (("expected_case_hits", s["expected"]), ("null_variance", s["variance"]),
                                     ("null_count_lo", s["lower"]), ("null_count_hi", s["upper"])):
                        check_close(float(val), prod[col], (screen, direction, col))
        numbers = (screen, direction, r["a"], r["n"] - r["a"], r["k"], r["h"] - r["k"],
                   float(s["expected"]), float(s["variance"]), s["lower"], s["upper"], float(effect),
                   float(s["p"]) if p1 is not None else "", float(p1) if p1 is not None else "",
                   float(r["fisher"]["p"]), float(p2) if p2 is not None else "", *r["odds"], r["maxdiff"])
        rows.append(numbers + ("PASS" if production else "FROZEN_INDEPENDENT",))
        print(screen, direction, "observed", r["k"], "expected", numbers[6], "p_holm", numbers[12],
              "Fisher OR", r["odds"][0], "Fisher Holm", numbers[14])
    save(checks_dir / "INDEPENDENT_LOF_EXACT.tsv",
         ("screen", "direction", "case_n", "control_n", "case_hits", "control_hits",
          "expected", "variance", "null_low95", "null_high95", "O_over_E", "exact_P", "Holm_P",
          "Fisher_P", "Fisher_Holm_P", "OR", "CI_low", "CI_high", "largest_probability_difference", "status"),
         rows)
    print(n_candidates, "frozen candidate-screen direction records and", n_frozen, "frozen R1 strict margins matched.")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("data", type=Path)
    cli.add_argument("output", type=Path)
    cli.add_argument("--production", type=Path, help="Contrast a finished R9 results directory with exact recomputation")
    cli.add_argument("--branch", default="all", choices=("a", "b", "all"))
    args = cli.parse_args()
    toy_n, toy_count = toy_check()
    keys, score1 = validate_background(load(args.data / "data/r7/inputs/R7_GENE_BACKGROUND.tsv.gz"))
    print(f"Confirmed {toy_n} toy strata using {toy_count} distinct exhaustive allocations and baseline 2×2 counts.")
    if args.branch in ("a", "all"):
        run_a(args.data, keys, score1, args.production, args.output)
    if args.branch in ("b", "all"):
        run_b(args.data, keys, args.production, args.output)


if __name__ == "__main__":
    main()
