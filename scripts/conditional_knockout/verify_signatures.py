"""Recompute scores and t contrasts independently from raw counts and frozen members."""

import argparse
import csv
import gzip
from hashlib import file_digest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


FILES = {
    "GSE121279": "GSE121279_rnaseq-counts-raw.txt.gz",
    "GSE263830": "GSE263830_Sin3aTreg_gene_count.txt.gz",
    "GSE121251": "GSE121251_Counts.mRNA.foxp1.txt.gz",
}
N_MIN = {"GSE121279": 3, "GSE263830": 3, "GSE121251": 2}
CONTRASTS = {
    "F1N": "naive",
    "F1A": "activated",
    "S3": None,
    "F2R": "resting",
    "F2A": "activated",
}
FIELDS = ["estimate", "SE", "CI_low", "CI_high", "P", "df", "n_KO", "n_control"]


def parse_tab(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def sha256(path):
    with path.open("rb") as file:
        return file_digest(file, "sha256").hexdigest()


def raw_counts(root, dataset):
    path = root / dataset / FILES[dataset]
    # An absent Index heading is explicitly supplied for the 2019 Konopacki file.
    if dataset == "GSE121279":
        with gzip.open(path, "rt") as stream:
            fields = next(csv.reader(stream, delimiter="\t"))
        return pd.read_csv(path, sep="\t", names=["gene_id", *fields], header=0, index_col=0)
    table = parse_tab(path)
    return table.set_index(table.columns[0])


def compare_membership(members, data, stage, mgi_path):
    active = members.loc[members.included]
    hill = active.loc[active.program.eq("Hill_LN")]
    targets = {"Foxp3", "Foxp1", "Sin3a"}
    assert targets.isdisjoint(hill.mouse_gene)
    assert not hill.mouse_gene.duplicated().any()
    if mgi_path is not None:
        table = pd.read_csv(mgi_path, sep="\t", header=None, dtype=str)
        conflict = pd.Series(False, index=table.index)
        for column in table.columns:
            conflict |= table[column].duplicated(keep=False) | table[column].isna()
        mgi = table.loc[~conflict].set_index(0)
        for entry in hill.itertuples():
            row = mgi.loc[entry.mgi_id]
            assert row[1] == entry.mouse_gene
            # R11 retains matches if official human symbols or Entrez IDs agree.
            assert row[4] == entry.source_symbol or float(row[5]) == entry.source_gene_id

    source = data / "data/r10/results/R10_IMMGEN_CELLTYPE_CPM.tsv.gz"
    atlas = parse_tab(source)
    atlas = atlas.loc[~atlas.mouse_gene.isin(targets)].copy()
    tconv = atlas[["T.4.Nve.Sp", "T.4.Nve.Fem.Sp"]].mean(axis=1)
    treg = atlas["Treg.4.25hi.Sp"]
    atlas["log_ratio"] = np.log2((treg + 0.5) / (tconv + 0.5))
    atlas = atlas.loc[np.maximum(treg, tconv) >= 1]
    for arm, sign in [("Treg_high", 1), ("Tconv_high", -1)]:
        base = atlas.loc[sign * atlas.log_ratio > 0]
        computed = base.sort_values(["log_ratio", "mouse_gene"], ascending=[sign == -1, True]).head(100)
        shipped = active.query("program == 'ImmGen_spleen' and component == @arm").sort_values("reference_rank")
        assert shipped.mouse_gene.tolist() == computed.mouse_gene.tolist()
        assert np.allclose(shipped.reference_log2_ratio, computed.log_ratio, rtol=0, atol=2e-13)

    gmt = stage / "00_download/files/MSigDB/R11_MSIGDB_MOUSE_MH_2025.1_SELECTED.gmt"
    with gmt.open() as stream:
        audit = {row[0]: row[2:] for row in csv.reader(stream, delimiter="\t")}
    for program, symbols in active[active.program.str.startswith("HALLMARK_")].groupby("program"):
        assert symbols.mouse_gene.tolist() == audit[program]
    assert len(audit) == 6


def derive_scores(count_root, sample_map, id_map, members, coverage):
    programs = {}
    measured = {}
    for dataset in FILES:
        columns = sample_map.loc[sample_map.dataset.eq(dataset), "sample_id"].tolist()
        matrix = raw_counts(count_root, dataset).loc[:, columns]
        id_table = id_map.loc[id_map.dataset.eq(dataset)].set_index("gene_id")
        id_table = id_table.loc[matrix.index.intersection(id_table.index)].rename_axis("gene_id")
        matrix = matrix.loc[id_table.index]
        keep = matrix.ge(10).sum(axis=1).ge(N_MIN[dataset])
        assert keep.eq(id_table.pass_count_filter).all()
        values = matrix.loc[keep]
        midrank = pd.DataFrame(
            (stats.rankdata(values.to_numpy(), axis=0, method="average") - 0.5) / len(values),
            index=values.index,
            columns=values.columns,
        )
        for (program, component), rows in members[members.included].groupby(["program", "component"]):
            by_symbol = id_table.loc[id_table.symbol_unique & keep].reset_index().set_index("gene_symbol")
            assert by_symbol.index.is_unique
            feature_ids = by_symbol.loc[by_symbol.index.intersection(rows.mouse_gene), "gene_id"].tolist()
            info = coverage.query("dataset == @dataset and program == @program and component == @component").iloc[0]
            assert len(feature_ids) == info.measured_n
            assert len(rows) == info.frozen_member_n
            assert info.measured_n + info.count_filter_failed_n + info.absent_or_unmapped_n + info.ambiguous_symbol_n == len(rows)
            assert len(values) == info.rank_universe_n
            measured[dataset, program, component] = len(feature_ids)
            means = midrank.loc[feature_ids].mean(axis=0)
            programs[dataset, program, "score" if component == "member_mean" else component] = means
        for program in ["Hill_LN", "ImmGen_spleen"]:
            programs[dataset, program, "score"] = (
                programs[dataset, program, "Treg_high"] - programs[dataset, program, "Tconv_high"]
            )
    return programs, measured


def independent_interval(values, design, contrast):
    """The two cohorts with interactions reduce to independent state by genotype cell means."""
    selected_state = CONTRASTS[contrast]
    values = values.reindex(design.index)
    if selected_state is not None:
        groups = pd.DataFrame({"score": values, "state": design.state, "genotype": design.genotype})
        means = groups.groupby(["state", "genotype"]).score.transform("mean")
        sigma_num = np.square(values - means).sum()
        df = len(design) - 4
        case = groups.loc[groups.state.eq(selected_state) & groups.genotype.eq("KO"), "score"]
        control = groups.loc[groups.state.eq(selected_state) & groups.genotype.eq("control"), "score"]
        estimate = case.mean() - control.mean()
        se = np.sqrt(sigma_num / df * (1 / len(case) + 1 / len(control)))
        n_ko, n_control = len(case), len(control)
    else:
        # Frisch-Waugh-Lovell via within-sex residualization plus one treatment coefficient.
        centered = values - values.groupby(design.sex).transform("mean")
        exposure = design.genotype.eq("KO").astype(float)
        centered_exposure = exposure - exposure.groupby(design.sex).transform("mean")
        denom = np.square(centered_exposure).sum()
        estimate = (centered_exposure * centered).sum() / denom
        resid = centered - estimate * centered_exposure
        df = len(design) - 3
        se = np.sqrt(np.square(resid).sum() / df / denom)
        n_ko, n_control = int(exposure.sum()), len(exposure) - int(exposure.sum())
    lo, hi = estimate + np.array([-1, 1]) * stats.t.ppf(0.975, df) * se
    p = stats.t.sf(abs(estimate / se), df) * 2
    return [estimate, se, lo, hi, p, df, n_ko, n_control]


def holm(values):
    order = np.argsort(values)
    adjusted = np.maximum.accumulate(values[order] * np.arange(len(values), 0, -1))
    out = np.ones(len(values))
    out[order] = np.minimum(adjusted, 1)
    return out


def run(data, stage, count_root, mgi_path=None):
    result_dir = stage / "03_results"
    members_path = result_dir / "R11_SIGNATURE_MEMBERS.tsv.gz"
    reference = data / "data/r11/results/R11_SIGNATURE_MEMBERS.tsv.gz"
    assert sha256(members_path) == sha256(reference)
    members = parse_tab(members_path)
    compare_membership(members, data, stage, mgi_path)
    meta = parse_tab(stage / "01_audit/R11_SAMPLE_MAP.tsv")
    meta = meta.loc[meta.include]
    ids = parse_tab(stage / "01_audit/R11_ID_MAP.tsv.gz")
    coverage = parse_tab(result_dir / "R11_SIGNATURE_COVERAGE.tsv")
    scores, _ = derive_scores(count_root, meta, ids, members, coverage)
    supplied = parse_tab(result_dir / "R11_SIGNATURE_SAMPLE_SCORES.tsv")
    errors = []
    for row in supplied.itertuples():
        expected = scores[row.dataset, row.program, row.component].at[row.sample_id]
        errors.append(abs(expected - row.score))
    assert max(errors) < 2e-13
    effects = parse_tab(result_dir / "R11_SIGNATURE_CONTRASTS.tsv")
    diff = []
    for index, row in effects.iterrows():
        design = meta.loc[meta.dataset.eq(row.dataset)].set_index("sample_id")
        expected = independent_interval(scores[row.dataset, row.program, row.component], design, row.contrast)
        source = row[FIELDS].to_numpy(dtype=float)
        compare_p = row.component == "score"
        if not compare_p:
            assert np.isnan(source[4])
            source[4] = expected[4]
        error = np.max(np.abs(expected - source))
        diff.append(error)
    assert max(diff) < 2e-11
    for family, df in effects.loc[effects.component.eq("score")].groupby("adjustment_family"):
        values = df.P.to_numpy()
        corrected = holm(values) if family.startswith("Hill") else stats.false_discovery_control(values, method="bh")
        assert np.allclose(corrected, df.adjusted_P, rtol=0, atol=2e-14)
    camera = parse_tab(result_dir / "R11_SIGNATURE_ROBUSTNESS.tsv")
    for _, df in camera.groupby("adjustment_family"):
        assert list(df.adjustment_family_size.unique()) == [len(df)]
        corrected = stats.false_discovery_control(df.CAMERA_P, method="bh")
        assert np.allclose(corrected, df.adjusted_P, rtol=0, atol=2e-14)
    assert len(camera) == 50 and camera.status.eq("estimated").all()
    hill_detail = ("Independent exact MGI checks for all Hill members." if mgi_path else
                   "Hill membership is frozen and SHA-256 verified.")
    rows = [
        ("frozen_members", "PASS", "Analysis membership exactly matches the released frozen membership file."),
        ("membership_reconstruction", "PASS", hill_detail + " Full ImmGen top-100 reconstruction on both arms; six native Hallmark lists identical."),
        ("cohort_coverage", "PASS", f"{len(coverage)} measured arm/universe margins from the raw matrices."),
        ("sample_scores", "PASS", f"{len(supplied)} scores recomputed via raw-rank invariance, max absolute difference {max(errors):.3g}."),
        ("ols_effects_ci", "PASS", f"{len(effects)} effects/SE/CI checked by cell means or within-sex regression; 40 tested score P, max absolute difference {max(diff):.3g}."),
        ("multiple_testing", "PASS", "Four score families and three analytic CAMERA BH families verified."),
        ("camera_cardinality", "PASS", "All 50 planned set-by-contrast competitive sensitivity tests present."),
    ]
    path = stage / "07_qa/R11_SIGNATURE_NUMERIC_QA.tsv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "status", "evidence"])
        writer.writerows(rows)
    for row in rows:
        print(*row, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("stage", type=Path)
    parser.add_argument("counts", type=Path)
    parser.add_argument("--mgi-report", type=Path)
    args = parser.parse_args()
    run(args.data, args.stage, args.counts, args.mgi_report)
