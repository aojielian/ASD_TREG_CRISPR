"""Audit and summarize Chen et al. Supplementary Tables 2, 3 and 5."""
import argparse
from collections import defaultdict
from io import BytesIO
from pathlib import Path
import re
from zipfile import ZipFile

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
from scipy.stats import false_discovery_control

PANEL = "CREBBP CUL3 DYRK1A EP300 FOXP1 KDM3B KDM6B KMT5B PHIP SIN3A SMARCA4 TBL1XR1 TLE3 ZMYND8".split()
SEED = 20261125
DRAWS = 20_000
PROGRAMS = [
    "Chromatin regulation", "Program 2", "T cell effector response", "Program 4",
    "IL-2 signalling", "Oxidative phosphorylation", "Aerobic glycolysis",
    "Cytoskeleton remodelling", "Active translation", "T cell effector response",
]
REF = "(mean protein in target minus mean protein in NTC) / mean protein in NTC"
TRANSCRIPT_LIMIT = ("No tabular target-by-program coefficients in the audited Supplementary Tables "
                    "or authors' code repository. Published RNA programs are gene memberships, "
                    "and the target-module network summarizes modules. The source coefficient "
                    "matrix referenced by the analysis code was unavailable to R10.")


def book(archive, number):
    members = [n for n in archive.namelist() if f"Table {number}." in n]
    assert len(members) == 1
    return load_workbook(BytesIO(archive.read(members[0])), read_only=True, data_only=True)


def membership(sheet):
    headers, *rows = list(sheet.values)
    return pd.DataFrame([
        dict(gene=str(value), module=header, source_sheet=sheet.title,
             source_row=row_number, source_column=get_column_letter(column_number))
        for row_number, row in enumerate(rows, 2)
        for column_number, (header, value) in enumerate(zip(headers, row), 1)
        if value is not None
    ])


def measured_sets(xlsx):
    table = book(xlsx, 2)
    proteins = table["ProteinScore_Targets"]
    headers, *rows = list(proteins.values)
    targets = list(headers[1:])
    assert len(targets) == len(set(targets)) and None not in targets
    labels = [row[0] for row in rows]
    assert len(labels) == len(set(labels))
    score = pd.DataFrame([row[1:] for row in rows], index=labels, columns=targets).T
    assert np.isfinite(score.to_numpy()).all()
    table3 = book(xlsx, 3)
    modules = membership(table3["Target Modules"])
    programs = membership(table3["Gene Programs"])
    assert not modules.gene.duplicated().any()
    assert set(modules.gene) == set(score.index) | {"NTC"}
    guides = book(xlsx, 5)["Fig. 2 Perturb-icCITE"]
    library = defaultdict(list)
    for i, row in enumerate(guides.iter_rows(min_row=2, max_col=1, values_only=True), 2):
        name = row[0]
        if name is not None:
            match = re.fullmatch(r"(.+)_\d+", name)
            assert match is not None, name
            library[match.group(1)].append(i)
    assert set(score.index) <= set(library)
    return score, modules, programs, library


def target_covariates(score, bg):
    attr = pd.read_csv(bg, sep="\t", low_memory=False).query("screen == 'CHEN'").copy()
    assert not attr.approved_symbol.duplicated().any()
    attr = attr.set_index("approved_symbol", drop=False)
    key = {}
    for g in score.index:
        if g in attr.index:
            key[g] = g
        else:
            aliases = attr.loc[attr.raw_primary_id.eq(g), "approved_symbol"].dropna().unique()
            if len(aliases) == 1:
                key[g] = aliases[0]
    assert len(set(key.values())) == len(key)
    cls, lib = "class_label", "library_member"
    exp, lof = "HEALTHY_CD4_MEAN_LOGCPM", "LOEUF_v2_1_1"
    attrs = attr.reindex([key.get(g) for g in score.index])[
        ["approved_symbol", "sfari_score1", cls, lib, exp, lof,
         "CLASS_ANNOTATION_UNCERTAIN", "annotation_status"]].copy()
    attrs.index = score.index
    for column in [lib, "sfari_score1", "CLASS_ANNOTATION_UNCERTAIN"]:
        attrs[column] = attrs[column].astype("string").str.upper().map(
            {"TRUE": 1, "FALSE": 0, "1": 1, "0": 0})
    attrs["candidate"] = attrs.index.isin(PANEL)
    attrs["exclusion"] = ""
    attrs.loc[attrs.approved_symbol.isna(), "exclusion"] = "No unique frozen approved-symbol link"
    attrs.loc[attrs.sfari_score1.eq(1) & ~attrs.candidate, "exclusion"] = "Other SFARI score-1 target"
    known = (attrs[cls].isin(["NEITHER_LISTED", "TF_ONLY", "EPIFACTOR_ONLY", "BOTH_LISTED"])
             & attrs.CLASS_ANNOTATION_UNCERTAIN.eq(0)
             & attrs.annotation_status.eq("CLASS_ANNOTATION_USABLE"))
    numeric = np.isfinite(attrs[exp]) & np.isfinite(attrs[lof])
    available = known & attrs[lib].isin([0, 1]) & attrs.sfari_score1.isin([0, 1]) & numeric
    attrs.loc[attrs.exclusion.eq("") & ~available, "exclusion"] = "Missing or ineligible frozen covariate"
    cc = attrs.loc[attrs.exclusion.eq("")].copy()
    for column in [exp, lof]:
        cc[column + "_bin"] = np.ceil(3 * cc[column].rank(method="average") / len(cc)).astype(int)
    cc["stratum"] = cc[[cls, lib, exp + "_bin", lof + "_bin"]].astype(str).agg("|".join, axis=1)
    strata_counts = cc.groupby("stratum")["candidate"].agg(["sum", "size"])
    mixed = strata_counts.index[strata_counts["sum"].gt(0) & strata_counts["sum"].lt(strata_counts["size"])]
    attrs["stratum"] = cc.stratum
    attrs.loc[attrs.exclusion.eq("") & ~attrs.stratum.isin(mixed), "exclusion"] = "No within-stratum comparison"
    return attrs, exp, lof


def conditional_protein(score, attrs, exp, lof):
    supported = attrs.loc[attrs.exclusion.eq("")]
    group_a = supported[supported.candidate]
    group_b = supported[~supported.candidate]
    n_a, n_b = len(group_a), len(group_b)
    group_n = supported.groupby("stratum").candidate.agg(["sum", "size"])
    group_n["other"] = group_n["size"] - group_n["sum"]
    weights = group_n["sum"] / group_n["other"]
    balance = []
    for name in [exp, lof]:
        a, b = group_a[name], group_b[name]
        if n_a >= 2 and n_b >= 2:
            denom = np.sqrt((a.var() + b.var()) / 2)
            if denom:
                unweighted = (a.mean() - b.mean()) / denom
                weighted = (a.mean() - np.average(b, weights=group_b.stratum.map(weights))) / denom
            else:
                unweighted = weighted = np.nan
        else:
            unweighted = weighted = np.nan
        balance.append(dict(attribute=name, n_candidate=n_a, n_nonASD=n_b,
                            smd_within_supported_unweighted=unweighted, smd_conditional=weighted))
    if n_a < 5 or n_b < 5:
        row = dict(analysis="protein", status="DESCRIPTIVE_ONLY", n_candidate=n_a, n_nonASD=n_b,
                   reason="Fewer than five fixed candidates or fewer than five non-ASD targets in informative strata.")
        return pd.DataFrame([row]), pd.DataFrame(balance), group_n
    rng = np.random.default_rng(SEED)
    expectation = np.zeros(score.shape[1])
    null = np.zeros((DRAWS, score.shape[1]))
    for _, part in supported.groupby("stratum"):
        n = int(part.candidate.sum())
        values = score.loc[part.index].to_numpy()
        expectation += n * values.mean(axis=0)
        indices = np.argsort(rng.random((DRAWS, len(part))), axis=1)[:, :n]
        for i in range(n):
            null += values[indices[:, i]]
    expectation /= n_a
    null /= n_a
    observed = score.loc[group_a.index].mean(axis=0).to_numpy()
    effect = observed - expectation
    extreme = (np.abs(null - expectation) >= np.abs(effect) - 1e-14).sum(axis=0)
    p = (extreme + 1) / (DRAWS + 1)
    bounds = np.quantile(null, [0.025, 0.975], axis=0)
    output = pd.DataFrame(dict(
        analysis="protein", status="EXPLORATORY", assay_label=score.columns, n_candidate=n_a,
        n_nonASD=n_b, n_strata=len(group_n), observed_mean=observed, expected_mean=expectation,
        observed_minus_expected=effect, null_low=bounds[0], null_high=bounds[1], permutations=DRAWS,
        seed=SEED, extreme=extreme, p=p, mc_se=np.sqrt(p * (1 - p) / (DRAWS + 1)),
        q_BH_88=false_discovery_control(p), reason="Fixed selected-target background; gene-level exchangeability"))
    return output, pd.DataFrame(balance), group_n


def main(download, data, out):
    path = download / "Chen2025" / "41586_2025_8795_MOESM5_ESM.zip"
    with ZipFile(path) as xlsx:
        score, modules, programs, library = measured_sets(xlsx)
    score_rows = []
    for gene in PANEL:
        for i, assay in enumerate(score.columns, 2):
            measured = gene in score.index
            score_rows.append(dict(
                gene=gene, protein=assay.split("-", 2)[2], assay_label=assay,
                compartment=assay.split("-", 1)[0], score=score.at[gene, assay] if measured else np.nan,
                reference_score=0 if measured else np.nan,
                score_minus_reference=score.at[gene, assay] if measured else np.nan,
                source_sheet="ProteinScore_Targets", source_row=i if measured else np.nan,
                source_column=get_column_letter(score.index.get_loc(gene) + 2) if measured else "",
                reference="Zero is the formula-defined NTC baseline; no NTC score column",
                missing_reason="" if measured else "No published Table 2 score",
            ))
    protein_long = pd.DataFrame(score_rows)
    coverage = []
    for gene in PANEL:
        entry = modules[modules.gene.eq(gene)]
        program = programs[programs.gene.eq(gene)]
        coverage.append(dict(
            gene=gene, targeted_in_input_library=gene in library, initial_guide_count=len(library.get(gene, [])),
            library_source_sheet="Fig. 2 Perturb-icCITE", library_source_rows=";".join(map(str, library.get(gene, []))),
            protein_score_available=gene in score.index, n_protein_scores=score.shape[1] if gene in score.index else 0,
            protein_source_sheet="ProteinScore_Targets",
            protein_source_column=get_column_letter(score.index.get_loc(gene) + 2) if gene in score.index else "",
            target_module=";".join(entry.module), target_module_source_sheet="Target Modules",
            target_module_source_rows=";".join(entry.source_row.astype(str)),
            transcript_program_member=";".join(program.module), program_source_sheet="Gene Programs",
            transcript_program_source_rows=";".join(program.source_row.astype(str)),
            perturbation_program_coefficient_available=np.nan,
            program_coefficient_status="NOT_AVAILABLE",
            status=("PROTEIN_SCORED" if gene in score.index else "LIBRARY_ONLY" if gene in library else "NOT_IN_LIBRARY"),
        ))
    program_summary = pd.DataFrame([
        dict(program=i, source_sheet="Gene Programs", source_column=get_column_letter(i),
             author_description=label, n_member_rows=int(programs.module.eq(f"Module {i}").sum()),
             panel_members=";".join(programs.loc[programs.module.eq(f"Module {i}") & programs.gene.isin(PANEL), "gene"]),
             program_coefficient_status="NOT_AVAILABLE", missing_reason=TRANSCRIPT_LIMIT)
        for i, label in enumerate(PROGRAMS, 1)
    ])
    modules["panel_candidate"] = modules.gene.isin(PANEL)
    modules["is_NTC"] = modules.gene.eq("NTC")
    modules["n_nontarget_control_rows"] = len(library.get("NTC", []))
    programs["panel_member"] = programs.gene.isin(PANEL)
    audits = pd.DataFrame([
        [2, "ProteinScore_Targets", 89, 228, "88 antibody/epitope rows; 227 targets; complete numeric matrix", REF],
        [3, "Gene Programs", 744, 10, "10 program member columns; 2,192 rows-of-memberships in total", "RNA program membership"],
        [3, "Target Modules", 75, 11, "11 target module columns; 227 scored targets plus NTC", "Perturbation-target membership"],
        [3, "GO Enrichment", 3378, 10, "3,377 table rows excluding header", "Author enrichment terms, not individual target responses"],
        [5, "Fig. 2 Perturb-icCITE", 908, 2, "907 guide IDs including 20 NTC, 295 gene targets before model QC", "Initial target library; sequences excluded"],
    ], columns=["table", "sheet", "rows_with_header", "columns_examined", "scope", "definition"])
    attrs, exp, lof = target_covariates(score, data / "data/r7/inputs/R7_GENE_BACKGROUND.tsv.gz")
    exploration, balance, strata = conditional_protein(score, attrs, exp, lof)
    exploration = pd.concat([exploration, pd.DataFrame([
        dict(analysis="transcript_program", status="NOT_AVAILABLE", reason=TRANSCRIPT_LIMIT)
    ])], ignore_index=True)
    assert len(coverage) == 14 and len(protein_long) == 14 * 88
    assert int(protein_long.score.notna().sum()) == sum(v["n_protein_scores"] for v in coverage)
    assert len(programs) == 2192
    outputs = {
        "SOURCE_SHEET_AUDIT": audits, "COVERAGE": pd.DataFrame(coverage),
        "PROTEIN": protein_long, "PROGRAM": program_summary, "MODULE": modules,
        "PROGRAM_MEMBERS": programs, "EXPLORATORY": exploration, "PROTEIN_BALANCE": balance,
        "PROTEIN_STRATA": strata.reset_index(),
        "TARGET_BACKGROUND": attrs.rename_axis("source_target").reset_index(),
    }
    for name, df in outputs.items():
        df.to_csv(out / f"R10_PERTURB_{name}.tsv", sep="\t", index=False, na_rep="NA")
    print(f"Library gene targets={len(library)-1}, model scored={len(score)}, epitopes={score.shape[1]}")
    print(pd.DataFrame(coverage).groupby("status").size().to_string())
    print(exploration.groupby(["analysis", "status"]).size().to_string())
    print(balance.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("download_files", type=Path)
    parser.add_argument("frozen_data", type=Path)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    main(args.download_files, args.frozen_data, args.results)
