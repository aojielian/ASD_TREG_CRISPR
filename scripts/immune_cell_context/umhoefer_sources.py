"""GSE286470 and frozen R9/R7 input readers for the R10 direct comparison."""

import gzip
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

from strict_xlsx import read_strict_xlsx

CONTEXTS = {"Treg": "UMHOEFER_TREG", "Tconv": "UMHOEFER_TCONV"}
GENE_SHEET = "Combined_gene_summary"
GUIDE_SHEET = "Combined_grna_summary"
SOURCE_STUDY = "GSE286470; PMID 41237776; DOI 10.1016/j.immuni.2025.10.020"
SOURCE_LINK = "https://doi.org/10.1016/j.immuni.2025.10.020"
CONTROL_IDS = {"Non-Targeting": "non-targeting guide pool", "EGFP": "non-human EGFP-targeting control"}
PANEL = (
    "CREBBP", "CUL3", "DYRK1A", "EP300", "FOXP1", "KDM3B", "KDM6B",
    "KMT5B", "PHIP", "SIN3A", "SMARCA4", "TBL1XR1", "TLE3", "ZMYND8",
)
INPUT_TABLES = {
    "background": "data/r7/inputs/R7_GENE_BACKGROUND.tsv.gz",
    "chen": "data/processed/chen_genomewide_analysis_ready.tsv",
    "shifrut_common": "data/processed/chen_shifrut_common_universe.tsv",
    "Treg": "data/processed/umhoefer_primary_treg_analysis_ready.tsv",
    "Tconv": "data/processed/umhoefer_primary_tconv_analysis_ready.tsv",
}
PROPERTIES = {
    "class_label": "class_label",
    "annotation_status": "annotation_status",
    "CLASS_ANNOTATION_UNCERTAIN": "class_uncertain",
    "library_member": "targeted_library",
    "TF_LISTED": "TF",
    "EPIFACTORS_GENES_LISTED": "EpiFactors",
    "HEALTHY_CD4_MEAN_LOGCPM": "expression",
    "LOEUF_v2_1_1": "LOEUF",
    "GENE_SPAN_BP": "gene_span_bp",
}
BIN_COUNT = {"Treg": 2, "Tconv": 3}


def write_tsv(frame, path):
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    pd.DataFrame(frame).to_csv(path, sep="\t", index=False, na_rep="NA",
                               compression=compression, float_format="%.17g")


def verify(condition, test, observation, checks):
    checks.append(dict(check=test, status="PASS" if condition else "FAIL", observation=observation))
    if not condition:
        raise ValueError(f"{test}: {observation}")


def workbook(path, sheet):
    # Reuse the strict-OOXML reader previously validated against the same
    # workbook hashes. The ordinary openpyxl loader misses these Strict sheets.

    with ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise ValueError(f"Workbook failed CRC: {path.name}")
        root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        source_rows = [int(row.attrib["r"]) for row in root.findall(".//{*}sheetData/{*}row")]
    sheets = read_strict_xlsx(path)
    if set(sheets) != {sheet} or source_rows != list(range(1, len(source_rows) + 1)):
        raise ValueError(f"Unexpected sheets or discontinuous physical row IDs: {path.name}")
    header, *rows = sheets[sheet]
    if len(rows) != len(source_rows) - 1 or len(set(header)) != len(header):
        raise ValueError(f"Invalid table structure for {path.name}")
    return pd.DataFrame(rows, columns=header).assign(source_row=source_rows[1:])


def release_tables(data, checks):
    # Each original v1.1.0 file is checked against the preserved manifest before replay.
    from replay_inputs import verify_baseline_file
    out = {}
    for key, relative in INPUT_TABLES.items():
        source = verify_baseline_file(data, relative)
        verify(True, f"frozen_hash_{key}", f"{relative} matches the v1.1.0 release manifest", checks)
        out[key] = pd.read_csv(source, sep="\t", low_memory=False)
    return out


def typed_numeric(data, cols):
    for col in cols:
        data[col] = pd.to_numeric(data[col].replace("", np.nan), errors="raise")


def soft_samples(path, checks):
    samples, current = [], None
    series_pmid, last_update, study_samples = [], None, []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            key, _, value = line.rstrip("\n").partition(" = ")
            if key == "!Series_pubmed_id":
                series_pmid.append(value)
            if key == "!Series_last_update_date":
                last_update = value
            if key == "!Series_sample_id":
                study_samples.append(value)
            if key == "^SAMPLE":
                current = {"accession": value}
                samples.append(current)
            if current is None:
                continue
            if key == "!Sample_characteristics_ch1":
                name, val = value.split(": ", maxsplit=1)
                current[name] = val
            elif key in ("!Sample_title", "!Sample_description", "!Sample_organism_ch1",
                         "!Sample_geo_accession"):
                current[key] = value
    if not series_pmid or "41237776" not in series_pmid:
        raise ValueError("SOFT series does not identify the experimental report PMID 41237776")
    rows = []
    for raw in samples:
        if raw["accession"] != raw["!Sample_geo_accession"]:
            raise ValueError("Conflicting SOFT GSM keys")
        rows.append(dict(
            dataset="GSE286470", GSM=raw["accession"], original_label=raw["!Sample_title"],
            matrix_sample_label=raw["!Sample_description"].removeprefix("Library name: "),
            donor=raw["donor"], replicate=raw["group"], species=raw["!Sample_organism_ch1"],
            cell_type=raw["cell type"], lineage="CD4 T cell", tissue=raw["tissue"],
            state="Stimulated / restimulated, followed by FOXP3 high/low sort",
            included=True, reason="Contributes to the author pooled gene-level context summary",
            donor_role="Within-context donor label shared across high and low bins; cross-context linkage unverified",
            source=path.name, series_last_update=last_update, GEO_pubmed_ids="|".join(series_pmid),
        ))
    table = pd.DataFrame(rows)
    verify(not table["GSM"].duplicated().any() and set(table["GSM"]) == set(study_samples),
           "geo_samples_unique", f"{len(table)} GSM records match the accession series list", checks)
    groups = table.groupby(["cell_type", "donor"])["replicate"].agg(frozenset)
    verify(groups.eq(frozenset({"FOXP3 high", "FOXP3 low"})).all()
           and groups.groupby(level=0).size().to_dict() == BIN_COUNT,
           "geo_bin_and_donor_resolution", "2 Treg and 3 Tconv donor labels, each with high/low FOXP3 bins", checks)
    return table


def workbook_genes(path, tables, checks):
    raw = workbook(path, GENE_SHEET)
    if set(raw["celltype"]) != set(CONTEXTS) or raw.duplicated(["celltype", "id"]).any():
        raise ValueError("Duplicate gene-sheet keys or unexpected dispersion of cell contexts")
    numeric = [col for col in raw if col not in ("id", "celltype", "source_row")]
    typed_numeric(raw, numeric)
    ferr = ["neg|fdr", "pos|fdr"]
    p_cols = ["neg|p-value", "pos|p-value", *ferr]
    verify(raw[p_cols].apply(lambda s: s.between(0, 1)).all(axis=None)
           and raw[["neg|lfc", "pos|lfc", "num"]].apply(np.isfinite).all(axis=None)
           and raw["num"].gt(0).all(),
           "gene_workbook_numeric_fields", f"{len(raw)} combined gene rows; p and FDR in [0,1]; finite signed LFC", checks)
    verify(raw["neg|lfc"].eq(raw["pos|lfc"]).all(),
           "gene_workbook_lfc_enumeration", "Both gene-level MAGeCK LFC columns agree exactly", checks)
    controls = raw[raw["id"].isin(CONTROL_IDS)].copy()
    controls["exclusion_reason"] = controls["id"].map(CONTROL_IDS)
    verify(len(controls) == len(CONTROL_IDS) * len(CONTEXTS), "workbook_controls",
           "Pooled Non-Targeting (593 guides) and EGFP control (13 guides) excluded in both cell contexts",
           checks)
    if not controls["num"].eq(controls["id"].map({"Non-Targeting": 593, "EGFP": 13})).all():
        raise ValueError("Unexpected guide counts among named controls")
    raw = raw[~raw["id"].isin(CONTROL_IDS)].copy()
    raw["hit"] = raw[ferr].le(.05).any(axis=1)
    raw["signed_foxp3_bin_lfc"] = -raw["pos|lfc"]
    background = tables["background"]
    all_cases = background.loc[background["sfari_score1"]].groupby("screen")["approved_symbol"].agg(frozenset)
    verify(all_cases.map(len).eq(245).all() and len(set(all_cases)) == 1,
           "frozen_full_asd_set", "245 unique resolved frozen SFARI score-1 genes in each v1.1.0 background",
           checks)
    score1 = all_cases.iloc[0]
    out = {}
    for context, screen in CONTEXTS.items():
        release = background[background["screen"].eq(screen) & background["measured"]].copy()
        release["source_row"] = release["source_row"].astype(int)
        frozen_names = {
            "neg|lfc": "frozen_neg_lfc", "pos|lfc": "frozen_pos_lfc",
            "neg|fdr": "frozen_neg_fdr", "pos|fdr": "frozen_pos_fdr",
            "num": "frozen_num", "raw_primary_id": "frozen_raw_id",
        }
        release = release.rename(columns=frozen_names | PROPERTIES)
        left = raw[raw["celltype"].eq(context)]
        pair = left.merge(release, on="source_row", how="outer", validate="one_to_one", indicator=True)
        verify(pair["_merge"].eq("both").all() and pair["id"].eq(pair["frozen_raw_id"]).all(),
               "georow_join_" + context, f"{len(left)} gene rows map uniquely by physical sheet row and raw ID; two named controls excluded",
               checks)
        verify(not pair["hgnc_id"].duplicated().any() and not pair["approved_symbol"].duplicated().any()
               and pair["hgnc_id"].str.fullmatch(r"HGNC:\d+").all() and pair["canonical_status"].eq("RETAIN_UNIQUE").all(),
               "approved_ids_" + context,
               f"{pair['mapping_route'].eq('approved_symbol_exact').sum()} exact and "
               f"{pair['mapping_route'].eq('previous_symbol_unique').sum()} frozen prior-symbol assignments",
               checks)
        fields = list(frozen_names)[:-1]
        verify(all(np.allclose(pair[col], pair[frozen_names[col]], atol=1e-12, rtol=1e-12) for col in fields)
               and pair["hit"].eq(pair["author_hit"]).all() and pair["sfari_score1"].eq(pair["approved_symbol"].isin(score1)).all(),
               "geo_release_parity_" + context, "All raw FDR, LFC, guide count, hit and frozen case labels reproduce",
               checks)
        mature = tables[context].set_index("gene")
        pair = pair.rename(columns={"approved_symbol": "gene"}).set_index("gene", verify_integrity=True)
        m = mature.reindex(pair.index)
        verify(len(pair) == len(mature) and m["source_row"].eq(pair["source_row"]).all()
               and m["screen_hit"].eq(pair["hit"]).all() and m["sfari_score1"].eq(pair["sfari_score1"]).all()
               and np.allclose(m["author_lfc"], pair["pos|lfc"], rtol=1e-12, atol=1e-12),
               "reproduce_frozen_individual_universe_" + context,
               "Physical Excel rows, cases, hits and effect values match the R9-verified v1.1.0 analysis tables", checks)
        verify(pair.loc["FOXP3", "signed_foxp3_bin_lfc"] < 0 and pair.loc["FOXP3", "hit"],
               "source_direction_control_" + context,
               "FOXP3 coding-target positive native LFC corresponds to operational lower FOXP3 after KO", checks)
        pair["direction"] = pd.Series(np.sign(pair["signed_foxp3_bin_lfc"]), index=pair.index).map(
            {-1.0: "KO-reduces-FOXP3", 0.0: "zero-estimate", 1.0: "KO-increases-FOXP3"})
        release_direction = pair["direction"].map(
            {"KO-reduces-FOXP3": "SUPPORTS_FOXP3", "KO-increases-FOXP3": "RESTRAINS_FOXP3",
             "zero-estimate": "SIGN_UNRESOLVED"})
        verify(release_direction[pair["hit"]].eq(pair.loc[pair["hit"], "FOXP3_REGULATOR_DIRECTION"]).all(),
               "direction_parity_" + context, "Every significant gene agrees with the frozen direction classification",
               checks)
        out[context] = pair
    return out, controls, raw, score1


def candidate_guides(path, screens, checks):
    guides = workbook(path, GUIDE_SHEET)
    numeric = ["control_mean", "treat_mean", "LFC", "control_var", "adj_var", "score", "p.low",
               "p.high", "p.twosided", "FDR", "high_in_treatment"]
    typed_numeric(guides, numeric)
    verify(not guides.duplicated(["celltype", "sgrna"]).any()
           and set(guides["celltype"]) == set(CONTEXTS)
           and guides[["p.low", "p.high", "p.twosided", "FDR"]].apply(lambda x: x.between(0, 1)).all(axis=None)
           and np.isfinite(guides["LFC"]).all(),
           "guide_workbook_numeric_fields", "Unique guide ID per cell context; guide-level MAGeCK FDR and P in [0,1]",
           checks)
    for field in ["control_count", "treatment_count"]:
        values = guides[field].str.split("/")
        arity = guides["celltype"].map(BIN_COUNT)
        verify(values.str.len().eq(arity).all() and values.explode().astype(float).ge(0).all(),
               "guide_vector_" + field,
               "Nonnegative normalized count vectors contain 2 Treg/3 Tconv unlabelled positions", checks)
    rows, stats, discrepancies = [], [], []
    for context, measured in screens.items():
        selected = guides[guides["celltype"].eq(context)]
        all_sizes = selected.groupby("Gene").size()
        observed = measured["id"].map(all_sizes)
        disagreements = ~observed.eq(measured["num"])
        verify(selected[selected["Gene"].eq("FOXP3")]["LFC"].gt(0).all(),
               "guide_positive_control_" + context,
               "All FOXP3 coding-control guides enrich the low FOXP3 bin",
               checks)
        for gene, rec in measured[disagreements].iterrows():
            discrepancies.append(dict(gene=gene, raw_id=rec["id"], cell_type=context,
                                      gene_sheet_num=rec["num"], guide_rows_in_guide_sheet=observed[gene],
                                      gene_sheet_source_row=rec["source_row"],
                                      interpretation="Retained the author's evaluable gene-level FDR and effect; guide discrepancy is unresolved"))
        checks.append(dict(check="guide_gene_coverage_" + context,
                           status="NOTE_SOURCE_DISCREPANCY" if disagreements.any() else "PASS",
                           observation=f"{int(disagreements.sum())} of {len(measured)} measured genes differ between pooled n and guide-table rows"))
        for gene in PANEL:
            if gene not in measured.index:
                rows.append(dict(gene=gene, cell_type=context, status="NOT_MEASURED",
                                 reason="Absent from this targeted human gene workbook and guide summary"))
                continue
            pool = selected[selected["Gene"].eq(measured.loc[gene, "id"])]
            gene_sign = np.sign(measured.loc[gene, "signed_foxp3_bin_lfc"])
            n_agree = np.sign(-pool["LFC"]).eq(gene_sign).sum()
            rows.append(dict(gene=gene, cell_type=context, status="DESCRIPTIVE_ONLY",
                             context_n_donor_labels=BIN_COUNT[context], guide_n=len(pool), gene_concordant_n=n_agree,
                             ko_lower_n=pool["LFC"].gt(0).sum(), ko_higher_n=pool["LFC"].lt(0).sum(),
                             zero_estimate_n=pool["LFC"].eq(0).sum(),
                             min_guide_FDR=pool["FDR"].min(), source_sheet=GUIDE_SHEET,
                             source_rows="|".join(pool["source_row"].astype(str)),
                             reason="Pooled guide summary; individual donor labels/order not given for slash-count positions"))
            for rec in pool.to_dict("records"):
                stats.append(dict(gene=gene, cell_type=context, source_guide_id=rec["sgrna"],
                                  source_sheet=GUIDE_SHEET, source_row=rec["source_row"],
                                  source_LFC=rec["LFC"], signed_foxp3_bin_lfc=-rec["LFC"],
                                  guide_FDR=rec["FDR"], count_vector_positions=BIN_COUNT[context],
                                  normalization="Author MAGeCK; no donor-level gene effect or donor crosswalk"))
    return pd.DataFrame(rows), pd.DataFrame(stats), guides, pd.DataFrame(discrepancies)


def shifrut_panel(shifrut_table, tables, checks):
    data = pd.read_csv(shifrut_table, sep="\t")
    verify(len(data) == 19094 and not data["approved_symbol"].duplicated().any()
           and data["SHIFRUT_HIT"].eq(data[["neg_fdr", "pos_fdr"]].le(.05).any(axis=1)).all(),
           "shifrut_canonical_background", "Historical frozen Shifrut gene-level data match size, unique symbols and FDR rule",
           checks)
    data = data.set_index("approved_symbol")
    common = tables["shifrut_common"].set_index("approved_symbol")
    for gene in PANEL:
        if gene in common.index and gene in data.index:
            legacy = common.loc[gene]
            original = data.loc[gene]
            if (legacy["SHIFRUT_HIT"], legacy["SHIFRUT_DIRECTION"]) != (
                    original["SHIFRUT_HIT"], original["SHIFRUT_DIRECTION"]):
                raise ValueError(f"Shifrut comparator provenance disagrees for a frozen gene: {gene}")
    return data


def functional_panel(screens, tables, shifrut, common, checks):
    chen = tables["chen"].set_index("gene")
    verify(set(chen.loc[chen["screen_hit"] & chen["sfari_score1"]].index) == set(PANEL),
           "all_14_discovery_hit_labels", "The fixed panel exactly matches all frozen Chen ASD-score-1 hits", checks)
    rows = []
    for gene in PANEL:
        row = {"gene": gene, "chen_measured": 1, "chen_hit": 1}
        src = chen.loc[gene]
        row.update(chen_author_lfc=src["author_lfc"], chen_signed_foxp3_bin_lfc=-src["author_lfc"],
                   chen_direction=src["regulator_direction"], chen_pos_fdr=src["chen_pos_fdr"],
                   chen_neg_fdr=src["chen_neg_fdr"], chen_source_row=src["source_row"],
                   chen_source=INPUT_TABLES["chen"], chen_screen_status="measured hit")
        row["treg_tconv_category"] = common.loc[gene, "category"] if gene in common.index else None
        for context, screen in screens.items():
            prefix = "treg" if context == "Treg" else "tconv"
            row[f"{prefix}_measured"] = int(gene in screen.index)
            row[f"{prefix}_source"] = f"GSE286470_Combined_gene_summary.xlsx / {GENE_SHEET}"
            if gene not in screen.index:
                row[f"{prefix}_screen_status"] = "not measured in targeted library"
                continue
            hit = screen.loc[gene]
            row.update({f"{prefix}_{new}": hit[old] for old, new in {
                "pos|lfc": "author_lfc", "signed_foxp3_bin_lfc": "signed_foxp3_bin_lfc",
                "direction": "effect_sign", "pos|fdr": "pos_fdr", "neg|fdr": "neg_fdr",
                "source_row": "source_row", "num": "guide_n",
                "FOXP3_REGULATOR_DIRECTION": "frozen_direction",
            }.items()})
            row[f"{prefix}_hit"] = int(hit["hit"])
            row[f"{prefix}_screen_status"] = "measured hit" if hit["hit"] else "measured non-hit"
        row["shifrut_measured"] = int(gene in shifrut.index)
        row["shifrut_source"] = "Shifrut 2018 integrated Table S4 / MAGeCK_Gene_GWscreen_rep1and2; frozen stage1d_r1"
        row["shifrut_screen_status"] = "not measured"
        if gene in shifrut.index:
            src = shifrut.loc[gene]
            row.update(shifrut_hit=int(src["SHIFRUT_HIT"]), shifrut_pos_fdr=src["pos_fdr"],
                       shifrut_neg_fdr=src["neg_fdr"], shifrut_pos_lfc=src["pos_lfc"],
                       shifrut_neg_lfc=src["neg_lfc"], shifrut_direction=src["SHIFRUT_DIRECTION"],
                       shifrut_screen_status="measured hit" if src["SHIFRUT_HIT"] else "measured non-hit")
        rows.append(row)
    return pd.DataFrame(rows)
