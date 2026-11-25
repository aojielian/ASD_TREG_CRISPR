#!/usr/bin/env python3
"""Audit GSE107011 and evaluate the fixed Monaco expression panel."""
import argparse
import csv
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


CANDIDATES = (
    "CREBBP", "CUL3", "DYRK1A", "EP300", "FOXP1", "KDM3B", "KDM6B",
    "KMT5B", "PHIP", "SIN3A", "SMARCA4", "TBL1XR1", "TLE3", "ZMYND8",
)
REFERENCES = ("FOXP3", "IL2RA", "CTLA4", "IKZF2")
PANEL = CANDIDATES + REFERENCES
SEED = 20261125
DRAWS = 20_000
TPM_THRESHOLDS = ((0.5, "ge0_5"), (1, "ge1"), (2, "ge2"))
CLASSES = {
    (0, 0): "NEITHER_LISTED", (1, 0): "TF_ONLY",
    (0, 1): "EPIFACTOR_ONLY", (1, 1): "BOTH_LISTED",
}
CONTINUOUS = ("healthy_cd4_mean_logcpm", "loeuf_v2_1_1", "atlas_median_log2_tpm_plus_1")
ENDPOINTS = ("breadth_ge1", "treg_rank_percentile")
MATCH_ATTRS = ("regulator_class", "library_member", *CONTINUOUS)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write(table, filename, index=False):
    table = pd.DataFrame(table)
    compression = {"method": "gzip", "mtime": 0} if filename.suffix == ".gz" else None
    table.to_csv(filename, sep="\t", index=index, na_rep="NA",
                 float_format="%.15g", compression=compression, lineterminator="\n")


def parse_soft(path):
    accession = None
    series = {}
    samples = {}
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            key, sep, value = line.rstrip("\r\n").partition(" = ")
            if not sep:
                continue
            if key == "^SAMPLE":
                require(value not in samples, f"SOFT duplicate {value}")
                accession = value
                samples[value] = {}
            elif key.startswith("!Series_"):
                series.setdefault(key, []).append(value)
            elif accession is not None and key.startswith("!Sample_"):
                samples[accession].setdefault(key.removeprefix("!Sample_"), []).append(value)
    require(series.get("!Series_geo_accession") == ["GSE107011"], "Expected GEO accession GSE107011")
    require("30726743" in series["!Series_pubmed_id"], "Expected Monaco PMID")
    require(set(series["!Series_sample_id"]) == set(samples), "Series and SOFT GSM lists disagree")
    return samples


def sample_crosswalk(samples, matrix_columns, groups):
    metadata = []
    for accession, record in samples.items():
        title, = record["title"]
        organism, = record["organism_ch1"]
        source, = record["source_name_ch1"]
        geo_accession, = record["geo_accession"]
        matrix_label, suffix, replicate = title.rpartition("_rep")
        require(accession == geo_accession and suffix and replicate.isdigit(), f"GSM or title mismatch: {accession}")
        require(matrix_label in matrix_columns and matrix_label in record["description"],
                f"Matrix key is absent from sample description: {accession}")
        donor, cell_type = matrix_label.split("_", 1)
        included = cell_type != "PBMC"
        expected_source = groups.loc[cell_type, "cell_type_label"] if included else "PBMCs"
        require(source == expected_source, f"Unexpected source label at {accession}")
        require(organism == "Homo sapiens", f"Unexpected organism at {accession}")
        chars = dict(value.split(": ", 1) for value in record["characteristics_ch1"])
        require(chars["cell type"] == source and chars["disease status"] == "healthy", f"SOFT metadata conflict: {accession}")
        require("Genome_build: GENCODE 26" in record["data_processing"], f"GENCODE provenance differs at {accession}")
        require(any("gene level" in value and "summarizing transcript" in value
                    for value in record["data_processing"]), "Expected gene-level author aggregation")
        metadata.append(dict(dataset="GSE107011", accession=accession, original_title=title,
                             author_description=matrix_label, author_source=source, matrix_label=matrix_label,
                             donor=donor, original_biological_replicate=int(replicate), species=organism,
                             cell_type=cell_type, lineage=groups.loc[cell_type, "lineage"] if included else "Mixed PBMC",
                             tissue="Peripheral blood, PBMC fraction", state="Healthy; no ex vivo stimulation reported",
                             included_primary=included,
                             reason="FACS-sorted S4 population" if included else "S13 mixed PBMC library"))
    meta = pd.DataFrame(metadata)
    require(meta["matrix_label"].is_unique and set(meta.matrix_label) == set(matrix_columns),
            "TPM sample columns and SOFT do not match one to one")
    require(meta.groupby("donor")["original_biological_replicate"].nunique().eq(1).all(),
            "Donor code and GEO suffix disagree")
    require(meta.groupby("original_biological_replicate")["donor"].nunique().eq(1).all(),
            "Donor and replicate suffix are not in one-to-one correspondence")
    return meta


def read_gtf(path):
    header = []
    records = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                header.append(line.strip())
                continue
            parts = line.split("\t", maxsplit=8)
            if parts[2] != "gene":
                continue
            fields = dict(re.findall(r'(\w+) "([^"]*)"', parts[8]))
            records.append((fields["gene_id"], fields["gene_name"], fields["gene_type"], parts[0]))
    require(any("GRCh38" in line and "version 26 (Ensembl 88)" in line for line in header),
            "GENCODE annotation header version mismatch")
    gtf = pd.DataFrame(records, columns=["ensembl_gene_id", "gencode_symbol", "gene_type", "chromosome"])
    require(gtf.ensembl_gene_id.is_unique, "Nonunique GENCODE gene feature IDs")
    return gtf


def annotate(tpm, gtf, hgnc_path):
    x = tpm.index.to_frame(index=False).merge(gtf, how="left", on="ensembl_gene_id", validate="one_to_one")
    hg = pd.read_csv(hgnc_path, sep="\t", dtype="string", usecols=["symbol", "hgnc_id", "ensembl_gene_id", "status"])
    hg = hg.loc[hg.status.eq("Approved")]
    require(hg.symbol.is_unique and hg.hgnc_id.is_unique, "Frozen HGNC approved IDs/symbols duplicate")
    hg_id = hg.loc[hg.ensembl_gene_id.notna()]
    hg_unique = hg_id.loc[~hg_id.ensembl_gene_id.duplicated(keep=False)].set_index("ensembl_gene_id")
    hg_by_symbol = hg.set_index("symbol")
    x["ensembl_stable_id"] = x.ensembl_gene_id.str.extract(r"^(ENSG[0-9]+)\.[0-9]+(?:_PAR_Y)?$")[0]
    x["approved_symbol"] = x.ensembl_stable_id.map(hg_unique.symbol)
    x["hgnc_id"] = x.ensembl_stable_id.map(hg_unique.hgnc_id)
    x["mapping_route"] = np.where(x.approved_symbol.notna(), "unique_frozen_HGNC_stable_id", "GENCODE_only")
    no_stable_frozen_id = set(hg_by_symbol.index[hg_by_symbol.ensembl_gene_id.isna()])
    fallback = (x.approved_symbol.isna() & x.gencode_symbol.isin(no_stable_frozen_id) &
                ~x.ensembl_stable_id.isin(hg_id.ensembl_gene_id))
    x.loc[fallback, "approved_symbol"] = x.loc[fallback, "gencode_symbol"]
    x.loc[fallback, "hgnc_id"] = x.loc[fallback, "gencode_symbol"].map(hg_by_symbol.hgnc_id)
    x.loc[fallback, "mapping_route"] = "frozen_HGNC_exact_symbol_no_stable_id"
    past_symbol = x.gencode_symbol.isin(hg_by_symbol.index)
    x.loc[x.approved_symbol.isna() & past_symbol, "mapping_route"] = "frozen_HGNC_link_disagreement_or_nonunique"
    x["gene"] = x.approved_symbol.fillna(x.gencode_symbol)
    y_par = x.ensembl_gene_id.str.endswith("_PAR_Y") & x.gene_type.notna()
    require(x.loc[y_par, "ensembl_gene_id"].str.removesuffix("_PAR_Y").isin(
        gtf.loc[gtf.chromosome.eq("chrX"), "ensembl_gene_id"]).all(), "Unpaired PAR_Y gene feature")
    x["included_protein_coding"] = x.gene_type.eq("protein_coding") & ~y_par
    duplicates = x.loc[x.included_protein_coding, "gene"].duplicated(keep=False)
    duplicated_rows = duplicates.index[duplicates]
    finite_rows = pd.Series(np.isfinite(tpm.to_numpy()).all(axis=1), x.index)
    x["missing_reason"] = ""
    x.loc[x.gene_type.isna(), "missing_reason"] = "No exact versioned GENCODE v26 gene feature"
    x.loc[x.gene_type.notna() & ~x.gene_type.eq("protein_coding"), "missing_reason"] = "GENCODE biotype not protein_coding"
    x.loc[y_par, "missing_reason"] = "Y pseudoautosomal duplicate of X feature"
    x.loc[duplicated_rows, "missing_reason"] = "Nonunique atlas symbol across multiple protein-coding loci"
    x.loc[x.included_protein_coding & ~finite_rows, "missing_reason"] = "Non-finite expression value"
    x["included_protein_coding"] &= x.missing_reason.eq("")
    require(x.loc[x.included_protein_coding, "gene"].is_unique, "Atlas symbol nonunique after annotation")
    return x


def atlas_summary(tpm, meta, order):
    selected = meta.loc[meta.included_primary]
    indexed = tpm.loc[:, selected.matrix_label].copy()
    indexed.columns = pd.MultiIndex.from_frame(selected[["donor", "cell_type"]])
    donor_type = indexed.T.groupby(level=["donor", "cell_type"], sort=False).median().T
    type_medians = donor_type.T.groupby(level="cell_type").median().T.loc[:, order]
    return donor_type, type_medians


def metrics(matrix):
    result = pd.DataFrame(index=matrix.index)
    n = matrix.notna().sum(axis=1)
    result["n_cell_types"] = n
    result["median_type_tpm"] = matrix.median(axis=1)
    for threshold, name in TPM_THRESHOLDS:
        result[f"detected_{name}"] = matrix.ge(threshold).sum(axis=1).where(n.gt(0))
        result[f"breadth_{name}"] = result[f"detected_{name}"] / n.where(n.gt(0))
    max_tpm = matrix.max(axis=1)
    denom = (n - 1).where(n.gt(1) & max_tpm.gt(0))
    result["tau"] = (1 - matrix.div(max_tpm.where(max_tpm.gt(0)), axis=0)).sum(axis=1, min_count=1) / denom
    result["treg_tpm"] = matrix["Treg"]
    result["treg_rank_desc"] = matrix.rank(axis=1, ascending=False, method="average")["Treg"]
    result["treg_rank_percentile"] = (n - result["treg_rank_desc"]) / (n - 1).where(n.gt(1))
    return result


def paired_contrast(donor_type, primary_types, all_types):
    results = []
    for donor in sorted(donor_type.columns.get_level_values("donor").unique()):
        matrix = donor_type[donor]
        for name, types in (("primary_balanced_six", primary_types), ("all_seven_complete_donors", all_types)):
            if not set(types).issubset(matrix.columns):
                continue
            kept = matrix.reindex(PANEL)
            subset_mean = kept.loc[:, types].mean(axis=1)
            result = pd.DataFrame({
                "gene": kept.index, "donor": donor, "contrast_set": name,
                "treg_tpm": kept["Treg"].to_numpy(), "tconv_equal_type_tpm": subset_mean.to_numpy(),
                "n_tconv_types": len(types), "tconv_types": "|".join(types),
                "log2_ratio": np.log2((kept["Treg"] + 1) / (subset_mean + 1)).to_numpy(),
            })
            results.append(result)
    details = pd.concat(results, ignore_index=True)
    details.insert(1, "role", details.gene.map(role))
    return details


def role(gene):
    return "candidate" if gene in CANDIDATES else "reference"


def heatmap_long(matrix, groups, meta=None):
    selected = matrix.reindex(PANEL)
    data = np.log2(selected + 1)
    z = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    records = []
    for order, kind in enumerate(selected.columns, 1):
        if meta is None:
            label, lineage = kind, kind
            n = np.nan
            type_n = int(groups.lineage.eq(kind).sum())
        else:
            label, lineage = groups.loc[kind, ["cell_type_label", "lineage"]]
            n = int(meta.loc[meta.cell_type.eq(kind), "donor"].nunique())
            type_n = 1
        for gene in selected.index:
            records.append(dict(gene=gene, role=role(gene), cell_type=kind, lineage=lineage,
                                cell_type_label=label, cell_type_order=order, donor_n=n,
                                n_subtypes=type_n, median_tpm=selected.loc[gene, kind],
                                log2_tpm_plus_1=data.loc[gene, kind], row_z=z.loc[gene, kind]))
    return pd.DataFrame(records)


def boolean(raw):
    values = raw.astype("string").str.strip().str.lower()
    encoded = values.map({"true": 1, "false": 0, "1": 1, "0": 0, "1.0": 1, "0.0": 0,
                          "yes": 1, "no": 0})
    require(encoded.notna().eq(raw.notna()).all(), "Unexpected frozen annotation Boolean")
    return encoded


def frozen_background(data_path, eligible, atlas_stats):
    source = pd.read_csv(data_path / "data/r7/inputs/R7_GENE_BACKGROUND.tsv.gz", sep="\t",
                         dtype="string", usecols=[
        "approved_symbol", "hgnc_id", "screen", "author_hit", "measured", "sfari_score1", "library_member",
        "class_label", "CLASS_ANNOTATION_UNCERTAIN", "annotation_status",
        "TF_LISTED", "EPIFACTORS_GENES_LISTED", "HEALTHY_CD4_MEAN_LOGCPM", "LOEUF_v2_1_1",
    ])
    cases = pd.read_csv(data_path / "data/r7/mechanism/ALLELIC_MECHANISM_MASTER.tsv", sep="\t",
                        dtype="string", usecols=["approved_symbol", "ensembl_id"])
    require(cases.approved_symbol.is_unique and len(cases) == 245, "Frozen 245-gene set not reproduced")
    chen = source.loc[source.screen.eq("CHEN") & boolean(source.measured).eq(1)].copy()
    require(chen.approved_symbol.is_unique, "Frozen Chen gene keys duplicate")
    chen["asd"] = chen.approved_symbol.isin(cases.approved_symbol)
    require(chen.asd.eq(boolean(chen.sfari_score1).eq(1)).all(), "ASD label parity failed")
    chen["chen_hit"] = boolean(chen.author_hit).eq(1)
    require(len(chen) == 18966 and chen.asd.sum() == 242 and chen.chen_hit.sum() == 257,
            "Frozen Chen baseline labels not reproduced")
    require(set(chen.loc[chen.asd & chen.chen_hit, "approved_symbol"]) == set(CANDIDATES),
            "Frozen Chen case hits differ from fixed panel")
    class_expected = pd.Series([CLASSES.get(pair) for pair in zip(boolean(chen.TF_LISTED),
                              boolean(chen.EPIFACTORS_GENES_LISTED))], index=chen.index, dtype="string")
    class_usable = (boolean(chen.CLASS_ANNOTATION_UNCERTAIN).eq(0) &
                    chen.annotation_status.eq("CLASS_ANNOTATION_USABLE") &
                    chen.class_label.isin(CLASSES.values()))
    require(chen.loc[class_usable, "class_label"].eq(class_expected.loc[class_usable]).all(),
            "Frozen gene classes disagree with TF/EpiFactors flags")
    chen["regulator_class"] = chen.class_label.where(class_usable)
    chen["library_member"] = boolean(chen.library_member)
    chen["healthy_cd4_mean_logcpm"] = pd.to_numeric(chen.HEALTHY_CD4_MEAN_LOGCPM, errors="raise")
    chen["loeuf_v2_1_1"] = pd.to_numeric(chen.LOEUF_v2_1_1, errors="raise")
    atlas = eligible.loc[eligible.approved_symbol.notna(), ["gene", "approved_symbol", "hgnc_id", "ensembl_stable_id",
                                                         "ensembl_gene_id", "mapping_route"]]
    atlas = atlas.join(atlas_stats.loc[:, ["median_type_tpm", *ENDPOINTS]], on="gene", validate="one_to_one")
    chen = chen[["approved_symbol", "hgnc_id", "asd", "chen_hit",
                 "regulator_class", "library_member", *CONTINUOUS[:2]]].merge(
        atlas, on="approved_symbol", how="left", validate="one_to_one", suffixes=("_r9", "_atlas"))
    chen["atlas_mapping"] = np.where(chen.gene.notna(), "included_unique_link", "no_unique_atlas_protein_coding_link")
    conflicts = chen.hgnc_id_r9.notna() & chen.hgnc_id_atlas.notna() & chen.hgnc_id_r9.ne(chen.hgnc_id_atlas)
    require(not conflicts.any(), "Frozen R9 against HGNC ID disagreement")
    stable = chen[["approved_symbol", "ensembl_stable_id"]].merge(cases, on="approved_symbol", validate="one_to_one")
    case_diff = stable.ensembl_stable_id.notna() & stable.ensembl_id.notna() & stable.ensembl_stable_id.ne(stable.ensembl_id)
    require(not case_diff.any(), "Frozen ASD against GENCODE locus disagreement")
    chen["atlas_median_log2_tpm_plus_1"] = np.log2(chen["median_type_tpm"] + 1)
    return chen


def balance(phase, subset, weights, scenario):
    flags = pd.get_dummies(subset.regulator_class, dtype=float).reindex(columns=list(CLASSES.values()), fill_value=0)
    numeric = pd.concat([subset.loc[:, ["library_member", *CONTINUOUS]], flags], axis=1)
    case = subset.asd.to_numpy()
    rows = []
    for name in numeric.columns:
        left = numeric.loc[case, name].to_numpy(float)
        right = numeric.loc[~case, name].to_numpy(float)
        wl, wr = weights[case], weights[~case]
        mean_a = np.average(left, weights=wl) if len(left) else np.nan
        mean_c = np.average(right, weights=wr) if len(right) else np.nan
        variance = (left.var(ddof=1) + right.var(ddof=1))/2 if min(len(left), len(right)) > 1 else np.nan
        smd = (mean_a-mean_c) / np.sqrt(variance) if variance > 0 else (0 if mean_a == mean_c else np.nan)
        rows.append(dict(contrast=scenario, phase=phase, covariate=name, n_asd=len(left), n_non_asd=len(right),
                         mean_asd=mean_a, mean_non_asd=mean_c, smd=smd, imbalance_flag=abs(smd) > 0.10))
    return rows


def weighted_median(values, weights):
    order = np.argsort(values, kind="stable")
    ordered = values[order]
    cdf = np.cumsum(weights[order]) / weights.sum()
    at = np.searchsorted(cdf, 0.5)
    if at + 1 < len(cdf) and np.isclose(cdf[at], 0.5, atol=1e-12, rtol=0):
        return (ordered[at] + ordered[at + 1])/2
    return ordered[at]


def matched_analysis(region, contrast, frozen_n, rng):
    region = region.copy()
    atlas = region.gene.notna()
    region["complete_case"] = atlas & region.loc[:, MATCH_ATTRS].notna().all(axis=1)
    complete = region.loc[region.complete_case].copy()
    bits = ["regulator_class", "library_member"]
    for name in CONTINUOUS:
        target = name + "_tertile"
        ranks = complete[name].rank(method="average")
        complete[target] = np.minimum(3, np.ceil(3 * ranks / len(complete))).astype(int)
        require(complete.groupby(name)[target].nunique().le(1).all(), f"Split tertile ties for {name}")
        bits.append(target)
    complete["stratum"] = complete.loc[:, bits].astype(str).agg("|".join, axis=1)
    counts = complete.groupby("stratum", sort=True).asd.agg(["sum", "size"]).rename(columns={"sum": "asd_n", "size": "n"})
    counts["control_n"] = counts["n"] - counts["asd_n"]
    supported_keys = counts.index[counts.asd_n.gt(0) & counts.control_n.gt(0)]
    supported = complete.loc[complete.stratum.isin(supported_keys)].copy()
    supported["weight"] = supported.stratum.map(counts.asd_n / counts.control_n.replace(0, np.nan))
    supported.loc[supported.asd, "weight"] = 1.0
    region["stratum"] = complete.stratum
    region = region.join(complete.loc[:, bits[2:]])
    region["common_support"] = region.index.isin(supported.index)
    region["weight"] = supported.weight
    region["contrast"] = contrast
    region["exclusion_reason"] = [
        "|".join([reason for reason, missing in zip(("atlas", *MATCH_ATTRS), flags) if missing])
        for flags in zip(~atlas, *[region[name].isna() for name in MATCH_ATTRS])
    ]
    region.loc[region.complete_case & ~region.common_support, "exclusion_reason"] = "one group absent in stratum"

    real_a = supported.loc[supported.asd]
    real_c = supported.loc[~supported.asd]
    n_asd = len(real_a)
    require(np.isclose(real_c.weight.sum(), n_asd), "Matched reference weights do not sum to case count")
    diagnostics = balance("complete_before_support", complete, np.ones(len(complete)), contrast)
    diagnostics.extend(balance("common_support_weighted", supported, supported.weight.to_numpy(), contrast))
    counts["contrast"] = contrast
    strata = list(supported.groupby("stratum", sort=True))
    full_n_asd = int(region.asd.sum())
    full_n_use = int((region.asd & atlas).sum())
    complete_n_asd = int(complete.asd.sum())

    estimates = []
    draw_out = []
    if n_asd:
        null = np.empty((DRAWS, n_asd, len(ENDPOINTS)), dtype=float)
        boot = np.empty_like(null)
        start = 0
        for name, group in strata:
            values = group.loc[:, ENDPOINTS].to_numpy(float)
            cases = values[group.asd.to_numpy()]
            n = len(cases)
            if n == 1:
                labels = rng.integers(len(group), size=(DRAWS, n))
            else:
                labels = np.fromiter(
                    (item for _ in range(DRAWS) for item in
                     rng.choice(len(group), size=n, replace=False, shuffle=False)),
                    dtype=np.int64, count=DRAWS*n).reshape(DRAWS, n)
            null[:, start:start+n, :] = values[labels]
            boot[:, start:start+n, :] = cases[rng.integers(n, size=(DRAWS, n))]
            start += n
        require(start == n_asd, "Randomized draw size differed from observed common-support case count")
        perm_stats = np.median(null, axis=1)
        boot_stats = np.median(boot, axis=1)

    for i, endpoint in enumerate(ENDPOINTS):
        measurement = "proportion of 29 types at median TPM >= 1" if i == 0 else "Treg midrank percentile (1 highest)"
        obs_group = real_a[endpoint]
        control_group = real_c[endpoint]
        observed = float(obs_group.median())
        reference = weighted_median(control_group.to_numpy(), real_c.weight.to_numpy()) if n_asd else np.nan
        info = {name: group[endpoint].nunique() > 1 for name, group in strata}
        informative_n = int(sum(counts.loc[name, "asd_n"] for name in info if info[name]))
        coverage80 = n_asd >= 0.8 * frozen_n
        coverage70 = informative_n >= 0.7 * complete_n_asd if complete_n_asd else False
        counts[f"{endpoint}_varying_common_support"] = [
            info.get(key, False) for key in counts.index]
        status = "ESTIMABLE_EXPLORATORY"
        null_center = lo = hi = effect = ci_lo = ci_hi = p = mcse = np.nan
        n_draw = DRAWS if n_asd else 0
        if n_asd:
            values = perm_stats[:, i]
            null_center = float(np.median(values))
            lo, hi = np.quantile(values, [0.025, 0.975])
            effect = observed - null_center
            ci_lo, ci_hi = np.quantile(boot_stats[:, i] - null_center, [0.025, 0.975])
            draw_out.append(pd.DataFrame({"contrast": contrast, "endpoint": endpoint,
                                          "draw": np.arange(1, DRAWS+1), "null_median": values}))
            if np.ptp(values) == 0:
                status = "DEGENERATE_NULL"
            elif not (coverage80 or coverage70):
                status = "INSUFFICIENT_MATCHED_COVERAGE"
            else:
                extreme = np.abs(values - null_center) >= abs(effect) - 1e-14
                p = (int(extreme.sum()) + 1)/(DRAWS + 1)
                mcse = np.sqrt(p * (1-p)/(DRAWS+1))
        else:
            status = "NO_COMPARATOR_SUPPORT"
        estimates.append(dict(contrast=contrast, endpoint=endpoint, unit=measurement,
                              frozen_asd_n=frozen_n, chen_measured_asd_n=full_n_asd,
                              chen_measured_non_asd_n=int((~region.asd).sum()), atlas_asd_n=full_n_use,
                              atlas_non_asd_n=int((~region.asd & atlas).sum()),
                              complete_case_asd_n=complete_n_asd, complete_case_non_asd_n=int((~complete.asd).sum()),
                              supported_asd_n=n_asd, supported_non_asd_n=len(real_c), supported_strata_n=len(strata),
                              unsupported_case_fraction=1-n_asd/frozen_n,
                              informative_asd_n=informative_n, coverage_80_pct_of_frozen=coverage80,
                              informative_70_pct_of_complete=coverage70,
                              unmatched_observed_asd_median=region.loc[region.asd, endpoint].median(),
                              unmatched_non_asd_median=region.loc[~region.asd, endpoint].median(),
                              observed_asd_median=observed, weighted_non_asd_median=reference,
                              null_median=null_center, effect_vs_null_median=effect,
                              null_lo=lo, null_hi=hi, ci_lo=ci_lo, ci_hi=ci_hi,
                              null_draws=n_draw, p=p, mcse=mcse, status=status,
                              residual_smd_flag=any(row["imbalance_flag"] for row in diagnostics
                                                    if row["phase"] == "common_support_weighted"),
                              method="Two-sided within-stratum gene-label permutation; bootstrap CI conditional on fixed null reference"))
    p = np.array([item["p"] for item in estimates])
    p_for_bh = np.where(np.isfinite(p), p, 1)
    rank_order = np.argsort(p_for_bh)
    adjusted = np.ones(len(p))
    adjusted[rank_order] = np.minimum.accumulate((p_for_bh[rank_order]*len(p) /
                                                   np.arange(1, len(p)+1))[::-1])[::-1]
    for item, q in zip(estimates, adjusted):
        item["q_bh"] = min(float(q), 1.0) if np.isfinite(item["p"]) else np.nan
    return estimates, diagnostics, draw_out, counts.reset_index(), region


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="Extracted processed-data deposit")
    parser.add_argument("--stage", type=Path, required=True, help="Replay output directory")
    parser.add_argument("--inputs", type=Path, required=True, help="Verified source file directory")
    parser.add_argument("--frozen-mappings", type=Path, help="Use the deposited sample and gene crosswalks")
    args = parser.parse_args()
    stage, raw = args.stage, args.inputs
    audit, results = stage / "01_audit", stage / "03_results"
    files = {
        "tpm": raw / "GSE107011/GSE107011_Processed_data_TPM.txt.gz",
        "soft": raw / "GSE107011/GSE107011_family.soft.gz",
        "gtf": raw / "annotation/gencode.v26.annotation.gtf.gz",
        "hgnc": raw / "annotation/hgnc_complete_set_snapshot.txt",
    }
    groups = pd.read_csv(stage / "07_config/monaco_cell_types.tsv", sep="\t", index_col="cell_type")
    require(groups.index.is_unique and len(groups) == 29 and groups.lineage.nunique() == 10,
            "Unexpected locked Monaco cell type dictionary")
    primary_types = list(groups.index[groups.conventional_primary.eq(1)])
    all_types = list(groups.index[groups.conventional_sensitivity.eq(1)])
    require(len(primary_types) == 6 and len(all_types) == 7, "Unexpected locked Tconv types")
    with gzip.open(files["tpm"], "rt", newline="") as stream:
        header = next(csv.reader(stream, delimiter="\t"))
    require(header[0] == "" and len(set(header[1:])) == len(header)-1, "TPM column header invalid")
    tpm = pd.read_csv(files["tpm"], sep="\t", index_col=0)
    tpm.index.name = "ensembl_gene_id"
    require(tpm.index.is_unique and list(tpm.columns) == header[1:], "TPM index/column IDs not unique")
    require(all(pd.api.types.is_numeric_dtype(col) for col in tpm.dtypes), "Nonnumeric TPM column")
    numeric = tpm.to_numpy()
    require(not (numeric < 0).any(), "Negative TPM")
    if args.frozen_mappings:
        files["soft"] = args.frozen_mappings / "R10_MONACO_SAMPLE_CROSSWALK.tsv"
        meta = pd.read_csv(files["soft"], sep="\t", low_memory=False)
        require(set(meta.matrix_label) == set(tpm.columns) and meta.matrix_label.is_unique,
                "Frozen Monaco sample crosswalk differs from the TPM matrix")
        require(meta.included_primary.dtype == bool, "Expected boolean Monaco inclusion flags")
    else:
        meta = sample_crosswalk(parse_soft(files["soft"]), list(tpm.columns), groups)
    sorted_meta = meta.loc[meta.included_primary]
    require(len(meta) == 127 and len(sorted_meta) == 114, "Observed libraries differ from GEO")
    require(sorted_meta.donor.nunique() == 4 and meta.loc[~meta.included_primary, "donor"].nunique() == 13,
            "Observed donors differ from GEO")
    by_type = sorted_meta.groupby("cell_type").donor.nunique().reindex(groups.index)
    require(by_type.drop("CD4_TE").eq(4).all() and by_type["CD4_TE"] == 2,
            "Expected 2 terminal effector missing donor×type combinations")
    missing = pd.MultiIndex.from_product([sorted_meta.donor.unique(), groups.index], names=["donor", "cell_type"]).difference(
        pd.MultiIndex.from_frame(sorted_meta[["donor", "cell_type"]])).to_frame(index=False)
    missing["reason"] = "Unavailable in author atlas; article reports insufficient RNA for these CD4 TE sorts"
    write(meta, audit / "R10_MONACO_SAMPLE_CROSSWALK.tsv")
    write(missing, audit / "R10_MONACO_MISSING_DONOR_TYPE.tsv")
    if args.frozen_mappings:
        files["gtf"] = args.frozen_mappings / "R10_MONACO_GENE_ID_CROSSWALK.tsv.gz"
        x = pd.read_csv(files["gtf"], sep="\t", dtype="string")
        require(x.ensembl_gene_id.tolist() == tpm.index.tolist(),
                "Frozen Monaco crosswalk differs from the TPM locus keys")
        require(x.included_protein_coding.isin(("True", "False")).all(),
                "Expected boolean Monaco locus inclusion flags")
        x["included_protein_coding"] = x.included_protein_coding.eq("True")
        x["missing_reason"] = x.missing_reason.fillna("")
        n_annotation_rows = len(x)
    else:
        gtf = read_gtf(files["gtf"])
        x = annotate(tpm, gtf, files["hgnc"])
        n_annotation_rows = len(gtf)
    write(x, audit / "R10_MONACO_GENE_ID_CROSSWALK.tsv.gz")
    eligible = x.loc[x.included_protein_coding]
    atlas = tpm.loc[eligible.ensembl_gene_id].copy()
    atlas.index = pd.Index(eligible.gene, name="gene")
    donor_type, median = atlas_summary(atlas, meta, groups.index)
    med_out = median.join(eligible.set_index("gene")[["ensembl_gene_id"]])
    write(med_out.loc[:, ["ensembl_gene_id", *groups.index]], results / "R10_MONACO_CELLTYPE_MEDIAN_TPM.tsv.gz", index=True)
    lineages = {}
    for lineage, group in groups.groupby("lineage", sort=False):
        lineages[lineage] = median.loc[:, group.index].mean(axis=1)
    lineage_matrix = pd.DataFrame(lineages, index=median.index)
    stat = metrics(median)
    lineage_stat = metrics(lineage_matrix)
    stat["tconv_equal_type_tpm"] = median.loc[:, primary_types].mean(axis=1)
    stat["treg_tconv_log2_ratio"] = np.log2((stat["treg_tpm"] + 1)/(stat["tconv_equal_type_tpm"] + 1))
    stat["tconv_all_seven_equal_type_tpm"] = median.loc[:, all_types].mean(axis=1)
    stat["treg_tconv_all_seven_log2_ratio"] = np.log2((stat["treg_tpm"] + 1)/(stat["tconv_all_seven_equal_type_tpm"] + 1))
    paired = paired_contrast(donor_type, primary_types, all_types)
    main_pairs = paired.loc[paired.contrast_set.eq("primary_balanced_six")]
    panel = stat.reindex(PANEL).join(eligible.set_index("gene")[["ensembl_gene_id", "gencode_symbol",
                                                              "hgnc_id", "mapping_route"]])
    paired_summary = main_pairs.groupby("gene").log2_ratio.agg(["count", "median", "min", "max"]).rename(columns={
        "count": "paired_donor_n", "median": "paired_median_log2_ratio", "min": "paired_min_log2_ratio",
        "max": "paired_max_log2_ratio"})
    panel = panel.join(paired_summary)
    panel.insert(0, "role", panel.index.map(role))
    panel["units"] = "TPM"
    panel["missing_reason"] = np.where(panel.ensembl_gene_id.notna(), "", "No uniquely evaluable protein-coding atlas locus")
    write(panel, results / "R10_HUMAN_EXPRESSION_PANEL.tsv", index=True)
    write(paired, results / "R10_EXPRESSION_DONOR_CONTRASTS.tsv")
    write(heatmap_long(median, groups, sorted_meta), results / "R10_MONACO_EXPRESSION_HEATMAP_SOURCE.tsv")
    write(heatmap_long(lineage_matrix, groups), results / "R10_MONACO_LINEAGE_HEATMAP_SOURCE.tsv")
    sensitivity = []
    for aggregation, m in (("sorted_types_equal", stat), ("ten_lineages_equal", lineage_stat)):
        for threshold, suffix in TPM_THRESHOLDS:
            rows = m.reindex(PANEL).loc[:, [f"detected_{suffix}", f"breadth_{suffix}", "n_cell_types",
                                          "tau", "treg_rank_percentile"]].copy()
            rows.columns = ["n_detected", "breadth", "n_evaluable", "tau", "treg_rank_percentile"]
            rows.insert(0, "role", rows.index.map(role))
            rows["dataset"], rows["species"], rows["expression_unit"] = "GSE107011", "Homo sapiens", "TPM"
            rows["aggregation"], rows["threshold"] = aggregation, threshold
            rows["missing_reason"] = np.where(rows.n_evaluable.notna(), "", "No uniquely evaluable protein-coding atlas locus")
            sensitivity.append(rows)
    write(pd.concat(sensitivity), results / "R10_HUMAN_EXPRESSION_SENSITIVITY.tsv", index=True)

    chen = frozen_background(args.data, eligible, stat)
    matching = [(chen.loc[chen.chen_hit], "fixed_14_vs_non_asd_chen_hits", 14),
                (chen, "full_asd_vs_non_asd_chen_measured", 245)]
    analyses = [matched_analysis(frame, label, n, np.random.default_rng(np.random.SeedSequence([SEED, i])))
                for i, (frame, label, n) in enumerate(matching)]
    combined = list(zip(*analyses))
    write([row for item in combined[0] for row in item], results / "R10_EXPRESSION_MATCHED_NULL.tsv")
    write([row for item in combined[1] for row in item], results / "R10_EXPRESSION_MATCH_BALANCE.tsv")
    draws = [table for item in combined[2] for table in item]
    if draws:
        write(pd.concat(draws, ignore_index=True), results / "R10_MONACO_MATCH_NULL_DRAWS.tsv.gz")
    write(pd.concat(combined[3], ignore_index=True), results / "R10_MONACO_MATCH_STRATA.tsv")
    write(pd.concat(combined[4], ignore_index=True), results / "R10_MONACO_MATCH_MEMBERS.tsv.gz")

    outcomes = pd.DataFrame([row for item in combined[0] for row in item])
    def check(label, passed, detail):
        return {"check": label, "status": "PASS" if passed else "FAIL", "detail": detail}
    qa = [
        check("TPM_shape_and_unique_keys", len(tpm) == 58311 and len(meta) == 127,
              f"{len(tpm)} gene-level source rows; {len(meta)} unique GEO GSM/TPM columns"),
        check("author_sorting_and_PBMC", len(sorted_meta) == 114 and len(missing) == 2,
              f"{len(sorted_meta)} sorted RNA libraries; 29 types; four donors; two absent CD4 TE; 13 bulk PBMC excluded"),
        check("annotation_finite_positive_range", np.isfinite(atlas.to_numpy()).all() and np.min(numeric) >= 0,
              f"PC atlas rows {len(eligible)}; exact GTF features {int(x.gene_type.notna().sum())}; source TPM {np.nanmin(numeric):g}–{np.nanmax(numeric):g}; nonfinite source cells {int((~np.isfinite(numeric)).sum())}"),
        check("frozen_panel_and_reference", panel.ensembl_gene_id.notna().all() and len(panel) == 18,
              f"{int(panel.ensembl_gene_id.notna().sum())}/18 uniquely mapped panel/reference genes; all 14 equal frozen Chen hits"),
        check("paired_donor_scope", len(main_pairs) == 18*4 and len(paired) == 18*(4+2),
              "4 donors with the fixed six-type primary contrast; 2 with all seven conventional CD4 types"),
        check("breadth_threshold_order", stat.breadth_ge0_5.ge(stat.breadth_ge1).all() and
              stat.breadth_ge1.ge(stat.breadth_ge2).all(), "All atlas rows obey 0.5/1/2 TPM detection monotonicity"),
        check("breadth_tau_rank_bounds", stat.breadth_ge1.between(0, 1).all() and
              stat.tau.dropna().between(0, 1).all() and stat.treg_rank_percentile.between(0, 1).all(),
              "Breadth, linear Tau and Treg type rank percentile within [0,1]"),
        check("matching_support_and_p_bounds", outcomes.supported_asd_n.le(outcomes.complete_case_asd_n).all() and
              outcomes.p.dropna().between(0, 1).all() and outcomes.q_bh.dropna().between(0, 1).all(),
              "Complete-case margin bound and legal P and BH q; unsupported or degenerate results explicitly coded"),
    ]
    write(qa, audit / "R10_MONACO_NUMERIC_QA.tsv")
    schemas = [
        dict(dataset="GSE107011", input=files["tpm"].name, rows=len(tpm), sample_columns=len(tpm.columns),
             row_key="versioned Ensembl gene ID; unmapped and spike rows in gene crosswalk",
             sample_key="exact SOFT !Sample_description and title without replicate suffix", value="author gene-level TPM",
             min_value=np.nanmin(numeric), max_value=np.nanmax(numeric), nonfinite=int((~np.isfinite(numeric)).sum()),
             output_scope=f"{len(eligible)} unique protein-coding atlas loci; 29 donor-median cell types"),
        dict(dataset="GSE107011", input=files["soft"].name, rows=len(meta), sample_columns=np.nan,
             row_key="unique GSM", sample_key="unique author label", value="sample metadata",
             output_scope="114 sorted, 13 mixed PBMC; 4 unique sorted donors and 13 PBMC donors (overlapping cohorts)"),
        dict(dataset="Frozen locus crosswalk" if args.frozen_mappings else "GENCODE human v26",
             input=files["gtf"].name, rows=n_annotation_rows, sample_columns=np.nan,
             row_key="unique versioned Ensembl gene feature", value="GRCh38 Ensembl 88",
             output_scope=f"{int(x.gene_type.notna().sum())} exact TPM row matches"),
    ]
    write(schemas, audit / "R10_MONACO_SCHEMA.tsv")
    print(pd.DataFrame(qa).to_string(index=False))
    print(panel[["detected_ge1", "n_cell_types", "treg_tpm", "treg_rank_desc",
                 "treg_tconv_log2_ratio"]].to_string())
    print(outcomes[["contrast", "endpoint", "supported_asd_n", "supported_non_asd_n", "observed_asd_median",
                    "null_median", "p", "q_bh", "status", "residual_smd_flag"]].to_string(index=False))
    require(all(row["status"] == "PASS" for row in qa), "Monaco QA reported an integrity error")


if __name__ == "__main__":
    main()
