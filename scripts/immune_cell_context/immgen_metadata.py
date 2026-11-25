"""Audit ImmGen library identities, sample states, and MGI orthology."""

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

SERIES = "GSE109125"
COUNTS = "GSE109125_Gene_count_table_GENCODE_vM25.csv.gz"
SORT_SOURCE = "https://www.immgen.org/ImmGenATAC1219Sorts.S1A.pdf"
SOFT = "GSE109125_family.soft.gz"
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
LINEAGES = (
    "Treg",
    "CD4 T",
    "CD8 T",
    "Unconventional T",
    "B",
    "NK",
    "ILC",
    "Dendritic",
    "Monocyte",
    "Macrophage",
    "Granulocyte",
    "Mast",
    "Progenitor",
    "Stromal",
)
EXTRA_B = {"B.GC.CB.Sp", "B.GC.CC.Sp", "B.PB.Sp"}
MAIN_TREG = ("Treg.4.25hi.Sp", "Treg.4.FP3+.Nrplo.Co")
TCONV = "T.4.Nve.Sp"

# Only remedy misleading/incomplete author descriptions or specify key gates.
DISPLAY_LABELS = {
    "Treg.4.25hi.Sp": "Splenic CD25hi Treg",
    "Treg.4.FP3+.Nrplo.Co": "Colonic Foxp3+ Nrp1lo Treg",
    "T.4.Nve.Sp": "Splenic naive CD4",
    "T.4.Nve.Fem.Sp": "Splenic naive CD4 (female)",
    "T.4.19-8-TCRb+CD4+.Sp": "Splenic bulk CD4",
    "T.8.Nve.Sp": "Splenic naive CD8",
    "DC.4+.Sp": "Splenic CD4+ DC",
    "DC.8+.Sp": "Splenic CD8+ DC",
    "DC.pDC.Sp": "Splenic plasmacytoid DC",
    "MF.PC": "Peritoneal F4/80+ ICAM2+ macrophage",
    "MF.102+480+.PC": "Peritoneal CD102+ F4/80+ MHCIIlo macrophage",
    "MF.226+II+480lo.PC": "Peritoneal CD226+ MHCII+ F4/80lo macrophage",
    "Mo.6C+II-.Bl": "Blood Ly6C-high monocyte",
    "Mo.6C-II-.Bl": "Blood Ly6C-low monocyte",
}

SOURCE_LINEAGE = {
    "B cells": "B",
    "Basophils": "Granulocyte",
    "Dendritic cells": "Dendritic",
    "Eosinophil": "Granulocyte",
    "Granulocytes": "Granulocyte",
    "Macrophages": "Macrophage",
    "Mast Cells": "Mast",
    "Monocytes": "Monocyte",
    "Stem Cells": "Progenitor",
    "Stromal Cells": "Stromal",
    "gd T cells": "Unconventional T",
}

DESCRIPTION_FIXES = {"Spleinc": "Splenic", "Spleenic": "Splenic", "Na<ef>ve": "naive"}


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No data to write: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def soft_samples(path):
    """Read each complete GEO sample block. The contact fields are not exported."""
    records = {}
    declared_ids = []
    declared_updates = []
    current = None
    with gzip.open(path, mode="rt", encoding="utf-8") as handle:
        for line in handle:
            key, _, value = line.rstrip().partition(" = ")
            if key == "!Series_last_update_date":
                declared_updates.append(value)
            if key == "!Series_sample_id":
                declared_ids.append(value)
            if key == "^SAMPLE":
                if value in records:
                    raise ValueError(f"Duplicate SOFT record: {value}")
                current = records[value] = defaultdict(list)
            elif current is not None and key.startswith("!Sample_"):
                current[key.removeprefix("!Sample_")].append(value)
    if set(records) != set(declared_ids) or len(records) != len(declared_ids):
        raise ValueError("The series sample index and all SOFT sample blocks disagree")
    if len(declared_updates) != 1:
        raise ValueError("Expected a single GEO series last update date")
    return records, declared_updates[0]


def counts_library_keys(path):
    with gzip.open(path, mode="rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        key, *libraries = next(reader)
        if key != "GeneSymbol" or len(libraries) != len(set(libraries)):
            raise ValueError("The matrix must have one GeneSymbol column and unique libraries")
        genes = set()
        for row in reader:
            if len(row) != len(libraries) + 1 or not row[0] or row[0] in genes:
                raise ValueError(f"Invalid count-row shape or duplicate symbol at row {reader.line_num}")
            genes.add(row[0])
    return libraries, genes


def metadata_value(record, key):
    (value,) = record[key]
    return value


def aliases_by_accession(path):
    result = {}
    for rule in read_tsv(path):
        for gsm in rule["gsm_pairs"].split(";"):
            if gsm in result:
                raise ValueError(f"GSM repeated in alias rules: {gsm}")
            result[gsm] = rule
    return result


def sample_identity(gsm, sample_label, ch, counts_libraries, aliases):
    group, suffix = sample_label.rsplit("#", 1)
    if not suffix.isdecimal():
        raise ValueError(f"Cannot identify author replicate suffix: {gsm}")
    if sample_label in counts_libraries:
        if gsm in aliases:
            raise ValueError(f"Alias {gsm} now has a native count-column match")
        return group, suffix, sample_label, "exact_title", "SOFT title stripped of _RNA-seq"
    if gsm in aliases:
        rule = aliases[gsm]
        count_group = rule["matrix_group"]
        count_label = f"{count_group}#{suffix}"
        if (group != rule["geo_group"]
                or count_label not in counts_libraries
                or rule["required_gate_or_characteristic"] not in ch["sorting markers"]):
            raise ValueError(f"The prespecified sample alias fails for {gsm}")
        return count_group, suffix, count_label, "source_documented_alias", rule["evidence"]
    return group, suffix, "NA", "not_in_count_matrix", "No verified vM25 count column"


def classify_lineage(group, soft_cell_type):
    if group.startswith("Treg."):
        return "Treg"
    if group.startswith(("NKT.", "MAIT.")):
        return "Unconventional T"
    if group.startswith(("T.4.", "T.4")):
        return "CD4 T"
    if group.startswith(("T.8.", "T8.")):
        return "CD8 T"
    if group.startswith("NK."):
        return "NK"
    if group.startswith("ILC"):
        return "ILC"
    if group.startswith(("preT.", "proB.", "T.D", "T.ISP", "B.FrE.", "B.T1.", "B.T2.", "B.T3.")):
        return "Progenitor"
    if soft_cell_type not in SOURCE_LINEAGE:
        raise ValueError(f"Unclassified author cell type: {group}, {soft_cell_type}")
    return SOURCE_LINEAGE[soft_cell_type]


def eligibility(group, lineage, age, strain):
    if re.search(r"\.(?:LPS|LCMV|pIC|Thio|aCD3|SeT)", group, flags=re.IGNORECASE):
        return "challenged_or_induced", "Explicit challenge/induction in source label"
    if ".Nb." in group:
        return "neonatal_label", "Neonatal-labelled group; SOFT generic age field cannot resolve"
    if group == "T8.TN.P14.Sp":
        return "antigen_specific_transgenic_control", "Author P14 antigen-specific CD8 group"
    if group == "MC.heparinase.PC":
        return "manipulated_preparation", "Heparinase-labelled preparation"
    if lineage == "Stromal":
        return "nonimmune", "Nonimmune stromal or epithelial population"
    if lineage == "Progenitor" or group in {"T.4.Th", "T.8.Th"} or (
        group.startswith("Tgd.") and group.endswith(".Th")
    ):
        return "developing_immune_cells", "Progenitor, transitional B cell, or developing thymocyte"
    if group.startswith(("NKT.", "MAIT.")) and group.endswith(".Th"):
        return "thymic_maturity_unresolved", "Thymic sort without mature-stage resolution"
    if strain != "C57BL/6J" or age != "6 weeks old":
        return "outside_main_age_or_strain", "Outside locked 6-week C57BL/6J atlas comparison"
    if group in EXTRA_B:
        return "activated_b_context_unconfirmed", "Germinal-center B or plasmablast; activation context kept separate"
    return "baseline_eligible", "Mature sorted immune group with no challenge indicated in metadata"


def clean_group_label(group, sample_label, descriptions):
    if group in DISPLAY_LABELS:
        return DISPLAY_LABELS[group]
    exact_descriptions = [value for value in descriptions if
                          value not in (sample_label, "NA")
                          and not value.startswith("processed data file:")
                          and not value.startswith("SK-")]
    if not exact_descriptions:
        return group
    result = re.sub(r"\s+", " ", exact_descriptions[0])
    for old, new in DESCRIPTION_FIXES.items():
        result = result.replace(old, new)
    return result


def sample_and_group_tables(records, libraries, alias_rules):
    data = []
    already_used = set()
    already_biosamples = set()
    reference_libraries = set(libraries)
    for gsm, fields in records.items():
        if metadata_value(fields, "geo_accession") != gsm:
            raise ValueError(f"Sample block/header mismatch: {gsm}")
        title = metadata_value(fields, "title")
        if not title.endswith("_RNA-seq"):
            raise ValueError(f"Unexpected RNA-seq title: {title}")
        sample_label = title.removesuffix("_RNA-seq")
        ch = dict(value.split(": ", 1) for value in fields["characteristics_ch1"])
        group, replicate, matrix_id, match, match_evidence = sample_identity(
            gsm, sample_label, ch, reference_libraries, alias_rules
        )
        lineage = classify_lineage(group, ch["cell type"])
        state, state_evidence = eligibility(group, lineage, ch["age"], ch["strain"])
        (biosample,) = [v.rsplit("/", 1)[-1] for v in fields["relation"] if v.startswith("BioSample: ")]
        if biosample in already_biosamples:
            raise ValueError(f"BioSample duplication requires a technical-replicate policy: {biosample}")
        already_biosamples.add(biosample)
        if matrix_id != "NA":
            if matrix_id in already_used:
                raise ValueError(f"Matrix column assigned to multiple GSMs: {matrix_id}")
            already_used.add(matrix_id)

        label = clean_group_label(group, sample_label, fields["description"])
        inclusion = state == "baseline_eligible" and matrix_id != "NA"
        sex_in_title = "female" if group.split(".")[-2:-1] == ["Fem"] or group.endswith(".F") else "NA"
        if group.endswith(".M"):
            sex_in_title = "male"
        phase = "primary" if inclusion else "excluded"
        reason = state_evidence if matrix_id != "NA" else f"No verified vM25 count column. {state_evidence}"
        data.append(
            {
                "dataset": SERIES,
                "gsm": gsm,
                "species": metadata_value(fields, "organism_ch1"),
                "source_title": title,
                "source_group": sample_label.rsplit("#", 1)[0],
                "source_description": " | ".join(fields["description"]) or "NA",
                "source_name": metadata_value(fields, "source_name_ch1"),
                "source_cell_type": ch["cell type"],
                "cell_type": group,
                "lineage": lineage,
                "label": label,
                "tissue": "Spleen" if ch["tissue"] == "Splenic" else ch["tissue"],
                "source_tissue": ch["tissue"],
                "sorting_markers": ch["sorting markers"],
                "strain": ch["strain"],
                "source_age": ch["age"],
                "sex_in_title": sex_in_title,
                "donor_id": "NA",
                "author_replicate": replicate,
                "replicate_id": biosample,
                "replicate_unit": "GEO biological library/BioSample; animal identities not provided",
                "state": state,
                "include_primary": "yes" if inclusion else "no",
                "inclusion": phase,
                "reason": reason,
                "matrix_library": matrix_id,
                "matrix_match": match,
                "matrix_match_evidence": match_evidence,
                "alias_evidence_url": alias_rules.get(gsm, {}).get("evidence_url", "NA"),
                "source_url": GEO + gsm,
            }
        )
    if reference_libraries != already_used:
        raise ValueError(f"Unlinked count columns: {reference_libraries - already_used}")
    by_group = defaultdict(list)
    for row in data:
        by_group[row["cell_type"]].append(row)

    def group_sort_key(key):
        order = LINEAGES.index(by_group[key][0]["lineage"])
        if key in MAIN_TREG:
            return order, MAIN_TREG.index(key), ""
        return order, 0, key.lower()

    groups = []
    for order, group in enumerate(sorted(by_group, key=group_sort_key), 1):
        rows = by_group[group]
        row = rows[0]
        if any(r["state"] != row["state"] or r["lineage"] != row["lineage"] or
               r["tissue"] != row["tissue"] for r in rows):
            raise ValueError(f"Inconsistent state, lineage or tissue within group {group}")
        matrix_n = sum(r["matrix_library"] != "NA" for r in rows)
        primary_n = sum(r["include_primary"] == "yes" for r in rows)
        groups.append(
            {
                "dataset": SERIES,
                "group_order": order,
                "cell_type": group,
                "label": row["label"],
                "lineage": row["lineage"],
                "tissue": row["tissue"],
                "state": row["state"],
                "source_age": row["source_age"],
                "source_titles": " | ".join(r["source_title"] for r in rows),
                "n_soft_gsm": len(rows),
                "n_count_libraries": matrix_n,
                "n_libraries": primary_n,
                "include_primary": "yes" if primary_n else "no",
                "source_sort": " | ".join(sorted({r["sorting_markers"] for r in rows})),
                "reason": row["reason"] if matrix_n == 0 else eligibility(
                    group, row["lineage"], row["source_age"], row["strain"]
                )[1],
                "n_exact_title": sum(r["matrix_match"] == "exact_title" for r in rows),
                "n_source_documented_alias": sum(r["matrix_match"] == "source_documented_alias" for r in rows),
            }
        )
    return data, groups


def orthologs(mgi_path, panel_path, matrix_symbols):
    fields = (
        "mouse_mgi", "mouse_gene", "mouse_ncbi", "human_hgnc", "human_gene", "human_ncbi",
    )
    with mgi_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, fieldnames=fields, delimiter="\t"))
    if any(None in row or any(not row[k] for k in fields) for row in rows):
        raise ValueError("Unexpected MGI HOM_ProteinCoding.rpt format")
    id_counts = {key: Counter(row[key] for row in rows) for key in fields}
    raw_pairs = defaultdict(list)
    for line_number, row in enumerate(rows, 1):
        raw_pairs[row["human_gene"]].append({"mgi_line": line_number, **row})
    mappings = []
    for panel_gene in read_tsv(panel_path):
        human = panel_gene["human_gene"]
        hits = raw_pairs[human]
        resolved = (len(hits) == 1 and
                    all(id_counts[key][hits[0][key]] == 1 for key in fields))
        mapping = {key: "NA" for key in ("mgi_line", *fields)}
        if resolved:
            mapping = hits[0].copy()
        status = "mapped_unique" if resolved else ("ambiguous" if hits else "no_mgi_ortholog")
        present = resolved and mapping["mouse_gene"] in matrix_symbols
        missing = "NA"
        if not resolved:
            missing = "No uniquely supported MGI one-to-one ortholog"
        elif not present:
            missing = "Ortholog symbol absent from vM25 gene-count matrix"
        mappings.append(
            {
                "dataset": SERIES,
                **mapping,
                **panel_gene,
                "orthology_status": status,
                "measured": "yes" if present else "no",
                "missing_reason": missing,
                "annotation_version": "GENCODE vM25 matrix GeneSymbol; MGI snapshot 2026-09-07",
                "mapping_basis": "Unique MGI/HGNC symbols and IDs and exact matrix GeneSymbol",
                "mgi_source_url": "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt",
            }
        )
    ambiguous = sum(any(id_counts[key][row[key]] != 1 for key in fields) for row in rows)
    return mappings, len(rows), ambiguous


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    config = args.out / "07_config"
    audit = args.out / "01_audit"
    source = args.download_dir / SERIES
    libraries, matrix_symbols = counts_library_keys(source / COUNTS)
    records, update_date = soft_samples(source / SOFT)
    sample_rows, group_rows = sample_and_group_tables(
        records, libraries, aliases_by_accession(config / "immgen_aliases.tsv")
    )
    mapping, n_mgi, n_ambiguous = orthologs(
        args.download_dir / "annotation/HOM_ProteinCoding.rpt",
        config / "immgen_panel.tsv",
        matrix_symbols,
    )
    if sum(row["panel_role"] == "candidate" for row in mapping) != 14:
        raise ValueError("The fixed mouse candidate panel requires 14 rows")
    write_tsv(audit / "R10_IMMGEN_SAMPLE_CROSSWALK.tsv", sample_rows)
    write_tsv(audit / "R10_IMMGEN_GROUP_CROSSWALK.tsv", group_rows)
    write_tsv(audit / "R10_IMMGEN_GENE_CROSSWALK.tsv", mapping)
    primary = [row for row in group_rows if row["include_primary"] == "yes"]
    print(
        f"{SERIES} SOFT updated {update_date}: {len(sample_rows)} GSM, "
        f"{len({row['source_group'] for row in sample_rows})} literal title groups, "
        f"{len(group_rows)} groups after resolving distinct sorting groups"
    )
    print(
        f"vM25 matrix: {len(matrix_symbols)} unique gene symbols x {len(libraries)} count libraries, "
        f"{sum(row['matrix_match'] == 'exact_title' for row in sample_rows)} exact title, "
        f"{sum(row['matrix_match'] == 'source_documented_alias' for row in sample_rows)} aliases"
    )
    print(
        f"Baseline main groups: {len(primary)} groups, "
        f"{sum(row['n_libraries'] for row in primary)} libraries; "
        f"lineages={dict(Counter(row['lineage'] for row in primary))}"
    )
    print(
        f"MGI report: {n_mgi} rows, {n_ambiguous} rows fail global ID/symbol uniqueness; "
        f"{sum(row['measured'] == 'yes' for row in mapping)}/{len(mapping)} panel genes map and are counted"
    )


if __name__ == "__main__":
    main()
