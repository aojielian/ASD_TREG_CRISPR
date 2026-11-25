# Prepare independent signatures without opening any knockout counts or effect files.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
    stop("usage: Rscript prepare_signature_membership.R DATA OUT HUMAN_MSIGDB_RDS MOUSE_HALLMARK_GMT MGI_REPORT")
}
data_root <- normalizePath(args[1], mustWork = TRUE)
out <- normalizePath(args[2])
human_cache <- normalizePath(args[3])
mouse_gmt <- normalizePath(args[4])
r10 <- file.path(data_root, "data/r10")
mgi_file <- normalizePath(args[5], mustWork = TRUE)
group_file <- file.path(r10, "reference/R10_IMMGEN_GROUP_CROSSWALK.tsv")
sample_file <- file.path(r10, "reference/R10_IMMGEN_SAMPLE_CROSSWALK.tsv")
atlas_file <- file.path(r10, "results/R10_IMMGEN_CELLTYPE_CPM.tsv.gz")
file_sha <- function(path) digest::digest(file = path, algo = "sha256")
write_tsv <- function(x, path) {
    write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}
hallmarks <- paste0("HALLMARK_", c("IL2_STAT5_SIGNALING", "E2F_TARGETS", "G2M_CHECKPOINT", "APOPTOSIS",
                                  "TNFA_SIGNALING_VIA_NFKB", "INTERFERON_GAMMA_RESPONSE"))
targets <- c("Foxp3", "Foxp1", "Sin3a")
hill_sets <- c("GSE7460_TCONV_VS_TREG_LN_DN", "GSE7460_TCONV_VS_TREG_LN_UP")
hill_urls <- paste0("https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/", hill_sets, ".html")
mh_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/mh.all.v2025.1.Mm.symbols.gmt"
atlas_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109125"

mgi <- read.delim(mgi_file, header = FALSE, sep = "\t", colClasses = "character", quote = "")
names(mgi) <- c("MGI", "mouse_symbol", "mouse_Entrez", "HGNC", "human_symbol", "human_Entrez")
non_unique <- function(x) is.na(x) | !nzchar(x) | duplicated(x) | duplicated(x, fromLast = TRUE)
mgi <- mgi[!Reduce("|", lapply(mgi, non_unique)), ]
stopifnot(all(vapply(mgi, function(x) !anyDuplicated(x), logical(1))))

human <- readRDS(human_cache)
hill <- human[match(human$gs_name, hill_sets, nomatch = 0L) > 0L, ]
rm(human)
stopifnot(
    identical(sort(unique(hill$db_version)), "2025.1.Hs"),
    all(hill$gs_collection == "C7"),
    all(hill$gs_subcollection == "IMMUNESIGDB"),
    all(hill$gs_pmid == "18024188"),
    all(hill$gs_geoid == "GSE7460"),
    all(table(hill$gs_name) == 200),
    !anyDuplicated(hill[c("gs_name", "db_gene_symbol")]),
    unname(identical(hill$gs_exact_source[match(hill_sets, hill$gs_name)], c("GSE7460_1327_200_DN", "GSE7460_1327_200_UP")))
)
match_id <- match(hill$db_ncbi_gene, mgi$human_Entrez)
match_symbol <- match(hill$db_gene_symbol, mgi$human_symbol)
mapping_conflict <- !is.na(match_id) & !is.na(match_symbol) & match_id != match_symbol
chosen <- ifelse(!is.na(match_id), match_id, match_symbol)
chosen[mapping_conflict] <- NA_integer_
mouse_symbol <- mgi$mouse_symbol[chosen]
map_route <- ifelse(mapping_conflict, "MGI_id_symbol_conflict",
                    ifelse(is.na(chosen), "no_unique_MGI_match",
                           ifelse(!is.na(match_id) & !is.na(match_symbol), "MGI_unique_id_symbol",
                                  ifelse(!is.na(match_id), "MGI_unique_Entrez_id", "MGI_unique_symbol"))))
arm <- ifelse(hill$gs_name == hill_sets[1], "Treg_high", "Tconv_high")
dual_member <- intersect(na.omit(mouse_symbol[arm == "Treg_high"]), na.omit(mouse_symbol[arm == "Tconv_high"]))
within_key <- paste(hill$gs_name, mouse_symbol)
duplicate_within_arm <- !is.na(mouse_symbol) & (duplicated(within_key) | duplicated(within_key, fromLast = TRUE))
eligible <- !is.na(mouse_symbol) & !mouse_symbol %in% c(targets, dual_member) & !duplicate_within_arm
hill_out <- data.frame(
    program = "Hill_LN", component = arm, source_set = hill$gs_name,
    source_symbol = hill$db_gene_symbol, source_original_symbol = hill$source_gene,
    source_gene_id = hill$db_ncbi_gene, mouse_gene = mouse_symbol, mouse_gene_id = mgi$mouse_Entrez[chosen],
    mgi_id = mgi$MGI[chosen], source_version = hill$db_version,
    included = eligible, mapping_route = map_route,
    reason = ifelse(is.na(mouse_symbol), map_route, ifelse(mouse_symbol %in% targets, "excluded_fixed_endpoints",
        ifelse(mouse_symbol %in% dual_member, "overlap_both_arms", ifelse(duplicate_within_arm, "duplicate_mouse_gene", "included")))),
    reference_log2_ratio = NA_real_, reference_rank = NA_integer_, reference_Treg_CPM = NA_real_, reference_Tconv_CPM = NA_real_
)

groups <- read.delim(group_file, sep = "\t", check.names = FALSE, quote = "")
compare <- c("Treg.4.25hi.Sp", "T.4.Nve.Sp", "T.4.Nve.Fem.Sp")
g <- groups[match(compare, groups$cell_type), ]
stopifnot(
    identical(g$cell_type, compare), all(g$tissue == "Spleen"), all(g$state == "baseline_eligible"),
    all(g$include_primary == "yes"), all(g$source_age == "6 weeks old"),
    grepl("CD25hi", g$source_sort[1], fixed = TRUE),
    all(grepl("CD25-", g$source_sort[2:3], fixed = TRUE))
)
eligible_lns <- groups$tissue %in% c("Lymph node", "Lymph Node", "Lymph nodes") & groups$include_primary == "yes"
stopifnot(!any(eligible_lns & groups$lineage == "Treg"))
atlas <- read.delim(gzfile(atlas_file), sep = "\t", check.names = FALSE, quote = "")
stopifnot("mouse_gene" %in% names(atlas), all(compare %in% names(atlas)), !anyDuplicated(atlas$mouse_gene))
treg_cpm <- atlas[[compare[1]]]
tconv_cpm <- rowMeans(atlas[compare[2:3]])
ratio <- log2((treg_cpm + 0.5) / (tconv_cpm + 0.5))
eligible_atlas <- nzchar(atlas$mouse_gene) & !atlas$mouse_gene %in% targets &
    is.finite(ratio) & pmax(treg_cpm, tconv_cpm) >= 1 & ratio != 0
positive <- which(eligible_atlas & ratio > 0)
negative <- which(eligible_atlas & ratio < 0)
by_magnitude <- function(i) i[order(-abs(ratio[i]), atlas$mouse_gene[i], method = "radix")]
pos_frozen <- head(by_magnitude(positive), 100)
neg_frozen <- head(by_magnitude(negative), 100)
frozen <- c(pos_frozen, neg_frozen)
tier <- rep(c("Treg_high", "Tconv_high"), c(length(pos_frozen), length(neg_frozen)))
imm_out <- data.frame(
    program = "ImmGen_spleen", component = tier,
    source_set = paste0("R10_ImmGen_spleen_", tier),
    source_symbol = atlas$mouse_gene[frozen], source_original_symbol = atlas$mouse_gene[frozen],
    source_gene_id = NA_character_, mouse_gene = atlas$mouse_gene[frozen], mouse_gene_id = NA_character_,
    mgi_id = NA_character_, source_version = "R10_GENCODE_vM25",
    included = TRUE, mapping_route = "unique_frozen_GENCODE_vM25_symbol", reason = "included",
    reference_log2_ratio = ratio[frozen],
    reference_rank = c(seq_along(pos_frozen), seq_along(neg_frozen)),
    reference_Treg_CPM = treg_cpm[frozen], reference_Tconv_CPM = tconv_cpm[frozen]
)
atlas_locks <- data.frame(
    dataset = "GSE109125", group = g$cell_type, comparator_role = c("Treg", "Tconv", "Tconv"),
    source_gate = g$source_sort, tissue = g$tissue, source_age = g$source_age,
    n_libraries = g$n_libraries, group_weight = c(1, 0.5, 0.5),
    atlas_eligible_positive = length(positive), atlas_eligible_negative = length(negative)
)
write_tsv(atlas_locks, file.path(out, "01_audit/R11_SIGNATURE_IMMGEN_GROUPS.tsv"))

gmt <- strsplit(readLines(mouse_gmt), "\t", fixed = TRUE)
set_names <- vapply(gmt, "[", character(1L), 1L)
stopifnot(length(gmt) == 50L, !anyDuplicated(set_names), all(hallmarks %in% set_names))
gmt <- gmt[match(hallmarks, set_names)]
mh_rows <- lapply(gmt, function(x) {
    members <- x[-c(1:2)]
    stopifnot(length(members) > 20, all(nzchar(members)), !anyDuplicated(members))
    data.frame(
        program = x[1], component = "member_mean", source_set = x[1],
        source_symbol = members, source_original_symbol = members,
        source_gene_id = NA_character_, mouse_gene = members, mouse_gene_id = NA_character_,
        mgi_id = NA_character_, source_version = "2025.1.Mm",
        included = TRUE, mapping_route = "official_mouse_MSigDB_symbol", reason = "included",
        reference_log2_ratio = NA_real_, reference_rank = NA_integer_, reference_Treg_CPM = NA_real_,
        reference_Tconv_CPM = NA_real_
    )
})
all_members <- rbind(hill_out, imm_out, do.call(rbind, mh_rows))
out_path <- file.path(out, "03_results/R11_SIGNATURE_MEMBERS.tsv.gz")
gz <- gzfile(out_path, "wt", compression = 9L)
write_tsv(all_members, gz)
close(gz)
subset_gmt <- file.path(out, "00_download/files/MSigDB/R11_MSIGDB_MOUSE_MH_2025.1_SELECTED.gmt")
writeLines(vapply(gmt, paste, collapse = "\t", FUN.VALUE = character(1L)), subset_gmt)
manifest <- data.frame(
    reference = c("Hill human C7 ImmuneSigDB", "R10 MGI Mouse-Human one-to-one",
                  "R10 ImmGen median TMM CPM", "R10 ImmGen group gates", "R10 ImmGen sample labels",
                  "Mouse Hallmark GMT", "Mouse Hallmark selected six sets"),
    version = c("2025.1.Hs", "R10 snapshot", "R10 GENCODE vM25", "R10", "R10", "2025.1.Mm", "2025.1.Mm"),
    location = c("external/MSigDB/msigdb.2025.1.Hs.rds", "external/annotation/HOM_ProteinCoding.rpt",
                 "data/r10/results/R10_IMMGEN_CELLTYPE_CPM.tsv.gz",
                 "data/r10/reference/R10_IMMGEN_GROUP_CROSSWALK.tsv",
                 "data/r10/reference/R10_IMMGEN_SAMPLE_CROSSWALK.tsv",
                 "external/MSigDB/mh.all.v2025.1.Mm.symbols.gmt",
                 "00_download/files/MSigDB/R11_MSIGDB_MOUSE_MH_2025.1_SELECTED.gmt"),
    sha256 = vapply(c(human_cache, mgi_file, atlas_file, group_file, sample_file, mouse_gmt, subset_gmt), file_sha, character(1L)),
    access_date = c("Official 2025.1.Hs release", rep("2026-11-25 verified against local frozen R10", 4),
                    "2026-11-25 official download", "2026-11-25 preserved exact subset"),
    source_url = c("https://zenodo.org/records/15800824/files/msigdb.2025.1.Hs.rds?download=1",
                   "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt",
                   rep(atlas_url, 3), mh_url, mh_url),
    pmid = c("18024188", NA, rep("30686579", 3), rep("26771021", 2)),
    doi = c("10.1016/j.immuni.2007.09.010", NA, rep("10.1016/j.cell.2018.12.036", 3), rep("10.1016/j.cels.2015.12.004", 2))
)
write_tsv(manifest, file.path(out, "01_audit/R11_SIGNATURE_SOURCE.tsv"))
defs <- unique(all_members[c("program", "component", "source_set", "source_version")])
defs$direction <- c("Tconv_high", "Treg_high")[1 + (defs$component == "Treg_high")]
defs$direction[defs$component == "member_mean"] <- "member_mean"
defs$source_url <- c(hill_urls[match(defs$source_set[defs$program == "Hill_LN"], hill_sets)],
                     rep(atlas_url, sum(defs$program == "ImmGen_spleen")),
                     vapply(gmt, "[", character(1L), 2L))
defs$starting_member_n <- as.integer(table(factor(all_members$source_set, levels = defs$source_set)))
defs$frozen_member_n <- as.integer(table(factor(all_members$source_set[all_members$included], levels = defs$source_set)))
write_tsv(defs, file.path(out, "01_audit/R11_SIGNATURE_DEFINITIONS.tsv"))
stopifnot(nrow(defs) == 10L, all(defs$frozen_member_n >= 20))
print(defs[c("program", "component", "starting_member_n", "frozen_member_n")], row.names = FALSE)
print(table(hill_out$component, hill_out$reason))
