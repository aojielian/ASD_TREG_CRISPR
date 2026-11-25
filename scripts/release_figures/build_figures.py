"""Generate four main figures and one paginated supplementary PDF from verified project data.

The released placement uses the direct phenotype and breadth comparison in Figure 3,
and the Foxp3 RNA and Hill identity comparison in Figure 4.
"""
from figure_common import *
from knockout_figures import (endpoint_figure,marker_figure,program_effects,qc_overview,qc_distances,count_distribution,evidence_graph,HALLMARKS)

def original_captions(r10):
    text=(r10/"R10_FIGURE_LEGENDS_EN.md").read_text()
    return {int(n):(title,body.strip()) for n,title,body in re.findall(r"## Figure (\d+)\. ([^\n]+)\n\n(.*?)(?=\n## Figure |\Z)",text,re.S)}

def build(stage,data):
    r10=data/"figure_sources/frozen_r10"
    src=data/"figure_sources/R10"
    (stage/"04_figures/source_data").mkdir(parents=True,exist_ok=True)
    sup=stage/"05_supplementary"; figs=stage/"04_figures"
    sup.mkdir(parents=True,exist_ok=True); figs.mkdir(parents=True,exist_ok=True)
    original=original_captions(r10);main=original.copy()
    for n in (1,2):
        name=f"Figure_{n}.pdf"
        shutil.copyfile(r10/name,figs/name)
        tiff_from_pdf(figs/name)
    endpoint_title="Foxp3 RNA and a lymph-node Treg identity signature after conditional deletion"
    endpoint_body=("Public mouse RNA-sequencing data comprise Foxp1 conditional deletion in naïve and activated Tregs from GSE121279 (three biological libraries per group), "
      "Sin3a conditional deletion in GSE263830 (six pooled libraries per genotype, three from each sex), and an independent Foxp1 study of resting and activated Tregs in GSE121251 (two biological libraries per group). "
      "(A) Each point represents one library on its study's log₂(DESeq2-normalized Foxp3 count + 1) scale. Lines indicate within-group medians. Open symbols indicate controls and filled symbols conditional deletions. "
      "In the Sin3a panel, circles indicate female pools and squares male pools. The y-axis limits are specific to each cohort. "
      "(B) Points show unshrunken knockout-minus-control log₂ fold changes with 95% Wald intervals. The Foxp1 models include genotype, state and their interaction. The Sin3a model adjusts for sex. "
      "(C) The signed identity score is the difference between mean within-sample expression percentile ranks of the two Hill lymph-node Treg-high and Tconv-high gene sets. "
      "Target genes and Foxp3 were excluded from the reference lists. (D) Points show model-based score differences with 95% t intervals. "
      "Within each endpoint, the three core contrasts above the dotted line and the two independent-study contrasts below it form separate Holm families. "
      "The RNA abundance and identity-score assays were reanalyzed from published data. Sources include Konopacki et al. (PMID 30643266), Christensen et al. (PMID 39156895), Ren et al. (PMID 31125332) and Hill et al. (PMID 18024188).")
    main_export(revised_context_figure(src),stage,3)
    main_export(endpoint_figure(stage),stage,4)
    old3=original[3][1];p3=old3.split("(B)")[0].strip()+" "
    rest=old3.split("(C)")[1].split(" Sources:")[0]
    primarybody=rest.split("(D)")[0].strip()
    continuousbody=rest.split("(D)")[1].strip()
    a=adjusted_continuous(src)
    br=original[4][1].split("(B) ")[1].split("(C)")[0].replace("The two panels of source data retain their original species-specific units.","").strip()
    main[3]=("Direct FOXP3 phenotype comparisons and immune-cell expression breadth",
      p3+"(B) "+primarybody+" (C) "+continuousbody+f" The property-adjusted continuous contrast had two-sided P = {a.P:.5f}. (D) "+br+
      " Supplementary Figures S15 and S16 retain the discordant Firth intervals, the full human expression map and the Chen protein readouts. Sources include Chen et al. (PMID 40140585), Shifrut et al. (PMID 30449619), Umhoefer et al. (PMID 41237776), Monaco et al. (PMID 30726743) and ImmGen ULI (GSE109125).")
    main[4]=(endpoint_title,endpoint_body)
    (figs/"Figure_Rendered_Legends.md").write_text("# Figure legends\n\n"+"\n\n".join(f"## Figure {i}. {title}\n\n{body}" for i,(title,body) in sorted(main.items()))+"\n")
    doc=fitz.open()
    with fitz.open(r10/"R10_Supplementary_Figures.pdf") as legacy:
        assert len(legacy)==14
        doc.insert_pdf(legacy)
    legends=(r10/"R10_SUPPLEMENTARY_FIGURE_LEGENDS_EN.md").read_text().rstrip()
    rows=[]
    def add(graphic,title,body,sources,fields,history):
        nonlocal legends
        i=doc.page_count+1
        with graphic:
            cap_page(doc,graphic,i,title,body)
        legends+=f"\n\n## Supplementary Figure S{i}. {title}\n\n{body}\n"
        rows.append([f"S{i}",history,title,sources,fields])
    for n in (3,4):
        add(fitz.open(r10/f"Figure_{n}.pdf"),*original[n],
            "figure_sources/R10/F3* and F4*",
            "Original vector page preserved with its complete legend",f"R10 Figure {n}")
    key=("F1N and F1A indicate naïve and activated Foxp1 deletion in GSE121279, S3 indicates sex-adjusted Sin3a deletion in GSE263830, and F2R and F2A indicate resting and activated Foxp1 deletion in GSE121251. "
         "Genotype groups contain three libraries per state in GSE121279, six pooled libraries in GSE263830 and two libraries per state in GSE121251.")
    add(figpdf(qc_overview(stage)),"Sample-level principal components and raw library sizes",
        "(A) Within each cohort, the principal components use centered log₂(TMM CPM + 0.5) over all genes passing the prespecified expression filter, without per-gene scaling. Open and filled symbols represent control and conditional knockout libraries. "
        "In the Foxp1 studies, circles denote naïve or resting and squares activated samples. In the Sin3a study, circles denote female pools and squares male pools. Axes give the fraction of within-cohort variance explained. "
        "(B) Bars show summed raw gene counts for each input library. Navy denotes control and teal conditional knockout. The final three digits of each GEO sample accession identify the library within its cohort. "+key,
        "data/r11/results/R11_QC_SAMPLE_PCA.tsv and R11_QC_SAMPLES.tsv","32 libraries; PC1, PC2, variance fractions and raw library size","R11")
    add(figpdf(qc_distances(stage)),"Within-study RNA expression similarity",
        "(A) Colors show pairwise Pearson correlations of library expression profiles across count-filtered genes. (B) Colors show Euclidean distances between profiles, using the same within-study log₂(TMM CPM + 0.5) expression space. "
        "Axes list GEO sample accession suffixes in matched order. Each cohort has its own observed scale, displayed beside the matrix. All biological libraries are included. "+key,
        "data/r11/results/R11_QC_SAMPLE_DISTANCES.tsv.gz and data/r11/reference/R11_SAMPLE_MAP.tsv","Pairwise Pearson correlation and Euclidean distance","R11")
    add(figpdf(count_distribution(stage)),"Raw gene-count distributions in every biological library",
        "Each vertical mark summarizes one input RNA-sequencing library before gene filtering. The thick segment spans the 25th to 75th percentiles of raw gene counts. A horizontal mark indicates the median, a dot the 95th percentile and the upper tip the 99th percentile. "
        "All percentiles include zero counts and are displayed after a log₂(count + 1) transformation. Navy denotes controls and teal conditional knockouts. The x-axis gives the final three digits of each GEO sample accession. "+key,
        "data/r11/results/R11_QC_SAMPLES.tsv","q25, median, q75, q95 and q99 for all genes (raw counts per library)","R11")
    add(figpdf(marker_figure(stage)),"Fixed marker effects and targeted gene RNA abundance",
        "(A) Colors show unshrunken DESeq2 knockout-minus-control log₂ fold changes for all 13 prespecified auxiliary genes. Large dots indicate Benjamini–Hochberg q < 0.05 across measurable marker-by-contrast tests. Small dots indicate q ≥ 0.05. A cross denotes an unestimable effect. "
        "(B) Points and lines show the targeted gene's RNA effect and its 95% Wald interval. Gene-level counts can include transcripts retained after conditional deletion. Targeted gene counts therefore describe RNA abundance. "+key,
        "data/r11/results/R11_MARKER_EFFECTS.tsv and R11_TARGET_QC.tsv","All 13 marker rows by 5 contrasts; BH and target log₂FC with 95% CI","R11")
    cov=read(stage/"03_results"/"R11_SIGNATURE_CONTRASTS.tsv")
    hcov=cov[cov.program.eq("Hill_LN")&cov.component.eq("score")].set_index("contrast")
    imm=cov[cov.program.eq("ImmGen_spleen")&cov.component.eq("score")].set_index("contrast")
    hallmark=cov[cov.program.str.startswith("HALLMARK_")&cov.component.eq("score")]
    nf=int(imm.adjustment_family_size.iloc[0]);hf=int(hallmark.adjustment_family_size.iloc[0])
    add(figpdf(program_effects(stage,["Hill_LN","ImmGen_spleen"],["Hill lymph node","ImmGen spleen"],"Signed Treg identity scores")),"Two independent Treg expression references",
        "The same library-level expression ranks were scored against the prespecified Hill lymph-node reference (GSE7460, PMID 18024188) and an independent baseline splenic ImmGen reference (GSE109125). For each reference, the score is the mean rank of the Treg-high genes minus the mean rank of the Tconv-high genes. "
        "Dots and intervals give knockout-minus-control changes and 95% t intervals from study-specific linear models. Filled points mark adjusted P < 0.05. Open points show estimates above that threshold. "
        f"Hill uses Holm adjustment separately across the three core contrasts and the two independent-study contrasts. ImmGen uses BH across {nf} auxiliary contrasts. "+key,
        "data/r11/results/R11_SIGNATURE_CONTRASTS.tsv","Hill_LN and ImmGen_spleen; component=score; five contrasts","R11")
    add(figpdf(program_effects(stage,[g for g,n in HALLMARKS],[n for g,n in HALLMARKS],"Prespecified Hallmark programs")),"IL-2/STAT5 and broad cellular program effects",
        "Each panel shows one of the six prespecified mouse MSigDB Hallmark programs. The score is the mean within-library expression percentile rank of measurable member genes. Points show knockout-minus-control score changes from separate cohort models, with 95% t intervals. "
        f"Filled symbols indicate BH q < 0.05 across the frozen family of {hf} Hallmark-by-contrast tests. Open symbols indicate q ≥ 0.05. Each panel spans its observed intervals. The IL-2/STAT5 panel measures a transcriptional program. "+key,
        "data/r11/results/R11_SIGNATURE_CONTRASTS.tsv","All six fixed mouse MH programs, component score, all five contrasts","R11")
    specs=read(stage/"01_audit"/"R11_FIGURE_FUNCTION_ASSAYS.tsv").to_dict("records")
    add(figpdf(evidence_graph(stage,src,specs)),"Cross-study evidence for FOXP1 and SIN3A",
        "(A) The human induced-Treg genome-wide screen and primary Treg/Tconv targeted screens measured FOXP3 protein-bin changes. The fourth column gives the independent stimulated T-cell proliferation screen. The original hit calls and tested library coverage are retained. "
        "(B) Mouse conditional-deletion rows show unshrunken Foxp3 RNA log₂ fold changes and signed Hill reference score differences on separate scales. Holm support uses the family for each endpoint and study role. "
        "(C) Published mouse Treg suppression experiments are identified by source figure. Each entry describes its experimental context. The new computational measurements in B use RNA from those public studies. "+
        " ".join(dict.fromkeys(str(v["detail"]) for v in specs))+" "+key,
        "figure_sources/R10/R10_14_GENE_FUNCTION_PROFILE.tsv, data/r11/results/R11_FOXP3_ENDPOINT.tsv, R11_SIGNATURE_CONTRASTS.tsv and data/r11/reference/R11_FIGURE_FUNCTION_ASSAYS.tsv","Source-specific screen calls, RNA/score estimates and verified assay outcomes","R11")
    fn=sup/"Supplementary_Figures_S1-S23.pdf"
    doc.save(fn,deflate=True,garbage=4,no_new_id=True);doc.close()
    (sup/"Supplementary_Figure_Rendered_Legends.md").write_text(legends.rstrip()+"\n")
    maps=read(r10/"R10_FIGURE_PANEL_MAP.tsv")
    changes={"3A":"3A; S15A","3B":"S15B","3C":"3B; S15C","3D":"3C; S15D","4A and S8":"S16A; S8","4B":"3D; S16B","4C":"S16C"}
    panelrows=[]
    for rec in maps.itertuples():
        identity=rec.revision_panel.split(" and legend")[0]
        new=changes.get(rec.revision_panel,identity)
        if new == "S1–S6":
            public_source = "figure_sources/supp/SOURCE_MAP.tsv"
        else:
            names = rec.source.replace("04_figures/source_data/","").replace("03_results/","").split(" + ")
            public_source = " + ".join("figure_sources/R10/" + name for name in names)
        panelrows.append([new,rec.revision_panel,rec.readout,public_source,rec.selection_and_fields])
    panelrows += [[f"4{panel}","New R11",title,sources,fields] for panel,title,sources,fields in [
        ("A","Foxp3 RNA samples","data/r11/results/R11_SAMPLE_EXPRESSION.tsv","gene_symbol=Foxp3, log2_norm_count_plus1"),
        ("B","Foxp3 RNA effects","data/r11/results/R11_FOXP3_ENDPOINT.tsv","unshrunken log2FC, CI_low/high, adjusted_p"),
        ("C","Hill scores per sample","data/r11/results/R11_SIGNATURE_SAMPLE_SCORES.tsv","Hill_LN, component=score for every library"),
        ("D","Hill score contrasts","data/r11/results/R11_SIGNATURE_CONTRASTS.tsv","Hill_LN, component=score; estimate, CI_low/high, adjusted_P")]]
    panelrows.append(["3C","R10","Property-adjusted continuous contrast",
                      "figure_sources/R10/R10_TREG_TCONV_ADJUSTED.tsv",
                      "analysis=continuous_delta_abs_rank_conditioned; P"])
    write(pd.DataFrame(panelrows+rows,columns=["R11_panel","earlier_panel_or_origin","readout","source","fields"]),figs/"SOURCE_PANEL_LINEAGE.tsv")
    with fitz.open(fn) as rendered:
        print(f"MAIN_FIGURES=4; SUPPLEMENT_PAGES={len(rendered)}")
