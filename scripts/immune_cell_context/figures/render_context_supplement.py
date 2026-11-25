"""Compose vector supplementary panels with an English justified legend under every figure."""
from context_figure_common import *
import fitz
import html
import io
import re

def figure_pdf(fig):
    buf=io.BytesIO();fig.savefig(buf,format="pdf",metadata={"CreationDate":None,"ModDate":None})
    plt.close(fig)
    return fitz.open(stream=buf.getvalue(),filetype="pdf")
def legends(data):
    path = data / "figure_sources/frozen_r10/R10_SUPPLEMENTARY_FIGURE_LEGENDS_EN.md"
    result = {}
    for section in re.split(r"\n## ", path.read_text(encoding="utf-8"))[1:]:
        heading, body = section.strip().split("\n", 1)
        number = int(re.match(r"Supplementary Figure S(\d+)", heading).group(1))
        if number <= 7:
            result[number] = (heading.split(". ", 1)[1], body.strip())
    return result


def ascii_sup(s):
    # Superscript digits in the R9 captions are retained as semantic superscripts in PDF.
    table=str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺","0123456789−+")
    return re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+",lambda m:"<sup>"+m.group().translate(table)+"</sup>",s)
def para_html(text):
    s=html.escape(text)
    s=re.sub(r"\*\*(.*?)\*\*",r"<b>\1</b>",s)
    s=re.sub(r"₁₀","<sub>10</sub>",s)
    return ascii_sup(s)
def caption_page(out, graphic, page_num, n, title, body):
    left=23; width=595; right=width-left
    available=right-left
    orig=graphic[page_num].rect
    art_width=available; art_h=orig.height*art_width/orig.width
    css="""* {font-family:sans-serif; color:#202d39;} body {margin:0; padding:0;}
    p {font-size:9pt; line-height:1.37; margin:0 0 5pt; text-align:justify;}
    h3 {font-size:10pt; line-height:1.33; margin:0 0 5pt; font-weight:bold;}"""
    paragraphs="".join("<p>"+para_html(p)+"</p>" for p in body.split("\n\n") if p.strip())
    label=f"Supplementary Figure S{n}. {title.rstrip('.')}."
    content="<h3>"+html.escape(label)+"</h3>"+paragraphs
    story=fitz.Story(content,user_css=css)
    more, occupied=story.place(fitz.Rect(0,0,available,2000));assert not more
    cap_h=occupied[3]-occupied[1]+9
    top=16; chart_end=top+art_h;cap_y=chart_end+12
    p=out.new_page(width=width,height=cap_y+cap_h+18)
    p.show_pdf_page(fitz.Rect(left,top,right,chart_end),graphic,page_num)
    p.draw_line((left,cap_y-6),(right,cap_y-6),color=(.83,.86,.86),width=.5)
    slack,scale=p.insert_htmlbox(fitz.Rect(left,cap_y,right,cap_y+cap_h),content,css=css,scale_low=1)
    assert slack>=0 and scale==1,(n,slack,scale)
    text=p.get_text(); assert f"Figure S{n}." in text and len(text)>len(body)*.75,(n,len(text),len(body))

def heatmap(ax,frame,cmap,bounds,role_axis="y",font=7):
    im=vector_cells(ax,frame,cmap,*bounds)
    ax.set_xticks(np.arange(frame.shape[1]),frame.columns,rotation=63,ha="right",rotation_mode="anchor",size=font)
    ax.set_yticks(np.arange(frame.shape[0]),frame.index,size=font,style="italic" if role_axis=="y" else "normal")
    ax.tick_params(length=0,pad=3)
    if role_axis=="y":ax.axhline(13.5,color="white",lw=1.8)
    if role_axis=="x":
        ax.axvline(13.5,color="white",lw=1.8)
        for t in ax.get_xticklabels():t.set_style("italic")
    ax.spines[["bottom","left"]].set_visible(False)
    return im

def human_heat(stage):
    d=read(stage/"03_results"/"R10_MONACO_EXPRESSION_HEATMAP_SOURCE.tsv")
    types=d.sort_values("cell_type_order").cell_type.unique()
    z=d.pivot(index="gene",columns="cell_type",values="row_z").reindex(index=PANEL,columns=types)
    raw=d.pivot(index="gene",columns="cell_type",values="log2_tpm_plus_1").reindex(index=PANEL,columns=types)
    fig=plt.figure(figsize=(8.15,8.9))
    for vals,cmap,row,label,y in [(raw,SEQ,"A","Absolute expression",.565),(z,DIV,"B","Within-gene relative expression",.115)]:
        header(fig,.03,y+.37,row,label)
        a=fig.add_axes([.12,y,.80,.337])
        vals=vals.rename(columns=lambda c:str(c).replace("_"," "))
        bound=(0,np.nanmax(vals)) if row=="A" else (-np.nanmax(abs(vals)),np.nanmax(abs(vals)))
        im=heatmap(a,vals,cmap,bound,font=6.95)
        cb=fig.colorbar(im,cax=fig.add_axes([.941,y+.04,.013,.21]));cb.outline.set_visible(False);cb.ax.tick_params(labelsize=7)
        for t in a.get_xticklabels():
            if t.get_text()=="Treg":t.set_weight("bold");t.set_color(TREG)
    return figure_pdf(fig)

def mouse_heat(stage,field):
    d=read(stage/"03_results"/"R10_MOUSE_EXPRESSION_PLOT.tsv").sort_values("group_order")
    groups=d.cell_type.unique()
    symbols=d.drop_duplicates("human_gene").set_index("human_gene").loc[PANEL,"mouse_gene"].tolist()
    mat=d.pivot(index="cell_type",columns="mouse_gene",values=field).reindex(index=groups,columns=symbols)
    lin=d.drop_duplicates("cell_type").set_index("cell_type").loc[groups,"lineage"]
    fig=plt.figure(figsize=(7.52,10.3))
    a=fig.add_axes([.379,.085,.524,.855])
    maxabs=np.nanmax(abs(mat)); signed=field=="zscore"
    bound=(-maxabs if signed else 0,maxabs)
    im=heatmap(a,mat,DIV if signed else SEQ,bound,role_axis="x",font=7.15)
    for i in range(1,len(lin)):
        if lin.iloc[i]!=lin.iloc[i-1]: a.axhline(i-.5,color="white",lw=1.05)
    c=fig.colorbar(im,cax=fig.add_axes([.935,.42,.013,.20]));c.outline.set_visible(False);c.ax.tick_params(labelsize=7)
    fig.text(.06,.970,"Mouse immune populations",fontsize=10,weight="bold")
    fig.text(.655,.970,"Gene-wise z score" if signed else "log₂(TMM CPM + 1)",fontsize=9)
    return figure_pdf(fig)

def donor_panel(stage):
    d=read(stage/"03_results"/"R10_EXPRESSION_DONOR_CONTRASTS.tsv")
    d=d[d.contrast_set.eq("primary_balanced_six") & d.gene.isin(PANEL)].copy()
    donors=sorted(d.donor.unique()); assert len(donors)==4 and len(d)==len(PANEL)*4
    fig=plt.figure(figsize=(8.0,8.8))
    header(fig,.03,.955,"A","Paired expression in four donors")
    header(fig,.53,.955,"B","Within-donor expression contrast")
    a=fig.add_axes([.15,.095,.32,.802]);b=fig.add_axes([.582,.095,.356,.802])
    for ax in (a,b):
        ax.set_ylim(len(PANEL)-.43,-.55);grid(ax);ax.axhline(13.5,lw=.75,color="#A7AFAF")
        ax.spines["left"].set_visible(False);ax.tick_params(axis="y",length=0)
    a.set_yticks(range(len(PANEL)),PANEL,style="italic"); b.set_yticks(range(len(PANEL)),[""]*len(PANEL))
    colors=["#183D60","#49A693","#B86A36","#755686"]; markers=["o","s","^","D"]
    for j,donor in enumerate(donors):
        part=d[d.donor.eq(donor)].set_index("gene").loc[PANEL]
        ys=np.arange(len(PANEL))+(j-1.5)*.155
        xt=np.log2(1+part.treg_tpm);xc=np.log2(1+part.tconv_equal_type_tpm)
        a.hlines(ys,np.minimum(xt,xc),np.maximum(xt,xc),lw=.6,color=GRAY)
        a.scatter(xc,ys,edgecolor=GRAY,color="white",s=9,lw=.7,zorder=3)
        a.scatter(xt,ys,color=TREG,s=8,lw=0,zorder=4)
        b.scatter(part.log2_ratio,ys,color=colors[j],marker=markers[j],s=15,lw=0,label=str(donor))
    a.set_xlabel("log₂(TPM + 1)");b.set_xlabel("log₂[(Treg TPM + 1) / (Tconv TPM + 1)]")
    b.axvline(0,color=GRAY,lw=.9,ls="--")
    a.legend(handles=[mark(TREG,"Treg"),mark(GRAY,"Tconv","o","white")],ncol=2,loc="lower left",bbox_to_anchor=(0,1.005),handletextpad=.15)
    b.legend(ncol=4,loc="lower left",bbox_to_anchor=(-.14,1.005),handletextpad=.05,columnspacing=.35,fontsize=7.2)
    return figure_pdf(fig)

def protein_panel(stage,compartment):
    d=read(stage/"03_results"/"R10_PERTURB_PROTEIN.tsv")
    all_=d[d.gene.isin(GENES)&d.score.notna()].copy()
    g=[g for g in GENES if g in set(all_.gene)]
    sel=all_[all_.compartment.eq(compartment)]
    order=list(sel.drop_duplicates("assay_label").sort_values("source_row").protein)
    m=sel.pivot(index="protein",columns="gene",values="score").reindex(index=order,columns=g)
    assert m.notna().all().all()
    limits=np.max(np.abs(all_.score))
    fig=plt.figure(figsize=(7.85,6.95))
    split=np.array_split(np.arange(len(m)),2) if compartment=="surface" else [np.arange(len(m))]
    for i,index in enumerate(split):
        place=([.15+i*.44,.10,.27,.76] if compartment=="surface" else [.37,.10,.43,.76])
        a=fig.add_axes(place); mat=m.iloc[index]
        im=heatmap(a,mat,DIV,(-limits,limits),role_axis=None,font=7.2)
        for t in a.get_xticklabels():t.set_style("italic")
    fig.text(.055,.935,"Intracellular proteins" if compartment=="intra" else "Surface proteins",weight="bold",size=11)
    c=fig.colorbar(im,cax=fig.add_axes([.905,.35,.014,.35]));c.outline.set_visible(False);c.ax.tick_params(labelsize=7)
    fig.text(.873,.73,"Score",size=8.2)
    return figure_pdf(fig)

def treg_ranks(stage):
    out=stage/"03_results"
    human=read(out/"R10_HUMAN_EXPRESSION_PANEL.tsv").set_index("gene").loc[PANEL]
    mouse=read(out/"R10_MOUSE_EXPRESSION_PANEL.tsv").set_index("human_gene").loc[PANEL]
    fig=plt.figure(figsize=(8,5.65))
    header(fig,.03,.95,"A","Human blood Tregs")
    header(fig,.515,.95,"B","Mouse splenic and colonic Tregs")
    a=fig.add_axes([.140,.13,.32,.71]); b=fig.add_axes([.630,.13,.319,.71])
    for ax,labels in [(a,PANEL),(b,mouse.mouse_gene)]:
        ax.set_yticks(range(len(PANEL)),labels,style="italic",size=7.2)
        ax.set_ylim(len(PANEL)-.5,-.55); ax.tick_params(axis="y",length=0)
        ax.axhline(13.5,color=GRID,lw=.8); ax.spines["left"].set_visible(False);grid(ax)
    ys=np.arange(len(PANEL)); a.scatter(human.treg_rank_desc,ys,s=18,color=TREG,lw=0,zorder=3)
    a.set_xlim(0,30);a.set_xticks([1,7,14,21,29]);a.set_xlabel("Rank among 29 cell types")
    b.scatter(mouse.splenic_treg_rank,ys-.13,color=ASD,s=16,lw=0,label="Spleen",zorder=3)
    b.scatter(mouse.colonic_treg_rank,ys+.13,color="white",edgecolor=SPECIES,marker="s",lw=.8,s=17,label="Colon",zorder=3)
    b.set_xlim(-1,58);b.set_xticks([1,14,28,42,56]);b.set_xlabel("Rank among 56 cell groups")
    b.legend(loc="lower left",bbox_to_anchor=(0,1.005),ncol=2,handletextpad=.15,columnspacing=.65)
    return figure_pdf(fig)

def added_legends(stage):
    out=stage/"03_results"
    h=read(out/"R10_HUMAN_EXPRESSION_PANEL.tsv"); m=read(out/"R10_MOUSE_EXPRESSION_PANEL.tsv")
    mp=read(out/"R10_MOUSE_EXPRESSION_PLOT.tsv")
    ncell=int(m.n_groups.dropna().iloc[0]); nl=int(mp.drop_duplicates("cell_type").n_libraries.sum())
    pro=read(out/"R10_PERTURB_PROTEIN.tsv"); n=pro.loc[pro.score.notna(),"gene"].nunique()
    return {
     14:("Within-atlas ranks of Treg expression",
         "(A) Each point gives the rank of Treg donor-median expression among the 29 human blood immune cell types in Monaco et al. (GSE107011). (B) Points show separate ranks for splenic and colonic Tregs across the 56 baseline cell groups in ImmGen (GSE109125), using the median TMM-normalized CPM in each group. Rank 1 indicates the highest measured expression of that gene within the respective atlas. Tied values receive average ranks. The reference genes follow the horizontal rule. Mouse symbols reflect one-to-one orthology. The two panels keep each species on its own observed rank scale."),
     8:("Human immune cell expression of the 14 discovery hits and four reference genes",
        "Monaco et al. (GSE107011) profiled sorted peripheral blood cells from four donors. Columns show all 29 cell types and rows show the 14 ASD discovery hits followed by FOXP3, IL2RA, CTLA4 and IKZF2. A white line separates candidates from reference genes. (A) Colors represent log₂(TPM + 1), calculated after taking the donor median within each cell type. Each terminal effector CD4 estimate uses two donors. The other 28 types each have four donors. (B) Expression is standardized across the 29 cell types within each gene on the log₂(TPM + 1) scale. The color bars give the complete observed ranges. PBMC mixtures are excluded. The accompanying cell-type crosswalk gives the original labels, donor coverage and broader lineage."),
     9:("Expression of mouse orthologs across baseline immune populations",
        f"Each row represents one of {ncell} mature cell populations from {nl} RNA-sequencing libraries in ImmGen ULI (GSE109125). GEO labels retain the tissue and gate distinctions, including separate splenic and colonic Treg populations. White horizontal lines mark lineage boundaries. Mouse gene symbols follow verified one-to-one human orthology. The first 14 columns are candidates and the last four are reference genes. Colors show log₂(CPM + 1) after TMM normalization across the count libraries and within-group median aggregation. The count-based units reflect relative RNA abundance and span different tissues. The source crosswalk expands all group labels and identifies replicate libraries."),
     10:("Relative expression of mouse orthologs across baseline immune populations",
         f"The same {ncell} cell populations, {nl} libraries, candidate orthologs and reference genes as Supplementary Figure S9 are shown. For each gene, log₂(TMM-normalized CPM + 1) was centered and scaled across eligible groups. The continuous color scale spans the observed gene-wise z scores. White horizontal lines separate lineages. The splenic and colonic Treg groups retain their original tissue labels. Between-population differences in gate, tissue and replicate count constrain the comparison."),
     11:("Paired Treg and conventional CD4 expression contrasts in human donors",
         "Four donors in GSE107011 each had RNA measurements for a Treg population and six prespecified conventional CD4 subsets. Each donor's conventional CD4 value gives equal weight to naive CD4, follicular helper, Th1, Th1/Th17, Th17 and Th2 TPM. (A) For every gene and donor, the line joins conventional CD4 (open circle) to Treg (filled circle) on the log₂(TPM + 1) scale. The four positions within each row correspond to individual donors. (B) The symbols show the matched Treg-minus-conventional-CD4 difference on this log scale. Positive values indicate higher Treg expression. The horizontal rule separates the 14 discovery hits from the four reference genes."),
     12:("Intracellular protein scores for discovery hits in the Chen follow-up",
         f"All 27 intracellular antibody measurements are shown for the {n} discovery genes with protein scores in Chen et al., Supplementary Table 2, ProteinScore_Targets. For each target and antibody, the authors calculated the mean protein measurement in targeted cells minus the mean in non-targeting controls, divided by the control mean. Positive values indicate a higher protein score after perturbation. The continuous scale is centered at zero and shared with Supplementary Figure S13 across the entire candidate-by-antibody panel. The original feature labels, antibody identifiers, target coverage and full numerical values accompany the figure. Unscored targets are documented in the coverage table."),
     13:("Surface protein scores for discovery hits in the Chen follow-up",
         f"All 61 surface antibody measurements are shown for the same {n} discovery genes as in Supplementary Figure S12. Rows retain the source-table order across both panels. Values are author-reported scores relative to the mean non-targeting control measurement, scaled by the control mean. The zero-centered scale spans all candidate intracellular and surface scores in this experiment. The target-by-antibody values are pooled gene-level summaries. No donor-level intervals are available for these scores.")
    }

def assemble(stage, data):
    titles = {**legends(data), **added_legends(stage)}
    if set(titles) != set(range(1, 15)):
        raise ValueError("Missing historical R10 supplementary legends")
    factories = {
        8: lambda: human_heat(stage),
        9: lambda: mouse_heat(stage, "log2_cpm1"),
        10: lambda: mouse_heat(stage, "zscore"),
        11: lambda: donor_panel(stage),
        12: lambda: protein_panel(stage, "intra"),
        13: lambda: protein_panel(stage, "surface"),
        14: lambda: treg_ranks(stage),
    }
    directory = stage / "05_supplementary"
    directory.mkdir(parents=True, exist_ok=True)
    fn = directory / "R10_Supplementary_Figures.pdf"
    frozen = data / "figure_sources/frozen_r10/R10_Supplementary_Figures.pdf"
    with fitz.open() as out:
        with fitz.open(frozen) as previous:
            out.insert_pdf(previous, from_page=0, to_page=6)
        for n, render in factories.items():
            title, body = titles[n]
            with render() as pdf:
                caption_page(out, pdf, 0, n, title, body)
        out.save(fn, deflate=True, garbage=4, no_new_id=True)
    legend_export = [f"## Supplementary Figure S{n}. {title}\n\n{body}"
                     for n, (title, body) in sorted(titles.items())]
    (directory / "R10_SUPPLEMENTARY_FIGURE_LEGENDS_EN.md").write_text(
        "# Supplementary figure legends\n\n" + "\n\n".join(legend_export) + "\n", encoding="utf-8")
    return fn


if __name__ == "__main__":
    opt = args()
    assemble(opt.stage, opt.data)
