"""Independently recompute key plotted point/interval values and test plotted arrays."""
from scipy import stats
from matplotlib.collections import PathCollection
from figure_common import *
from knockout_figures import endpoint_figure,prepared_samples,HALLMARKS,MARKERS

def audit(stage,data):
    res=stage/"03_results";src=data/"figure_sources/R10"
    qa=[]
    def check(panel,description,test,count):
        if not test: raise ValueError(panel+" "+description)
        qa.append([panel,description,"PASS",count])
    c=read(src/"F3C_hit_classes.tsv").set_index("category")
    p=read(src/"F3C_discordant_fisher.tsv").iloc[0]
    cells=c.loc[["Treg-only","Tconv-only"],["asd_n","non_asd_n"]].to_numpy().T
    ft=stats.fisher_exact(cells); odd=stats.contingency.odds_ratio(cells);ci=odd.confidence_interval()
    check("3B","Two-sided Fisher P, conditional OR and interval recalculated from discordant counts",
        np.allclose([ft.pvalue,odd.statistic,ci.low,ci.high],[p.P,p.OR,p.CI_low,p.CI_high],rtol=5e-6),int(cells.sum()))
    r=read(src/"F3D_rank_difference.tsv.gz")
    e=read(src/"F3D_continuous_contrast.tsv").iloc[0]
    ranks=r.groupby("asd").delta_abs_rank.mean()
    check("3C","ECDF source indices reproduce the ASD-minus-other mean rank difference",np.isclose(ranks[1]-ranks[0],e.effect),len(r))
    ref=read(src/"F4B_breadth.tsv").set_index("gene").loc[PANEL]
    human=read(src/"F4A_human_expression.tsv");human=human[human.gene.isin(PANEL)]
    calc=human.groupby("gene").log2_tpm_plus_1.apply(lambda x:(x>=1).sum()).reindex(PANEL)
    mouse=read(src/"R10_MOUSE_EXPRESSION_PLOT.tsv")
    mouse=mouse[mouse.human_gene.isin(PANEL)] if "human_gene" in mouse else mouse[mouse.gene.isin(PANEL)]
    gene="human_gene" if "human_gene" in mouse else "gene"
    mcalc=mouse.groupby(gene).log2_cpm1.apply(lambda x:(x>=1).sum()).reindex(PANEL)
    check("3D","Human ≥1 TPM/mouse ≥1 TMM CPM detection recalculated from per-cell/gene abundance",np.array_equal(calc,ref.human_k) and np.array_equal(mcalc,ref.mouse_k),len(ref))
    ex=read(res/"R11_SAMPLE_EXPRESSION.tsv")
    ex=prepared_samples(ex[ex.gene_symbol.eq("Foxp3")])
    library=read(res/"R11_QC_SAMPLES.tsv").set_index("GSM")
    norm=ex.raw_count.to_numpy()/library.loc[ex.GSM,"size_factor"].to_numpy()
    check("4A","Library-specific raw counts/DESeq2 size factors independently reproduce sample y coordinates",np.allclose(np.log2(1+norm),ex.log2_norm_count_plus1,rtol=1e-12),len(ex))
    de=read(res/"R11_FOXP3_ENDPOINT.tsv").set_index("contrast").loc[CONTRASTS]
    z=stats.norm.ppf(.975)
    check("4B","Unshrunken estimates and standard errors reproduce every plotted Wald interval",
        np.allclose(de.log2FC-z*de.SE,de.CI_low) and np.allclose(de.log2FC+z*de.SE,de.CI_high),len(de))
    sig=read(res/"R11_SIGNATURE_CONTRASTS.tsv");score=read(res/"R11_SIGNATURE_SAMPLE_SCORES.tsv")
    values=sig[sig.component.eq("score")]
    groups=values.groupby(["program","dataset"])
    effects=[]
    for (program,dataset),subset in groups:
        v=score[score.program.eq(program)&score.dataset.eq(dataset)&score.component.eq("score")]
        y=v.score.to_numpy()
        geno=v.genotype.eq("KO").astype(float).to_numpy()
        if dataset=="GSE263830":
            X=np.column_stack([np.ones(len(v)),v.sex.eq("male").astype(float),geno])
            vec={"S3":[0,0,1]}
        else:
            act=v.state.eq("activated").astype(float).to_numpy()
            X=np.column_stack([np.ones(len(v)),act,geno,act*geno])
            keys=["F1N","F1A"] if dataset=="GSE121279" else ["F2R","F2A"]
            vec={keys[0]:[0,0,1,0], keys[1]:[0,0,1,1]}
        beta=np.linalg.lstsq(X,y,rcond=None)[0]
        df=len(y)-np.linalg.matrix_rank(X); cov=np.linalg.inv(X.T@X)*np.sum((y-X@beta)**2)/df
        for cont,vector in vec.items():
            t=stats.t.ppf(.975,df);v=np.array(vector);b=v@beta;se=np.sqrt(v@cov@v)
            row=subset.set_index("contrast").loc[cont]
            effects.append(np.allclose([b,b-t*se,b+t*se],[row.estimate,row.CI_low,row.CI_high],atol=2e-12,rtol=1e-12))
    check("4C/D; S21/S22","All eight programs, model contrasts and t intervals independently refitted from library scores",all(effects),len(effects))
    fig=endpoint_figure(stage)
    # Source order is explicit in the renderer. Check the actual Matplotlib artists before export.
    for ax_range,frame,value,lab in [(range(3),ex,"log2_norm_count_plus1","4A"),(range(3,6),score[score.program.eq("Hill_LN")&score.component.eq("score")],"score","4C")]:
        for aix,dataset in zip(ax_range,STUDIES):
            observed=np.sort(np.concatenate([np.asarray(line.get_offsets()[:,1],dtype=float) for line in fig.axes[aix].collections if isinstance(line,PathCollection)]))
            target=np.sort(frame[frame.dataset.eq(dataset)][value])
            check(lab,f"{dataset} plotted scatter coordinates matched to every biological library",np.allclose(observed,target,atol=1e-12),len(target))
    for axidx,value,data,name in [(6,"log2FC",de,"4B"),(7,"estimate",values[values.program.eq("Hill_LN")].set_index("contrast").loc[CONTRASTS],"4D")]:
        points=np.array([o.get_offsets()[0] for o in fig.axes[axidx].collections if isinstance(o,PathCollection)],dtype=float)
        ci=np.array([l.get_xdata() for l in fig.axes[axidx].lines[-5:]])
        check(name,"Actual plotted five forest coordinates, order and 95% interval endpoints",np.allclose(points[:,0],data[value]) and np.allclose(points[:,1],range(5)) and np.allclose(ci,data[["CI_low","CI_high"]]),len(data))
    plt.close(fig)
    marks=read(res/"R11_MARKER_EFFECTS.tsv").dropna(subset=["pvalue"])
    adjusted=stats.false_discovery_control(marks.pvalue.to_numpy())
    check("S20","Row families and BH dot classifications rederived from the raw marker P values",set(marks.gene_symbol)==set(MARKERS) and np.allclose(adjusted,marks.adjusted_p),len(marks))
    plain=read(stage/"04_figures"/"source_data"/"figure_4_candidate_evidence_plot.tsv")
    check("S23A","Eight targeted gene/readout combinations, including two unmeasured assays and one measured nonhit",
           len(plain[plain.system.eq("Human screen")])==8 and sum(plain.plotted_value_or_label.eq("Unmeasured"))==2 and sum(plain.plotted_value_or_label.eq("Measured nonhit"))==1,8)
    render = read(stage/"07_qa/FIGURE_RENDER_QA.tsv")
    check("1/2; S1–S14","Frozen SHA-256 and fresh page pixel preservation",
          render.iloc[:14].identical_to_fixed_vector_source.eq(True).all() and
          render.iloc[23:27].identical_to_fixed_vector_source.eq(True).all(),16)
    qa_path=stage/"07_qa"/"FIGURE_NUMERIC_QA.tsv"
    write(pd.DataFrame(qa,columns=["panels","independent_check","status","units_or_values_checked"]),qa_path)
    print(f"{len(qa)} independent numeric checks passed.")
