"""Plot public conditional-deletion RNA data and independently defined program results."""
from figure_common import *
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

STATE={"naive":"F1N","resting":"F2R","bulk":"S3"}
PCMARK={"F1N":"o","F1A":"s","S3":"o","F2R":"o","F2A":"s"}
COHORT_COLOR=[NAVY,TEAL,PURPLE]
HALLMARKS=[
    ("HALLMARK_IL2_STAT5_SIGNALING","IL-2 / STAT5 signaling"),("HALLMARK_E2F_TARGETS","E2F targets"),
    ("HALLMARK_G2M_CHECKPOINT","G2/M checkpoint"),("HALLMARK_APOPTOSIS","Apoptosis"),
    ("HALLMARK_TNFA_SIGNALING_VIA_NFKB","TNFα signaling via NF-κB"),
    ("HALLMARK_INTERFERON_GAMMA_RESPONSE","Interferon-γ response")]
MARKERS="Il2ra Ctla4 Ikzf2 Tnfrsf18 Tigit Entpd1 Lrrc32 Il10 Tgfb1 Satb1 Tcf7 Tet1 Tet3".split()

def sample_contrast(row):
    if row.state=="activated":
        return {"GSE121251":"F2A","GSE121279":"F1A"}[row.dataset]
    return STATE[row.state]

def prepared_samples(df):
    df=df.copy()
    df["contrast"]=df.apply(sample_contrast,axis=1)
    assert not df.duplicated(["GSM"]).any()
    return df

def sample_grid(fig,data,value,top,ytitle,shade=False):
    entries=[(.101,.136),(.287,.112),(.450,.136)]
    for i,((dataset,contrasts,title),(left,width)) in enumerate(zip(DESIGN,entries)):
        a=fig.add_axes([left,top-.299,width,.243])
        sub=data[data.dataset.eq(dataset)].copy()
        a.set_title(title,fontsize=6.95,pad=5.5,fontweight="bold",color=COHORT_COLOR[i])
        a.set_xlim(-.53,len(contrasts)-.47)
        a.set_xticks(np.arange(len(contrasts)),[SHORT[c].replace("Pooled Tregs","Tregs") for c in contrasts],fontsize=6.8)
        a.tick_params(axis="y",labelsize=6.65,pad=2)
        a.yaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=3))
        grid(a,"y")
        if i==0:a.set_ylabel(ytitle,labelpad=2,fontsize=7.4)
        finite=sub[value].dropna()
        if finite.empty:
            a.text(.5,.5,"Unavailable",transform=a.transAxes,ha="center",size=7,color=GREY)
            continue
        for j,contrast in enumerate(contrasts):
            group=sub[sub.contrast.eq(contrast)]
            for g,(gt,pos,color) in enumerate([("control",-.16,NAVY),("KO",.16,TEAL)]):
                samples=group[group.genotype.eq(gt)].sort_values("GSM")
                vals=samples[value].to_numpy()
                if not len(vals):continue
                xs=j+pos+np.linspace(-.058,.058,len(vals))
                for x,(_,row) in zip(xs,samples.iterrows()):
                    marker="s" if dataset=="GSE263830" and row.sex=="male" else "o"
                    a.scatter(x,row[value],s=17,color="white" if g==0 else color,edgecolor=color,lw=.85,marker=marker,zorder=3)
                a.hlines(np.nanmedian(vals),j+pos-.095,j+pos+.095,lw=1.1,color=color,zorder=2)
        pad=(finite.max()-finite.min())*.17
        pad=pad if pad>0 else abs(finite.max())*.05+.03
        a.set_ylim(finite.min()-pad,finite.max()+pad)
    if shade:
        handles=[Line2D([],[],marker="o",linestyle="",color=NAVY,mfc="white",ms=4,label="Control"),
                 Line2D([],[],marker="o",linestyle="",color=TEAL,mfc=TEAL,ms=4,label="Conditional KO")]
        fig.legend(handles=handles,ncol=2,loc="center",bbox_to_anchor=(.351,top+.036),fontsize=7,handletextpad=.2,columnspacing=1.0)

def forest_axis(fig,data,top,xlabel,col,adjust_col,left=.723,width=.156):
    ax=fig.add_axes([left,top-.305,width,.266])
    vals=data.set_index("contrast").reindex(CONTRASTS)
    ax.set_yticks(range(5),ROW_NAMES,fontsize=6.5);ax.set_ylim(4.55,-.52)
    ax.axhline(2.5,color=GREY,lw=.5,ls=(0,(1,3)))
    ax.axvline(0,color=GREY,lw=.75,ls="--")
    ax.spines["left"].set_visible(False);ax.tick_params(axis="y",length=0,pad=3)
    ax.tick_params(axis="x",labelsize=6.9);ax.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=3));grid(ax)
    for j,(key,row) in enumerate(vals.iterrows()):
        color=COHORT_COLOR[0 if j<2 else 1 if j==2 else 2]
        if pd.notna(row[col]) and pd.notna(row.CI_low):
            ax.plot([row.CI_low,row.CI_high],[j,j],lw=1.25,color=color,solid_capstyle="round")
            ax.scatter(row[col],j,color=color,marker="o" if j<3 else "D",s=22,zorder=3)
        y=top-.305+.266*(4.55-j)/(4.55+.52)
        fig.text(.943,y,ptext(row[adjust_col]),va="center",ha="center",size=7)
    fig.text(.943,top-.016,"Holm P",va="center",ha="center",weight="bold",fontsize=7)
    ax.set_xlabel(xlabel,fontsize=7.5)
    return ax

def endpoint_figure(root):
    """All five prespecified primary and auxiliary contrasts, separate units and adjustment families."""
    result=root/"03_results"
    expr=read(result/"R11_SAMPLE_EXPRESSION.tsv")
    expr=prepared_samples(expr[expr.gene_symbol.eq("Foxp3")])
    de=read(result/"R11_FOXP3_ENDPOINT.tsv")
    sig=read(result/"R11_SIGNATURE_CONTRASTS.tsv")
    scores=read(result/"R11_SIGNATURE_SAMPLE_SCORES.tsv")
    primary=sig[sig.program.eq("Hill_LN")&sig.component.eq("score")]
    hscores=prepared_samples(scores[scores.program.eq("Hill_LN")&scores.component.eq("score")])
    assert set(de.contrast)==set(primary.contrast)==set(CONTRASTS)
    assert len(expr)==32
    for rec in de.itertuples():
        sub=expr[expr.contrast.eq(rec.contrast)]
        assert (sum(sub.genotype.eq("KO")),sum(sub.genotype.eq("control")))==(rec.n_KO,rec.n_control)
    fig=plt.figure(figsize=(8.95,5.9))
    header(fig,.024,.948,"A","Foxp3 RNA per biological library")
    header(fig,.618,.948,"B","Foxp3 RNA contrasts")
    header(fig,.024,.455,"C","Hill Treg identity score per library")
    header(fig,.618,.455,"D","Hill Treg identity contrasts")
    sample_grid(fig,expr,"log2_norm_count_plus1",.867,"log₂(normalized count + 1)",shade=True)
    sample_grid(fig,hscores,"score",.380,"Signed rank-score difference")
    forest_axis(fig,de,.867,"KO − control (log₂ fold change)","log2FC","adjusted_p")
    forest_axis(fig,primary,.380,"KO − control (rank score)","estimate","adjusted_P")
    return fig

def effect_matrix(ax,data,index,columns,value,q,title,fmt=None):
    matrix=data.pivot(index=index,columns="contrast",values=value).reindex(index=columns,columns=CONTRASTS)
    adjusted=data.pivot(index=index,columns="contrast",values=q).reindex(index=columns,columns=CONTRASTS)
    finite=np.isfinite(matrix)
    bound=np.nanmax(np.abs(matrix)) if finite.any().any() else 1
    bound=max(bound,.001)
    image=ax.pcolormesh(np.arange(6)-.5,np.arange(len(columns)+1)-.5,np.ma.masked_invalid(matrix.to_numpy()),cmap=SIGNED,vmin=-bound,vmax=bound,edgecolors="white",linewidth=.35,shading="flat")
    labels=["F1N","F1A","S3","F2R","F2A"]
    ax.set_xticks(range(5),labels)
    ax.set_yticks(range(len(columns)),columns,style="italic")
    ax.set_ylim(len(columns)-.5,-.5); ax.set_xlim(-.5,4.5)
    ax.vlines(2.5,-.5,len(columns)-.5,color=DARK,lw=.8)
    ax.tick_params(length=0);ax.spines[["left","bottom"]].set_visible(False)
    for j in range(5):
        for i in range(len(columns)):
            v=matrix.iloc[i,j]; qval=adjusted.iloc[i,j]
            if pd.isna(v): text="×"
            elif pd.isna(qval): text="?"
            elif qval<.05: text="●"
            else: text="·"
            txcolor="white" if np.isfinite(v) and abs(v)>.67*bound else DARK
            ax.text(j,i,text,fontsize=6.5 if text=="●" else 9,ha="center",va="center",color=txcolor)
    ax.set_title(title,loc="left",weight="bold",pad=8)
    return image

def marker_figure(root):
    r=root/"03_results"; df=read(r/"R11_MARKER_EFFECTS.tsv")
    tary=read(r/"R11_TARGET_QC.tsv")
    fig=plt.figure(figsize=(7.4,4.6))
    a=fig.add_axes([.14,.15,.43,.72])
    im=effect_matrix(a,df,"gene_symbol",MARKERS,"log2FC","adjusted_p","A    Fixed marker effects")
    bar=fig.colorbar(im,cax=fig.add_axes([.598,.33,.018,.33]));bar.outline.set_visible(False);bar.ax.tick_params(labelsize=7)
    bar.ax.set_title("log₂FC",size=7.5,pad=7)
    b=fig.add_axes([.747,.26,.216,.48])
    tab=tary.set_index("contrast").loc[CONTRASTS]
    for j,(_,row) in enumerate(tab.iterrows()):
        color=COHORT_COLOR[0 if j<2 else 1 if j==2 else 2]
        if pd.notna(row.log2FC):
            b.plot([row.CI_low,row.CI_high],[j,j],color=color,lw=1.2)
            b.scatter(row.log2FC,j,color=color,s=23,zorder=3)
    b.axvline(0,lw=.7,color=GREY,ls="--")
    b.set_ylim(4.55,-.5);b.set_yticks(range(5),CONTRASTS);b.spines["left"].set_visible(False);b.tick_params(axis="y",length=0)
    b.set_xlabel("KO − control (log₂FC)");grid(b);b.xaxis.set_major_locator(MaxNLocator(4))
    fig.text(.696,.89,"B    Target RNA",weight="bold",size=9)
    fig.text(.696,.834,"Foxp1 or Sin3a",style="italic",size=8)
    fig.text(.14,.04,"●   BH q < 0.05      ·   BH q ≥ 0.05      ×   Unavailable      ?   q unavailable",size=7.5)
    return fig

def program_effects(root,programs,labels,title):
    data=read(root/"03_results"/"R11_SIGNATURE_CONTRASTS.tsv")
    data=data[data.component.eq("score") & data.program.isin(programs)]
    n=len(programs);cols=2 if n<=2 else 3;rows=(n+cols-1)//cols
    fig,axs=plt.subplots(rows,cols,figsize=(8,2.6*rows+.75),squeeze=False)
    fig.subplots_adjust(left=.11,right=.97,bottom=.16 if rows==1 else .12,top=.75 if rows==1 else .85,hspace=.75,wspace=.7)
    fig.text(.04,.94,title,weight="bold",size=11)
    for j,(prog,label) in enumerate(zip(programs,labels)):
        a=axs.flat[j]; subset=data[data.program.eq(prog)].set_index("contrast").reindex(CONTRASTS)
        a.set_title(label,weight="bold",loc="left",size=8.15,pad=7)
        a.set_yticks(range(5),CONTRASTS);a.set_ylim(4.5,-.5)
        a.axvline(0,color=GREY,lw=.75,ls="--");a.axhline(2.5,color=GREY,lw=.5,ls=(0,(1,3)))
        for i,(_,row) in enumerate(subset.iterrows()):
            color=COHORT_COLOR[0 if i<2 else 1 if i==2 else 2]
            if pd.notna(row.estimate):
                a.plot([row.CI_low,row.CI_high],[i,i],color=color,lw=1.3)
                a.scatter(row.estimate,i,s=23,facecolor=color if row.adjusted_P<.05 else "white",edgecolor=color,lw=1.1,zorder=3)
            else:
                a.annotate("NA",(0,i),va="center",size=7,color=GREY)
        grid(a);a.tick_params(axis="y",length=0);a.spines["left"].set_visible(False)
        a.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=3)); a.set_xlabel("KO − control (rank score)",fontsize=7.5)
    return fig

def qc_overview(root):
    p=read(root/"03_results"/"R11_QC_SAMPLE_PCA.tsv")
    size=read(root/"03_results"/"R11_QC_SAMPLES.tsv")
    overlap=[x for x in ["state","genotype","sex","GSM"] if x not in size.columns]
    size=size.merge(p[["dataset","sample_id"]+overlap],on=["dataset","sample_id"],validate="one_to_one")
    fig,axes=plt.subplots(2,3,figsize=(8,5.5));fig.subplots_adjust(left=.09,right=.978,bottom=.13,top=.81,wspace=.43,hspace=.65)
    for j,(study,contrasts,title) in enumerate(DESIGN):
        s=prepared_samples(p[p.dataset.eq(study)]); a=axes[0,j];b=axes[1,j]
        a.set_title(title,weight="bold",size=8,loc="left",pad=7)
        for row in s.itertuples():
            color=NAVY if row.genotype=="control" else TEAL
            shape=("s" if row.sex=="male" else "o") if study=="GSE263830" else PCMARK[row.contrast]
            a.scatter(row.PC1,row.PC2,facecolor="white" if row.genotype=="control" else color,edgecolor=color,s=30,marker=shape,lw=.95,zorder=3)
        a.set_xlabel(f"PC1 ({s.pc1_variance_fraction.iloc[0]:.1%})")
        a.set_ylabel(f"PC2 ({s.pc2_variance_fraction.iloc[0]:.1%})");grid(a,"both")
        lib=prepared_samples(size[size.dataset.eq(study)])
        lib=lib.sort_values(["contrast","genotype","GSM"],ascending=[True,False,True])
        for i,row in enumerate(lib.itertuples()):
            color=NAVY if row.genotype=="control" else TEAL
            b.bar(i,row.raw_library_size/1e6,color=color,alpha=1,width=.7)
        b.set_xticks(range(len(lib)),[x[-3:] for x in lib.GSM],rotation=70,size=6.2)
        b.set_ylabel("Read counts (millions)")
        b.set_xlabel("GSM suffix");grid(b,"y")
    fig.text(.028,.944,"A    Sample principal components",weight="bold",size=11)
    fig.text(.028,.470,"B    Library sizes",weight="bold",size=10)
    fig.legend(handles=[Line2D([],[],ls="",marker="o",mfc="white",mec=NAVY,ms=5,label="Control"),
                        Line2D([],[],ls="",marker="o",color=TEAL,ms=5,label="Conditional KO")],ncol=2,loc="upper right",bbox_to_anchor=(.98,.986),handletextpad=.3)
    return fig

def qc_distances(root):
    data=read(root/"03_results"/"R11_QC_SAMPLE_DISTANCES.tsv.gz")
    annot=read(root/"01_audit"/"R11_SAMPLE_MAP.tsv")
    fig,axes=plt.subplots(2,3,figsize=(8.55,6.45))
    fig.subplots_adjust(left=.075,right=.966,bottom=.095,top=.872,hspace=.41,wspace=.33)
    for j,(study,contrasts,title) in enumerate(DESIGN):
        df=data[data.dataset.eq(study)]
        meta=annot[(annot.dataset==study)&annot.sample_id.isin(df.sample_i)].sort_values(["state","genotype","GSM"],ascending=[True,False,True])
        ids=meta.sample_id.tolist();assert len(ids) in [8,12] and not meta.sample_id.duplicated().any()
        labs=[x[-3:] for x in meta.GSM]
        for i,(key,cmap) in enumerate([("pearson_correlation","viridis"),("euclidean_distance","magma_r")]):
            a=axes[i,j]
            matrix=df.pivot(index="sample_i",columns="sample_j",values=key).reindex(index=ids,columns=ids)
            matrix=matrix.combine_first(matrix.T)
            values=matrix.to_numpy(copy=True);np.fill_diagonal(values,1.0 if i==0 else 0.0)
            assert np.isfinite(values).all() and np.allclose(values,values.T)
            im=a.pcolormesh(np.arange(len(ids)+1)-.5,np.arange(len(ids)+1)-.5,values,cmap=cmap,shading="flat")
            a.set_aspect("equal");a.set_ylim(len(ids)-.5,-.5)
            a.set_xticks(range(len(ids)),labs,rotation=90,fontsize=5.9)
            a.set_yticks(range(len(ids)),labs,fontsize=5.9);a.tick_params(length=0,pad=2);a.spines[["left","bottom"]].set_visible(False)
            if i==0:a.set_title(title,fontsize=7.8,weight="bold",pad=5)
            cb=fig.colorbar(im,ax=a,shrink=.65,pad=.035,aspect=19);cb.ax.tick_params(labelsize=6.1);cb.outline.set_visible(False)
    fig.text(.025,.954,"A    Pearson correlation",weight="bold",size=10.5)
    fig.text(.025,.495,"B    Euclidean distance",weight="bold",size=10.5)
    return fig

def count_distribution(root):
    data=read(root/"03_results"/"R11_QC_SAMPLES.tsv")
    fig,axes=plt.subplots(1,3,figsize=(8.3,3.3));fig.subplots_adjust(left=.087,right=.98,top=.77,bottom=.23,wspace=.38)
    fig.text(.03,.93,"Count distributions across input genes",weight="bold",size=11)
    for j,(study,_,title) in enumerate(DESIGN):
        a=axes[j];df=data[data.dataset.eq(study)].sort_values(["state","genotype","GSM"],ascending=[True,False,True])
        for i,row in enumerate(df.itertuples()):
            color=NAVY if row.genotype=="control" else TEAL
            low,med,up,p95,p99=np.log2(np.array([row.raw_count_q25,row.raw_count_median,row.raw_count_q75,row.raw_count_q95,row.raw_count_q99])+1)
            a.vlines(i,low,p99,lw=.75,color=color)
            a.vlines(i,low,up,lw=3.5,color=color)
            a.plot(i,med,"_",ms=5,color=DARK,mew=.85)
            a.plot(i,p95,"o",ms=2.5,color=color)
        a.set_title(title,loc="left",weight="bold",size=8.1)
        a.set_xticks(range(len(df)),[x[-3:] for x in df.GSM],rotation=75,size=6.1)
        a.set_xlabel("GSM suffix")
        if j==0: a.set_ylabel("log₂(raw count + 1)")
        grid(a,"y")
    return fig

def evidence_graph(root,source,lit_rows):
    """Each block retains the evidence's own units and study-specific assay context."""
    profile=read(source/"R10_14_GENE_FUNCTION_PROFILE.tsv").set_index("gene").loc[["FOXP1","SIN3A"]]
    f=read(root/"03_results"/"R11_FOXP3_ENDPOINT.tsv").set_index("contrast").loc[CONTRASTS]
    allscores=read(root/"03_results"/"R11_SIGNATURE_CONTRASTS.tsv")
    h=allscores[allscores.program.eq("Hill_LN")&allscores.component.eq("score")].set_index("contrast").loc[CONTRASTS]
    fig=plt.figure(figsize=(8.15,5.8))
    header(fig,.026,.946,"A","Human CRISPR screen readouts")
    a=fig.add_axes([.176,.723,.781,.161]); a.axis("off"); a.set(xlim=(-.5,3.5),ylim=(1.55,-.87))
    columns=[("chen","iTreg FOXP3"),("treg","Primary Treg FOXP3"),("tconv","Tconv FOXP3"),("shifrut","TCR proliferation")]
    srcrows=[]
    for j,(key,title) in enumerate(columns):
        a.text(j,-.70,title,ha="center",va="center",fontsize=7.5,weight="bold")
        for i,(gene,row) in enumerate(profile.iterrows()):
            status=row[f"{key}_screen_status"]
            if not row[f"{key}_measured"]:label="Unmeasured";color=GREY
            elif not row[f"{key}_hit"]:label="Measured nonhit";color=GREY
            elif key=="shifrut":
                label="Lower · hit" if row.shifrut_direction=="NEG_TCR_PHENOTYPE" else "Higher · hit";color=TEAL
            else:
                val=row[f"{key}_signed_foxp3_bin_lfc"]
                label=("Lower" if val<0 else "Higher")+" · hit";color=TEAL
            a.text(j,i,label,ha="center",va="center",size=7.4,color=color)
            a.axhline(i+.45,color=GRID,lw=.5)
            srcrows.append([gene,"Human screen",title,label,status])
    for i,g in enumerate(profile.index):a.text(-.58,i,g,ha="right",va="center",style="italic",weight="bold")
    header(fig,.026,.634,"B","Mouse conditional deletions")
    b=fig.add_axes([.277,.314,.66,.252]);b.axis("off");b.set(xlim=(-.5,2.5),ylim=(4.45,-1.03))
    for j,lab in enumerate(["Foxp3 RNA, log₂FC","Hill score difference","Holm support"]):
        b.text(j,-.80,lab,ha="center",va="center",weight="bold",fontsize=7.6)
    for i,c in enumerate(CONTRASTS):
        row=f.loc[c];hr=h.loc[c]; rsign=row.adjusted_p<.05;hsign=hr.adjusted_P<.05
        norm="RNA & Hill" if rsign and hsign else "RNA only" if rsign else "Hill only" if hsign else "Neither"
        b.text(-.64,i,f"{c}   {ROW_NAMES[i]}",ha="right",va="center",fontsize=7.5)
        b.text(0,i,f"{row.log2FC:+.2f}",ha="center",va="center",color=NAVY,size=8)
        b.text(1,i,f"{hr.estimate:+.3f}",ha="center",va="center",color=TEAL,size=8)
        b.text(2,i,norm,ha="center",va="center",fontsize=7.5)
        b.axhline(i+.47,color=GRID,lw=.6)
        srcrows.append([row.target,"Mouse RNA-seq",c,row.log2FC,norm])
        srcrows.append([row.target,"Mouse Hill score",c,hr.estimate,norm])
    header(fig,.026,.222,"C","Published mouse Treg suppression assays")
    c=fig.add_axes([.055,.018,.90,.160]);c.axis("off");c.set(xlim=(0,1),ylim=(len(lit_rows)-.60,-.60))
    for i,rec in enumerate(lit_rows):
        c.text(.01,i,rec["target"],va="center",weight="bold",style="italic",size=7.7)
        c.text(.174,i,rec["figure"],va="center",size=7.2)
        c.text(.53,i,rec["readout"],va="center",size=7.2,color=DARK)
        c.axhline(i+.42,color=GRID,lw=.6)
        srcrows.append([rec["target"],"Published functional assay",rec["figure"],rec["readout"],rec["PMID"]])
    write(pd.DataFrame(srcrows,columns=["gene","system","assay_or_contrast","plotted_value_or_label","source_status_or_PMID"]),root/"04_figures"/"source_data"/"figure_4_candidate_evidence_plot.tsv")
    return fig
