"""Render four journal figures from the frozen archive and R10 result tables."""
from context_figure_common import *
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import Normalize

def proportions(ax, df, primary_color, xmax):
    vals=df[df.row_type.eq("proportion")]
    for i,r in enumerate(vals.itertuples()):
        frac=100*r.numerator/r.denominator
        y=1-i
        ax.barh(y,frac,height=.42,color=[primary_color,GRAY][i],zorder=3)
        ax.text(frac+.025*xmax,y,f"{int(r.numerator):,}/{int(r.denominator):,}",va="center",size=7.5)
    ax.set_yticks([1,0],["ASD","Other"]); ax.tick_params(axis="y",length=0)
    ax.set_xlim(0,xmax);ax.set_ylim(-.58,1.58);ax.set_xlabel("Genes called as hits (%)")
    ax.spines["left"].set_visible(False)
    grid(ax)

def f1(stage,src):
    data=read(src/"F1_chen.tsv.gz")
    f=coords(data,"screen_gene")
    asd,hit=truth(f.sfari_score1),truth(f.screen_hit)
    found=f[asd & hit].sort_values("author_lfc")
    assert len(f)==18966 and asd.sum()==242 and set(found.gene)==set(GENES)
    d=coords(data,"enrichment"); effect=d[d.row_type.eq("effect")].iloc[0]
    fig=plt.figure(figsize=(7.8,4.85))
    header(fig,.035,.94,"A","Genome-wide iTreg FOXP3 screen")
    a=fig.add_axes([.096,.445,.43,.425]); grid(a,"y")
    a.scatter(f.loc[~asd,"author_lfc"],f.loc[~asd,"screen_strength"],s=4,linewidths=0,color="#A3ACB5",alpha=.45)
    a.scatter(f.loc[asd & ~hit,"author_lfc"],f.loc[asd & ~hit,"screen_strength"],s=11,linewidths=.55,edgecolors=ASD,facecolors="none")
    a.axvline(0,color="#919AA2",lw=.6,ls="--")
    for direction,color,marker in [("SUPPORTS_FOXP3",TREG,"v"),("RESTRAINS_FOXP3",ORANGE,"^")]:
        g=found[found.regulator_direction.eq(direction)]
        a.scatter(g.author_lfc,g.screen_strength,marker=marker,s=27,color=color,edgecolor="white",lw=.35,zorder=6)
    a.set_ylabel("−log₁₀ (smaller raw tail P)");a.set_xlabel("Native log fold change, FOXP3-low vs high")
    a.legend(handles=[mark(GRAY,"Other genes"),mark(ASD,"ASD nonhit",mfc="none")],loc="upper right",handletextpad=.2,borderaxespad=.1)
    header(fig,.035,.284,"B","Hit rates in the measured universe")
    b=fig.add_axes([.112,.108,.40,.125]); proportions(b,d,ASD,8.3); b.set_xticks([0,2,4,6,8])
    fig.text(.113,.011,f"OR {effect.effect:.2f}  (95% CI {effect.effect_low:.2f}–{effect.effect_high:.2f})     P = {pval(effect.p)}",size=7.5)
    header(fig,.566,.94,"C","Effects of the 14 ASD hits")
    c=fig.add_axes([.698,.158,.275,.705]);grid(c)
    for i,r in enumerate(found.itertuples()):
        color,marker=(TREG,"v") if r.regulator_direction=="SUPPORTS_FOXP3" else (ORANGE,"^")
        c.hlines(i,0,r.author_lfc,color=color,lw=1.3)
        c.scatter(r.author_lfc,i,s=31,marker=marker,color=color,zorder=3)
    c.axvline(0,color="#728087",lw=.75)
    c.set_yticks(range(len(found)),found.gene,style="italic"); c.tick_params(axis="y",length=0,pad=4)
    c.set_ylim(len(found)-.35,-.7);c.spines["left"].set_visible(False);c.set_xlabel("Native log fold change")
    fig.legend(handles=[mark(TREG,"Lower FOXP3 after knockout","v"),mark(ORANGE,"Higher FOXP3 after knockout","^")],
               loc="lower left",bbox_to_anchor=(.604,.009),ncol=1,handletextpad=.25,labelspacing=.3)
    save_main(fig,stage,1)

def f2(stage,src):
    f=read(src/"F2_cross_screen.tsv")
    cohort=coords(f,"common_asd_gene")
    assert len(cohort)==33
    genes=cohort.loc[truth(cohort.chen_hit)|truth(cohort.umhoefer_treg_hit)].sort_values(["cross_screen_category","gene"])
    assert len(genes)==5
    fig=plt.figure(figsize=(7.85,6.05))
    header(fig,.035,.956,"A","Primary Treg FOXP3 hits")
    a=fig.add_axes([.122,.707,.327,.155])
    props=coords(f,"enrichment"); proportions(a,props,TREG,23)
    a.set_xticks([0,5,10,15,20])
    e=props[props.row_type.eq("effect")].iloc[0]
    fig.text(.107,.614,f"OR {e.effect:.2f}  (95% CI {e.effect_low:.2f}–{e.effect_high:.2f})",size=8)
    fig.text(.107,.585,f"One-sided Fisher P = {pval(e.p)}",size=7.6)
    header(fig,.501,.956,"B","Hits among 33 common ASD genes")
    b=fig.add_axes([.655,.65,.255,.219]); b.set_ylim(4.65,-.65); b.set_xlim(-.45,1.45)
    for i,r in enumerate(genes.itertuples()):
        for j,(hit,direction) in enumerate([(r.chen_hit,r.chen_regulator_direction),(r.umhoefer_treg_hit,r.umhoefer_treg_direction)]):
            marker,col=("v",TREG) if direction=="SUPPORTS_FOXP3" else ("^",ORANGE)
            if hit=="NO":marker,col="o",GRAY
            b.scatter(j,i,marker=marker,s=35,edgecolor=col,linewidth=.9,facecolor="white" if hit=="NO" else col,zorder=3)
    b.set_xticks([0,1],["iTreg","Primary Treg"]);b.xaxis.tick_top(); b.set_yticks(range(len(genes)),genes.gene,style="italic")
    for sp in b.spines.values():sp.set_visible(False)
    b.tick_params(length=0);grid(b,"y")
    fig.legend(handles=[mark(TREG,"Lower","v"),mark(ORANGE,"Higher","^"),mark(GRAY,"Nonhit","o","white")],
               loc="lower left",bbox_to_anchor=(.563,.566),ncol=3,columnspacing=.7,handletextpad=.1)
    header(fig,.035,.484,"C","ASD hit odds after covariate adjustment")
    c=fig.add_axes([.125,.144,.313,.262]);grid(c)
    m=read(src/"F2_models.tsv")
    for i,model in enumerate(["M0","M1","M2","M3"]):
        for screen,y,col,marker in [("CHEN",i-.13,ASD,"o"),("UMHOEFER_TREG",i+.13,TREG,"s")]:
            r=m[m.screen.eq(screen)&m.model.eq(model)].iloc[0]
            c.plot([r.CI_low,r.CI_high],[y,y],color=col,lw=1.1)
            c.scatter(r.OR,y,facecolor=col if model=="M2" else "white",edgecolor=col,lw=1,s=25,marker=marker,zorder=4)
    c.axvline(1,color="#737D83",lw=.75,ls="--")
    c.set_xscale("log");c.set_xlim(.47,42);c.set_xticks([.5,1,2,5,10,30],["0.5","1","2","5","10","30"])
    c.set_yticks(range(4),["M0","M1","M2","M3"]);c.tick_params(axis="y",length=0)
    c.set_ylim(3.55,-.55);c.spines["left"].set_visible(False);c.set_xlabel("Odds ratio (95% profile CI)")
    fig.legend(handles=[mark(ASD,"iTreg"),mark(TREG,"Primary Treg","s")],ncol=2,loc="lower left",bbox_to_anchor=(.138,.023),handletextpad=.15,columnspacing=.85)
    header(fig,.501,.484,"D","Loss-of-function subset by direction")
    l=read(src/"F2_lof.tsv"); l=l[l.analysis.eq("STRICT_CONDITIONAL_EXACT")]
    assert len(l)==4
    d=fig.add_axes([.697,.144,.155,.262]);grid(d)
    for i,r in enumerate(l.itertuples()):
        col=TREG if r.direction=="KO_REDUCES_FOXP3" else ORANGE
        d.hlines(i,0,r.effect,color=col,lw=1)
        d.scatter(r.effect,i,marker="v" if r.direction=="KO_REDUCES_FOXP3" else "^",s=32,color=col,zorder=4)
        fig.text(.911,.144+.262*(1-(i+.55)/4.1),pval(r.p_holm),ha="center",va="center",size=7.5)
    d.axvline(1,color=GRAY,lw=.75,ls="--")
    d.set_yticks(range(4),["iTreg · lower","iTreg · higher","Primary · lower","Primary · higher"])
    d.set_ylim(3.55,-.55);d.tick_params(axis="y",length=0);d.spines["left"].set_visible(False)
    d.set_xlim(-.1,3.35);d.set_xticks([0,1,2,3]);d.set_xlabel("Observed / expected")
    fig.text(.911,.425,"Holm P",size=7.4,weight="bold",ha="center")
    save_main(fig,stage,2)

def f3(stage,src):
    out=stage/"03_results"
    r=read(src/"F3_foxp3_tcr.tsv").set_index("analysis")
    null=read(src/"F3_foxp3_tcr_null.tsv.gz")
    groups=read(out/"R10_TREG_TCONV_COUNTS.tsv").set_index("category").loc[["Treg-only","Tconv-only","Both","Neither"]]
    primary=read(out/"R10_TREG_TCONV_PRIMARY.tsv").set_index("analysis").loc["discordant_primary"]
    u=read(out/"R10_TREG_TCONV_UNIVERSE.tsv.gz")
    continuous=read(out/"R10_TREG_TCONV_CONTINUOUS.tsv").set_index("analysis").loc["abs_rank_mean"]
    assert sum(groups.total_n)==len(u) and sum(groups.asd_n)==u.asd.sum()
    assert int(primary.n)==groups.iloc[:2].total_n.sum() and len(u)==continuous.n
    write(groups.reset_index(),src/"F3C_hit_classes.tsv")
    write(primary.to_frame().T.reset_index(),src/"F3C_discordant_fisher.tsv")
    write(continuous.to_frame().T.reset_index(),src/"F3D_continuous_contrast.tsv")
    write(u[["gene","asd","delta_abs_rank"]].sort_values(["asd","delta_abs_rank","gene"]),src/"F3D_rank_difference.tsv.gz")
    fig=plt.figure(figsize=(7.85,6.30))
    header(fig,.035,.951,"A","FOXP3 vs TCR proliferation")
    a=fig.add_axes([.11,.625,.354,.245]);grid(a,"y")
    for key,col,label in [("FULL_MEASURED_UNADJUSTED",ASD,"Full measured"),("COMPLETE_CASE_CONDITIONED",TREG,"Attribute-conditioned")]:
        curves=null[null.analysis.eq(key)].sort_values("S")
        a.plot(curves.S,curves.probability,color=col,lw=1.5,marker="o",ms=2.3,
               label=f"{label}   P = {pval(r.loc[key,'raw_P'])}")
    observed=r.loc["FULL_MEASURED_UNADJUSTED","observed_S"]
    a.axvline(observed,color=ORANGE,ls=(0,(3,2)),lw=1)
    a.text(observed+.8,.061,f"Observed\n{observed:+.0f}",color=ORANGE,size=7)
    a.set_xlim(-28,9);a.set_ylim(bottom=0);a.set_ylabel("Exact null probability")
    a.set_xlabel("ASD FOXP3-only minus TCR-only hit count")
    a.legend(loc="upper left",bbox_to_anchor=(-.08,1.29),fontsize=6.8,handlelength=1.6,labelspacing=.24)
    header(fig,.523,.951,"B","FOXP3-only vs TCR-only hit odds")
    b=fig.add_axes([.681,.643,.256,.197]);grid(b)
    for i,(key,col) in enumerate([("DISCORDANT_UNADJUSTED_FIRTH",ASD),("DISCORDANT_ADJUSTED_FIRTH",TREG)]):
        row=r.loc[key]
        b.plot([row.CI_low,row.CI_high],[i,i],color=col,lw=1.5)
        b.scatter(row.OR,i,color=col,marker="s" if i else "o",s=28,zorder=3)
        b.text(row.OR,i-.20,f"{row.OR:.2f}",size=7.3,ha="center",color=col)
    b.axvline(1,color=GRAY,lw=.8,ls="--"); b.set_xscale("log")
    b.set_xlim(.65,20); b.set_xticks([1,2,5,10,20],["1","2","5","10","20"])
    b.set_yticks([0,1],["Unadjusted","Adjusted"]);b.tick_params(axis="y",length=0);b.spines["left"].set_visible(False)
    b.set_ylim(1.62,-.55);b.set_xlabel("Firth odds ratio (95% profile CI)")
    row=r.loc["DISCORDANT_ADJUSTED_FIRTH"]
    fig.text(.555,.543,f"Adjusted P = {pval(row.raw_P)}  ·  {int(row.N):,} discordant genes",size=7.5)
    header(fig,.035,.443,"C","Primary Treg and Tconv FOXP3 hits")
    c=fig.add_axes([.067,.138,.391,.24]);c.axis("off")
    c.set_xlim(0,1);c.set_ylim(-.52,4.38)
    for x,head in zip([.02,.61,.91],["Hit class","ASD","Other"]):
        c.text(x,4.09,head,ha="left" if x==.02 else "right",size=7.5,weight="bold")
    for i,(name,row) in enumerate(groups.iterrows()):
        y=3.3-i*.91
        c.axhline(y-.42,color=GRID,lw=.63)
        c.text(.02,y,name,va="center",size=8.2)
        c.text(.61,y,f"{int(row.asd_n):,}",va="center",ha="right",weight="bold",size=8.2,color=ASD)
        c.text(.91,y,f"{int(row.non_asd_n):,}",va="center",ha="right",size=8.2)
    fig.text(.075,.094,f"Discordant OR {primary.OR:.2f}  (95% CI {primary.CI_low:.2f}–{primary.CI_high:.1f})",size=7.6)
    fig.text(.075,.064,f"Two-sided Fisher P = {pval(primary.P)}  ·  n = {int(primary.n)}",size=7.5)
    header(fig,.523,.443,"D","Relative FOXP3 effect strength")
    d=fig.add_axes([.617,.159,.32,.218]);grid(d,"y")
    for asd,color,label in [(0,"#737F89","Other genes"),(1,ASD,"ASD genes")]:
        vals=u.loc[u.asd.eq(asd),"delta_abs_rank"].dropna().to_numpy()
        vals=np.sort(vals); d.step(np.r_[-1,vals,1],np.r_[0,np.arange(1,len(vals)+1)/len(vals),1],where="post",color=color,lw=1.25,label=f"{label} ({len(vals):,})")
    d.axvline(0,color=GRAY,lw=.6,ls="--")
    d.set_xlim(-1,1);d.set_xticks([-1,-.5,0,.5,1]);d.set_ylim(0,1.03)
    d.set_xlabel("Absolute-effect percentile, Treg − Tconv")
    d.set_ylabel("Cumulative fraction"); d.legend(loc="upper left",fontsize=6.8,handlelength=1.25,handletextpad=.4,labelspacing=.28)
    fig.text(.555,.064,f"Mean difference {continuous.effect:+.3f}  (95% CI {continuous.CI_low:+.3f} to {continuous.CI_high:+.3f})",size=7.3)
    fig.text(.555,.034,f"Label-permutation P = {continuous.P:.5f}",size=7.3)
    save_main(fig,stage,3)

def f4(stage,src):
    out=stage/"03_results"
    h=read(out/"R10_MONACO_EXPRESSION_HEATMAP_SOURCE.tsv")
    h=h[h.gene.isin(PANEL)]
    types=h.sort_values("cell_type_order").cell_type.unique().tolist()
    palette_data=h.pivot(index="gene",columns="cell_type",values="log2_tpm_plus_1").reindex(index=PANEL,columns=types)
    humans=read(out/"R10_HUMAN_EXPRESSION_PANEL.tsv").set_index("gene").loc[PANEL]
    mice=read(out/"R10_MOUSE_EXPRESSION_PANEL.tsv").set_index("human_gene").loc[PANEL]
    proteins=read(out/"R10_PERTURB_PROTEIN.tsv")
    proteins=proteins[proteins.gene.isin(GENES)]
    tags=["FOXP3","CD25","CTLA4","Helios"]
    panel=proteins[proteins.protein.isin(tags)].pivot(index="gene",columns="protein",values="score").reindex(columns=tags)
    panel=panel.dropna(how="all")
    panel=panel.reindex([g for g in GENES if g in panel.index])
    assert len(types)==29 and panel.notna().all().all() and len(panel)>0
    write(h,src/"F4A_human_expression.tsv")
    both=pd.DataFrame({"gene":PANEL,"human_breadth":humans.breadth_ge1.values,
                       "human_n":humans.n_cell_types.values,"human_k":humans.detected_ge1.values,
                       "mouse_breadth":mice.breadth_1.values,"mouse_n":mice.n_groups.values,
                       "mouse_k":mice.detected_groups_1.values})
    write(both,src/"F4B_breadth.tsv");write(panel.reset_index(),src/"F4C_protein_scores.tsv")
    assert np.allclose(both.human_breadth,both.human_k/both.human_n,equal_nan=True)
    assert np.allclose(both.mouse_breadth,both.mouse_k/both.mouse_n,equal_nan=True)
    fig=plt.figure(figsize=(8.15,6.85))
    header(fig,.035,.954,"A","Expression across 29 human immune cell types")
    a=fig.add_axes([.115,.432,.65,.46]); im=vector_cells(a,palette_data,SEQ,vmin=0)
    a.set_yticks(range(len(PANEL)),PANEL,style="italic",size=7.2);a.tick_params(length=0,pad=3)
    a.set_xticks(range(len(types)),[str(c).replace("_"," ") for c in types],rotation=61,ha="right",rotation_mode="anchor",size=6.9)
    a.axhline(13.5,color="white",lw=2);a.spines[["left","bottom"]].set_visible(False)
    cb=fig.colorbar(im,cax=fig.add_axes([.79,.48,.012,.20]));cb.outline.set_visible(False);cb.ax.tick_params(labelsize=7)
    fig.text(.788,.716,"log₂(TPM + 1)",size=7.4)
    for t in a.get_xticklabels():
        if t.get_text()=="Treg": t.set_weight("bold"); t.set_color(TREG)
    header(fig,.035,.270,"B","Detection breadth by species")
    b=fig.add_axes([.197,.065,.247,.155])
    labels=[]; yi=0
    for species, color in [("human",TREG),("mouse",SPECIES)]:
        counts=(both[both.gene.isin(GENES)].groupby([species+"_k",species+"_n"]).size())
        for (k,n),number in counts.items():
            b.barh(yi,number,height=.50,color=color)
            b.text(number+.25,yi,str(int(number)),va="center",fontsize=7.4,weight="bold")
            labels.append(f"{species.title()}   {int(k)}/{int(n)}")
            yi+=1
    b.set_yticks(range(yi),labels); b.set_ylim(yi-.4,-.60)
    b.tick_params(axis="y",length=0);b.spines["left"].set_visible(False)
    b.set_xlim(0,len(GENES)+1.8); b.set_xticks([0,7,14]);grid(b)
    b.set_xlabel("Number of discovery genes")
    header(fig,.512,.270,"C","Protein scores after perturbation")
    c=fig.add_axes([.65,.055,.199,.151])
    maxabs=np.nanmax(np.abs(panel)); img=vector_cells(c,panel,DIV,-maxabs,maxabs)
    c.set_yticks(range(len(panel)),panel.index,style="italic",size=7.2)
    c.set_xticks(range(len(tags)),tags,fontsize=7);c.tick_params(length=0,pad=3)
    c.spines[["left","bottom"]].set_visible(False)
    bar=fig.colorbar(img,cax=fig.add_axes([.874,.06,.012,.125]));bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=7);fig.text(.867,.208,"Score",size=7.1)
    save_main(fig,stage,4)

if __name__=="__main__":
    opt=args()
    s=legacy_inputs(opt.stage,opt.data)
    for func in (f1, f2, f3, f4):
        func(opt.stage,s)
