#!/usr/bin/env python3
"""Draw the six supplementary figures from companion gene- and program-level tables.

Usage: python figures/render_supp.py DATA_ROOT FIGURES_OUTPUT_DIR
Individual PDFs and 600-dpi TIFFs appear in individual_supplementary/. The six
pages of Supplementary_Figures_FINAL.pdf use identical source geometry.
"""
import argparse
import datetime
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

INK = "#183043"
NAVY = "#244B66"
TEAL = "#087D78"
RUST = "#A84F31"
MID = "#465966"
GRAY = "#85929C"
PALE = "#D6DEE2"
PAPER = "white"
SCREENS = ["CHEN", "UMHOEFER_TREG"]
NAMES = ["Chen iTreg", "Primary Treg"]
COLORS = [NAVY, TEAL]
DATE = datetime.datetime(2000, 1, 1, tzinfo=datetime.timezone.utc)
META = {"Creator": "ASD_TREG_CRISPR", "CreationDate": DATE, "ModDate": DATE}
mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "text.color": INK,
    "axes.labelcolor": INK, "axes.edgecolor": GRAY, "xtick.color": MID,
    "ytick.color": MID, "axes.linewidth": .6, "xtick.major.width": .6,
    "ytick.major.width": .6, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": PAPER,
})


def read(source, file):
    return pd.read_csv(source / file, sep="\t", dtype=str, keep_default_na=False)


def subset(data, kind):
    s = data.loc[data.record_type.eq(kind)].copy()
    if s.empty:
        raise ValueError(f"Missing {kind} records")
    return s


def number(value):
    return int(float(value))


def ptext(value):
    value = float(value)
    if value == 1:
        return "1"
    if value < .001:
        mant, exp = f"{value:.2e}".split("e")
        return f"{mant}e{int(exp)}"
    return f"{value:.3f}"


def page(n, heading, height):
    fig = plt.figure(figsize=(7.1, height))
    fig.text(.04, .981, f"Supplementary Figure S{n}", fontsize=9, weight="bold",
             color=TEAL, va="top")
    fig.text(.04, .932, heading, fontsize=14, weight="bold", va="top")
    return fig


def header(fig, y, letter, heading, x=.04):
    fig.text(x, y, letter, fontsize=11, weight="bold", va="center")
    fig.text(x + .038, y, heading, fontsize=10, weight="bold", va="center")


def canvas(fig, bounds):
    ax = fig.add_axes(bounds)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    return ax


def rule(ax, x1, x2, y):
    ax.plot([x1, x2], [y, y], color=PALE, lw=.6, clip_on=False)


def cell(ax, x, y, value, align="center", size=8.4, color=INK, bold=False, **kwargs):
    ax.text(x, y, value, ha=align, va="center", fontsize=size, color=color,
            weight="bold" if bold else "normal", **kwargs)


def table(ax, xs, headers, rows, left=(0,), size=8.2, header_size=None, emphasis=()):
    """Use equal row centers; the first and final horizontal rules frame the table."""
    step = 1 / (len(rows) + 1.15)
    for j, label in enumerate(headers):
        cell(ax, xs[j], 1 - .5*step, label, "left" if j in left else "center",
             header_size or size, MID, True)
    rule(ax, 0, 1, 1-step)
    for i, row in enumerate(rows):
        y = 1-(i+1.5)*step
        for j, value in enumerate(row):
            cell(ax, xs[j], y, str(value), "left" if j in left else "center",
                 size, INK, i in emphasis)
    rule(ax, 0, 1, 1-(len(rows)+1)*step)


def clean(ax, grid="x"):
    ax.spines["left" if grid == "x" else "top"].set_visible(False)
    ax.tick_params(length=3, labelsize=8)
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color=PALE, lw=.45)


def intfrac(k, n):
    return f"{number(k):,}/{number(n):,}"


def program_order(data):
    d = data.copy()
    d["module"] = d.program_id.str.extract(r"(\d+)").astype(int)
    return d.sort_values("module")


def s1(src):
    data = read(src, "S1_source.tsv")
    genes = subset(data, "common_asd_gene").sort_values("gene")
    prefixes = ["chen", "umhoefer_treg", "umhoefer_tconv"]
    vals = genes[[p+"_author_lfc" for p in prefixes]].astype(float)
    limit = np.abs(vals).to_numpy().max()
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    cmap = LinearSegmentedColormap.from_list("signed", ["#EACBBC", PAPER, "#A0CFC9"])
    assert len(genes) == genes.gene.nunique() == 33
    assert all(genes[p+"_measured"].eq("TRUE").all() for p in prefixes)
    fig = page(1, "FOXP3 screen measurements across 33 common genes", 6.65)
    header(fig, .855, "a", "Source log-fold change and author hit calls")
    for j, frame in enumerate([genes.iloc[:17], genes.iloc[17:]]):
        ax = canvas(fig, [.04 + j*.48, .272, .447, .545])
        step = 1/18.35
        xs = [.377, .634, .882]
        for x, label in zip(xs, ["iTreg", "Treg", "Tconv"]):
            cell(ax, x, 1-.5*step, label, bold=True)
        cell(ax, 0, 1-.5*step, "Gene", "left", bold=True)
        rule(ax, 0, 1, 1-step)
        for i, (_, row) in enumerate(frame.iterrows()):
            y = 1-(i+1.5)*step
            cell(ax, 0, y, row.gene, "left", size=8.2, fontstyle="italic")
            for x, prefix in zip(xs, prefixes):
                value = float(row[prefix+"_author_lfc"])
                hit = row[prefix+"_hit"] == "YES"
                ax.add_patch(Rectangle((x-.106,y-.43*step), .208,.86*step,
                                       lw=1.05 if hit else 0,
                                       ec=INK, fc=cmap(norm(value))))
                display = f"{value:.2f}".replace("-0.00","0.00")
                cell(ax,x,y,display,size=8.1,bold=hit)
    fig.text(.04,.244,"Outlined bold cells indicate author hits.  Other cells are measured nonhits.",
             fontsize=8.4,color=MID)
    header(fig,.186,"b","Screen context")
    header(fig,.186,"c","Conventional T-cell hits",x=.53)
    a=canvas(fig,[.04,.035,.445,.115])
    for i,(name,desc) in enumerate([("iTreg","Genome-wide, induced Tregs"),("Treg","Targeted, primary Tregs"),("Tconv","Targeted, conventional T cells")]):
        cell(a,0,.84-i*.335,name,"left",bold=True,size=8.3)
        cell(a,.19,.84-i*.335,desc,"left",size=8.2)
    r=subset(data,"tconv_enrichment").iloc[0]
    a=canvas(fig,[.53,.032,.43,.12])
    rows=[]
    for label,k,n in [("SFARI score 1",r.ASD_hits_or_observed,r.ASD_n),("Other genes",r.comparator_hits_or_expected,r.comparator_n)]:
        rows.append((label,intfrac(k,n),f"{100*number(k)/number(n):.1f}%"))
    table(a,[0,.58,.87],["Genes","Hits / tested","Rate"],rows)
    return fig


def s2(src):
    d=read(src,"S2_source.tsv")
    recal=subset(d,"recalibration")
    omission=subset(d,"omission")
    fig=page(2,"Sensitivity to hit definition and candidate omission",6.2)
    header(fig,.845,"a","Screen hits after recalibration or removal of shared candidates")
    a=canvas(fig,[.045,.542,.91,.26])
    rows=[]
    for screen,screen_name in zip(SCREENS,NAMES):
        s=recal[recal.screen.eq(screen)].set_index("analysis_type_or_removal")
        for tag,label in [("FULL_UNADJUSTED_RECALIBRATED","Recalibrated"),("REMOVE_FOXP1;ZMYND8","Omit FOXP1 + ZMYND8")]:
            r=s.loc[tag]
            rows.append((screen_name,label,intfrac(r.K,r.A),f"{float(r.effect):.2f}",ptext(r.one_sided_P)))
    table(a,[0,.24,.64,.78,.922],["Screen","Sensitivity","Hits / ASD","OR","Fisher P"],rows,left=(0,1))
    header(fig,.498,"b","Regulator-class conditioning after recalibration")
    a=canvas(fig,[.045,.394,.91,.065])
    for j,(screen,label,color) in enumerate(zip(SCREENS,NAMES,COLORS)):
        r=recal[(recal.screen==screen)&recal.analysis_type_or_removal.eq("CLASS_CONDITIONED_RECALIBRATED")].iloc[0]
        cell(a,.01+j*.52,.48,label,"left",size=9,bold=True,color=color)
        cell(a,.18+j*.52,.48,f"O/E {float(r.observed_expected_ratio):.2f}","left",size=9)
    header(fig,.34,"c","Omission of each measured score-1 gene")
    for j,(screen,name,color,marker) in enumerate(zip(SCREENS,NAMES,COLORS,["o","D"])):
        a=fig.add_axes([.12+j*.47,.087,.345,.188]); clean(a,"y")
        s=omission[omission.screen.eq(screen)].sort_values("removed_gene")
        assert len(s)==[242,33][j] and s.removed_gene.is_unique
        x=np.arange(1,len(s)+1)
        y=-np.log10(s.one_sided_fisher_p.astype(float))
        hit=s.removed_was_author_hit.eq("TRUE")
        a.scatter(x[~hit],y[~hit],s=9,facecolors="white",edgecolors=GRAY,lw=.5)
        a.scatter(x[hit],y[hit],s=12,marker=marker,color=color,lw=0)
        a.axhline(-np.log10(.05),color=MID,lw=.7,ls="--")
        a.set(xlim=(0,len(s)+1),ylim=(0,6),xticks=[1,len(s)],yticks=[0,2,4,6])
        a.set_xlabel("Omitted gene (alphabetical rank)",fontsize=8)
        if j==0: a.set_ylabel("−log₁₀(Fisher P)",fontsize=8.4)
        a.set_title(name,color=color,fontsize=9,pad=5,weight="bold")
    return fig


def conditioning_rows(grid):
    """For a screen, the historical class and joint reference names differ."""
    out=[]
    for tag,threshold,name in [("old",None,"Original"),("CPM","0.5","CPM 0.5"),("CPM","1","CPM 1"),("CPM","2","CPM 2"),("HPA","1","HPA 1")]:
        row=[name]
        for screen,reference in zip(SCREENS,["CHEN_ITREG_NTC","CHEN_HTREG"]):
            for analysis in ["class_conditioned","strict_joint"]:
                ref = ("FROZEN_MEASURED" if analysis=="class_conditioned" else "FROZEN_GSE90600") if tag=="old" else ("HPA_IMMUNE" if tag=="HPA" else reference)
                s=grid[(grid.screen==screen)&(grid.reference==ref)&(grid.analysis==analysis)&grid["set"].eq("SFARI_SCORE1")]
                if threshold is not None: s=s[s.threshold.astype(float).eq(float(threshold))]
                assert len(s)==1
                r=s.iloc[0]
                row.append(f"{float(r.effect):.2f} / {float(r.adjusted_p):.4f}")
        out.append(row)
    return out


def s3(src):
    d=read(src,"S3_source.tsv")
    strata=subset(d,"class_stratum")
    g=read(src,"S3_conditioning.tsv")
    m=read(src,"S3_regression.tsv")
    fig=page(3,"Regulator classes, expression references and coverage",8.8)
    header(fig,.863,"a","Hit counts within measured regulator classes")
    a=canvas(fig,[.045,.704,.915,.125])
    rows=[]
    for cl,lab in [("BOTH_LISTED","TF + EpiFactors"),("EPIFACTOR_ONLY","EpiFactors only"),("TF_ONLY","TF only"),("NEITHER_LISTED","Neither")]:
        row=[lab]
        for screen,lib in [("CHEN","TRUE"),("CHEN","FALSE"),("UMHOEFER_TREG","TRUE")]:
            r=strata[(strata.screen==screen)&(strata.library_member==lib)&(strata.class_label==cl)].iloc[0]
            row.append(f"{intfrac(r.K,r.A)}   {intfrac(number(r.H)-number(r.K),number(r.N)-number(r.A))}")
        rows.append(row)
    table(a,[0,.384,.636,.891],["Catalog class","iTreg, library in","iTreg, library out","Primary Treg"],rows,size=8.1)
    fig.text(.045,.686,"Each pair gives hits/tested for score-1 genes (left) and other genes (right).",fontsize=8,color=MID)
    header(fig,.641,"b","Class and strict joint conditioning")
    a=canvas(fig,[.045,.431,.915,.177])
    cell(a,.36,.968,"Chen iTreg",size=8.8,color=NAVY,bold=True)
    cell(a,.787,.968,"Primary Treg",size=8.8,color=TEAL,bold=True)
    b=canvas(fig,[.045,.428,.915,.152])
    table(b,[0,.257,.46,.688,.9],["Reference","Class","Strict","Class","Strict"],conditioning_rows(g),size=8.2,emphasis=(2,))
    fig.text(.045,.398,"Entries give observed/expected hits and two-screen Holm P. HPA threshold is 1 nTPM.",fontsize=8,color=MID)
    header(fig,.352,"c","Coverage and adjusted estimates across cohorts")
    c=canvas(fig,[.045,.178,.915,.148])
    rr=[]
    for label,block,context in [("Strict cohort, original","enrichment","original"),("Strict cohort, CPM 1","enrichment","CPM1"),("Historical adjusted OR","historical_regression",None),("M2 adjusted OR, CPM 1","regression",None)]:
        row=[label]
        for screen in SCREENS:
            sub=m[(m.screen==screen)&(m.block==block)]
            if block=="enrichment":
                r=sub[sub.analysis.eq("strict_joint")&sub.display_context.eq(context)].iloc[0]
                row.append(f"{intfrac(r.asd_n,r.original_asd_n)}    {intfrac(r.informative_asd_n,r.asd_n)}")
            elif block=="historical_regression":
                r=sub.iloc[0];row.append(f"{float(r.effect):.2f}    P {ptext(r.p)}")
            else:
                r=sub[sub.model.eq("M2")].iloc[0];row.append(f"{float(r.OR):.2f}    Holm {ptext(r.holm_p_m2)}")
        rr.append(row)
    table(c,[0,.536,.83],["Cohort / estimate","Chen iTreg","Primary Treg"],rr,size=8.3)
    fig.text(.045,.151,"Cohort cells give complete/original, then informative/complete score-1 genes.",fontsize=8,color=MID)
    header(fig,.11,"d","Mouse Foxp3 screen")
    q=subset(d,"mouse_sensitivity").set_index("row_type")
    prim=q.loc["PRIMARY"]; fit=q.loc["FITNESS_SENSITIVITY"]; clas=q.loc["CLASS_CONDITIONED"]
    a=canvas(fig,[.045,.023,.915,.059])
    for j,(name,metric,est,p) in enumerate([("Primary","OR",prim.odds_ratio,prim.one_sided_fisher_p),("Fitness sensitivity","OR",fit.fitness_odds_ratio,fit.fitness_one_sided_p),("Class-conditioned","O/E",clas.observed_expected_ratio,clas.class_conditioned_p)]):
        x=[.0,.34,.67][j]
        cell(a,x,.78,name,"left",size=8.2,bold=True)
        cell(a,x,.22,f"{metric} {float(est):.2f}    P {ptext(p)}","left",size=8.2)
    return fig


def s4(src):
    d=read(src,"S4_source.tsv")
    blood=program_order(subset(d,"blood_program"))
    pair=subset(d,"blood_pair")
    yan=program_order(subset(d,"yan_program"))
    assert blood.module.tolist()==yan.module.tolist()==list(range(1,11))
    assert not pair.duplicated(["program_id","anonymous_pair_index"]).any()
    fig=page(4,"Family-paired whole-blood scores and T-cell list overlap",6.9)
    header(fig,.845,"a","Whole-blood differences across 20 discordant sibling pairs")
    y=np.arange(10)
    ax=fig.add_axes([.107,.447,.31,.306]); clean(ax)
    mu=blood.mean_difference.astype(float).to_numpy()
    lo=blood.ci95_low.astype(float).to_numpy(); hi=blood.ci95_high.astype(float).to_numpy()
    ax.errorbar(mu,y,xerr=[mu-lo,hi-mu],fmt="o",mfc="white",mec=NAVY,ecolor=NAVY,
                ms=4,lw=.85,capsize=2)
    ax.axvline(0,color=GRAY,lw=.75,ls="--")
    ax.set_ylim(9.55,-.55)
    lim=max(np.max(-lo),np.max(hi))*1.08
    ax.set_xlim(-lim,lim)
    ax.set_xticks([-.01, 0, .01],["−0.01","0","0.01"])
    ax.set_yticks(y,[f"M{i}" for i in blood.module])
    ax.tick_params(axis="y",length=0)
    ax.set_xlabel("Mean ASD − sibling score",fontsize=8.3)
    fig.text(.107,.782,"Paired mean (95% CI)",fontsize=8.6,weight="bold")
    aux=fig.add_axes([.45,.447,.504,.306])
    aux.set(xlim=(0,1),ylim=(9.55,-.55));aux.axis("off")
    for x,t in zip([.12,.41,.74],["Exact P","BH10 q","Pairs (+ / −)"]):
        aux.text(x,1.09,t,transform=aux.transAxes,ha="center",fontsize=8.6,weight="bold")
    for i,(_,r) in enumerate(blood.iterrows()):
        these=pair[pair.program_id.eq(r.program_id)].ASD_minus_sibling.astype(float)
        assert len(these)==number(r.n_pairs)==20
        assert np.isclose(these.mean(),float(r.mean_difference),atol=1e-14)
        assert [(these>0).sum(),(these<0).sum(),(these==0).sum()]==[number(r.n_positive),number(r.n_negative),number(r.n_zero)]
        for x,val in zip([.12,.41,.74],[ptext(r.p_exact_two_sided),ptext(r.q_BH_10),f"{r.n_positive} / {r.n_negative}"]):
            aux.text(x,i,val,va="center",ha="center",fontsize=8.4)
    header(fig,.345,"b","Program genes in the Treg-like cluster reported by Yan et al.")
    for j,g in enumerate([yan.iloc[:5],yan.iloc[5:]]):
        a=fig.add_axes([.106+j*.489,.087,.324,.197]); clean(a)
        y=np.arange(5)
        a.barh(y,100*g.overlap_fraction.astype(float),height=.48,color=TEAL)
        for i,(_,r) in enumerate(g.iterrows()):
            a.text(10.6,i,intfrac(r.program_genes_listed,r.program_size),va="center",fontsize=8)
        a.set(xlim=(0,15.5),ylim=(4.5,-.5),xticks=[0,5,10],yticks=y)
        a.set_yticklabels([f"M{i}" for i in g.module]);a.tick_params(axis="y",length=0)
        a.set_xlabel("Program genes listed (%)",fontsize=8.3)
    return fig


def s5(src):
    d=read(src,"S5_source.tsv")
    programs=program_order(subset(d,"program"))
    excl=subset(d,"exact_overlap").set_index("dataset")
    hits=subset(d,"asd14_feature")
    assert len(hits)==hits.mapped_symbol.nunique()==14
    assert hits.tested_valid.eq("TRUE").all() and hits.deg_fdr05.eq("FALSE").all()
    fig=page(5,"Program-level shifts and exact-gene overlap",6.55)
    header(fig,.844,"a","Purified Treg gene-level program tests")
    a=canvas(fig,[.045,.418,.91,.38])
    rows=[]
    for _,r in programs.iterrows():
        eligible=r.formal_eligibility=="ELIGIBLE"
        status=("Higher" if r.effect_direction=="ASD_HIGHER" else "Lower") if r.interpretation_status=="FORMAL_SUPPORTED_Q_LT_005" else ("No BH10 support" if eligible else "Coverage-ineligible")
        rows.append([f"M{r.module}",intfrac(r.in_valid_universe,r.source_size),f"{float(r.coverage)*100:.1f}%",ptext(r.formal_q_BH10) if eligible else "—",status])
    table(a,[0,.245,.425,.615,.785],["Program","Measured / total","Coverage","BH10 q","Result"],rows,left=(0,4),size=8.3)
    header(fig,.357,"b","Differential expression among the 14 discovery hits")
    r=excl.loc["Moreno2026"]; y=excl.loc["Yan2025"]
    a=canvas(fig,[.045,.195,.91,.121])
    rows=[["Discovery hits",intfrac(r.asd14_deg,r.asd14_in_universe),f"{100*number(r.asd14_deg)/number(r.asd14_in_universe):.1f}%"],
          ["Other features",intfrac(r.other_deg,number(r.other_deg)+number(r.other_non_deg)),f"{100*number(r.other_deg)/(number(r.other_deg)+number(r.other_non_deg)):.1f}%"]]
    table(a,[0,.532,.89],["Moreno purified Tregs","DEGs / measured","Rate"],rows)
    fig.text(.045,.158,f"One-sided Fisher OR = {r.odds_ratio}, P = {r.p_value}.",fontsize=8.5,color=MID)
    header(fig,.099,"c","Coverage in the reported Treg-like gene list")
    fig.text(.078,.045,f"Yan et al.: {intfrac(y.asd14_deg,y.asd14_size)} queried genes listed; measured-gene coverage unavailable.",fontsize=9)
    return fig


def ko_direction(row, is_chen):
    direction=row.Chen_FOXP3_regulator_direction if is_chen else row.Umhoefer_Treg_direction
    measured=is_chen or row.Umhoefer_Treg_measured=="TRUE"
    if not measured:
        return "NM"
    return {"SUPPORTS_FOXP3":"↓", "RESTRAINS_FOXP3":"↑", "NOT_SIGNIFICANT":"○"}[direction]


def s6(src):
    mem=read(src,"S6_membership.tsv");cands=read(src,"S6_candidates.tsv").sort_values("gene")
    assert len(mem)==mem.approved_symbol.nunique()==245 and len(cands)==cands.gene.nunique()==17
    elig=mem[mem.lof_eligible.eq("1")]
    g2p=mem.mechanism_class.value_counts()
    support=elig.lof_candidate_supports.value_counts()
    assert len(elig)==104 and g2p["LoF-compatible"]==130
    assert cands.loc[~cands.display_sup.eq("–"),"gene"].tolist()==["SATB1"]
    fig=page(6,"Disease mechanisms and candidate-level functional assays",8.25)
    header(fig,.863,"a","Gene2Phenotype classifications across 245 score-1 genes")
    cat=["LoF-compatible","GoF","dominant negative","mixed","unknown"]
    ax=fig.add_axes([.238,.643,.344,.162]);clean(ax)
    vals=g2p.reindex(cat).to_numpy()
    ax.barh(range(5),vals,height=.55,color=[NAVY,NAVY,NAVY,NAVY,GRAY])
    ax.set(ylim=(4.55,-.55),xlim=(0,151),xticks=[0,50,100],yticks=range(5))
    ax.set_yticklabels(["LoF-compatible","Gain of function","Dominant negative","Mixed","Unknown"],fontsize=8.3)
    ax.tick_params(axis="y",length=0)
    for i,v in enumerate(vals):ax.text(v+3,i,str(v),fontsize=8.4,va="center",weight="bold")
    ax.set_xlabel("Score-1 genes",fontsize=8.3)
    a=canvas(fig,[.646,.647,.32,.162])
    cell(a,0,.89,"LoF sensitivity subset","left",size=9.1,bold=True)
    cell(a,0,.61,f"{len(elig)} / {len(mem)}","left",size=18,color=TEAL,bold=True)
    cell(a,0,.34,f"{support['G2P_inferred_absent_product']} inferred absent product","left",size=8.5)
    cell(a,0,.13,f"{support['G2P_evidence']} evidence-tagged","left",size=8.5)
    header(fig,.547,"b","Screen directions and additional published assays for 17 candidates")
    a=canvas(fig,[.04,.096,.92,.423])
    xs=[.003,.14,.35,.472,.566,.70,.81,.945]
    step=1/19.6
    for xx,txt in [(.52,"FOXP3 after KO"),(.822,"Additional readouts")]:
        cell(a,xx,.983,txt,bold=True,size=8.6)
    rule(a,.423,.619,.951);rule(a,.65,1,.951)
    for j,txt in enumerate(["Gene","G2P class","Eligible","iTreg","Treg","FOXP3","Phenotype","Suppression"]):
        cell(a,xs[j],.922,txt,"left" if j<2 else "center",size=7.9,color=MID,bold=True)
    rule(a,0,1,.891)
    for i,(_,r) in enumerate(cands.iterrows()):
        y=.863-i*step
        vals=[r.gene,{"LoF-compatible":"LoF-compatible","unknown":"Unknown","mixed":"Mixed"}[r.mechanism_class],
              "Yes" if r.lof_eligible=="1" else "No",ko_direction(r,True),ko_direction(r,False)]
        for name in ["display_expr","display_phen","display_sup"]:
            vals.append("●" if r[name]!="–" else "—")
        for j,(x,val) in enumerate(zip(xs,vals)):
            col=TEAL if val=="↓" else RUST if val=="↑" else MID if val in ("NM","—","○") else INK
            cell(a,x,y,val,"left" if j<2 else "center",size=8.2 if j<3 else 9.3,color=col,
                 fontstyle="italic" if j==0 else "normal")
        if (i+1)%4==0 and i<16:rule(a,0,1,y-.50*step)
    fig.text(.042,.063,"↑ higher and ↓ lower FOXP3 after knockout     ○ measured nonhit     NM unmeasured",fontsize=8.2,color=MID)
    fig.text(.042,.035,"● additional assay reported (including null results)     — no qualifying assay in these sources",fontsize=8.2,color=MID)
    return fig


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("data_root",type=Path)
    ap.add_argument("output_dir",type=Path)
    args=ap.parse_args()
    src=args.data_root/"figure_sources"/"supp"
    individual=args.output_dir/"individual_supplementary"
    individual.mkdir(parents=True,exist_ok=True)
    with PdfPages(args.output_dir/"Supplementary_Figures_FINAL.pdf",metadata=META) as combined:
        for i,draw in enumerate([s1,s2,s3,s4,s5,s6],start=1):
            fig=draw(src)
            stem=individual/f"Supplementary_Figure_S{i}"
            fig.savefig(stem.with_suffix(".pdf"),metadata=META)
            fig.savefig(stem.with_suffix(".tif"),dpi=600,pil_kwargs={"compression":"tiff_lzw"})
            combined.savefig(fig)
            plt.close(fig)
            print(f"Rendered Supplementary Figure S{i}",flush=True)


if __name__=="__main__":
    main()
