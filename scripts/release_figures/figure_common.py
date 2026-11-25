"""R11 vector plotting, frozen-panel preservation and deterministic PDF/TIFF export."""
from pathlib import Path
import hashlib
import html
import io
import re
import shutil
import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image

DARK="#20303F"; NAVY="#243F65"; TEAL="#15776F"; ORANGE="#B86128"; PURPLE="#76568A"
GREY="#798791"; GRID="#DCE1E3"; MISSING="#D5D8DB"
SIGNED=LinearSegmentedColormap.from_list("signed",["#924F32","#D89871","#FCFAF3","#88C2B2","#19665F"])
SIGNED.set_bad(MISSING)
PANEL="CREBBP CUL3 DYRK1A EP300 FOXP1 KDM3B KDM6B KMT5B PHIP SIN3A SMARCA4 TBL1XR1 TLE3 ZMYND8".split()
STUDIES=["GSE121279","GSE263830","GSE121251"]
CONTRASTS=["F1N","F1A","S3","F2R","F2A"]
SHORT={"F1N":"Naïve","F1A":"Activated","S3":"Pooled Tregs","F2R":"Resting","F2A":"Activated"}
DESIGN=[("GSE121279",["F1N","F1A"],"Foxp1 · GSE121279"),("GSE263830",["S3"],"Sin3a · GSE263830"),("GSE121251",["F2R","F2A"],"Foxp1 · GSE121251")]
ROW_NAMES=["Foxp1 · naïve","Foxp1 · activated","Sin3a · sex-adjusted","Foxp1 · resting","Foxp1 · activated"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.titlesize":9,
    "axes.labelsize":8,"axes.linewidth":.65,"axes.spines.top":False,"axes.spines.right":False,
    "axes.edgecolor":GREY,"text.color":DARK,"axes.labelcolor":DARK,"xtick.color":DARK,"ytick.color":DARK,
    "xtick.major.size":3,"ytick.major.size":3,"xtick.major.width":.6,"ytick.major.width":.6,
    "xtick.labelsize":7.5,"ytick.labelsize":7.5,"pdf.fonttype":42,"ps.fonttype":42,
    "pdf.compression":9,"svg.fonttype":"none","legend.frameon":False,"legend.fontsize":7,
    "savefig.facecolor":"white","axes.facecolor":"white"})

def read(path):
    return pd.read_csv(path,sep="\t",low_memory=False)

def write(frame,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    frame.to_csv(path,index=False,sep="\t",na_rep="NA",float_format="%.17g",
                 compression={"method":"gzip","mtime":0} if path.suffix==".gz" else None)

def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()

def bytehash(data):
    return hashlib.sha256(data).hexdigest()

def grid(ax,axis="x"):
    ax.set_axisbelow(True)
    ax.grid(axis=axis,color=GRID,lw=.55)

def header(fig,x,y,letter,title):
    fig.text(x,y,letter,fontsize=13,weight="bold",va="baseline")
    fig.text(x+.035,y,title,fontsize=9.1,weight="bold",va="baseline")

def ptext(value):
    if pd.isna(value):
        return "NA"
    if value<.0001:
        return f"{value:.1e}".replace("e-0","e−").replace("e-","e−")
    return f"{value:.4g}"

def figpdf(fig):
    out=io.BytesIO()
    fig.savefig(out,format="pdf",metadata={"CreationDate":None,"ModDate":None})
    plt.close(fig)
    return fitz.open(stream=out.getvalue(),filetype="pdf")

def main_export(fig,stage,number):
    path=stage/"04_figures"/f"Figure_{number}.pdf"
    doc=figpdf(fig)
    doc.save(path,deflate=True,garbage=4,no_new_id=True)
    doc.close()
    tiff_from_pdf(path)

def tiff_from_pdf(path):
    with fitz.open(path) as doc:
        pix=doc[0].get_pixmap(dpi=600,colorspace=fitz.csRGB,alpha=False)
        im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
    im.save(path.with_suffix(".tif"),compression="tiff_lzw",dpi=(600,600))

def cap_page(target,graphic,number,title,body):
    """Full vector art and separate, unscaled justified caption on the same page."""
    margin=23; width=595; inner=width-margin*2
    css="""*{font-family:sans-serif;color:#20303f} body{margin:0;padding:0}
    p{font-size:9pt;line-height:1.38;margin:0 0 5pt;text-align:justify}
    h3{font-size:10pt;line-height:1.33;margin:0 0 5pt;font-weight:bold}"""
    heading=f"Supplementary Figure S{number}. {title.rstrip('.')}."
    cap="<h3>"+html.escape(heading)+"</h3>"+"".join("<p>"+html.escape(p)+"</p>" for p in body.split("\n\n"))
    story=fitz.Story(cap,user_css=css)
    more,rect=story.place(fitz.Rect(0,0,inner,2200)); assert not more
    art_height=inner*graphic[0].rect.height/graphic[0].rect.width
    cap_y=16+art_height+12; cap_h=rect[3]-rect[1]+9
    page=target.new_page(width=width,height=cap_y+cap_h+18)
    page.show_pdf_page(fitz.Rect(margin,16,width-margin,16+art_height),graphic,0)
    page.draw_line((margin,cap_y-6),(width-margin,cap_y-6),color=(.83,.86,.86),width=.5)
    slack,scale=page.insert_htmlbox(fitz.Rect(margin,cap_y,width-margin,cap_y+cap_h),cap,css=css,scale_low=1)
    assert scale==1 and slack>=0 and f"Figure S{number}." in page.get_text()
    return heading

def adjusted_continuous(src):
    data=read(src/"R10_TREG_TCONV_ADJUSTED.tsv")
    subset=data[data.analysis.eq("continuous_delta_abs_rank_conditioned")]
    if len(subset)!=1:
        raise ValueError("Expected one property-adjusted continuous contrast: "+str(data.analysis.tolist()))
    return subset.iloc[0]

def revised_context_figure(src):
    """Reuse the frozen estimates and show expression breadth with the two direct comparisons."""
    r=read(src/"F3_foxp3_tcr.tsv").set_index("analysis")
    null=read(src/"F3_foxp3_tcr_null.tsv.gz")
    margins=read(src/"F3C_hit_classes.tsv").set_index("category").loc[["Treg-only","Tconv-only","Both","Neither"]]
    primary=read(src/"F3C_discordant_fisher.tsv").iloc[0]
    continuous=read(src/"F3D_continuous_contrast.tsv").iloc[0]
    adj=adjusted_continuous(src)
    ranks=read(src/"F3D_rank_difference.tsv.gz")
    breadth=read(src/"F4B_breadth.tsv").set_index("gene").loc[PANEL]
    assert margins.asd_n.sum()==ranks.asd.sum() and margins.total_n.sum()==len(ranks)==continuous.n
    assert primary.n==margins.iloc[:2].total_n.sum()
    fig=plt.figure(figsize=(8.1,6.35))
    header(fig,.028,.944,"A","FOXP3 vs TCR proliferation")
    a=fig.add_axes([.107,.602,.335,.254]); grid(a,"y")
    for key,color,name in [("FULL_MEASURED_UNADJUSTED",NAVY,"Measured"),("COMPLETE_CASE_CONDITIONED",TEAL,"Conditioned")]:
        subset=null[null.analysis.eq(key)].sort_values("S")
        a.plot(subset.S,subset.probability,color=color,lw=1.4,marker="o",ms=2.1,label=f"{name}   P = {ptext(r.loc[key,'raw_P'])}")
    observed=r.loc["FULL_MEASURED_UNADJUSTED","observed_S"]
    a.axvline(observed,color=ORANGE,lw=1,ls="--")
    a.annotate(f"Observed\n{observed:+.0f}",(observed+.9,.062),fontsize=7,color=ORANGE)
    a.set(xlim=(-28,10),ylim=(0,None),xlabel="ASD FOXP3-only minus TCR-only hit count",ylabel="Exact null probability")
    a.legend(loc="upper left",bbox_to_anchor=(0,1.25),labelspacing=.23,handlelength=1.1)
    header(fig,.521,.944,"B","FOXP3 hits in primary Tregs and Tconvs")
    b=fig.add_axes([.566,.624,.358,.245]); b.axis("off"); b.set_xlim(0,1);b.set_ylim(-.25,4.4)
    for j,(x,label) in enumerate([(.015,"Hit class"),(.65,"ASD"),(.96,"Other")]):
        b.text(x,4.12,label,ha="left" if not j else "right",weight="bold")
    for j,(name,row) in enumerate(margins.iterrows()):
        y=3.26-j*.88
        b.axhline(y-.32,color=GRID,lw=.65)
        b.text(.015,y,name,va="center")
        b.text(.65,y,str(row.asd_n),va="center",ha="right",weight="bold",color=NAVY)
        b.text(.96,y,f"{row.non_asd_n:,}",va="center",ha="right")
    fig.text(.567,.580,f"Discordant OR {primary.OR:.2f}  (95% CI {primary.CI_low:.2f}–{primary.CI_high:.1f})",size=7.6)
    fig.text(.567,.552,f"Two-sided Fisher P = {ptext(primary.P)}  ·  n = {int(primary.n)}",size=7.5)
    header(fig,.028,.441,"C","FOXP3 effect strength in Tregs vs Tconvs")
    c=fig.add_axes([.107,.166,.335,.210]); grid(c,"y")
    for asd,color,name in [(0,GREY,"Other"),(1,NAVY,"ASD")]:
        vals=np.sort(ranks.loc[ranks.asd.eq(asd),"delta_abs_rank"].dropna())
        c.step(np.r_[-1,vals,1],np.r_[0,np.arange(1,len(vals)+1)/len(vals),1],where="post",color=color,lw=1.3,label=f"{name} ({len(vals):,})")
    c.set(xlim=(-1,1),ylim=(0,1.03),xlabel="Absolute-effect percentile, Treg − Tconv",ylabel="Cumulative fraction")
    c.axvline(0,color=GREY,lw=.6,ls="--")
    c.legend(loc="upper left",labelspacing=.24)
    fig.text(.071,.057,f"Mean difference {continuous.effect:+.3f}   P = {ptext(continuous.P)}",size=7.45)
    fig.text(.071,.029,f"Property-adjusted P = {ptext(adj.P)}",size=7.45)
    header(fig,.521,.441,"D","Detection breadth across immune cells")
    d=fig.add_axes([.713,.150,.218,.219]);grid(d)
    rows=[]
    for species,color in [("human",TEAL),("mouse",PURPLE)]:
        grouped=breadth.groupby([f"{species}_k",f"{species}_n"]).size()
        for (k,n),value in grouped.items():
            rows.append((f"{species.title()}   {k:.0f}/{n:.0f}",int(value),color))
    for i,(name,value,color) in enumerate(rows):
        d.barh(i,value,color=color,height=.49)
        d.text(value+.25,i,str(value),va="center",weight="bold",fontsize=7.5)
    d.set_yticks(range(len(rows)),[v[0] for v in rows]);d.invert_yaxis()
    d.tick_params(axis="y",length=0);d.spines["left"].set_visible(False)
    d.set_xlim(0,16);d.set_xticks([0,7,14]);d.set_xlabel("Number of candidate genes")
    return fig
