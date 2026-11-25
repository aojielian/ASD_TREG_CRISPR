"""Render the historical gene-property panel layout from published source tables.

A separately supplied historical vector page is required only to compose the full
Figure 3 with its unchanged earlier panels. The bare overlay retains the same
coordinates, and S6 is rebuilt solely from the released source tables.
"""
import argparse
import csv
import hashlib
import math
import os
from pathlib import Path
from types import SimpleNamespace
import tempfile
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "release"))
from layout import REPLAY_MAP, replay_view

os.environ.setdefault("SOURCE_DATE_EPOCH","1795392000")
os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"asd-treg-public-mpl"))
import fitz
import matplotlib
matplotlib.use("Agg")
from matplotlib import get_data_path
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import matplotlib.pyplot as plt
import pandas as pd

fontdir = Path(get_data_path()) / "fonts/ttf"
fontpaths = [fontdir / "DejaVuSans.ttf", fontdir / "DejaVuSans-Bold.ttf"]
fonts = [fitz.Font(fontfile=str(path)) for path in fontpaths]
colors = {"navy":"#183B56","teal":"#007F86","orange":"#BD603D","dark":"#344553",
          "gray":"#6B7781","line":"#CCD4D8","soft":"#EBF0F2","black":"#24313A","white":"#FFFFFF"}
CLASSES = ["LoF-compatible", "GoF", "dominant negative", "mixed", "unknown"]

def color(name):
    s = colors.get(name, name).lstrip("#")
    return tuple(int(s[i:i+2], 16)/255 for i in (0, 2, 4))


def textat(page, x, y, value, size=8, bold=False, ink="dark", align="left"):
    value = str(value)
    font = fonts[bool(bold)]
    width = font.text_length(value, fontsize=size)
    if align == "center":
        x -= width / 2
    elif align == "right":
        x -= width
    assert x >= -0.01 and x + width <= page.rect.width + 0.1, (x, y, value, width)
    assert y > size and y <= page.rect.height, (x, y, value)
    page.insert_text((x, y), value, fontname="FigBold" if bold else "FigSans",
                     fontfile=str(fontpaths[bool(bold)]), fontsize=size, color=color(ink))


def paragraph(page, x, y, width, value, size=8.6, leading=12):
    line = ""
    for word in value.split():
        candidate = f"{line} {word}".strip()
        if fonts[0].text_length(candidate, fontsize=size) > width and line:
            textat(page, x, y, line, size)
            y += leading
            line = word
        else:
            line = candidate
    if line:
        textat(page, x, y, line, size)
        y += leading
    return y


def rule(page, x0, x1, y, ink="line", width=0.55):
    page.draw_line((x0, y), (x1, y), color=color(ink), width=width)


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def savepdf(doc, path):
    doc.subset_fonts()
    doc.save(path, garbage=4, deflate=True, no_new_id=True)


SCREENS = ("CHEN", "UMHOEFER_TREG")
REFS = dict(CHEN="CHEN_ITREG_NTC", UMHOEFER_TREG="CHEN_HTREG")
OLDREFS = {"unadjusted":"FROZEN_MEASURED","class_conditioned":"FROZEN_MEASURED","strict_joint":"FROZEN_GSE90600"}
MODEL_LABEL = {"M0":"M0  Unadjusted","M1":"M1  Regulator class","M2":"M2  + expression, LOEUF","M3":"M3  + span, guides"}


def floatval(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def bound(x):
    # Preserve which side of the null a nearly unit confidence limit falls on.
    return f"{x:.3f}" if x!=1 and abs(x-1)<.005 else f"{x:.2f}"


NAVY = '#183B56'
TEAL = '#007F86'
ORANGE = '#BD603D'
GRAY = '#BDC4CA'
DARK = '#344553'
LIGHT = '#E8EDF0'
PALE = '#F5F7F8'

def configure():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 8,
        'axes.labelsize': 8, 'axes.titlesize': 9, 'axes.titleweight': 'bold',
        'axes.edgecolor': DARK, 'axes.labelcolor': DARK, 'text.color': NAVY,
        'xtick.color': DARK, 'ytick.color': DARK,
        'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'axes.linewidth': .6, 'xtick.major.width': .6, 'ytick.major.width': .6,
        'xtick.major.size': 3, 'ytick.major.size': 3,
        'lines.linewidth': 1, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'svg.fonttype': 'none', 'savefig.facecolor': 'white',
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
    })

def new_figure(title, height=6.6):
    fig = plt.figure(figsize=(7.2, height))
    fig.text(.055, .973, title, ha='left', va='top', fontsize=11, fontweight='bold')
    return fig


STYLE = SimpleNamespace(NAVY=NAVY,TEAL=TEAL,ORANGE=ORANGE,GRAY=GRAY,DARK=DARK,LIGHT=LIGHT,PALE=PALE,new_figure=new_figure)

def add_text(fig, x, y, text, size=8, color=None, **kwargs):
    return fig.text(x, y, text, ha=kwargs.pop("ha", "left"), va="center", fontsize=size,
                    color=color or STYLE.DARK, **kwargs)


def title(fig, letter, y, text):
    add_text(fig, .028, y, letter, 11, STYLE.NAVY, weight="bold")
    add_text(fig, .066, y, text, 9.3, STYLE.NAVY, weight="bold")


def panel_rule(fig, y, left=.065, right=.967, lw=.7):
    fig.add_artist(plt.Line2D([left, right], [y, y], transform=fig.transFigure,
                              color=STYLE.LIGHT, linewidth=lw, zorder=0))


def export(fig, folder, stem):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for obj in fig.findobj(Text):
        if not obj.get_visible() or not obj.get_text():
            continue
        extent = obj.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        if min(extent.x0, extent.y0) < -.0005 or max(extent.x1, extent.y1) > 1.0005:
            raise ValueError(f"{stem} text outside page: {obj.get_text()}")
    fig.savefig(folder / (stem + ".pdf"), dpi=300, metadata={"Creator": "Matplotlib"})
    fig.savefig(folder / (stem + ".png"), dpi=300)
    plt.close(fig)


def render_fig3(source, out, underlay=None):
    legacy, models, bridge = {}, {}, {}
    for row in rows(source/"Figure_3d_source.tsv"):
        screen = row["screen"]
        if row["block"] == "enrichment":
            legacy[screen,row["analysis"],row["display_context"]] = row
        elif row["block"] == "historical_regression":
            legacy[screen,"preflight","original"] = row
        elif row["block"] == "regression":
            models[screen,row["model"]] = row
        else:
            category = "OLD" if "_OLD_" in row["analysis"] else "NEW"
            bridge[screen,category] = row
    final = fitz.open()
    page = final.new_page(width=518.4,height=693)
    if underlay:
        expected = "a934ad8aed7e54cac4d35b059c517efb529be220b197050e3472dd416f111782"
        if hashlib.sha256(underlay.read_bytes()).hexdigest()!=expected:
            raise ValueError("Historical Figure 3 vector does not match the frozen source SHA-256.")
        with fitz.open(underlay) as doc:
            pg = doc[0]
            pg.add_redact_annot(fitz.Rect(290,253,518.4,464.4),fill=False)
            pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED)
            page.show_pdf_page(doc[0].rect, doc, 0, keep_proportion=False)
    textat(page,307,273,"d",10.7,True,"navy")
    textat(page,325,273,"Screen enrichment and conditioning",8.4,True,"navy")
    textat(page,301,293,"Background / analysis",7.3,True)
    textat(page,420,293,"iTreg",7.3,True,"navy","center")
    textat(page,484,293,"Treg",7.3,True,"teal","center")
    rule(page,301,508,299)
    def summary_line(y, title, analysis, cohort):
        textat(page,301,y,title,7.6)
        for screen, x in zip(SCREENS,(420,484)):
            r=legacy[screen,analysis,cohort]
            if analysis=="unadjusted":
                cell=f"{float(r['effect']):.2f} ({float(r['ci_low']):.2f}–{float(r['ci_high']):.2f})"
            else:
                cell=f"{float(r['effect']):.2f} ({float(r['adjusted_p']):.4f})"
            textat(page,x,y,cell,6.55,False,"black","center")
    textat(page,301,312,"Odds ratio (95% CI)",7.4,True,"navy")
    summary_line(326,"All measured", "unadjusted","original")
    summary_line(339,"Detectable (CPM ≥1)", "unadjusted","CPM1")
    rule(page,301,508,345)
    textat(page,301,358,"Observed / expected hits (Holm P)",7.4,True,"navy")
    summary_line(372,"Class, all measured", "class_conditioned","original")
    summary_line(385,"Class, detectable", "class_conditioned","CPM1")
    summary_line(398,"Strict, original", "strict_joint","original")
    summary_line(411,"Strict, detectable", "strict_joint","CPM1")
    rule(page,301,508,418)
    for x,sc in zip((420,484),SCREENS):
        before = 100*float(legacy[sc,"strict_joint","original"]["asd_original_coverage"])
        after = 100*float(legacy[sc,"strict_joint","CPM1"]["asd_original_coverage"])
        textat(page,x,431,f"{before:.1f} → {after:.1f}",7.1,True,"navy","center")
    textat(page,301,431,"Score 1 coverage (%)",7.05)
    rule(page,14,508,444,"gray",0.65)
    textat(page,16,458,"Gene-property adjustment within the same expression-detectable cohort",9.4,True,"navy")
    textat(page,16,470,"Firth odds ratios with 95% profile intervals. Each screen uses one complete-case cohort for M0–M3.",7.3)
    lower = [floatval(r["CI_low"]) for r in models.values()]
    upper = [floatval(r["CI_high"]) for r in models.values()]
    finite = [v for v in lower+upper if v is not None and v>0]
    if not finite:
        finite=[1.]
    xmin = min(0.5, 2**math.floor(math.log2(min(finite))))
    xmax = max(4, 2**math.ceil(math.log2(max(finite))))
    def project(v, left, right):
        return left + (math.log(v/xmin)/math.log(xmax/xmin))*(right-left)
    positions={"CHEN":(143,214,220,"navy"),"UMHOEFER_TREG":(325,398,404,"teal")}
    for sc,title in zip(SCREENS,("Chen iTreg","Umhoefer primary Treg")):
        left,right,numeric,ink=positions[sc]
        cohort=models[(sc,"M0")]
        textat(page,left,486,title,8.7,True,ink)
        textat(page,left,498,f"{int(cohort['n_genes']):,} genes  ·  {int(cohort['asd_n']):,} score 1",7)
        page.draw_line((project(1,left,right),504),(project(1,left,right),562),color=color("gray"),width=0.6,dashes="[2 2] 0")
        ticks=[t for t in (.125,.25,.5,1,2,4,8,16,32,64,128,256) if xmin<=t<=xmax]
        keep=[t for i,t in enumerate(ticks) if t==1 or t==xmax or (t>=2 and (round(math.log2(t))%2==0)) or (t<1 and t==ticks[0])]
        rule(page,left,right,564)
        for tv in keep:
            xx=project(tv,left,right)
            page.draw_line((xx,564),(xx,567),color=color("line"),width=.5)
            textat(page,xx,575,f"{tv:g}",6,False,"gray","center")
        for k,(mod,label) in enumerate(MODEL_LABEL.items()):
            y=512+15*k
            if sc==SCREENS[0]:
                textat(page,16,y+2,label,7.3,mod=="M2")
            r=models[(sc,mod)]
            or_,lo,hi=(floatval(r[v]) for v in ("OR","CI_low","CI_high"))
            if or_ and lo and hi and lo>0 and hi>=lo:
                assert lo<=or_<=hi
                px, x0, x1 = [project(v,left,right) for v in (or_,lo,hi)]
                rule(page,x0,x1,y,ink,1.15)
                for xcap in (x0,x1):
                    page.draw_line((xcap,y-1.6),(xcap,y+1.6),color=color(ink),width=.8)
                page.draw_circle((px,y),2.5,color=color(ink),fill=color(ink),width=.5)
                detail = f"{or_:.2f}  ({bound(lo)}–{bound(hi)})"
            elif or_:
                detail = f"{or_:.2f}  (CI unstable)"
            else:
                detail = "Inestimable"
            textat(page,numeric,y+2,detail,7.25,mod=="M2",ink if mod=="M2" else "dark")
    for sc,x in zip(SCREENS,(143,325)):
        textat(page,x,589,f"M2 Holm P = {float(models[sc,'M2']['holm_p_m2']):.3f}",7.3,True,"navy" if sc=="CHEN" else "teal")
    rule(page,16,508,597)
    pre = [legacy[sc,"preflight","original"] for sc in SCREENS]
    textat(page,16,610,"Historical adjusted model (different cohorts and expression reference)",7.8,True,"navy")
    p0,p1=pre
    textat(page,16,624,f"Chen: OR {float(p0['effect']):.2f} ({float(p0['ci_low']):.2f}–{float(p0['ci_high']):.2f}), n = {int(p0['total_n']):,}.   Primary Treg: OR {float(p1['effect']):.2f}, CI unstable, n = {int(p1['total_n']):,}.",7.6)
    rule(page,16,508,632)
    textat(page,16,646,"Expression-reference comparison within the same cohort and model",7.8,True,"navy")
    textat(page,16,660,"Expression reference",7.2,True)
    for sc,x in zip(SCREENS,(162,342)):
        textat(page,x,660,"Chen iTreg" if sc=="CHEN" else "Umhoefer primary Treg",7.2,True,"navy" if sc=="CHEN" else "teal")
    for which,y,label in (("OLD",672,"Original CD4 reference"),("NEW",684,"Treg-state reference")):
        textat(page,16,y,label,7.2)
        for sc,x in zip(SCREENS,(162,342)):
            r=bridge[sc,which]
            assert r["interval_status"]=="ESTIMABLE" and int(r["n_genes"])==int(models[sc,"M0"]["n_genes"])
            lo_,hi_=map(float,(r["CI_low"],r["CI_high"]))
            textat(page,x,y,f"{float(r['OR']):.2f} ({bound(lo_)}–{bound(hi_)})",7.4)

    final_path = out/("Figure_3_gene_property.pdf" if underlay else "Figure_3d_overlay.pdf")
    savepdf(final,final_path)
    final[0].get_pixmap(dpi=300,alpha=False).save(final_path.with_suffix(".png"))
    final.close()


def render_s6(source,out):
    summary = pd.read_csv(source/"Supplementary_Figure_S6_membership_source.tsv",
                          sep="\t",dtype=str,keep_default_na=False)
    genes = pd.read_csv(source/"Supplementary_Figure_S6_candidate_source.tsv",
                        sep="\t",dtype=str,keep_default_na=False).sort_values("gene").reset_index(drop=True)
    eligible = summary.lof_eligible.eq("1")
    assert len(summary) == 245 and eligible.sum() == 104 and len(genes) == 17
    configure()
    fig=STYLE.new_figure("Disease mechanisms and published functional readouts",8.45)
    title(fig,"a",.910,"G2P disease-mechanism labels across 245 SFARI score 1 genes")
    counts=summary.mechanism_class.value_counts().reindex(CLASSES,fill_value=0)
    scale=counts.max()
    for i,(cl,count) in enumerate(counts.items()):
        y=.865-i*.028
        add_text(fig,.066,y,{"GoF":"Gain of function","LoF-compatible":"Loss-of-function compatible",
                             "dominant negative":"Dominant negative","mixed":"Mixed","unknown":"Unknown"}[cl],7.5)
        fig.add_artist(Rectangle((.365,y-.0064),.175*count/scale,.0128,
                                 transform=fig.transFigure,fc=STYLE.NAVY if cl!="unknown" else STYLE.GRAY,
                                 ec="none"))
        add_text(fig,.55,y,str(count),7.8,ha="left")
    add_text(fig,.659,.855,"Eligible for LoF sensitivity",8.0,STYLE.TEAL,weight="bold")
    add_text(fig,.659,.808,f"{int(eligible.sum())} / {len(summary)}",14.0,STYLE.NAVY,weight="bold")
    add_text(fig,.659,.760,"Gene–disease evidence and\nallelic requirement applied.",7.8,linespacing=1.4)
    panel_rule(fig,.713)
    title(fig,"b",.682,"Screen hits and additional published assays for 17 candidates")
    positions={"gene":.066,"mech":.199,"i":.413,"t":.512,"expr":.644,"phen":.774,"sup":.911}
    add_text(fig,.461,.642,"FOXP3 after KO",8.1,STYLE.NAVY,ha="center",weight="bold")
    add_text(fig,.777,.642,"Additional readouts",8.1,STYLE.NAVY,ha="center",weight="bold")
    panel_rule(fig,.628,.366,.562)
    panel_rule(fig,.628,.578,.967)
    for key,label in zip(positions,["Gene","G2P label","iTreg","Treg","FOXP3","Phenotype","Suppression"]):
        add_text(fig,positions[key],.608,label,7.7,STYLE.NAVY,
                 ha="left" if key in ["gene","mech"] else "center",weight="bold")
    for i,r in genes.iterrows():
        y=.578 - i*.0231
        if i%2==0:
            fig.add_artist(Rectangle((.060,y-.01155),.906,.0231,transform=fig.transFigure,
                                     fc=STYLE.PALE,ec="none",zorder=-1))
        add_text(fig,positions["gene"],y,r.gene,8,STYLE.NAVY,fontstyle="italic")
        mech={"GoF":"GoF","LoF-compatible":"LoF compatible","dominant negative":"Dom. negative",
              "mixed":"Mixed","unknown":"Unknown"}[r.mechanism_class]
        eligible_mark=r.lof_eligible.lower() in ["1","yes","true"]
        add_text(fig,positions["mech"],y,mech+(" *" if eligible_mark else ""),7.15)
        for x,direction,hit in [(positions["i"],r.Chen_FOXP3_regulator_direction,r.Chen_hit),
                                (positions["t"],r.Umhoefer_Treg_direction,r.Umhoefer_Treg_hit)]:
            if hit=="YES":
                label,color={"SUPPORTS_FOXP3":("↓",STYLE.TEAL),"RESTRAINS_FOXP3":("↑",STYLE.ORANGE)}[direction]
            else:
                label,color={"NO":("○",STYLE.DARK),"NOT_MEASURED":("–",STYLE.GRAY)}[hit]
            add_text(fig,x,y,label,9.0,color,ha="center",weight="bold")
        for key in ("expr", "phen", "sup"):
            disp = r["display_"+key]
            add_text(fig,positions[key],y,disp,7.6,STYLE.NAVY if disp!="–" else STYLE.GRAY,ha="center")
    add_text(fig,.066,.170,"↓ / ↑  Lower / higher FOXP3 after KO.  ○  Measured non-hit.  –  Unmeasured or no additional assay found.",7.1)
    add_text(fig,.066,.145,"* Eligible for LoF sensitivity. Letters identify a study reporting the assay, including null results.",7.1)
    add_text(fig,.066,.120,"Screens use Chen iTreg and Umhoefer Treg. Additional studies   C  Chen     S  Schumann     K  Kolb",7.1)
    panel_rule(fig,.097)
    add_text(fig,.066,.073,"SATB1",8.3,STYLE.TEAL,weight="bold",fontstyle="italic")
    add_text(fig,.155,.073,"KO reduced Treg suppression at all three tested ratios (Kolb Fig. 5C, n = 4 in the legend).",7.45)
    add_text(fig,.066,.044,"Across the four studies, variant-specific and ASD patient Treg functional assays were unavailable.",7.3)
    caption1 = ("Gene2Phenotype (G2P) records classified 245 SFARI score 1 genes. The 104 genes eligible for the gene-loss sensitivity "
                "comprised 100 with an inferred absent product and four with direct database evidence for the mechanism. "
                "For seven genes, the disorder name specified ASD. For 16, ASD was a disease feature. The other 81 genes had broader neurodevelopmental associations. "
                "PTEN was excluded because some variants exhibit dominant-negative activity (PMID 32350270).")
    caption2 = ("Arrows show FOXP3 protein change after gene knockout (KO) in the source screen. Circles indicate measured non-hits. A dash "
                "denotes an unmeasured gene in a screen column or absence of a qualifying readout in the four reviewed sources in an assay column. "
                "Letters mark assay availability, including null observations. Only SATB1 had a direct suppression assay in these sources. "
                "Kolb et al. reported reduced suppression after SATB1 KO and preserved Treg recovery (Fig. 5C,D). The Fig. 5C legend specifies "
                "four biological replicates; the mapping to nine workbook columns per arm is unresolved. Variant dosage and patient Treg function "
                "remain untested by these KO results. See Supplementary Table S12 for variant scope, source figures and readout-specific directions.")
    short = "Supplementary Figure S6. Disease mechanisms and published functional readouts."

    with tempfile.TemporaryDirectory(prefix="s6-vector-") as tmp:
        export(fig,Path(tmp),"S6_vector")
        with fitz.open(Path(tmp)/"S6_vector.pdf") as fig, fitz.open() as outdoc:
            size = fitz.Rect(0,0,595.2755737304688,841.8897705078125)
            page = outdoc.new_page(width=size.width,height=size.height)
            clip = fig[0].rect
            x = (size.width-clip.width)/2
            page.show_pdf_page(fitz.Rect(x,13,x+clip.width,13+clip.height),fig,0)
            textat(page, x, 640, short, size=10.2, bold=True, ink="navy")
            yy = paragraph(page, x, 659, clip.width, "(a) " + caption1, size=8.65, leading=12.15) + 2.5
            yy = paragraph(page, x, yy, clip.width, "(b) " + caption2, size=8.65, leading=12.15) + 2.5
            yy = paragraph(page, x, yy, clip.width, "Source PMIDs: Chen, 40140585 · Umhoefer, 41237776 · Schumann, 32989329 · Kolb, 42286212.", size=8.25, leading=11)
            assert yy < 823, yy
            textat(page, size.width-33, 823, "6", size=8.5, align="right", ink="gray")
            single = out / "Supplementary_Figure_S6.pdf"
            savepdf(outdoc, single)
            outdoc[0].get_pixmap(dpi=300,alpha=False).save(single.with_suffix(".png"))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_root",type=Path)
    p.add_argument("out_dir",type=Path)
    p.add_argument("--historical-figure3",type=Path)
    a=p.parse_args()
    data = a.data_root.resolve()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    if (data / REPLAY_MAP).is_file():
        with replay_view(data) as compatibility:
            render(compatibility, a.out_dir, a.historical_figure3)
    else:
        render(data, a.out_dir, a.historical_figure3)


def render(data, output, underlay):
    source = data / "data/r7/figure_sources"
    render_fig3(source, output, underlay)
    render_s6(source, output)


if __name__ == "__main__":
    main()
