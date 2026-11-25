"""Rendering utilities for the editor revision. No statistical estimators are fitted here."""
from pathlib import Path
from hashlib import sha256
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

ASD = "#243F65"
TREG = "#15776F"
ORANGE = "#B86128"
GRAY = "#A6ADB5"
DARK = "#202D39"
GRID = "#DCE1E3"
LIGHT = "#F2F4F4"
SPECIES = "#76568A"
GENES = "CREBBP CUL3 DYRK1A EP300 FOXP1 KDM3B KDM6B KMT5B PHIP SIN3A SMARCA4 TBL1XR1 TLE3 ZMYND8".split()
REFS = "FOXP3 IL2RA CTLA4 IKZF2".split()
PANEL = GENES + REFS
SEQ = LinearSegmentedColormap.from_list("teal_amount", ["#F7FAF7","#D3E8DF","#88C2B2","#398F85","#125E61","#173F4D"])
DIV = LinearSegmentedColormap.from_list("signed", ["#97522D","#DAA06C","#F9F7F1","#7CBDB0","#135F65"])
for cm in (SEQ, DIV):
    cm.set_bad("#D1D3D5")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.titlesize":9,
    "axes.labelsize":8,"axes.linewidth":0.65,"axes.spines.top":False,"axes.spines.right":False,
    "axes.edgecolor":"#7A858B","text.color":DARK,"axes.labelcolor":DARK,
    "xtick.color":DARK,"ytick.color":DARK,"xtick.major.size":3,"ytick.major.size":3,
    "xtick.major.width":0.6,"ytick.major.width":0.6,"ytick.labelsize":7.5,"xtick.labelsize":7.5,
    "pdf.fonttype":42,"ps.fonttype":42,"pdf.compression":9,"svg.fonttype":"none",
    "legend.frameon":False,"legend.fontsize":7.2,"savefig.facecolor":"white","axes.facecolor":"white"})
def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    return parser.parse_args()


def read(path,**kwargs):
    return pd.read_csv(path,sep="\t",low_memory=False,**kwargs)
def write(frame,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    frame.to_csv(path,index=False,sep="\t",na_rep="NA",float_format="%.17g",
                 compression={"method":"gzip","mtime":0} if path.suffix==".gz" else None)
def hashfile(p):
    return sha256(p.read_bytes()).hexdigest()
def pval(x):
    return f"{x:.2g}" if x<0.001 else f"{x:.3f}".rstrip("0").rstrip(".")
def header(fig,x,y,letter,label):
    fig.text(x,y,letter,weight="bold",size=13,va="baseline",color=DARK)
    fig.text(x+.034,y,label,weight="bold",size=9,va="baseline",color=DARK)
def grid(ax,axis="x"):
    ax.set_axisbelow(True)
    ax.grid(axis=axis,color=GRID,lw=.5)
def mark(color,label,marker="o",mfc=None,ls="none"):
    return Line2D([],[],ls=ls,marker=marker,ms=5,color=color,
                  markerfacecolor=color if mfc is None else mfc,label=label)
def coords(data,category):
    return data.loc[data["record_type"].eq(category)].copy()
def save_main(fig,stage,num):
    out=stage/"04_figures"
    out.mkdir(exist_ok=True)
    title=f"Figure_{num}"
    fig.savefig(out/f"{title}.pdf",metadata={"Title":title,"CreationDate":None,"ModDate":None})
    plt.close(fig)
    tiff_from_pdf(out/f"{title}.pdf")
def truth(series):
    return series.astype(str).str.lower().isin(["true","yes","1"])
def legacy_inputs(stage, data):
    src = stage / "04_figures/source_data"
    index = []
    for name in ("F1_chen.tsv.gz", "F2_cross_screen.tsv", "F2_models.tsv", "S7_moreno.tsv.gz"):
        file = data / "figure_sources/R10" / name
        frame = read(file)
        write(frame, src / name)
        index.append((str(file.relative_to(data)), hashfile(file), len(frame), name))
    for name, output in [("PHENOTYPE_PREFERENCE_R9.tsv", "F3_foxp3_tcr.tsv"),
                         ("PHENOTYPE_PREFERENCE_NULL.tsv", "F3_foxp3_tcr_null.tsv.gz"),
                         ("LOF_DIRECTIONAL_ENRICHMENT.tsv", "F2_lof.tsv")]:
        file = data / "data/r9/results" / name
        frame = read(file)
        write(frame, src / output)
        index.append((str(file.relative_to(data)), hashfile(file), len(frame), output))
    frozen = data / "figure_sources/frozen_r10/R10_Supplementary_Figures.pdf"
    index.append((str(frozen.relative_to(data)), hashfile(frozen), 7, "Frozen S1–S7 vector pages"))
    write(pd.DataFrame(index, columns=["source_path", "source_sha256", "rows_used", "render_source"]),
          src / "FROZEN_FIGURE_INPUTS.tsv")
    return src


def vector_cells(ax,values,cmap,vmin=None,vmax=None):
    array=np.asarray(values,dtype=float)
    rows,cols=array.shape
    artist=ax.pcolormesh(np.arange(cols+1)-.5,np.arange(rows+1)-.5,array,
                         cmap=cmap,vmin=vmin,vmax=vmax,shading="flat",linewidth=0,
                         edgecolors="face",antialiased=False,rasterized=False)
    ax.set_xlim(-.5,cols-.5)
    ax.set_ylim(rows-.5,-.5)
    return artist

def tiff_from_pdf(path):
    # Generate the companion TIFF from the delivered vector page itself.
    import fitz
    from PIL import Image
    with fitz.open(path) as doc:
        pix=doc[0].get_pixmap(dpi=600,colorspace=fitz.csRGB,alpha=False)
        image=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
        image.save(path.with_suffix(".tif"),compression="tiff_lzw",dpi=(600,600))
