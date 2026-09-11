"""Shared rendering primitives for the final main figures."""
from pathlib import Path
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

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

def panel(fig, rect, letter, title, title_y=1.055, header_left=None):
    ax = fig.add_axes(rect)
    hx=rect[0] if header_left is None else header_left
    hy=rect[1]+title_y*rect[3]
    fig.text(max(.026,hx-.032),hy,letter,fontsize=12,
             fontweight='bold',ha='right',va='bottom')
    fig.text(hx,hy,title,fontsize=9,fontweight='bold',ha='left',va='bottom')
    return ax

def clean(ax, grid='y'):
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(axis=grid, color=LIGHT, lw=.5)
    return ax

def blank(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    return ax

def save(fig, output_root, number):
    target = Path(output_root); target.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    texts=[]
    from matplotlib.text import Text
    for artist in fig.findobj(match=Text):
        value=artist.get_text()
        if not value or not artist.get_visible():continue
        box=artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        texts.append({'text':value,'font_size_pt':artist.get_fontsize(),
                      'x0':box.x0,'y0':box.y0,'x1':box.x1,'y1':box.y1,
                      'outside_canvas':bool(box.x0<0 or box.y0<0 or box.x1>1 or box.y1>1)})
    qc=target.parent/'04_qc';qc.mkdir(parents=True,exist_ok=True)
    with (qc/f'Figure{number}_text_artists.tsv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=texts[0].keys(),delimiter='\t')
        writer.writeheader();writer.writerows(texts)
    with (qc/f'Figure{number}_render_metadata.json').open('w') as f:
        json.dump({'figure':number,'size_inches':fig.get_size_inches().tolist(),
                   'png_dpi':600,'text_artist_count':len(texts),
                   'outside_canvas':[x for x in texts if x['outside_canvas']],
                   'font_family':'DejaVu Sans','pdf_fonttype':42,'svg_fonttype':'none',
                   'renderer':'public figure renderer'},f,indent=2)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(target / f'Figure{number}_FINAL.{ext}', dpi=600)
    plt.close(fig)

def pformat(value, digits=3):
    x = float(value)
    if x == 0:
        return '0'
    if x >= .001:
        return f'{x:.{digits}g}'
    exponent = int(np.floor(np.log10(x)))
    mantissa = x / 10 ** exponent
    return rf'${mantissa:.{digits-1}f}\times10^{{{exponent}}}$'

def marker_legend(ax, entries, **kwargs):
    handles = [Line2D([], [], marker=marker, ls='none', markersize=4,
                      color=color, label=label) for label,color,marker in entries]
    return ax.legend(handles=handles, frameon=False, fontsize=7,
                     handletextpad=.35, borderpad=0, **kwargs)

def ecdf(ax, values, color, label, lw=1.2):
    values = np.sort(np.asarray(values, dtype=float))
    values = values[np.isfinite(values)]
    if len(values) == 0:
        raise ValueError(f'Empty observed distribution: {label}')
    ax.step(values, np.arange(1,len(values)+1)/len(values), where='post',
            color=color, label=label, lw=lw)

def table_rule(ax, y, xmin=0, xmax=1, color=LIGHT, lw=.55):
    ax.plot([xmin,xmax], [y,y], color=color, lw=lw, clip_on=False)
