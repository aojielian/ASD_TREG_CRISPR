#!/usr/bin/env python3
"""Assemble the five supplementary PDFs with their public legends."""
import argparse,io,json
from pathlib import Path
from pypdf import PdfReader,PdfWriter,PageObject,Transformation
from reportlab.pdfgen import canvas
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import Paragraph
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from matplotlib.font_manager import findfont

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--legends',type=Path,default=Path(__file__).with_name('figure_legends.json'))
    ap.add_argument('--font-path',type=Path,default=None)
    args=ap.parse_args();payload=json.loads(args.legends.read_text())
    font_path=args.font_path or Path(findfont('DejaVu Sans'))
    bold_path=font_path.with_name('DejaVuSans-Bold.ttf')
    if not bold_path.exists(): bold_path=Path(findfont('DejaVu Sans:style=bold'))
    pdfmetrics.registerFont(TTFont('SuppDejaVu',str(font_path)))
    pdfmetrics.registerFont(TTFont('SuppDejaVu-Bold',str(bold_path)))
    pdfmetrics.registerFontFamily('SuppDejaVu',normal='SuppDejaVu',bold='SuppDejaVu-Bold')
    style=ParagraphStyle('legend',fontName='SuppDejaVu',fontSize=9,leading=12,
                         alignment=TA_JUSTIFY,textColor='#183B56',spaceAfter=0)
    W,H=595.2756,841.8898;margin=36;figure_width=W-2*margin
    writer=PdfWriter();geometry=[]
    for i in range(1,6):
        reader=PdfReader(args.output_dir/f'Supplementary_Figure_S{i}_FINAL.pdf')
        source=reader.pages[0];sw,sh=float(source.mediabox.width),float(source.mediabox.height)
        scale=figure_width/sw;figure_height=sh*scale
        top=25;legend_top=top+figure_height+12
        legend=Paragraph(payload[str(i)]['html'],style);_,lh=legend.wrap(figure_width,H)
        if legend_top+lh>H-35:
            raise ValueError(f'S{i} figure + legend cannot fit without resizing: {legend_top+lh:.1f} > {H-35:.1f}')
        page=PageObject.create_blank_page(width=W,height=H)
        page.merge_transformed_page(source,Transformation().scale(scale).translate(margin,H-top-figure_height))
        bio=io.BytesIO();c=canvas.Canvas(bio,pagesize=(W,H))
        legend.drawOn(c,margin,H-legend_top-lh)
        c.setFont('SuppDejaVu',7);c.setFillColorRGB(.35,.40,.44);c.drawRightString(W-margin,19,str(i))
        c.showPage();c.save();bio.seek(0)
        page.merge_page(PdfReader(bio).pages[0]);writer.add_page(page)
        geometry.append({'figure':f'S{i}','page':i,'figure_rect_top_coordinates_pt':[margin,top,W-margin,top+figure_height],
                         'legend_rect_top_coordinates_pt':[margin,legend_top,W-margin,legend_top+lh],
                         'legend_font_pt':9,'legend_leading_pt':12,'figure_scale':scale,'vector_composite':True})
    target=args.output_dir/'Supplementary_Figures_FINAL_WITH_LEGENDS.pdf'
    with target.open('wb') as f:writer.write(f)
    qc=args.output_dir.parent/'04_qc';qc.mkdir(parents=True,exist_ok=True)
    (qc/'SUPPLEMENT_ASSEMBLY_GEOMETRY.json').write_text(json.dumps(geometry,indent=2))
    print(f'ASSEMBLED={target}')
if __name__=='__main__':main()
