"""Check real rendered pages, hashes, inherited pixels, captions, font embedding and 600-dpi TIFF pixels."""
import fitz
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from figure_common import read,write,sha,bytehash
import pandas as pd

def rgba(page,dpi):
    p=page.get_pixmap(dpi=dpi,alpha=False,colorspace=fitz.csRGB)
    return Image.frombytes("RGB",(p.width,p.height),p.samples)

def contact(images,folder,prefix,columns=2,width=650,maxrows=2):
    font=ImageFont.truetype(findfont("DejaVu Sans"),19)
    for start in range(0,len(images),columns*maxrows):
        group=images[start:start+columns*maxrows]
        thumbs=[]
        for name,im in group:
            thumb=im.copy(); thumb.thumbnail((width,950))
            thumbs.append((name,thumb))
        rows=(len(thumbs)+columns-1)//columns
        heights=[max(v[1].height for v in thumbs[i*columns:(i+1)*columns])+37 for i in range(rows)]
        canvas=Image.new("RGB",((width+7)*columns,sum(heights)),(229,232,233))
        draw=ImageDraw.Draw(canvas);y=0
        for i,h in enumerate(heights):
            for j,(name,im) in enumerate(thumbs[i*columns:(i+1)*columns]):
                x=(width+7)*j
                draw.text((x+6,y+4),name,font=font,fill=(36,50,60))
                canvas.paste(im,(x+(width-im.width)//2,y+30))
            y+=h
        canvas.save(folder/f"{prefix}_{start+1:02}-{start+len(group):02}.png")

def verify(stage,r10,out):
    out.mkdir(parents=True,exist_ok=True)
    records=[]; mains=[]; supp=[]
    sf=stage/"05_supplementary"/"Supplementary_Figures_S1-S23.pdf"
    with fitz.open(sf) as doc,fitz.open(r10/"R10_Supplementary_Figures.pdf") as old:
        assert len(doc)==23 and len(old)==14
        all_fonts=set()
        for i,page in enumerate(doc):
            for xref,*_ in doc.get_page_fonts(i):
                if xref and xref not in all_fonts:
                    assert len(doc.extract_font(xref)[3])>0, (i,xref)
                    all_fonts.add(xref)
            num=i+1
            pix=rgba(page,92)
            oldmatch=bytehash(rgba(old[i],92).tobytes())==bytehash(pix.tobytes()) if i<14 else None
            assert oldmatch is not False
            text=page.get_text()
            assert f"Figure S{num}." in text
            blocks=page.get_text("blocks")
            cap=[b for b in blocks if f"Figure S{num}." in str(b[4])]
            assert len(cap)==1 and min(b[1] for b in cap)>page.rect.height*.40
            chars=[w for w in page.get_text("words") if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
            assert not chars, (num,chars[:2])
            path=out/f"S{num:02}.png";pix.save(path)
            records.append([str(sf.relative_to(stage)),num,sha(sf),sha(path),oldmatch,"PDF; full-page caption",None])
            supp.append((f"Supplementary Figure S{num}",pix))
    for n in range(1,5):
        pdf=stage/"04_figures"/f"Figure_{n}.pdf";tif=pdf.with_suffix(".tif")
        with fitz.open(pdf) as doc,Image.open(tif) as im:
            assert len(doc)==1 and abs(float(im.info["dpi"][0])-600)<.1 and abs(float(im.info["dpi"][1])-600)<.1
            page=doc[0]
            oldmatch=sha(r10/pdf.name)==sha(pdf) if n in (1,2) else None
            assert oldmatch is not False
            fonts=doc.get_page_fonts(0); assert fonts
            for xref,*_ in fonts:
                if xref: assert len(doc.extract_font(xref)[3])>0, (pdf.name,xref)
            raster=rgba(page,600)
            assert im.size==raster.size and bytehash(im.convert("RGB").tobytes())==bytehash(raster.tobytes()),n
            pix=rgba(page,120);path=out/f"F{n}.png";pix.save(path)
            records.append([str(pdf.relative_to(stage)),1,sha(pdf),sha(path),oldmatch,"PDF; embedded fonts",None])
            records.append([str(tif.relative_to(stage)),1,sha(tif),sha(path),oldmatch,"600 dpi RGB; full-resolution PDF pixel agreement",f"{im.width} x {im.height}"])
            mains.append((f"Figure {n}",pix))
    contact(mains,out,"main",2,750)
    contact(supp,out,"supp",2,710)
    write(pd.DataFrame(records,columns=["artifact","page","artifact_sha256","page_preview_sha256","identical_to_fixed_vector_source","export_check","TIFF_dimensions"]),stage/"07_qa"/"FIGURE_RENDER_QA.tsv")
