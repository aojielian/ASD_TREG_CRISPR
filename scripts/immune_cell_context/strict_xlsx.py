"""Historical strict OOXML reader for the two GSE286470 workbooks."""

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


def read_strict_xlsx(path: Path):
    """Minimal read-only OOXML reader for strict-conformance XLSX files."""
    ns = {'main':'http://purl.oclc.org/ooxml/spreadsheetml/main', 'rel':'http://purl.oclc.org/ooxml/officeDocument/relationships', 'pkgrel':'http://schemas.openxmlformats.org/package/2006/relationships'}
    with zipfile.ZipFile(path) as z:
        names=set(z.namelist())
        ss=[]
        if 'xl/sharedStrings.xml' in names:
            root=ET.fromstring(z.read('xl/sharedStrings.xml'))
            for si in root.findall('.//{*}si'):
                ss.append(''.join(t.text or '' for t in si.findall('.//{*}t')))
        wb=ET.fromstring(z.read('xl/workbook.xml'))
        relroot=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        rels={r.attrib.get('Id'):r.attrib.get('Target') for r in relroot.findall('.//{*}Relationship')}
        result={}
        for sh in wb.findall('.//{*}sheet'):
            title=sh.attrib.get('name','')
            rid=sh.attrib.get('{http://purl.oclc.org/ooxml/officeDocument/relationships}id') or sh.attrib.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
            target=rels.get(rid,'worksheets/sheet1.xml')
            if target.startswith('/'): member=target.lstrip('/')
            elif target.startswith('xl/'): member=target
            else: member='xl/'+target
            if member not in names: continue
            root=ET.fromstring(z.read(member)); rows=[]; maxcol=0
            for rr in root.findall('.//{*}sheetData/{*}row'):
                cells={}
                for c in rr.findall('{*}c'):
                    ref=c.attrib.get('r','A1'); m=re.match(r'([A-Z]+)([0-9]+)',ref)
                    if not m: continue
                    letters=m.group(1); col=0
                    for ch in letters: col=col*26+(ord(ch)-64)
                    maxcol=max(maxcol,col)
                    v=c.find('{*}v'); val=''
                    if c.attrib.get('t')=='s' and v is not None:
                        try: val=ss[int(v.text or 0)]
                        except Exception: val=v.text or ''
                    elif c.attrib.get('t')=='inlineStr': val=''.join(t.text or '' for t in c.findall('.//{*}t'))
                    elif v is not None: val=v.text or ''
                    cells[col-1]=val
                row=['']*(max(cells.keys(),default=-1)+1)
                for i,v in cells.items(): row[i]=v
                rows.append(row)
            result[title]=rows
        return result
