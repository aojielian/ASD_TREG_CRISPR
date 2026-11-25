#!/usr/bin/env python3
"""Assemble Supplementary Tables S1–S20 from the release's frozen OOXML sources.

Requires Python 3.9+ and only the standard library. Does not calculate statistics.
"""

import argparse
import csv
import hashlib
import io
import re
from copy import deepcopy
from pathlib import Path
from tempfile import NamedTemporaryFile
from xml.etree import ElementTree as ET
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

SPREADSHEETS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
OFFICE_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACKAGE_REL = "http://schemas.openxmlformats.org/package/2006/relationships"
CONTENT_TYPES = "http://schemas.openxmlformats.org/package/2006/content-types"
SOURCE_ARCHIVE = "SUPPLEMENTARY_TABLE_WORKBOOK_SOURCE.zip"
VOLUMES = ("S1-S12", "S13", "S14-S16", "S17-S20")
ET.register_namespace("x", SPREADSHEETS)
ET.register_namespace("r", OFFICE_REL)

def tag(local):
    return f"{{{SPREADSHEETS}}}{local}"

def xml(element):
    return ET.tostring(element, encoding="utf-8", xml_declaration=True)

def add(parent, element, **attributes):
    return ET.SubElement(parent, tag(element), {key: str(value) for key, value in attributes.items()})

def package_xml(element, namespace):
    ET.register_namespace("", namespace)
    return xml(element)

def column_letter(index):
    letters = ""
    while index:
        index, digit = divmod(index - 1, 26)
        letters = chr(65 + digit) + letters
    return letters

def digest(data):
    return hashlib.sha256(data).hexdigest()

def write_member(archive, name, content):
    member = ZipInfo(name, date_time=(2026, 11, 25, 0, 0, 0))
    member.external_attr = 0o100644 << 16
    archive.writestr(member, content, compress_type=ZIP_DEFLATED, compresslevel=9)

def read_tsv(data):
    return list(csv.reader(io.StringIO(data.decode("utf-8"), newline=""), delimiter="\t"))

def style_collection(root, name):
    part = root.find(tag(name))
    if part is None:
        raise ValueError(f"Missing stylesheet collection: {name}")
    return part

def extend_styles(target, incoming):
    """Append the volume's styles, remapping local resources and number formats."""
    if xml(style_collection(target, "cellStyleXfs")) != xml(style_collection(incoming, "cellStyleXfs")):
        raise ValueError("Unexpected base-style definitions across the frozen volumes")
    numeric_formats = style_collection(target, "numFmts")
    next_format = max(int(fmt.get("numFmtId")) for fmt in numeric_formats) + 1
    by_code = {fmt.get("formatCode"): int(fmt.get("numFmtId")) for fmt in numeric_formats}
    format_map = {}
    for fmt in style_collection(incoming, "numFmts"):
        format_id, code = int(fmt.get("numFmtId")), fmt.get("formatCode")
        if code not in by_code:
            add(numeric_formats, "numFmt", numFmtId=next_format, formatCode=code)
            by_code[code] = next_format
            next_format += 1
        format_map[format_id] = by_code[code]
    numeric_formats.set("count", str(len(numeric_formats)))

    offsets = {}
    for resource, collection in (("fontId", "fonts"), ("fillId", "fills"), ("borderId", "borders")):
        part = style_collection(target, collection)
        offsets[resource] = len(part)
        part.extend(deepcopy(item) for item in style_collection(incoming, collection))
        part.set("count", str(len(part)))

    styles = style_collection(target, "cellXfs")
    offset = len(styles)
    for source_style in style_collection(incoming, "cellXfs"):
        item = deepcopy(source_style)
        for resource, base in offsets.items():
            item.set(resource, str(base + int(item.get(resource, "0"))))
        format_id = int(item.get("numFmtId", "0"))
        if format_id >= 164:
            item.set("numFmtId", str(format_map[format_id]))
        styles.append(item)
    styles.set("count", str(len(styles)))
    return offset

def remap_sheet_styles(root, offset):
    for cell in root.iter(tag("c")):
        if cell.get("t") == "s":
            raise ValueError("Unexpected external shared-string reference")
        cell.set("s", str(int(cell.get("s", "0")) + offset))
    for node_name, style_attribute in (("row", "s"), ("col", "style")):
        for item in root.iter(tag(node_name)):
            if style_attribute in item.attrib:
                item.set(style_attribute, str(int(item.get(style_attribute)) + offset))

def meta_worksheet(rows, widths, text_style, header_style, integer_columns=()):
    """Build the unified public Index or Dictionary from its TSV."""
    if not rows or any(len(row) != len(widths) for row in rows):
        raise ValueError("Nonrectangular index or dictionary")
    root = ET.Element(tag("worksheet"))
    views = add(root, "sheetViews")
    view = add(views, "sheetView", workbookViewId=0, showGridLines=0)
    add(view, "pane", ySplit=1, topLeftCell="A2", activePane="bottomLeft", state="frozen")
    add(view, "selection", pane="bottomLeft", activeCell="A2", sqref="A2")
    add(root, "sheetFormatPr", defaultRowHeight=18)
    cols = add(root, "cols")
    for i, width in enumerate(widths, 1):
        add(cols, "col", min=i, max=i, width=width, customWidth=1)
    cells = add(root, "sheetData")
    for row_number, values in enumerate(rows, 1):
        max_chars = max((len(value) / max(1, width - 2) for value, width in zip(values, widths)), default=0)
        height = 28 if row_number == 1 else max(22, int(max_chars + 0.999) * 14 + 8)
        row = add(cells, "row", r=row_number, ht=height, customHeight=1)
        for index, value in enumerate(values, 1):
            address = f"{column_letter(index)}{row_number}"
            style = header_style if row_number == 1 else text_style
            if row_number > 1 and index in integer_columns and value:
                cell = add(row, "c", r=address, s=style, t="n")
                add(cell, "v").text = value
            else:
                cell = add(row, "c", r=address, s=style, t="str")
                add(cell, "v").text = value
    add(root, "autoFilter", ref=f"A1:{column_letter(len(widths))}{len(rows)}")
    return root

def metadata_styles(styles):
    fonts = style_collection(styles, "fonts")
    fills = style_collection(styles, "fills")
    records = style_collection(styles, "cellXfs")
    normal_font = len(fonts)
    font = add(fonts, "font")
    add(font, "name", val="Arial")
    add(font, "sz", val=10)
    add(font, "color", rgb="FF222222")
    heading_font = len(fonts)
    font = add(fonts, "font")
    add(font, "name", val="Arial")
    add(font, "sz", val=10)
    add(font, "b")
    add(font, "color", rgb="FFFFFFFF")
    fill = add(add(fills, "fill"), "patternFill", patternType="solid")
    add(fill, "fgColor", rgb="FF265661")
    add(fill, "bgColor", indexed=64)
    normal_style = len(records)
    xf = add(records, "xf", numFmtId=0, fontId=normal_font, fillId=0, borderId=0, xfId=0,
             applyFont=1, applyAlignment=1)
    add(xf, "alignment", vertical="top", wrapText=1)
    header_style = len(records)
    xf = add(records, "xf", numFmtId=0, fontId=heading_font, fillId=len(fills)-1, borderId=0, xfId=0,
             applyFont=1, applyFill=1, applyAlignment=1)
    add(xf, "alignment", vertical="center", wrapText=1)
    for part in (fonts, fills, records):
        part.set("count", str(len(part)))
    return normal_style, header_style

def relationship(parent, name, destination, index):
    ET.SubElement(parent, f"{{{PACKAGE_REL}}}Relationship",
                  {"Id": index, "Type": f"{OFFICE_REL}/{name}", "Target": destination})

def write_workbook(path, theme, styles, sheets):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    workbook = ET.Element(tag("workbook"))
    book_view = add(workbook, "bookViews")
    add(book_view, "workbookView", activeTab=0, firstSheet=0)
    sheet_list = add(workbook, "sheets")
    root_rels = ET.Element(f"{{{PACKAGE_REL}}}Relationships")
    relationship(root_rels, "officeDocument", "xl/workbook.xml", "rId1")
    book_rels = ET.Element(f"{{{PACKAGE_REL}}}Relationships")
    relationship(book_rels, "styles", "styles.xml", "rId1")
    relationship(book_rels, "theme", "theme/theme1.xml", "rId2")
    ct = ET.Element(f"{{{CONTENT_TYPES}}}Types")
    for ext, content_type in (("rels", "application/vnd.openxmlformats-package.relationships+xml"),
                              ("xml", "application/xml")):
        ET.SubElement(ct, f"{{{CONTENT_TYPES}}}Default", Extension=ext, ContentType=content_type)
    def override(part, content_type):
        ET.SubElement(ct, f"{{{CONTENT_TYPES}}}Override", PartName=f"/xl/{part}", ContentType=content_type)
    prefix = "application/vnd.openxmlformats-officedocument."
    override("workbook.xml", prefix + "spreadsheetml.sheet.main+xml")
    override("styles.xml", prefix + "spreadsheetml.styles+xml")
    override("theme/theme1.xml", prefix + "theme+xml")
    sheets_xml = []
    for number, (name, content) in enumerate(sheets, 1):
        if len(name) > 31 or re.search(r"[\\/*?:\[\]]", name):
            raise ValueError(f"Invalid tab name: {name}")
        relation_id = f"rId{number + 2}"
        sheet_item = add(sheet_list, "sheet", name=name, sheetId=number)
        sheet_item.set(f"{{{OFFICE_REL}}}id", relation_id)
        part = f"worksheets/sheet{number}.xml"
        relationship(book_rels, "worksheet", part, relation_id)
        override(part, prefix + "spreadsheetml.worksheet+xml")
        sheets_xml.append((f"xl/{part}", xml(content)))
    with NamedTemporaryFile(dir=path.parent, suffix=".xlsx", delete=False) as temp:
        temporary = Path(temp.name)
    try:
        with ZipFile(temporary, "w") as archive:
            for member, content in [
                ("[Content_Types].xml", package_xml(ct, CONTENT_TYPES)),
                ("_rels/.rels", package_xml(root_rels, PACKAGE_REL)),
                ("xl/workbook.xml", xml(workbook)),
                ("xl/_rels/workbook.xml.rels", package_xml(book_rels, PACKAGE_REL)),
                ("xl/styles.xml", xml(styles)), ("xl/theme/theme1.xml", theme), *sheets_xml,
            ]:
                write_member(archive, member, content)
        temporary.chmod(0o644)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)

def rebuild(source_dir, output):
    with ZipFile(Path(source_dir) / SOURCE_ARCHIVE) as sources:
        path_rows = read_tsv(sources.read("SHA256SUMS.tsv"))
        expected = {path: sha for path, sha in path_rows[1:]}
        if len(sources.namelist()) != len(expected) + 1 or set(sources.namelist()) != set(expected) | {"SHA256SUMS.tsv"}:
            raise ValueError("Workbook source archive membership does not match its checksums")

        def verified(member):
            data = sources.read(member)
            if digest(data) != expected[member]:
                raise ValueError(f"Workbook source checksum mismatch: {member}")
            return data

        sheet_rows = read_tsv(verified("SHEETS.tsv"))
        if sheet_rows[0] != ["table_id", "worksheet", "volume", "xml_path"]:
            raise ValueError("Unexpected worksheet source schema")
        index = read_tsv(verified("TABLE_INDEX.tsv"))
        dictionary = read_tsv(verified("TABLE_DICTIONARY.tsv"))
        if index[1:77] and [row[1] for row in index[1:77]] != [row[1] for row in sheet_rows[1:]]:
            raise ValueError("Unified index differs from workbook science-tab order")
        if len(sheet_rows) != 77 or len(set(row[1] for row in sheet_rows[1:])) != 76:
            raise ValueError("The complete release requires 76 unique scientific worksheets")

        theme = verified("theme.xml")
        styles = ET.fromstring(verified(f"styles/{VOLUMES[0]}.xml"))
        offsets = {VOLUMES[0]: 0}
        for volume in VOLUMES[1:]:
            incoming = ET.fromstring(verified(f"styles/{volume}.xml"))
            offsets[volume] = extend_styles(styles, incoming)
        text_format, header_format = metadata_styles(styles)
        shelves = [("Index", meta_worksheet(index, (14, 31, 14, 13, 47, 37, 48, 52, 46, 63, 60),
                                            text_format, header_format, integer_columns=(3, 4)))]
        for _, name, volume, member in sheet_rows[1:]:
            sheet = ET.fromstring(verified(member))
            remap_sheet_styles(sheet, offsets[volume])
            shelves.append((name, sheet))
        shelves.append(("Dictionary", meta_worksheet(dictionary, (14, 30, 15, 40, 34, 25, 70),
                                                      text_format, header_format)))
    write_workbook(output, theme, styles, shelves)
    print(f"Wrote {len(shelves) - 2} scientific sheets and two metadata sheets to {output}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True, help="Extracted data path tables/submission")
    parser.add_argument("--output", type=Path, required=True, help="Destination workbook path")
    args = parser.parse_args()
    rebuild(args.source_dir, args.output)

if __name__ == "__main__":
    main()
