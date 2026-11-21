#!/usr/bin/env python
"""Step 15: build Supplementary Table S1 in Word, Excel and CSV.

Core table, as requested: accession, virus name, abbreviation, strain, isolate,
country, region.  An extended sheet adds host, year, genome length and the
evidence behind each strain call, so a reader can see which assignments rest on
a complete genome and which are provisional.
"""
import csv, os, json
from collections import Counter, OrderedDict

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR

meta = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "03_curated_metadata.csv"), encoding="utf-8-sig"))}
sass = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "06_strain_assignment.csv"), encoding="utf-8-sig"))}
frag = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "13b_partial_assignment.csv"), encoding="utf-8-sig"))}
cfg = json.load(open(os.path.join(OUT, "10_analysis_set.json")))
BEYDV = set(cfg["excluded_beydv"])

def strain_cell(a):
    """final strain designation plus the evidence class behind it"""
    m = meta[a]
    ann = m["strain"]
    if a in BEYDV:
        return "n.a.", "not CpCDV: bean yellow dwarf virus lineage"
    if len(ann) == 1 and ann.isalpha():
        return ann, "as published"
    if a in sass:
        s = sass[a]
        if s["strain_final"] == "unassigned":
            return "unassigned", (f"below 94% to every strain reference "
                                  f"(closest {s['nearest_identity']}%)")
        return s["strain_final"], (f"assigned here, complete genome "
                                   f"({s['nearest_identity']}% to {s['nearest']})")
    if a in frag:
        f = frag[a]
        if f["assignment_confident"] == "yes":
            return f["nearest_strain"] + " (prov.)", (
                f"provisional, partial sequence ({f['identity_to_nearest']}%)")
        return "undetermined", (
            f"partial sequence; {f['nearest_strain']} and {f['second_strain']} "
            f"within {f['margin']}%")
    if ann:
        return ann, "as published"
    return "undetermined", "no strain information"

rows = []
for a in sorted(meta, key=lambda x: (meta[x]["region"], meta[x]["country"],
                                     meta[x]["strain"] or "~", x)):
    m = meta[a]
    s, ev = strain_cell(a)
    vname = ("bean yellow dwarf virus" if a in BEYDV
             else "chickpea chlorotic dwarf virus")
    abbr = "BeYDV" if a in BEYDV else "CpCDV"
    rows.append(OrderedDict([
        ("Accession", a),
        ("Virus name", vname),
        ("Abbreviation", abbr),
        ("Strain", s),
        ("Isolate / clone", m["isolate"]),
        ("Country", m["country"]),
        ("Region", m["subregion"]),
        # extended fields
        ("Broad region", m["region"]),
        ("Locality", m["locality"]),
        ("Host (species)", m["host_species"]),
        ("Host (common)", m["host_common"]),
        ("Host family", m["host_family"]),
        ("Collection year", m["collection_year"]),
        ("Length (nt)", m["length"]),
        ("Sequence type", m["molecule_class"]),
        ("Basis of strain assignment", ev),
        ("Reference", m["reference_title"]),
        ("Journal", m["reference_journal"]),
    ]))

CORE = ["Accession", "Virus name", "Abbreviation", "Strain", "Isolate / clone",
        "Country", "Region"]

# ---------------------------------------------------------------- CSV
p = os.path.join(OUT, "Table_S1_core.csv")
with open(p, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=CORE, extrasaction="ignore")
    w.writeheader(); w.writerows(rows)
p2 = os.path.join(OUT, "Table_S1_extended.csv")
with open(p2, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("wrote", p, "and", p2)

# ---------------------------------------------------------------- Excel
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

wb = Workbook()
HDR = PatternFill("solid", fgColor="1F3864")
HF = Font(color="FFFFFF", bold=True, size=10, name="Calibri")
thin = Side(style="thin", color="BFBFBF")

def sheet(ws, cols, data, widths):
    ws.append(cols)
    for c in range(1, len(cols) + 1):
        cell = ws.cell(1, c); cell.fill = HDR; cell.font = HF
        cell.alignment = Alignment(vertical="center", wrap_text=True)
    for r in data:
        ws.append([r[c] for c in cols])
    for i, wd in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = wd
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{len(data)+1}"
    for row in ws.iter_rows(min_row=1, max_row=len(data) + 1, max_col=len(cols)):
        for cell in row:
            cell.border = Border(bottom=thin)
            if cell.row > 1:
                cell.font = Font(size=9.5, name="Calibri")
                cell.alignment = Alignment(vertical="top", wrap_text=False)
    ws.row_dimensions[1].height = 30

ws = wb.active; ws.title = "Table S1 (core)"
sheet(ws, CORE, rows, [13, 32, 13, 17, 26, 14, 18])
ws2 = wb.create_sheet("Table S1 (extended)")
sheet(ws2, list(rows[0].keys()), rows,
      [13, 30, 12, 16, 24, 13, 17, 20, 22, 24, 22, 17, 13, 10, 20, 42, 52, 26])

# a summary sheet
ws3 = wb.create_sheet("Summary")
summ = [
    ["Supplementary Table S1 - summary of the sequence dataset", ""],
    ["", ""],
    ["Total GenBank records analysed", len(rows)],
    ["Complete genomes", sum(1 for r in rows if r["Sequence type"] == "full genome")],
    ["Partial / fragment sequences",
     sum(1 for r in rows if r["Sequence type"] != "full genome")],
    ["Countries represented", len({r["Country"] for r in rows})],
    ["Plant host species represented",
     len({r["Host (species)"] for r in rows
          if r["Host (species)"] and "insect" not in r["Host family"]})],
    ["Plant host families represented",
     len({r["Host family"] for r in rows
          if r["Host family"] and "insect" not in r["Host family"]})],
    ["Records derived from insects (Cicadellidae)",
     sum(1 for r in rows if "insect" in r["Host family"])],
    ["Collection years spanned (explicit collection dates only)",
     (lambda y: f"{min(y)}-{max(y)}")(
         [int(m["collection_year"]) for m in meta.values()
          if m["year_source"] == "collection_date"])],
    ["Records with an explicit collection date",
     sum(1 for m in meta.values() if m["year_source"] == "collection_date")],
    ["", ""],
    ["Records that are not CpCDV (bean yellow dwarf virus lineage)", len(BEYDV)],
    ["Genomes in the diversity analysis set", len(cfg["keep"])],
    ["", ""],
]
for k, v in Counter(r["Region"] for r in rows).most_common():
    summ.append([f"Records from {k}", v])
summ.append(["", ""])
for k, v in Counter(r["Strain"] for r in rows).most_common():
    summ.append([f"Strain {k}", v])
for r in summ:
    ws3.append(r)
ws3.column_dimensions["A"].width = 58
ws3.column_dimensions["B"].width = 16
ws3["A1"].font = Font(bold=True, size=12)
px = os.path.join(OUT, "Table_S1_CpCDV_sequences.xlsx")
wb.save(px)
print("wrote", px)

# ---------------------------------------------------------------- Word
import docx
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

doc = Document()
sec = doc.sections[0]
sec.orientation = WD_ORIENT.LANDSCAPE
sec.page_width, sec.page_height = sec.page_height, sec.page_width
for m in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
    setattr(sec, m, Cm(1.6))
st = doc.styles["Normal"]
st.font.name = "Times New Roman"; st.font.size = Pt(8)

h = doc.add_paragraph()
r = h.add_run("Supplementary Table S1.")
r.bold = True; r.font.size = Pt(9)
r2 = h.add_run(
    " Chickpea chlorotic dwarf virus (CpCDV; species Mastrevirus cicerparvi) "
    "sequences analysed in this review, with the strain designation, country "
    "and region of origin of each isolate. Strain designations given without "
    "qualification are those published with the sequence; designations marked "
    "\u201cprov.\u201d are provisional assignments inferred from a partial "
    "sequence and should not be relied upon (see Section 6.1). Isolates listed "
    "as \u201cundetermined\u201d could not be assigned because the two closest "
    "strain references differed by less than 1% over the available fragment. "
    "Accessions AM849096 and DQ458791 are retained in the table but belong to "
    "the bean yellow dwarf virus lineage within the same species and were "
    "excluded from the CpCDV diversity analyses.")
r2.font.size = Pt(8)
h.paragraph_format.space_after = Pt(8)

tbl = doc.add_table(rows=1, cols=len(CORE))
tbl.style = "Table Grid"
widths = [Cm(2.3), Cm(5.6), Cm(2.0), Cm(2.9), Cm(4.6), Cm(2.6), Cm(3.4)]
hdr = tbl.rows[0].cells
for i, c in enumerate(CORE):
    hdr[i].text = ""
    p = hdr[i].paragraphs[0]
    run = p.add_run(c); run.bold = True; run.font.size = Pt(8)
    shd = OxmlElement("w:shd"); shd.set(qn("w:fill"), "D9E2F3")
    hdr[i]._tc.get_or_add_tcPr().append(shd)
for row in rows:
    cells = tbl.add_row().cells
    for i, c in enumerate(CORE):
        cells[i].text = ""
        p = cells[i].paragraphs[0]
        p.paragraph_format.space_after = Pt(0)
        run = p.add_run(str(row[c]))
        run.font.size = Pt(7.5)
        if c == "Virus name":
            run.italic = False
for r_ in tbl.rows:
    for i, c in enumerate(r_.cells):
        c.width = widths[i]
# repeat the header row across pages
trPr = tbl.rows[0]._tr.get_or_add_trPr()
el = OxmlElement("w:tblHeader"); el.set(qn("w:val"), "true"); trPr.append(el)

pd = os.path.join(OUT, "Table_S1_CpCDV_sequences.docx")
doc.save(pd)
print("wrote", pd)
print(f"\nTable S1: {len(rows)} records, {len(CORE)} core columns")
print("strain column values:", dict(Counter(r["Strain"] for r in rows).most_common()))
