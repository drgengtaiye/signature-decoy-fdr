"""Assemble the Genome Biology submission package under the new numbering.

usage: build_package.py <submission_GM dir> <revision_GB dir> <out dir>

Old figures 1, 3, 4-6 become 1, 3, 5-7; Figure 2 is the redrawn fig2_gb.py; Fig_targetdecoy is
Figure 4; old 7-8 are dropped. Additional files 1-3 are the rebuilt versions in deliverables/
(build_addfile1/2/3.py), 4 the decoy workbook. STATUS.md lists what each file is and what is open.
"""
import shutil
import sys
import zipfile
from pathlib import Path

import fitz

GM, GB, OUT = (Path(a) for a in sys.argv[1:4])
DL = GB / "deliverables"
FIGS = GB.parent / "figs_pub"

# (published name, source, status note; empty note = ready)
FILES = [
    ("Cover_letter.docx", DL / "COVER_LETTER_GenomeBiology_v2.docx", ""),
    ("Manuscript.docx", DL / "MANUSCRIPT_GB_v3.docx",
     "title page gives e-mails for the three corresponding authors only; enter every author's e-mail in the submission system"),
    ("Figure1.pdf", GM / "Figure1.pdf", ""),
    ("Figure2.pdf", FIGS / "Figure2_GB.pdf", ""),
    ("Figure3.pdf", GM / "Figure3.pdf", ""),
    ("Figure4.pdf", DL / "Fig_targetdecoy.pdf", ""),
    ("Figure5.pdf", GM / "Figure4.pdf", ""),
    ("Figure6.pdf", GM / "Figure5.pdf", ""),
    ("Figure7.pdf", GM / "Figure6.pdf", ""),
    ("Additional_file_1_Notes_and_Figures.docx", DL / "AF1" / "Additional_file_1_Notes_and_Figures.docx", ""),
    ("Additional_file_2_derived_tables.xlsx", DL / "Additional_file_2_derived_tables.xlsx", ""),
    ("Additional_file_3_analysis_code.zip", DL / "Additional_file_3_analysis_code.zip", ""),
    ("Additional_file_4_target_decoy.xlsx", DL / "Additional_file_4_target_decoy.xlsx", ""),
]
PDFDIR = DL / "AF1" / "Additional_file_1_figures_pdf"

OUT.mkdir(parents=True, exist_ok=False)
for name, src, _ in FILES:
    shutil.copy2(src, OUT / name)
shutil.copytree(PDFDIR, OUT / PDFDIR.name)

# previews of the seven main and six supplementary figures, one page each at 110 dpi
prev = OUT / "_preview"
prev.mkdir()
pages = [(OUT / f"Figure{i}.pdf", f"Figure{i}.png") for i in range(1, 8)]
pages += [(p, p.stem.replace("Supplementary_Figure_", "Fig") + ".png") for p in sorted((OUT / PDFDIR.name).glob("*.pdf"))]
for src, png in pages:
    with fitz.open(src) as doc:
        doc[0].get_pixmap(dpi=110).save(prev / png)

lines = ["# Genome Biology submission package", "",
         "| File | Source | Status |", "|---|---|---|"]
for name, src, note in FILES:
    lines.append(f"| `{name}` | `{src.relative_to(src.parents[1])}` | {note or 'ready'} |")
lines.append(f"| `{PDFDIR.name}/` | `deliverables/AF1/{PDFDIR.name}` | ready (S1-S6) |")
(OUT / "STATUS.md").write_text("\n".join(lines) + "\n")

with zipfile.ZipFile(OUT.with_suffix(".zip"), "w", zipfile.ZIP_DEFLATED) as z:
    for name, _, _ in FILES:
        z.write(OUT / name, f"{OUT.name}/{name}")
    for p in sorted((OUT / PDFDIR.name).glob("*.pdf")):
        z.write(p, f"{OUT.name}/{PDFDIR.name}/{p.name}")
    z.write(OUT / "STATUS.md", f"{OUT.name}/STATUS.md")
print(OUT.with_suffix(".zip"), OUT.with_suffix(".zip").stat().st_size)
