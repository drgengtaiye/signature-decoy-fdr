"""Code availability for the Genome Biology manuscript.

usage: build_manuscript_GB_v3.py <MANUSCRIPT_GB_v2.docx> <out.docx>

The analysis code is released as a public repository holding exactly Additional file 3; the
earlier immunoTriage-paper repository (which also holds the separate paper's analyses) and the
immunoTriage R package (which this paper's code does not call) are no longer cited. Reference 99
was cited only here and is the last reference, so no renumbering is needed.
"""
import sys

import docx

SRC, OUT = sys.argv[1:3]
REPO = "https://github.com/drgengtaiye/signature-decoy-fdr"
ACCESSED = "20 Oct 2026"
d = docx.Document(SRC)
ps = d.paragraphs
find = lambda start: [p for p in ps if p.text.startswith(start)]


def set_text(p, new):
    """p has one run; replace its text and keep the run's format."""
    assert len(p.runs) == 1, p.text[:60]
    p.runs[0].text = new


(avail,) = find("All analysis code is provided as Additional file 3")
set_text(avail, "All analysis code is provided as Additional file 3 and in the signature-decoy-fdr repository [98], "
         "including the scripts that score each signature by its publication’s own method against its size- and "
         "engine-matched null (null_matched_v2.R, null_matched_prog.R) and those that grow and score the decoy "
         "modules and estimate the false discovery rate (calib_decoys_v3.R, calib_score_v3.R, null_natpool_v4.R, "
         "calib_v3_fdr.py, calib_v3_robust.py). The repository will be archived with a DOI on acceptance.")
(r98,) = find("98. immunoTriage-paper")
set_text(r98, f"98. signature-decoy-fdr: analysis code for this study. {REPO}. Accessed {ACCESSED}.")
(r99,) = find("99. immunoTriage")
r99._p.getparent().remove(r99._p)

d.save(OUT)
d2 = docx.Document(OUT)
txt = "\n".join(p.text for p in d2.paragraphs)
assert "[99]" not in txt and "dryichen" not in txt and "immunoTriage" not in txt, "residual"
print("saved", OUT)
