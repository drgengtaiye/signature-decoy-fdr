#!/usr/bin/env python3
from docx import Document
from docx.shared import Pt

OUT = "/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-self-shanghai-6gastirc/30d20b73-d629-4090-b43c-c8ca139220a1/scratchpad/GM_revision"

PARA = (
 "How a gene set should be compared with chance has been studied for two decades, and the "
 "answer turns on what the comparison is against: a self-contained test asks whether a set’s "
 "genes are associated with the outcome, a competitive test whether they are more associated "
 "than other genes, and only the second speaks to whether a published signature beats an "
 "arbitrary gene list [9, 10]. Competitive statistics are inflated by correlation between "
 "genes, which do not contribute independent evidence [11]; the established remedies estimate "
 "that inflation and divide it out, as camera does with a variance-inflation factor [12], or "
 "resample around it by rotation [13], structured permutation [14] or restandardization [15]. "
 "None applies directly here. Those methods test a per-gene differential-expression statistic "
 "inside a linear model, where the inflation has a closed form. A published signature is "
 "instead a fixed gene list carrying its own weights, reduced to one score per patient by its "
 "own engine, and what is asked of it is how well that score discriminates a clinical "
 "endpoint. No closed-form null exists for that quantity, so it has to be drawn, from "
 "comparison sets matched to the signature on every property that moves the statistic. Gene "
 "count and scoring engine are two such properties; co-expression is a third, and the one a "
 "size-matched null omits — in these cohorts the genes of a published signature have a mean "
 "pairwise correlation of 0.41 against 0.09 for random sets of the same size, and are the more "
 "correlated of the two in every cohort. What a signature must be compared with is therefore a "
 "measurement question, and it has to be settled before the literature can be audited."
)

REFS = [
 ("9",  "Goeman JJ, Bühlmann P. Analyzing gene expression data in terms of gene sets: "
        "methodological issues. Bioinformatics. 2007;23:980–987.", "PMID 17303618"),
 ("10", "Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, "
        "Pomeroy SL, Golub TR, Lander ES, Mesirov JP. Gene set enrichment analysis: a "
        "knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl "
        "Acad Sci U S A. 2005;102:15545–15550.", "PMID 16199517"),
 ("11", "Tamayo P, Steinhardt G, Liberzon A, Mesirov JP. The limitations of simple gene set "
        "enrichment analysis assuming gene independence. Stat Methods Med Res. "
        "2016;25:472–487.", "PMID 23070592"),
 ("12", "Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene "
        "correlation. Nucleic Acids Res. 2012;40:e133.", "PMID 22638577"),
 ("13", "Wu D, Lim E, Vaillant F, Asselin-Labat ML, Visvader JE, Smyth GK. ROAST: rotation gene "
        "set tests for complex microarray experiments. Bioinformatics. 2010;26:2176–2182.",
        "PMID 20610611"),
 ("14", "Barry WT, Nobel AB, Wright FA. Significance analysis of functional categories in gene "
        "expression studies: a structured permutation approach. Bioinformatics. "
        "2005;21:1943–1949.", "PMID 15647293"),
 ("15", "Efron B, Tibshirani R. On testing the significance of sets of genes. Ann Appl Stat. "
        "2007;1:107–129.", "doi:10.1214/07-AOAS101 (not indexed in PubMed)"),
]

d = Document()
st = d.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(11)

def head(t, sz=13):
    p = d.add_paragraph(); r = p.add_run(t); r.bold = True; r.font.size = Pt(sz)
    p.paragraph_format.space_before = Pt(10); return p

def note(t):
    p = d.add_paragraph(); r = p.add_run(t); r.italic = True; r.font.size = Pt(9.5); return p

head("Background: new positioning paragraph", 14)
note("Placement: insert as the third paragraph of Background, directly after the paragraph "
     "ending “…has not been tested systematically” and before the paragraph "
     "beginning “There is also a conceptual issue…”. It answers the question a "
     "methods reviewer asks first — how this differs from the existing correlation-aware "
     "gene-set tests — which the current Background does not address. 247 words.")
p = d.add_paragraph(PARA); p.paragraph_format.space_before = Pt(8)
p.paragraph_format.line_spacing = 2.0

head("References to add")
note("Citation numbers above assume these enter as [9]–[15]; every later reference shifts "
     "by seven. Because the author lists of 65 of the existing 96 entries belong to other "
     "papers (see REFERENCES_corrected.docx), the bibliography should be rebuilt in a "
     "reference manager from the identifiers rather than renumbered by hand. Identifiers are "
     "given for that purpose. All seven were verified against PubMed or Crossref.")
for n, t, pid in REFS:
    p = d.add_paragraph(f"{n}. {t}")
    p.paragraph_format.space_after = Pt(3)
    r = p.add_run(f"  [{pid}]"); r.font.size = Pt(9); r.italic = True

head("What the paragraph commits the paper to")
for t in [
 "The paper is a measurement study: it defines the comparison, then applies it. That is the "
 "methods-journal identity, and it is stated in the Background rather than left implicit.",
 "The co-expression gap is given as a number in the Background (0.41 against 0.09), so the "
 "new matched null is motivated before it is described.",
 "camera is named as the closest prior art and the difference is stated precisely — an "
 "arbitrary scoring engine and a discrimination endpoint admit no closed-form inflation "
 "factor — which converts the missing-citation problem into the paper’s rationale.",
]:
    p = d.add_paragraph(t, style="List Bullet"); p.paragraph_format.space_after = Pt(4)

d.save(f"{OUT}/REVISED_background_paragraph.docx")
print(f"saved {OUT}/REVISED_background_paragraph.docx")
print(f"paragraph words: {len(PARA.split())}")
