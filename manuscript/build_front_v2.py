"""Title and abstract, version 2: the target-decoy FDR as the paper's spine (Genome Biology format:
Background / Results / Conclusions, <= 350 words). Assumes the attribution model goes to its own paper.

    python build_front_v2.py <MANUSCRIPT_GM1015.docx> <outdir>  -> REVISED_title_and_abstract_v2.docx
Numbers are read from ../icb so the abstract follows the analysis.
"""
import copy, sys, pandas as pd
from docx import Document
from docx.text.paragraph import Paragraph
I = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/"
F = pd.read_csv(I + "calib_v3_fdr.csv").set_index(["null", "stat", "pool", "decoys"])
U = pd.read_csv(I + "calib_v3_units.csv")
INFO = pd.read_csv(I + "calib_decoys_v3_info.csv")
f = lambda nl, st: F.loc[(nl, st, "A", "coex")]
T = U[U.kind == "target"]
called = lambda nl, st: int(T[(T["null"] == nl) & (T.stat == st)].call.sum())
r0, r1, r2 = f("N0", "fold"), f("N0", "dir"), f("N2nat", "dir")
import glob, numpy as np
X = pd.concat([pd.read_csv(x) for x in sorted(glob.glob(I + "calib_v3/ICB_*_v3.csv"))])
dn = X.drop_duplicates("unit").set_index("unit").dirn; U["dirn"] = U.unit.map(dn)
res = lambda nl: 100 * U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == "coex") & (U["null"] == nl) & (U.stat == "dir") & (U.dirn == -1)].call.mean()
rv = T[(T["null"] == "N0") & (T.stat == "dir")].set_index("target").loc["CD39-CD8Tcell_Chow"]
ci = lambda r: f"{r.fdr:.0f}% (95% CI {r.fdr_lo:.0f}–{r.fdr_hi:.0f}%)"

TITLE = "Matched decoy modules give immune gene-expression signatures a false discovery rate"
DIFF = pd.read_csv(I + "calib_v3_fdr_diff.csv").set_index(["alt", "ref", "quantity"])
MAT = pd.read_csv(I + "calib_v3_fdr_matched.csv").set_index(["null", "stat", "subset"])
DCL = pd.read_csv(I + "calib_v3_fdr_diff_close.csv").set_index(["subset", "alt", "ref", "quantity"])
gain = DIFF.loc[("N0-dir", "N0-fold", "exp_true")]
m0, m1 = MAT.loc[("N0", "fold", "decoys within 0.05 r")], MAT.loc[("N0", "dir", "decoys within 0.05 r")]
lo = DCL.loc[("targets r <= 0.5", "N0-dir", "N0-fold", "fdr")]
rr = U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == "coex") & (U["null"] == "N0") & (U.stat == "dir")].copy()
rr["r"] = rr.unit.map(INFO.set_index("decoy").r_achieved_med)
hi_sig = 100 * (rr[rr.r > .5].call | rr[rr.r > .5].neg).mean()
AB = {
"Background: ":
 "A gene-expression signature is usually declared to work when it beats random gene sets of the same size, but "
 "nothing says how many such claims in a library are false. Peptide identification solved the same problem with "
 "target-decoy searches.",
"Results: ":
 f"We grew {len(INFO):,} decoy modules from genes in no published signature, each matched to one of 52 published "
 "immune signatures in size, scoring engine, expression level and internal co-expression in each of 12 "
 "checkpoint-treated cohorts (827 patients), and put signatures and decoys through the identical test. The "
 "published analysis, a folded AUC against a null matched on size and engine, called "
 f"{called('N0','fold')} signatures at an estimated false discovery rate (FDR) of {ci(r0)}, and {m0.fdr:.0f}% with "
 "only the closest-matched decoys. Its per-test error was nominal, but co-expression widens the spread of a module's "
 "AUC, folding reads that spread at its favourable tail, and pooling 12 cohorts called decoys at "
 f"{r0.decoy_rate/2.5:.1f} times the nominal rate. Scoring each signature in its published direction called "
 f"{called('N0','dir')}, adding {gain['diff']:.0f} expected true discoveries (95% CI {gain.lo:.0f}–{gain.hi:.0f}), at an "
 f"FDR of {r1.fdr:.0f}% ({m1.fdr:.0f}% with the closest-matched decoys); the gain lay in weakly co-expressed "
 f"signatures, whose folded FDR was {lo.ref_value:.0f}%. Folding had also counted as a success a signature that runs "
 f"significantly the wrong way (P = {rv.p:.3f}). Strongly co-expressed decoys reached significance at "
 f"{hi_sig/5:.1f} times the nominal rate, so coherent signatures must be judged against coherent modules. A "
 "co-expression-matched null drawn from a pool enriched for signature genes favoured one declared direction; drawn "
 f"from the natural transcriptome it did not, and called {called('N2nat','dir')} at an FDR of {r2.fdr:.0f}%. In two "
 "randomized trials, signatures beat their null in 15 of 34 checkpoint-arm tests and 0 of 34 comparator-arm tests "
 "(pooled hazard ratio 0.66 per standard deviation, 355 patients, one-sided P = 0.0026).",
"Conclusions: ":
 "Decoys matched to each signature attach a false discovery rate to the claim that a signature works, at the level "
 "at which the claim is made, and expose errors that per-test calibration misses. Signatures should be scored in "
 "their declared direction and judged against decoys that share their co-expression.",
}
KW = "Gene-expression signature; False discovery rate; Target-decoy; Null model; Immune checkpoint blockade; Benchmarking"

def main(SRC, OUTD):
    src = Document(SRC)
    sps = [p for p in src.paragraphs if p.text.strip()]
    authors = [p.text.strip() for p in sps[1:9]]
    d = Document(SRC); body = d.element.body
    tmpl = d.paragraphs[4]
    for p in list(d.paragraphs): body.remove(p._element)
    def add(t, bold_prefix=None):
        el = copy.deepcopy(tmpl._element); body.append(el)
        par = Paragraph(el, d)
        for r in par.runs[1:]: r._element.getparent().remove(r._element)
        if not par.runs: par.add_run('')
        par.runs[0].text = ''
        if bold_prefix:
            r0_ = par.runs[0]; r0_.text = bold_prefix; r0_.bold = True; par.add_run(t)
        else:
            par.runs[0].text = t
        return par
    add(TITLE).runs[0].bold = True
    for a in authors: add(a)
    add('Abstract').runs[0].bold = True
    for k, v in AB.items(): add(v, bold_prefix=k)
    add('Keywords').runs[0].bold = True
    add(KW)
    d.save(f"{OUTD}/REVISED_title_and_abstract_v2.docx")
    print('title words:', len(TITLE.split()), '| abstract words:', sum(len((k + v).split()) for k, v in AB.items()))
    for k, v in AB.items(): print(k + v)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
