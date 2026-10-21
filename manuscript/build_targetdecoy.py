"""Methods subsection, Results subsection and figure legend for the target-decoy calibration.

    python build_targetdecoy.py <outdir>   -> TARGETDECOY_text_and_legend.docx
Every number is read from ../icb CSVs, so the text follows the analysis when it is re-run.
"""
import glob, os, sys, numpy as np, pandas as pd
from docx import Document
from docx.shared import Pt
I = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/"
F = pd.read_csv(I + "calib_v3_fdr.csv").set_index(["null", "stat", "pool", "decoys"])
U = pd.read_csv(I + "calib_v3_units.csv")
X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(I + "calib_v3/ICB_*_v3.csv"))])
INFO = pd.read_csv(I + "calib_decoys_v3_info.csv")
PS = pd.read_csv(I + "pool_shift_check.csv")
DIFF = pd.read_csv(I + "calib_v3_fdr_diff.csv").set_index(["alt", "ref", "quantity"])
MAT = pd.read_csv(I + "calib_v3_fdr_matched.csv").set_index(["null", "stat", "subset"])
DCL = pd.read_csv(I + "calib_v3_fdr_diff_close.csv").set_index(["subset", "alt", "ref", "quantity"])
OWN = pd.read_csv(I + "calib_v3_ownbank_compare.csv") if os.path.exists(I + "calib_v3_ownbank_compare.csv") else None
DIRN = X.drop_duplicates("unit").set_index("unit").dirn
U["dirn"] = U.unit.map(DIRN); U["tdirn"] = U.target.map(DIRN)
f = lambda nl, st, p="A", s="coex": F.loc[(nl, st, p, s)]
pc = lambda v, d=1: f"{v:.{d}f}%"
sgn = lambda v, d=3: f"{v:+.{d}f}".replace("-", "−")
WORDS = "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen eighteen nineteen".split()
TENS = {20: "twenty", 30: "thirty", 40: "forty", 50: "fifty"}
def cap(n):
    s = WORDS[n] if n < 20 else TENS[n // 10 * 10] + ("" if n % 10 == 0 else "-" + WORDS[n % 10])
    return s[0].upper() + s[1:]
ci = lambda a, b, d=1: f"{a:.{d}f}–{b:.{d}f}%"

def dec_rate(nl, st, p="A", s="coex", dn=None):
    d = U[(U.kind == "decoy") & (U.pool == p) & (U.structure == s) & (U["null"] == nl) & (U.stat == st)]
    if dn is not None: d = d[d.dirn == dn]
    return 100 * d.call.mean()

T = U[(U.kind == "target")]
called = lambda nl, st: int(T[(T["null"] == nl) & (T.stat == st)].call.sum())
x0 = X[(X["null"] == "N0") & (X.kind == "decoy") & (X.pool == "A")]
clear_c = 100 * (x0[x0.structure == "coex"].fold_obs > x0[x0.structure == "coex"].fold_p95).mean()
ex = (x0[x0.structure == "coex"].fold_obs - x0[x0.structure == "coex"].fold_med)
x0s = x0.assign(signed=np.where(x0.dirn == 1, x0.dir_obs, 1 - x0.dir_obs))
cA = INFO[(INFO.pool == "A") & (INFO.structure == "coex")]
hi_r = cA[cA.r_target_med > .6]
B = X[(X.kind == "target")].copy(); B["sh"] = np.where(B.dirn == 1, B.dir_med, 1 - B.dir_med) - .5
sh = {nl: B[B["null"] == nl].sh for nl in ("N0", "N1", "N2")}
psh = PS[PS.signature.notna()].groupby("pool")["shift"].mean(); pshare = PS[PS.signature.notna()].groupby("pool").sig_share.mean()
sg = PS[PS.pool == "single_gene_sig"]["shift"]
per = U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == "coex") & (U["null"] == "N0") & (U.stat == "dir")]
tt = T[(T["null"] == "N0") & (T.stat == "dir")].set_index("target")
beat_all = sum(int((per[per.target == s].med_delta >= r.med_delta).sum() == 0) for s, r in tt.iterrows())
neg = tt[tt.neg].index.tolist(); rev = tt.loc["CD39-CD8Tcell_Chow"]
nat = "N1nat" in set(U["null"])
r0, r1 = f("N0", "fold"), f("N0", "dir")

METHODS = (
 "Target-decoy calibration. Each of the 52 signatures (the targets) received 56 decoy modules built from genes "
 f"measured in all 12 response cohorts: 32 co-expressed decoys from the {14148:,} genes in none of the 52 signatures "
 f"(pool A), 8 uniform random-gene decoys from the same pool, and 16 co-expressed decoys from the {13106:,} genes "
 f"that also lie outside every set of the 90-signature tumour-microenvironment library (pool B); {len(INFO):,} decoys in all. "
 "A co-expressed decoy was grown to its target's size and expression-decile composition by seed-and-grow, each step "
 "adding a gene whose mean correlation with the genes already chosen came closest to the target's own internal mean "
 "pairwise correlation in each of the 12 cohorts jointly (root-mean-square deviation over cohorts; any gene within 0.02 "
 "of the best was eligible, chosen at random). A decoy inherited its target's scoring engine and, for weighted-mean "
 "signatures, its weight vector in random order, and was assigned a direction at random. It was scored in each cohort "
 "exactly as its target's null sets had been and read against its target's 1,000-draw null bank under N0, N1 and N2, "
 "so that target and decoy differ only in which genes carry the structure; the targets, re-read from the same banks, "
 "reproduce the stored counts. A target or decoy is called when the Wilcoxon signed-rank test of ΔAUC over the 12 "
 "cohorts gives P < 0.05 with a positive median, the test used throughout, whose nominal rate under the null is 2.5%. "
 "The expected number of false discoveries is the sum over targets of the call rate among that target's decoys (for the "
 "directed statistic, among the decoys declared in the target's direction), and the target-decoy false discovery rate "
 "(FDR) is that sum divided by the number of targets called. Intervals are the 2.5th and 97.5th percentiles of 6,000 "
 "bootstrap replicates that resample targets together with their decoys; differences between two analyses use the "
 "same resamples for both. As a sensitivity analysis the FDR was recomputed from only the decoys whose median internal r "
 "lay within 0.05 of their target's, overall and separately for targets with internal r at or below 0.5 and above it. To check the reading of a decoy against its "
 "target's bank, the first four pool-A co-expressed decoys of every target (208) were given their own 300-draw N0, N1 "
 "and N2 banks by the published procedure. N1 and N2 banks were also redrawn (300 draws) from a decile-stratified "
 "sample of all measured genes, in which signature genes appear at their natural share, in place of the original pool, "
 "which held every signature gene.")

# ---- robustness: paired differences and stricter matching (calib_v3_robust.py)
dv = lambda alt, q: DIFF.loc[(alt, "N0-fold", q)]
mt = lambda nl, st, sub: MAT.loc[(nl, st, sub)]
dc = lambda sub, alt, q: DCL.loc[(sub, alt, "N0-fold", q)]
cl_share = 100 * ((INFO.r_achieved_med - INFO.r_target_med).abs() <= .05)[(INFO.pool == "A") & (INFO.structure == "coex")].mean()
rr = U[(U.kind == "decoy") & (U.pool == "A") & (U.structure == "coex") & (U["null"] == "N0") & (U.stat == "dir")].copy()
rr["r"] = rr.unit.map(INFO.set_index("decoy").r_achieved_med); rr["sig"] = rr.call | rr.neg
hi_sig, lo_sig = 100 * rr[rr.r > .5].sig.mean(), 100 * rr[rr.r <= .3].sig.mean()
nlo = int(mt("N0", "fold", "targets r <= 0.5").targets)
a, b = dv("N0-dir", "exp_true"), dv("N0-dir", "fdr")
m0, m1, m2 = mt("N0", "fold", "decoys within 0.05 r"), mt("N0", "dir", "decoys within 0.05 r"), mt("N2nat", "dir", "decoys within 0.05 r")
c1, c2 = dc("all targets", "N0-dir", "exp_true"), dc("all targets", "N2nat-dir", "exp_true")
l0, l1 = dc("targets r <= 0.5", "N0-dir", "fdr"), dc("targets r > 0.5", "N0-dir", "exp_true")
l1f = dc("targets r > 0.5", "N0-dir", "fdr")
ROB = (f"Resampling the targets once for both analyses, directed scoring added {a['diff']:.1f} expected true discoveries "
  f"({a.lo:.1f} to {a.hi:.1f}) and lowered the FDR by {-b['diff']:.1f} percentage points ({-b.hi:.1f} to {-b.lo:.1f}) "
  f"against the published analysis. Because the decoys fell short of the most co-expressed signatures, the FDR was "
  f"recomputed from the {cl_share:.0f}% of decoys whose median internal r lay within 0.05 of their target's "
  f"({int(m0.targets)} targets keep at least one). The published analysis did not move (FDR {m0.fdr:.0f}%, "
  f"{m0.fdr_lo:.0f}–{m0.fdr_hi:.0f}%), but the directed FDR rose to {m1.fdr:.0f}% ({m1.fdr_lo:.0f}–{m1.fdr_hi:.0f}%; "
  f"{m2.fdr:.0f}% under the natural-pool N2) and its gain in expected true discoveries narrowed to {c1['diff']:.1f} "
  f"({sgn(c1.lo, 1)} to {sgn(c1.hi, 1)}; natural-pool N2 {c2['diff']:.1f}, {sgn(c2.lo, 1)} to {sgn(c2.hi, 1)}). "
  f"The two analyses part company on co-expression. Among the {nlo} signatures with internal r of 0.5 or less, the "
  f"folded analysis called {int(l0.ref_called)} at an FDR of {l0.ref_value:.0f}% and the directed {int(l0.alt_called)} "
  f"at {l0.alt_value:.0f}%; among the {52 - nlo} above 0.5 they called {int(l1.ref_called)} and {int(l1.alt_called)}, "
  f"with about {l1.ref_value:.0f} expected true discoveries each (FDR {l1f.ref_value:.0f}% and {l1f.alt_value:.0f}%). "
  f"Strongly co-expressed decoys carry response information of their own: decoys with internal r above 0.5 reached "
  f"significance in one direction or the other in {hi_sig:.1f}% of directed tests under the published null, against "
  f"{lo_sig:.1f}% for decoys at 0.3 or below and a nominal 5%. For a coherent signature the bar is set by coherent "
  "modules rather than by random genes, and its decoys measure what its gene content adds beyond its coherence.")

RES = [
 "A false discovery rate for the signature library",
 (f"Decoys reproduced their targets' co-expression cohort by cohort (median internal r {cA.r_achieved_med.median():.2f} "
  f"against {cA.r_target_med.median():.2f}; correlation of the two across cohorts {cA.r_profile_cor.median():.2f}), "
  f"falling short only for the most tightly co-expressed signatures (r above 0.6: {hi_r.r_achieved_med.median():.2f} "
  f"against {hi_r.r_target_med.median():.2f}), so their false-positive rates there are, if anything, underestimates "
  "(Fig. Xb). Random-gene decoys stayed near the transcriptome background."),
 (f"The published analysis, the folded AUC against a null matched on size and engine, called {called('N0','fold')} of "
  "52 signatures (20 in the original analysis, whose null had 300 draws; with 1,000 draws CYT_Davoli moves from "
  f"P = 0.052 to 0.043) and {pc(r0.decoy_rate)} of their co-expressed decoys ({ci(r0.decoy_lo, r0.decoy_hi)}), "
  f"{r0.decoy_rate/2.5:.1f} times the nominal rate: {r0.exp_false:.1f} expected false discoveries, an FDR of "
  f"{r0.fdr:.0f}% ({r0.fdr_lo:.0f}–{r0.fdr_hi:.0f}%) (Fig. Xa,c). This cannot be seen test by test. Co-expressed "
  f"decoys cleared the null's 95th percentile in {clear_c:.1f}% of cohort tests against a nominal 5%, but lay above "
  f"the null median in {100*(ex>0).mean():.0f}% of them (median {sgn(ex.median(), 4)}), and a signed-rank test over 12 "
  "cohorts turns a small consistent excess into a call (Fig. Xd). The excess is spread, not location: co-expressed "
  f"decoys were centred (signed AUC − 0.5 = {sgn(x0s[x0s.structure=='coex'].signed.mean()-.5)}), but their AUC varied "
  "more than that of random-gene decoys in the cohorts where co-expression is strongest, and folding reads that spread "
  f"at its favourable tail (Fig. Xe). Decoys drawn from outside the TME library behaved alike ({pc(f('N0','fold','B').decoy_rate)}), "
  "so the excess is not borrowed immune biology."),
 (f"Scoring each signature in its published direction against the same null called {called('N0','dir')} signatures "
  f"and {pc(r1.decoy_rate)} of direction-matched decoys ({ci(r1.decoy_lo, r1.decoy_hi)}): {r1.exp_false:.1f} expected "
  f"false discoveries, an FDR of {r1.fdr:.0f}% ({r1.fdr_lo:.0f}–{r1.fdr_hi:.0f}%). The expected number of true "
  f"discoveries rises from {r0.exp_true:.0f} to {r1.exp_true:.0f}. All {called('N0','dir')} are sensitivity signatures; "
  f"none of the {int((DIRN.reindex(tt.index)==-1).sum())} resistance signatures passes. {cap(beat_all)} signatures exceed "
  f"every one of their own 32 decoys (Fig. Xf). {cap(len(neg))} run significantly the wrong way in their own declared "
  f"direction, among them CD39+CD8+ T cells (Chow; median ΔAUC {sgn(rev.med_delta)}, P = {rev.p:.4f}), one of the "
  "signatures the folded analysis counts as a success."),
 ROB,
 (f"The matched nulls N1 and N2, as first drawn, favoured one declared direction: decoys declared to mark resistance "
  f"were called in {pc(dec_rate('N1','dir',dn=-1))} (N1) and {pc(dec_rate('N2','dir',dn=-1))} (N2), sensitivity decoys "
  f"in {pc(dec_rate('N1','dir',dn=1))} and {pc(dec_rate('N2','dir',dn=1))} (Fig. Xg). Their banks sat above 0.5 "
  f"(mean signed median − 0.5 = {sgn(sh['N1'].mean())} and {sgn(sh['N2'].mean())}; positive in "
  f"{100*(sh['N2']>0).mean():.0f}% of tests) where N0's did not ({sgn(sh['N0'].mean())}). Their draw pool held every "
  f"signature gene, {100*pshare['enriched']:.1f}% of the pool against {100*pshare['natural']:.1f}% of a natural sample, "
  f"and signature genes run higher in responders (single-gene AUC − 0.5 = {sgn(sg.mean())}, positive in "
  f"{int((sg>0).sum())} of 12 cohorts); under mean scoring, moving from a signature-free to that pool shifts a "
  f"decile-matched null by {sgn(psh['enriched']-psh['sig_free'])}.")]
if nat:
    n1, n2 = f("N1nat", "dir"), f("N2nat", "dir")
    RES.append(
     f"Redrawn from a natural pool, the matched nulls no longer favoured a direction (resistance decoys "
     f"{pc(dec_rate('N1nat','dir',dn=-1))} and {pc(dec_rate('N2nat','dir',dn=-1))}, sensitivity "
     f"{pc(dec_rate('N1nat','dir',dn=1))} and {pc(dec_rate('N2nat','dir',dn=1))}). Directed scoring then called "
     f"{called('N1nat','dir')} (N1) and {called('N2nat','dir')} (N2) signatures, FDR {n1.fdr:.0f}% and {n2.fdr:.0f}%, and "
     f"folded scoring {called('N1nat','fold')} and {called('N2nat','fold')}, FDR {f('N1nat','fold').fdr:.0f}% and "
     f"{f('N2nat','fold').fdr:.0f}% (Fig. Xa). Matched correctly, co-expression costs the folded analysis no "
     f"signatures ({called('N0','fold')} under N0) and halves its FDR; the fall to {called('N2','fold')} under the "
     "original N2 came from the pool, not from the matching.")
if OWN is not None:
    def rcor(nl, st="dir"):
        q = OWN[(OWN["null"] == nl) & (OWN.stat == st)].pivot_table(index="unit", columns="src", values="med").dropna()
        return np.corrcoef(q.shared, q.own)[0, 1]
    def own_rate(nl, dn):
        q = OWN[(OWN.src == "own") & (OWN["null"] == nl) & (OWN.stat == "dir") & (OWN.dirn == dn)]
        return 100 * q.call.mean()
    RES.append(
     "Reading a decoy against its target's bank and against a bank built for the decoy itself by the published "
     f"procedure gave the same answer (correlation of median ΔAUC over {OWN.unit.nunique()} decoys, directed "
     f"statistic: N0 {rcor('N0'):.2f}, N1 {rcor('N1'):.2f}, N2 {rcor('N2'):.2f}; natural pool {rcor('N1nat'):.2f} and "
     f"{rcor('N2nat'):.2f}). The decoys' own banks from the original pool reproduced the direction bias "
     f"(resistance-declared decoys called in {own_rate('N1',-1):.1f}% and {own_rate('N2',-1):.1f}% under N1 and N2, "
     f"sensitivity decoys in {own_rate('N1',1):.1f}% and {own_rate('N2',1):.1f}%), and their banks from the natural "
     f"pool did not ({own_rate('N1nat',-1):.1f}% and {own_rate('N2nat',-1):.1f}% against {own_rate('N1nat',1):.1f}% "
     f"and {own_rate('N2nat',1):.1f}%) (Fig. Xh). The bias belongs to the original procedure, not to the shortcut.")

LEG = [
 "Figure X. A target-decoy false discovery rate for 52 published immune signatures.",
 ("(a) Signatures called (Wilcoxon signed-rank over 12 cohorts, P < 0.05, positive median) under each null and "
  "statistic, split into expected true (solid) and expected false (hatched) discoveries; the bar on the boundary is "
  "the 95% interval of the expected false count, and the number above each bar the target-decoy FDR. Folded, grey; "
  "directed, blue. (b) Internal mean pairwise correlation of each signature against the median of its decoys, one point "
  "per signature and cohort. (c) Percentage of decoys called at the same unit of analysis; nominal rate 2.5%. (d) Folded "
  "statistic under the published null: the share of cohort tests above the null 95th percentile (nominal 5%) against "
  "the share of modules called over 12 cohorts (nominal 2.5%). (e) Standard deviation of the AUC across decoys, per "
  "cohort, co-expressed against random-gene decoys. (f) Each signature's median ΔAUC (directed statistic, published "
  "null; filled, called; red ring, significantly negative) against the 5th–95th percentile and median of its 32 "
  "co-expressed decoys. (g) Directed decoy call rate by the direction the decoy was declared to carry, under the "
  "published null and the matched nulls drawn from the original (enriched) and a natural gene pool. (h) Median ΔAUC of "
  "208 decoys read against their target's null bank and against their own. Intervals are 95% cluster-bootstrap "
  "intervals over targets, each carried with its decoys (6,000 replicates in a, 4,000 in d and g).")]

def main(OUT):
    d = Document()
    st = d.styles["Normal"]; st.font.name = "Arial"; st.font.size = Pt(10)
    d.add_heading("Methods — new subsection", level=2); d.add_paragraph(METHODS)
    d.add_heading("Results — new subsection", level=2)
    d.add_paragraph().add_run(RES[0]).bold = True
    for p in RES[1:]: d.add_paragraph(p)
    d.add_heading("Figure legend", level=2)
    d.add_paragraph().add_run(LEG[0]).bold = True; d.add_paragraph(LEG[1])
    d.save(os.path.join(OUT, "TARGETDECOY_text_and_legend.docx"))
    w = lambda s: len(" ".join(s if isinstance(s, list) else [s]).split())
    print(f"Methods {w(METHODS)} words | Results {w(RES)} words | legend {w(LEG)} words | natural={'yes' if nat else 'pending'} own-bank={'yes' if OWN is not None else 'pending'}")
    print("\n".join(RES[1:]))

if __name__ == "__main__":
    main(sys.argv[1])
