"""Figure 8 (v18) in the look of pixl_nc.py: published prognostic risk models over the whole search (round 15 part A),
clinical multigene assays as positive controls (round 15 part D), and the published pipeline on theme pools (round 13
part B, v17 Fig. 9 n-r). Plans: validation/PLAN_round13_published_risk_models.md (amendments 1-2) and
validation/PLAN_round15_literature_cohorts_controls.md (amendments 1-4). Tables: icb/v20_round15 (decoys/v20_audit.py,
v20_citations.py, v20_assays_audit.py; applications/v20_published_models.R, v20_assays.R) and icb/v19_audit (round 13).
Shares, medians and quantiles read from those tables are descriptive; no test, interval or cut-off is made here."""
import glob
import json
import os
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd
from matplotlib.patches import Polygon

A13 = D + "/icb/v19_audit/"
A15 = D + "/icb/v20_round15/"
V = D + "/signature-decoy-fdr/validation/"
THEME = {"emt": "EMT", "m6a": "m6A"}
TYPES = ["BLCA", "BRCA", "CESC", "COAD", "HNSC", "LUAD", "LUSC", "OV", "SARC", "SKCM", "STAD"]

# ---- round 15 part A
U = pd.read_csv(A15 + "audit_models_units.csv")
W = U[U.n_coex.notna()]                                   # models with decoys
IND = W[~W.uses_external.astype(bool)]
USED = W[W.uses_external.astype(bool)]
R = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A15 + "models_[A-Z]*.csv.gz"))])
DEC = R[R.kind == "decoy"]
AM = pd.read_csv(A15 + "audit_models.csv")
AM13 = pd.read_csv(A13 + "audit_models.csv")
BYT = pd.read_csv(A15 + "audit_models_by_type_independent.csv").set_index("type")
BYF = pd.read_csv(A15 + "audit_models_by_family_independent.csv")
BYY = pd.read_csv(A15 + "audit_models_by_years_independent.csv").set_index("years")
BYS = pd.read_csv(A15 + "audit_models_by_size_independent.csv").set_index("size")
CIT = pd.read_csv(A15 + "citations.csv")
CITS = pd.read_csv(A15 + "citations_summary.csv").set_index("population")
PM = pd.read_csv(V + "published_models_r15.csv")
SEARCH = json.load(open(V + "round13_search/search_summary.json"))
SCR13 = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(V + "round13_search/screen_*.csv"))])
SCR15 = pd.read_csv(V + "round15_search/screen_r15_final.csv")


def am(T, pop, q, dec="coex"):
    return T[T.population.str.startswith(pop) & (T.quantity == q) & (T.decoys == dec)].iloc[0]


# the brief's reference numbers: stop if a table disagrees
_chk = [(am(AM, "primary", "fdr").estimate, 0.395), (am(AM, "primary", "fdr").lo, 0.357), (am(AM, "primary", "fdr").hi, 0.442),
        (am(AM, "primary", "fdr", "rand").estimate, 0.347), (am(AM, "primary", "fdr", "rand").lo, 0.308),
        (am(AM, "primary", "fdr", "rand").hi, 0.394), (am(AM, "all models", "fdr").estimate, 0.318), (am(AM, "all models", "fdr").lo, 0.293),
        (am(AM, "all models", "fdr").hi, 0.347), (am(AM, "models whose", "fdr").estimate, 0.209),
        (am(AM, "primary", "decoy_share").estimate, 0.167), (am(AM, "primary", "rho").estimate, 0.726)]
assert all(abs(a - b) < 0.0006 for a, b in _chk), _chk
assert len(IND) == 427 and IND.validated.sum() == 180

# ---- round 15 part D (clinical assays); drawn only once v20_assays_audit.py has written its tables
AS = A15
HAVE_ASSAY = all(os.path.exists(AS + f) for f in ("assay_units_all.csv", "assay_h7.csv", "assay_h7_sensitivity.csv", "assay_by_assay.csv"))
if HAVE_ASSAY:
    AU = pd.read_csv(AS + "assay_units_all.csv")
    AH = pd.concat([pd.read_csv(AS + "assay_h7.csv").assign(set="primary"), pd.read_csv(AS + "assay_h7_sensitivity.csv")]).set_index("set")
    ABA = pd.read_csv(AS + "assay_by_assay.csv").set_index("assay")
ASSAYS = [("RS21_breast", "21-gene RS"), ("MammaPrint70", "70-gene"), ("EndoPredict", "EndoPredict"), ("PAM50_ROR", "PAM50 ROR"),
          ("BCI", "BCI")]
COLON = ("RS12_colon", "12-gene colon RS")
BR_PRIM = ["TCGA", "BRCA_GSE20685", "BRCA_GSE42568", "BRCA_GSE58812", "BRCA_GSE21653", "BRCA_GSE20711", "BRCA_GSE96058", "BRCA_METABRIC"]
BR_ER = ["BRCA_GSE21653_ERpos", "BRCA_GSE96058_ERpos", "BRCA_METABRIC_ERpos"]
CO_PRIM = ["TCGA", "COAD_GSE39582", "COAD_GSE17538", "COAD_GSE38832", "COAD_GSE103479", "COAD_GSE72970"]
CO_OVL = ["COAD_GSE17536", "COAD_GSE14333"]

# ---- round 13 part B (pipeline), as fig9_pixl.py reads it
PU = pd.read_csv(A13 + "audit_pipeline_units.csv")
PA = pd.read_csv(A13 + "audit_pipeline.csv")
PD = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A13 + "pipeline_decoys_*.csv.gz"))])


def colhead(x, y, s, ha="center"):
    pg.txt(x, y, s, ha=ha, weight="bold")


def ymm(top, h, lo, hi, v, invert=False):
    """page y (mm) of data value v on an axis placed at top..top+h with limits lo..hi"""
    f = (v - lo) / (hi - lo)
    return top + h * (f if invert else 1 - f)


def sina(ax, x, vals, c, half=0.42, bw=0.05, step=0.045, ms=1.2, filled=True):
    """points of one group spread sideways in proportion to their local density, with a black median line"""
    vals = np.asarray(vals, float); vals = vals[np.isfinite(vals)]
    off = np.zeros(len(vals)); b = np.floor(vals / bw).astype(int)
    for k in np.unique(b):
        ii = np.where(b == k)[0]; n = len(ii)
        w = min(half, step * (n - 1) / 2)
        off[ii[np.argsort(vals[ii])]] = np.linspace(-w, w, n)[np.argsort(np.arange(n) % 2 * 2 - 1, kind="stable")] if n > 1 else 0
    point(ax, x + off, vals, c, ms=ms, filled=filled, mew=0.15 if filled else 0.35)
    q = float(np.median(vals))
    ax.plot([x - half - 0.05, x + half + 0.05], [q, q], color=INK, lw=0.8, zorder=6, solid_capstyle="butt")


def corner_dot(x, y, w, leg=1.05):
    """a red triangle in the top-right corner of a tile (validated); it stays clear of the printed value"""
    xr = x + w - 0.15
    pg.c.add_patch(Polygon([(xr - leg, y), (xr, y), (xr, y + leg)], closed=True, fc=SIG, ec="none", zorder=4))


pg = Page(240)

# ================================================================ a  the search, both passes
pg.panel(1, 3, "a", "Search to independent models")
r13, r15 = (lambda X: X[X["round"] == 13]), (lambda X: X[X["round"] == 15])
el13 = int(SCR13.eligible.astype(str).str.upper().eq("TRUE").sum()); el15 = int(SCR15.eligible.astype(str).str.upper().eq("TRUE").sum())
assert (el13, el15) == tuple(PM.groupby("round").size().loc[[13, 15]])
steps = [("found by the search", sum(v["count"] for v in SEARCH.values()), None, None),
         ("screened", len(SCR13) + len(SCR15), len(SCR13), len(SCR15)),
         ("eligible", el13 + el15, el13, el15),
         ("scored", int(PM.scored.sum()), int(r13(PM).scored.sum()), int(r15(PM).scored.sum())),
         ("with decoys (≥ 3 genes)", len(W), len(r13(W)), len(r15(W))),
         ("independent of the cohort", len(IND), len(r13(IND)), len(r15(IND))),
         ("paper used the cohort", len(USED), len(r13(USED)), len(r15(USED)))]
AY, AH_ = 12.5, 33
ax = pg.ax(23, AY, 10.5, AH_)
y = np.arange(len(steps))[::-1]
ax.plot([s[1] for s in steps[:6]], y[:6], color=OTHER, lw=0.5, zorder=2)
point(ax, [s[1] for s in steps], y, INK, ms=3.0)
ax.set_xscale("log"); ax.set_xlim(200, 2500); ax.set_xticks([300, 1000], ["300", "1,000"]); ax.minorticks_off()
ax.set_yticks(y, [s[0] for s in steps], fontsize=5); ax.set_ylim(-0.6, len(steps) - 0.4)
ax.set_xlabel("records or models")
colhead(39.3, AY - 1.6, "all", ha="right"); colhead(44.2, AY - 1.6, "1st", ha="right"); colhead(50.3, AY - 1.6, "2nd", ha="right")
pg.line(41.2, AY - 3.3, 50.4, AY - 3.3, lw=0.4); pg.txt(45.8, AY - 4.6, "pass", ha="center")
for yy, s in zip(y, steps):
    yp = ymm(AY, AH_, -0.6, len(steps) - 0.4, yy)
    for xc, v in ((39.3, s[1]), (44.2, s[2]), (50.3, s[3])):
        if v is not None:
            pg.txt(xc, yp, f"{v:,}", ha="right")
pg.src("a", pd.DataFrame(steps, columns=["step", "whole_search", "first_pass_round13", "second_pass_round15"]))

# ================================================================ b  validated, TCGA and independent
pg.panel(53, 3, "b", "Validated")
ax = pg.ax(59.5, AY, 14.5, AH_)
dt = DEC.assign(v=((DEC.tcga_p < 0.05) & (DEC.tcga_worse == 1)).astype(int)).groupby(["model", "class"]).v.mean().unstack()
dt = dt.loc[IND.model]
vals = {"TCGA": (IND.tcga_validated.mean(), dt["coex"].mean(), dt["rand"].mean()),
        "indep.": (IND.validated.mean(), IND.d_coex.mean(), IND.d_rand.mean())}
rows = []
for i, (k, v) in enumerate(vals.items()):
    for j, (val, u) in enumerate(zip(v, ("sig", "coex", "rand"))):
        xx = i + (j - 1) * 0.31
        point(ax, xx, 100 * val, UNIT[u], MARK[u], ms=3.2)
        ax.text(xx, 100 * val + 3.5, f"{100 * val:.0f}", ha="center", va="bottom", fontsize=5)
        rows.append((k, u, 100 * val))
ax.set_xticks([0, 1], list(vals), fontsize=5.5); ax.set_xlim(-0.55, 1.5); ax.set_ylim(0, 95); ax.set_yticks([0, 25, 50, 75])
ax.set_ylabel("validated (%)")
pg.src("b", pd.DataFrame(rows, columns=["cohort", "unit", "validated_pct"]).assign(n_models=len(IND)))

# ================================================================ c  one FDR forest: first pass, whole search, assays
pg.panel(77, 3, "c", "False discovery rate")
FX0, FY0, FH = 77.5, 11.6, 2.25                          # label x, first row y, row pitch (mm)
fr = [("h", "First pass, 294 records"),
      ("r", "independent", AM13, "primary", "coex"),
      ("r", "independent, random genes", AM13, "primary", "rand"),
      ("r", "all models", AM13, "all models", "coex"),
      ("r", "paper used the cohort", AM13, "models whose", "coex"),
      ("h", "Whole search, 1,892 records"),
      ("r", "independent", AM, "primary", "coex"),
      ("r", "independent, random genes", AM, "primary", "rand"),
      ("r", "all models", AM, "all models", "coex"),
      ("r", "paper used the cohort", AM, "models whose", "coex"),
      ("r", "without hazard-ratio weights", AM, "sensitivity", "coex"),
      ("r", "first-pass models, new decoys", AM, "independent, round-13", "coex"),
      ("h", "Clinical assays"),
      ("a", "primary cohorts", "primary"), ("a", "all built cohorts", "all_built"), ("a", "ER-positive", "erpos")]
nrow = len(fr)
CX = {"n": 107.5, "v": 114.5, "d": 122.5}
for k, h in (("n", "n"), ("v", "valid."), ("d", "decoys")):
    colhead(CX[k], FY0 - 2.9, h, ha="right")
colhead(159.5, FY0 - 2.9, "FDR, % (95% interval)", ha="left")
ax = pg.ax(125.5, FY0 - FH / 2, 32, FH * nrow)
ax.set_ylim(nrow - 0.5, -0.5); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlim(0, 100); ax.set_xticks([0, 25, 50, 75, 100], ["0", "", "50", "", "100"]); ax.set_xlabel("estimated FDR (%)", labelpad=1)
out = []
for i, row in enumerate(fr):
    yp = FY0 + FH * i
    if row[0] == "h":
        if row[1] != "Clinical assays" or HAVE_ASSAY:
            pg.txt(FX0 - 0.5, yp, row[1], weight="bold")
        continue
    if row[0] == "r":
        _, lab, T, pop, dec = row
        f = am(T, pop, "fdr", dec); M = am(T, pop, "M", dec).estimate; v = am(T, pop, "validated", dec).estimate
        ds = am(T, pop, "decoy_share", dec).estimate
        est, lo, hi, u, grp = f.estimate, f.lo, f.hi, dec, ("first pass" if T is AM13 else "whole search")
    else:
        if not HAVE_ASSAY:
            continue
        _, lab, s = row
        h7 = AH.loc[s]
        M, v, ds = h7.units, h7.validated, h7.decoy_expected / h7.units
        est, lo, hi, u, grp = h7.fdr_assays, h7.fdr_lo, h7.fdr_hi, "coex", "clinical assays"
    pg.txt(FX0 + 1.0, yp, lab)
    pg.txt(CX["n"], yp, f"{M:,.0f}", ha="right"); pg.txt(CX["v"], yp, f"{v:.0f}", ha="right")
    pg.txt(CX["d"], yp, f"{100 * ds:.1f}", ha="right")
    ci(ax, 100 * lo, min(100 * hi, 100), i, UNIT[u], lw=0.7)
    point(ax, 100 * est, i, UNIT[u], MARK[u], ms=3.0)
    pg.txt(159.5, yp, f"{100 * est:.1f} ({100 * lo:.1f}–{100 * hi:.1f})")
    out.append((grp, lab, u, M, v, ds, est, lo, hi))
if HAVE_ASSAY:
    pg.txt(159.5, FY0 + FH * 12, f"P = {AH.loc['primary'].p_one_sided:.2g}")
key(pg, 133, 4.1, [(COEX, "co-expressed decoys", "s"), (RAND, "random genes", "D")])
pg.src("c", pd.DataFrame(out, columns=["group", "population", "decoys", "n", "validated", "decoy_share", "fdr", "lo", "hi"]))

# ================================================================ d  427 independent models among their decoys
DY = 55.5
pg.panel(1, DY, "d", "Independent models among their co-expressed decoys")
Q = DEC[DEC["class"] == "coex"].groupby("model").ext_cdir.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
ax = pg.ax(11, DY + 5.5, 168, 25)
xpos, xc, out = 0, {}, []
for tp in TYPES:
    o = IND[IND.type == tp].sort_values("ext_cdir")
    slot = max(len(o), 13)
    x0 = xpos + (slot - len(o)) / 2
    for j, (_, r) in enumerate(o.iterrows()):
        q = Q.loc[r.model]; xx = x0 + j
        ax.plot([xx, xx], [q[0.05], q[0.95]], color=tint(COEX, 0.38), lw=0.35, solid_capstyle="butt", zorder=2)
        ax.plot([xx, xx], [q[0.25], q[0.75]], color=COEX, lw=0.55, solid_capstyle="butt", zorder=3)
        out.append((r.model, tp, r.ext_cdir, r.validated, r.rho, *q.values))
    for v_, f_ in ((1, True), (0, False)):
        s = o.validated.values == v_
        point(ax, (x0 + np.arange(len(o)))[s], o.ext_cdir.values[s], SIG, ms=1.9, filled=f_, mew=0.2 if f_ else 0.4, zorder=5)
    xc[tp] = xpos + slot / 2 - 0.5
    xpos += slot + 6
zero(ax, 0.5, axis="y")
ax.set_xlim(-3, xpos - 3); ax.set_ylim(0.37, 0.70); ax.set_yticks([0.4, 0.5, 0.6, 0.7])
ax.set_xticks([xc[t] for t in TYPES], TYPES, fontsize=5.5); ax.tick_params(axis="x", length=0, pad=2)
ax.spines["bottom"].set_visible(False)
ax.set_ylabel("C-index, independent")
key(pg, 103, DY + 0.1, [(SIG, "validated", "o"), (SIG, "not validated", "open-o"), (COEX, "decoys, 25–75%", "line"),
                        (tint(COEX, 0.38), "5–95%", "line")])
pg.src("d", pd.DataFrame(out, columns=["model", "type", "cdir", "validated", "rho", "d05", "d25", "d50", "d75", "d95"]))

# ================================================================ e, f  rho
RY = 93.0
for (lt, xa, S, name, pop) in (("e", 8, IND, "independent", "primary"), ("f", 41.5, USED, "paper used the cohort", "models whose")):
    pg.panel(xa - 7, RY, lt, f"ρ, {name}")
    ax = pg.ax(xa, RY + 6, 24.5, 28.5)
    hist(ax, S.rho, np.linspace(0, 1, 11), weights=100 * np.ones(len(S)) / len(S))
    zero(ax, 10, axis="y")
    ax.set_ylim(0, 52); ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"]); ax.set_yticks([0, 25, 50])
    r_ = am(AM, pop, "rho")
    ax.text(0.03, 51, f"n = {len(S)}\nmean ρ {S.rho.mean():.2f}\n({r_.lo:.2f}–{r_.hi:.2f})", va="top", fontsize=5, linespacing=1.15)
    ax.set_xlabel("ρ"); ax.set_ylabel("models (%)")
    pg.src(lt, S[["model", "type", "rho"]])

# ================================================================ g  by type: counts and the independent models' audit
pg.panel(68, RY, "g", "Models by type")
cols = [("all", "all", "{:.0f}"), ("ind", "indep.", "{:.0f}"), ("used", "used", "{:.0f}"), ("genes", "genes", "{:.0f}"),
        ("validated", "valid.", "{:.0f}"), ("decoy_exp", "exp.", "{:.1f}"), ("rho", "ρ", "{:.2f}"), ("fdr", "FDR", "{:.0f}")]
G = pd.DataFrame({"all": W.groupby("type").size(), "ind": IND.groupby("type").size(), "used": USED.groupby("type").size(),
                  "genes": W.groupby("type").n_genes.median()}).join(BYT[["validated", "decoy_exp", "rho", "fdr"]]).loc[TYPES]
G["fdr"] = 100 * G.fdr
x0, y0, cw, ch = 76.5, RY + 10.6, 5.75, 2.72
for j, (c, h, _) in enumerate(cols):
    colhead(x0 + cw * j + cw / 2, y0 - 1.6, h)
pg.line(x0 + cw * 4 + 0.3, y0 - 4.0, x0 + cw * 8 - 0.4, y0 - 4.0, lw=0.4)
pg.txt(x0 + cw * 6, y0 - 5.3, "independent", ha="center")
for i, (tp, r) in enumerate(G.iterrows()):
    pg.txt(x0 - 0.8, y0 + ch * i + ch / 2, tp, ha="right")
    for j, (c, _, f) in enumerate(cols):
        val = r[c]
        if not np.isfinite(val):
            tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, EMPTY, "–")
            continue
        fc = {"all": BLUE(val / 140), "ind": BLUE(val / 100), "used": BLUE(val / 100), "genes": BLUE((val - 2) / 10),
              "validated": BLUE(val / 80), "decoy_exp": PURPLE(val / 35), "rho": BLUE((val - 0.3) / 0.65), "fdr": PURPLE(min(val, 100) / 100)}[c]
        tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, fc, f.format(val))
yb = y0 + ch * len(G) + 0.25
pg.txt(x0 - 0.8, yb + ch / 2, "all", ha="right", weight="bold")
tot = {"all": len(W), "ind": len(IND), "used": len(USED), "genes": W.n_genes.median(), "validated": IND.validated.sum(),
       "decoy_exp": IND.d_coex.sum(), "rho": IND.rho.mean(), "fdr": 100 * am(AM, "primary", "fdr").estimate}
for j, (c, _, f) in enumerate(cols):
    pg.txt(x0 + cw * j + cw / 2 - 0.07, yb + ch / 2, f.format(tot[c]), ha="center", weight="bold")
pg.src("g", G.reset_index().rename(columns={"index": "type", "ind": "independent", "fdr": "fdr_pct"}))

# ================================================================ h  by theme family, independent models
pg.panel(125, RY, "h", "By theme")
key(pg, 149, RY + 0.1, [(SIG, "validated", "o"), (COEX, "decoys", "s")])
b = BYF.sort_values(["models", "family"]).reset_index(drop=True)
HY, HH = RY + 6.5, 34.0
ax = pg.ax(141.5, HY, 24, HH)
for i, r in b.iterrows():
    ax.plot([100 * r.decoy_share, 100 * r.validated_share], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, 100 * r.decoy_share, i, COEX, "s", ms=2.4)
    point(ax, 100 * r.validated_share, i, SIG, ms=2.6)
    yp = ymm(HY, HH, -0.7, len(b) - 0.3, i)
    pg.txt(172.3, yp, f"{int(r.models)}", ha="right")
    pg.txt(179.6, yp, "–" if not np.isfinite(r.fdr) else f"{100 * r.fdr:.0f}", ha="right")
colhead(172.3, HY - 1.2, "n", ha="right"); colhead(179.6, HY - 1.2, "FDR", ha="right")
ax.set_yticks(range(len(b)), [f.replace("immune-related", "immune") for f in b.family], fontsize=5)
ax.set_ylim(-0.7, len(b) - 0.3); ax.set_xlim(0, 70); ax.set_xticks([0, 25, 50]); ax.set_xlabel("models (%)")
pg.src("h", BYF)

# ================================================================ i  by years and by size
IY = 142.5
pg.panel(1, IY, "i", "By year and size")
rows_i = [("2019–2021", BYY.loc["2019-2021"]), ("2022–2025", BYY.loc["2022-2025"]),
          ("2–5 genes", BYS.loc["2-5"]), ("6–10 genes", BYS.loc["6-10"]), ("11+ genes", BYS.loc["11+"])]
YI = [5.4, 4.4, 2.2, 1.2, 0.2]
ax = pg.ax(15.5, IY + 7, 17, 21)
for yy, (lab, r) in zip(YI, rows_i):
    ax.plot([100 * r.decoy_share, 100 * r.validated_share], [yy, yy], color=OTHER, lw=0.5, zorder=2)
    point(ax, 100 * r.decoy_share, yy, COEX, "s", ms=2.4)
    point(ax, 100 * r.validated_share, yy, SIG, ms=2.6)
    yp = ymm(IY + 7, 21, -0.5, 6.0, yy)
    pg.txt(38.5, yp, f"{int(r.models)}", ha="right"); pg.txt(44.8, yp, f"{100 * r.fdr:.0f}", ha="right")
ax.set_ylim(-0.5, 6.0); ax.set_yticks(YI, [r[0] for r in rows_i], fontsize=5)
ax.set_xlim(0, 60); ax.set_xticks([0, 25, 50]); ax.set_xlabel("models (%)")
colhead(38.5, IY + 5.6, "n", ha="right"); colhead(44.8, IY + 5.6, "FDR", ha="right")
pg.src("i", pd.concat([BYY.reset_index().rename(columns={"years": "group"}).assign(by="years"),
                       BYS.reset_index().rename(columns={"size": "group"}).assign(by="genes")]))

# ================================================================ j  citations per year by rho (exploratory)
pg.panel(47, IY, "j", "Citations per year")
JX, JY, JW, JH = 55.5, IY + 7, 24, 24
ax = pg.ax(JX, JY, JW, JH)
tr = lambda v: np.log10(1 + np.asarray(v, float))  # noqa: E731
pos = {("independent", True): 0, ("independent", False): 1.3, ("all models with decoys", True): 3.0, ("all models with decoys", False): 4.3}
srcj = []
for (pop, hi), xx in pos.items():
    S = CIT[~CIT.uses_external.astype(bool)] if pop == "independent" else CIT
    v = S[(S.rho >= 0.95) == hi].cites_per_year.dropna()
    assert len(v) == CITS.loc[pop, "n_hi" if hi else "n_lo"]
    sina(ax, xx, tr(v), SIG, half=0.5, bw=0.06, step=0.05, ms=1.15, filled=hi)
for pop, xx in (("independent", 0.65), ("all models with decoys", 3.65)):
    ax.text(xx, tr(80), f"P = {CITS.loc[pop, 'p']:.2g}", ha="center", va="bottom", fontsize=5)
    ax.plot([xx - 0.65, xx + 0.65], [tr(70), tr(70)], color=INK, lw=0.4)
yt = [0, 1, 3, 10, 30]
ax.set_yticks(tr(yt), [str(t) for t in yt]); ax.set_ylim(0, tr(140))
ax.set_xticks(list(pos.values()), ["≥ 0.95", "< 0.95", "≥ 0.95", "< 0.95"], fontsize=5)
ax.set_xlim(-0.7, 5.0); ax.set_ylabel("citations per year")
pg.txt(JX - 1.6, JY + JH + 2.0, "ρ", ha="right")
for xx, lab in ((0.65, "independent"), (3.65, "all")):
    pg.txt(JX + JW * (xx + 0.7) / 5.7, JY + JH + 5.0, lab, ha="center")
pg.src("j", CIT)

# ================================================================ k, l  clinical assays in each cohort
KX, KY, kw, kh = 97.0, 156.4, 4.55, 3.1                    # tile origin, cell width, cell height (mm)
pg.panel(81, IY, "k", "Clinical assays, ρ per cohort")
pg.panel(150.5, IY, "l", "Cohorts validated")
KGAP = 1.2
CY = KY + kh * len(ASSAYS) + 1.2                            # the colon row


def kcol(j, n_first):
    return KX + kw * j + (KGAP if j >= n_first else 0)


def lab_of(c):
    return c.replace("_ERpos", "").split("_", 1)[-1].replace("GSE", "")


if HAVE_ASSAY:
    UA = AU.set_index(["assay", "cohort"])
    srck = []
    for cols_, rows_, ytop, first, above in ((BR_PRIM + BR_ER, ASSAYS, KY, len(BR_PRIM), True), (CO_PRIM + CO_OVL, [COLON], CY, len(CO_PRIM), False)):
        for j, c in enumerate(cols_):
            xm = kcol(j, first) + (kw - 0.15) / 2
            if above:
                pg.txt(xm, ytop - 0.6, lab_of(c), ha="left", va="center", rotation=90, rotation_mode="anchor")
            else:
                pg.txt(xm, ytop + kh + 0.45, lab_of(c), ha="right", va="center", rotation=90, rotation_mode="anchor")
        for i, (a, al) in enumerate(rows_):
            yr = ytop + kh * i
            pg.txt(KX - 0.8, yr + kh / 2, al, ha="right")
            for j, c in enumerate(cols_):
                xx = kcol(j, first)
                if (a, c) not in UA.index:
                    tile(pg, xx, yr, kw, kh, EMPTY)
                    continue
                r = UA.loc[(a, c)]
                tile(pg, xx, yr, kw, kh, BLUE(float(r.rho)), f"{r.rho:.2f}")
                if r.validated == 1:
                    corner_dot(xx, yr, kw)
                srck.append(dict(assay=a, cohort=c, set=r["set"], rho=r.rho, validated=r.validated, cdir=r.cdir, d_coex=r.d_coex, d_rand=r.d_rand))
        # group labels with hairlines: primary, then the set apart
        yg = ytop - 10.6 if above else ytop + kh + 7.8
        for (j0, j1, lab) in ((0, first, "primary"), (first, len(cols_), "ER-positive" if above else "overlapping")):
            xa_, xb_ = kcol(j0, first) + 0.2, kcol(j1 - 1, first) + kw - 0.35
            pg.line(xa_, yg + (0.9 if above else -0.9), xb_, yg + (0.9 if above else -0.9), lw=0.4)
            pg.txt((xa_ + xb_) / 2, yg - (0.35 if above else -0.35), lab, ha="center")
    pg.src("k", pd.DataFrame(srck))
    # l: per assay, cohorts validated against decoys, rows aligned with k
    LX, LW = 152.5, 16.5
    ax = pg.ax(LX, KY, LW, CY + kh - KY)
    ax.set_ylim(CY + kh, KY); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    srcl = []
    for (a, al), yc in [(x, KY + kh * i + kh / 2) for i, x in enumerate(ASSAYS)] + [(COLON, CY + kh / 2)]:
        r = ABA.loc[a]
        sv, sd, sr = 100 * r.validated / r.cohorts, 100 * r.decoy_expected / r.cohorts, 100 * r.decoy_rand / r.cohorts
        ax.plot([sd, sv], [yc, yc], color=OTHER, lw=0.5, zorder=2)
        point(ax, sr, yc, RAND, "D", ms=2.0, zorder=3)
        point(ax, sd, yc, COEX, "s", ms=2.4)
        point(ax, sv, yc, SIG, ms=2.6)
        pg.txt(173.2, yc, f"{int(r.cohorts)}", ha="right")
        pg.txt(179.6, yc, "–" if not np.isfinite(r.fdr) else f"{100 * r.fdr:.0f}", ha="right")
        srcl.append(dict(assay=a, cohorts=r.cohorts, validated=r.validated, decoy_expected=r.decoy_expected, decoy_rand=r.decoy_rand,
                         fdr=r.fdr, mean_rho=r.mean_rho))
    colhead(173.2, KY - 1.6, "n", ha="right"); colhead(179.6, KY - 1.6, "FDR", ha="right")
    ax.set_xlim(-3, 103); ax.set_xticks([0, 50, 100]); ax.set_xlabel("cohorts (%)")
    pg.src("l", pd.DataFrame(srcl))

# ================================================================ m  pipeline by cancer type
PY = 188.5
pg.panel(1, PY, "m", "Pipeline, by cancer type")
ax = pg.ax(8.5, PY + 6, 28, 34)
by = PU.groupby("type").agg(t=("validated", "mean"), d=("d_coex", "mean"), r=("d_rand", "mean")).sort_values("d")
for i, (tp, r) in enumerate(by.iterrows()):
    ax.plot([100 * r.d, 100 * r.t], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, 100 * r.d, i, COEX, "s", ms=2.6)
    point(ax, 100 * r.t, i, SIG, ms=2.8)
    point(ax, 100 * r.r, i, RAND, "D", ms=2.0, zorder=3)
ax.set_yticks(range(len(by)), by.index, fontsize=5); ax.set_ylim(-0.7, len(by) - 0.3)
ax.set_xlim(-2, 100); ax.set_xticks([0, 50, 100]); ax.set_xlabel("pools validated (%)")
fb = PA[(PA.population == "theme pools") & (PA.quantity == "fdr") & (PA.decoys == "coex")].iloc[0]
ax.text(100, 0.6, f"FDR {100 * fb.estimate:.0f}%\n({100 * fb.lo:.0f}–{100 * fb.hi:.0f}%)", ha="right", va="center", fontsize=5, linespacing=1.15)
for yy, (c, lab, k) in zip((PY + 20, PY + 23, PY + 26), ((SIG, "15 theme", "o"), (COEX, "20 co-expr.", "s"), (RAND, "5 random", "D"))):
    key(pg, 23, yy, [(c, lab, k)])
pg.src("m", by.reset_index().rename(columns={"t": "theme_pools", "d": "coex_decoy_pools", "r": "random_decoy_pools"}))

# ================================================================ n  pipeline by theme pool
pg.panel(39.5, PY, "n", "Pipeline, by theme")
ax = pg.ax(55, PY + 6, 19, 39)
bp = PU.groupby("pool").agg(v=("validated", "sum"), e=("d_coex", "sum")).sort_values("e")
for i, (p, r) in enumerate(bp.iterrows()):
    ax.plot([r.e, r.v], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, r.e, i, COEX, "s", ms=2.6)
    point(ax, r.v, i, SIG, ms=2.8)
ax.set_yticks(range(len(bp)), [THEME.get(p, p) for p in bp.index], fontsize=5); ax.set_ylim(-0.7, len(bp) - 0.3)
ax.set_xlim(0, 7); ax.set_xticks([0, 2, 4, 6]); ax.set_xlabel("of 11 cancer types")
pg.src("n", bp.reset_index().rename(columns={"v": "validated_units", "e": "expected_from_decoys"}))

# ================================================================ o  pipeline, type by theme
pg.panel(77, PY, "o", "Decoy pools validated, by type and theme")
tt = by.index[::-1]; pp = bp.index[::-1]
x0, y0, cw, ch = 83.5, PY + 14.6, 4.42, 2.8
for j, p in enumerate(pp):
    pg.txt(x0 + cw * j + (cw - 0.15) / 2, y0 - 0.6, THEME.get(p, p), ha="left", va="center", rotation=90, rotation_mode="anchor")
M_ = PU.set_index(["type", "pool"])
for i, tp in enumerate(tt):
    pg.txt(x0 - 0.8, y0 + ch * i + ch / 2, tp, ha="right")
    for j, p in enumerate(pp):
        r = M_.loc[(tp, p)]
        tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, BLUE(float(r.d_coex)), f"{100 * r.d_coex:.0f}")
        if r.validated == 1:
            corner_dot(x0 + cw * j, y0 + ch * i, cw)
pg.cbar(117, 235.4, 16, 1.4, BLUE, 0, 100, [0, 50, 100], "co-expr. pools (%)")
pg.c.add_patch(Polygon([(84.1, 235.5), (85.4, 235.5), (85.4, 236.8)], closed=True, fc=SIG, ec="none", zorder=4))
pg.txt(86.2, 236.1, "theme model validated")
pg.src("o", PU[["type", "pool", "d_coex", "validated"]])

# ================================================================ p  pipeline models, TCGA and independent
C_ = PD[(PD["class"] == "coex") & PD.ext_cdir.notna() & PD.tcga_cdir.notna()]
pg.panel(152, PY, "p", "Pipeline models")
ax = pg.ax(161, PY + 6, 16.5, 15)
ax.hexbin(C_.tcga_cdir, C_.ext_cdir, gridsize=13, cmap=GREYHEX, mincnt=1, linewidths=0, rasterized=True, extent=(0.5, 0.9, 0.3, 0.75))
T_ = PU[PU.ext_cdir.notna()]
for v, f in ((1, True), (0, False)):
    s = T_.validated == v
    point(ax, T_.tcga_cdir[s], T_.ext_cdir[s], SIG, ms=1.8, filled=f, mew=0.35 if not f else 0.2)
zero(ax, 0.5, axis="y")
ax.set_xlim(0.5, 0.9); ax.set_ylim(0.3, 0.75); ax.set_xticks([0.5, 0.7, 0.9]); ax.set_yticks([0.3, 0.5, 0.7])
ax.set_xlabel("C-index, TCGA"); ax.set_ylabel("C-index, ind.")
pg.src("p", pd.concat([C_[["type", "pool", "kind", "tcga_cdir", "ext_cdir", "validated"]],
                       T_.assign(kind="target")[["type", "pool", "kind", "tcga_cdir", "ext_cdir", "validated"]]]))

# ================================================================ q  pipeline model size
pg.panel(152, PY + 28.5, "q", "Model size")
ax = pg.ax(161, PY + 34, 16.5, 9.5)
bins = np.arange(-0.5, 38.5, 3)
nc = (PD["class"] == "coex").sum()
hist(ax, PD[PD["class"] == "coex"].n_model, bins, weights=100 * np.ones(nc) / nc)
h, _ = np.histogram(PU.n_model, bins)
ax.step(bins[:-1], 100 * h / len(PU), where="post", color=SIG, lw=0.8, zorder=4)
ax.set_xlabel("genes in model"); ax.set_ylabel("pools (%)"); ax.set_xlim(-0.5, 37); ax.set_xticks([0, 12, 24, 36]); ax.set_yticks([0, 20])
pg.src("q", pd.concat([PD[PD["class"] == "coex"][["type", "pool", "n_model"]].assign(kind="decoy"), PU[["type", "pool", "n_model"]].assign(kind="target")]))

pg.save("Fig8v18")
