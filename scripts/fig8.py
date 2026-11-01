"""Figure 8 (GB v14: panels f and g now show the independent cohorts of round 11,
validation/PLAN_round11_external_replication.md; the v13 panels f and g are Additional file 1 Fig. S8, figS8_v14.py)
Figure 8 (GB v13) - gene content and co-expression structure, by endpoint, family and cancer type
(validation/PLAN_round10_content_and_structure.md; icb/v15_content from decoys/v15_content_structure.py).

rho(s, c) = the share of a signature's 320 co-expressed decoys whose gain over the null in cohort c lies below its own
(0.5 when its genes add nothing to a module of its structure).
a  rho over every signature x cohort test, by family, for response and for overall survival
b  each signature's mean rho for response against overall survival; content calls
c  family means of rho with intervals over correlation clusters, per setting
d  the atlas: rho for every signature in every response cohort and TCGA type; calls and single-cohort estimates above
e  the five tumor types of the response cohorts: inflamed and confounder rho, response against survival
f  single-cohort estimates against the inflamed family's rho
g  libraries: structure (decoys' per-test clearance) against content (mean rho of called targets)
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
from matplotlib.colors import LinearSegmentedColormap, Normalize
K.apply()
V = K.D + "/icb/v15_content/"
FAMD = pd.read_csv(K.D + "/icb/v5_signature_clusters.csv").set_index("signature").family
FAMS = ["Inflamed/IFN", "Other", "Confounder"]
CELL = {k: pd.read_csv(V + f"cells_{k}.csv.gz") for k in ("response", "tcga_OS")}
SIG = {k: pd.read_csv(V + f"signature_{k}.csv", index_col=0) for k in ("response", "tcga_OS")}
FBS = pd.read_csv(V + "family_by_setting.csv")
PC = pd.read_csv(V + "primary_contrast.csv")
MT = pd.read_csv(V + "matched_tumor_types.csv")
AT = pd.read_csv(V + "atlas_by_type.csv")
LIB = pd.read_csv(V + "libraries.csv")
ICB_TUM = {"ICB_Gide": "melanoma", "ICB_Hugo": "melanoma", "ICB_Liu": "melanoma", "ICB_Riaz": "melanoma",
           "ICB_Van_Allen": "melanoma", "ICB_Puch": "melanoma", "ICB_Braun": "kidney", "ICB_Miao1": "kidney",
           "ICB_Jung": "lung", "ICB_Limagne1": "lung", "ICB_Mariathasan": "bladder", "ICB_Kim": "gastric"}
RHO = LinearSegmentedColormap.from_list("rho", [K.AQUA, "#9fd3c3", "#f4f4f2", "#9aa3ad", "#2b3440"])
SETL = {"response": "response", "tcga_OS": "overall survival"}

C = K.Canvas(240)

# ================================================================== a  rho over every test
ax = C.panel("a", 14, 12, 40, 38, head="Place among decoys, every test")
bins = np.linspace(0, 1, 21)
ys = {}
for j, k in enumerate(["tcga_OS", "response"]):
    cc = CELL[k].assign(family=lambda d: d.target.map(FAMD))
    for i, f in enumerate(FAMS):
        y0 = (j * 3.75) + i * 1.12
        ys[(k, f)] = y0
        h, _ = np.histogram(cc[cc.family == f].rho, bins=bins, density=True)
        hh = np.minimum(h / 6.0, 1.05)
        ax.fill_between(np.repeat(bins, 2)[1:-1], y0, y0 + np.repeat(hh, 2), color=K.tint(K.FAM[f], 0.55), lw=0,
                        zorder=3)
        ax.step(bins, y0 + np.r_[hh, hh[-1]], where="post", color=K.FAM[f], lw=0.6, zorder=4)
        ax.plot([0, 1], [y0 + 1 / 6.0] * 2, color=K.INK2, lw=0.45, ls=(0, (2, 1.5)), zorder=5)
        ax.text(-0.03, y0 + 0.12, K.FAMLAB[f], ha="right", va="bottom", fontsize=K.T_SMALL, color=K.FAM[f])
    ax.text(0.0, j * 3.75 + 3.3, SETL[k], ha="left", va="bottom", fontsize=K.T_NOTE, fontweight="bold", color=K.INK)
ax.set_ylim(-0.15, 7.75); ax.set_xlim(0, 1); ax.set_yticks([]); K.frame(ax, left=False)
ax.set_xticks([0, .5, 1]); ax.set_xlabel("ρ, signature's place among its decoys")
ax.text(1.0, 7.7, "dashed: genes add nothing", ha="right", va="top", fontsize=K.T_SMALL - 0.3, color=K.INK2)

# ================================================================== b  response vs survival, per signature
ax = C.panel("b", 76, 12, 42, 38, head="Response against survival")
r, o = SIG["response"], SIG["tcga_OS"].reindex(SIG["response"].index)
ax.fill_between([0.2, 1], [0.2, 1], [0.2, 0.2], color=K.tint(K.BLUE, 0.92), lw=0, zorder=0)
ax.fill_between([0.2, 1], [0.2, 1], [1, 1], color=K.tint(K.ORANGE, 0.92), lw=0, zorder=0)
K.ref(ax, 0.5); K.ref(ax, 0.5, axis="v"); K.identity(ax, 0.2, 1)
for f in ["Other", "Confounder", "Inflamed/IFN"]:
    sel = r.index[r.family == f]
    for s in sel:
        cc_ = bool(r.loc[s, "content_call"]) or bool(o.loc[s, "content_call"])
        K.dots(ax, [r.loc[s, "rho_bar"]], [o.loc[s, "rho_bar"]], K.FAM[f], s=13 if cc_ else 9, filled=cc_, z=5)
    ax.scatter([r.loc[sel, "rho_bar"].mean()], [o.loc[sel, "rho_bar"].mean()], marker="+", s=46, color=K.FAM[f],
               lw=1.2, zorder=6)
ax.set_xlim(0.2, 1); ax.set_ylim(0.2, 1); ax.set_xticks([.2, .4, .6, .8, 1]); ax.set_yticks([.2, .4, .6, .8, 1])
ax.set_xlabel("mean ρ, response (12 cohorts)"); ax.set_ylabel("mean ρ, overall survival (20 types)")
p0 = PC[PC.b == "tcga_OS"].iloc[0]
ax.text(0.98, 0.04, f"inflamed − confounder\n{p0.contrast_inflamed_minus_confounder:.2f} "
        f"({p0.lo:.2f}–{p0.hi:.2f}), P = {K.sci(p0.p_perm)}", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=K.T_SMALL - 0.2, color=K.INK, linespacing=1.05)
ax.text(0.98, 0.25, "more for response", transform=ax.transAxes, ha="right", va="center", fontsize=K.T_SMALL,
        color=K.BLUE, linespacing=1.0)
ax.text(0.03, 0.97, "more for survival", transform=ax.transAxes, ha="left", va="top", fontsize=K.T_SMALL,
        color=K.ORANGE, linespacing=1.0)
K.legend(ax, [(K.INK, "content call", "o", True), (K.INK, "none", "o", False), (K.INK, "family mean", "+")],
         x=0.06, y=0.86, dy=0.075, s=7, size=K.T_SMALL)

# ================================================================== c  family means per setting
ax = C.panel("c", 138, 12, 40, 38, head="Family means")
SETS = [("response", "response"), ("tcga_OS", "OS"), ("tcga_DSS", "DSS"), ("tcga_PFI", "PFI"),
        ("tcga_OS_matched", "OS, same tumor types")]
yt = []
for i, (k, lab_) in enumerate(SETS):
    for jf, f in enumerate(FAMS):
        y = i * 1.0 + (jf - 1) * 0.24
        q = FBS[(FBS.setting == k) & (FBS.family == f)].iloc[0]
        K.hbar_interval(ax, y, q.lo, q.hi, K.tint(K.FAM[f], 0.45), lw=1.3)
        K.dots(ax, [q.rho_bar], [y], K.FAM[f], s=11, z=5)
    yt.append(i)
K.ref(ax, 0.5, axis="v")
ax.set_yticks(yt); ax.set_yticklabels([s[1] for s in SETS], fontsize=K.T_SMALL); ax.set_ylim(len(SETS) - 0.5, -0.5)
ax.set_xlim(0.15, 1.0); ax.set_xticks([.25, .5, .75, 1]); K.frame(ax, left=False); ax.set_xlabel("mean ρ")
K.legend(ax, [(K.FAM[f], K.FAMLAB[f]) for f in FAMS], x=0.03, y=0.985, dy=0.07, s=7, size=K.T_SMALL - 0.3, xpad=0.03)

# ================================================================== d  the atlas
r_cols = AT[AT.setting == "response"].assign(tum=lambda d: d.cohort.map(ICB_TUM))
TUMS = ["melanoma", "lung", "bladder", "gastric", "kidney"]
r_cols = r_cols.sort_values(["tum", "rho_inflamed"], key=lambda s: s.map(TUMS.index) if s.name == "tum" else -s)
o_cols = AT[AT.setting == "tcga_OS"].sort_values("rho_inflamed", ascending=False)
Cs = SIG["response"].rho_bar - SIG["tcga_OS"].rho_bar.reindex(SIG["response"].index)
rows = []
for f in FAMS:
    rows += sorted([s for s in Cs.index if FAMD[s] == f], key=lambda s: -Cs[s])
MR = CELL["response"].pivot_table(index="target", columns="cohort", values="rho").reindex(index=rows)[list(r_cols.cohort)]
MO = CELL["tcga_OS"].pivot_table(index="target", columns="cohort", values="rho").reindex(index=rows)[list(o_cols.cohort)]
SUMC = pd.DataFrame({"R": SIG["response"].rho_bar, "OS": SIG["tcga_OS"].rho_bar}).reindex(rows)
nr = len(rows); g1, g2 = 1.0, 1.0
xr = np.arange(MR.shape[1]); xo = xr[-1] + 1 + g1 + np.arange(MO.shape[1]); xs = xo[-1] + 1 + g2 + np.arange(2)
HX, HY, HW, HH = 35, 80, 131, 92
ax = C.ax(HX, HY, HW, HH)
C.letter(HX - 28.5, 60.8, "d"); C.head(HX - 23.5, 61.0, "Where the genes matter: every signature, every response cohort and cancer type")
norm = Normalize(0, 1)
for M_, xx in ((MR, xr), (MO, xo), (SUMC, xs)):
    ax.imshow(M_.to_numpy(), cmap=RHO, norm=norm, aspect="auto", interpolation="nearest",
              extent=(xx[0] - .5, xx[-1] + .5, nr - .5, -.5))
for k_, xx in (("response", xs[0]), ("tcga_OS", xs[1])):
    cc_ = SIG[k_].reindex(rows).content_call.to_numpy()
    yy_ = np.where(cc_)[0]
    ax.scatter([xx] * len(yy_), yy_, s=2.2, color="white", lw=0, zorder=5)
cum = 0
for f in FAMS[:-1]:
    cum += int((FAMD.reindex(rows) == f).sum())
    ax.axhline(cum - .5, color="white", lw=1.0, zorder=4)
for b_ in np.cumsum(r_cols.tum.value_counts().reindex(TUMS).to_numpy())[:-1]:
    ax.axvline(b_ - .5, color="white", lw=0.8, zorder=4)
ax.set_xlim(-0.6, xs[-1] + 0.6); ax.set_ylim(nr - .5, -.5)
ax.set_yticks(range(nr)); ax.set_yticklabels([K.lab(s) for s in rows], fontsize=3.9)
ax.set_xticks(list(xr) + list(xo) + list(xs))
TAB = {"melanoma": "mel", "lung": "lung", "bladder": "blad", "gastric": "gast", "kidney": "kid"}
ax.set_xticklabels([f"{K.COHORT[c]} ({TAB[ICB_TUM[c]]})" for c in r_cols.cohort] + list(o_cols.cohort) +
                   ["response", "OS"], rotation=90, fontsize=4.3)
ax.tick_params(axis="x", length=0, pad=1); ax.tick_params(axis="y", length=0, pad=4.6); [sp.set_visible(False) for sp in ax.spines.values()]
# family strip
axf = C.ax(HX - 1.5, HY, 1.0, HH); K.clean(axf)
axf.imshow(np.array([[0 if FAMD[s] == "Inflamed/IFN" else 1 if FAMD[s] == "Other" else 2] for s in rows]),
           cmap=LinearSegmentedColormap.from_list("f", [K.FAM[f] for f in FAMS], N=3), aspect="auto",
           interpolation="nearest")
# column headers and calls track
axt = C.ax(HX, HY - 11.5, HW, 10); axt.set_xlim(ax.get_xlim())
for co, xx in list(zip(r_cols.cohort, xr)) + list(zip(o_cols.cohort, xo)):
    q = AT[AT.cohort == co].iloc[0]
    K.stacked_calls(axt, xx, q.calls, min(q.exp_false, q.calls), K.SLATE, width=0.72)
    lbl = "–" if not q.calls else (f"{q.fdr:.0f}" if q.fdr < 1000 else ">999")
    axt.text(xx, q.calls + 0.8, lbl, ha="center", va="bottom", fontsize=3.6, rotation=90, color=K.INK2)
axt.set_ylim(0, 40); axt.set_yticks([0, 20]); axt.tick_params(axis="y", labelsize=4.3, length=1.2, pad=1)
axt.set_xticks([]); K.frame(axt, bottom=False); axt.spines["bottom"].set_visible(True)
axt.set_ylabel("calls", fontsize=K.T_SMALL, labelpad=1)
axt.text(xr.mean(), 41, "response to checkpoint blockade (12 cohorts)", ha="center", va="bottom", fontsize=K.T_NOTE,
         fontweight="bold", color=K.INK)
axt.text(xo.mean(), 41, "overall survival, TCGA (20 cancer types)", ha="center", va="bottom", fontsize=K.T_NOTE,
         fontweight="bold", color=K.INK)
axt.text(xs.mean(), 41, "mean", ha="center", va="bottom", fontsize=K.T_NOTE, fontweight="bold", color=K.INK)
axt.text(xs[-1] + 0.7, 30, "bars: signatures above the null's\n95th percentile (light: expected\nfalse); "
         "numbers: estimated FDR (%)", ha="right", va="top", fontsize=3.9, color=K.INK2, linespacing=1.0)
K.cbar(C, HX + HW + 3.5, HY + 2, 1.6, 26, RHO, 0, 1, [0, .5, 1], label="ρ", horizontal=False)
C.text(HX + HW + 2.2, HY + 31, "white dot:\ncontent\ncall", size=3.9, color=K.INK2, linespacing=1.0)

# ================================================================== e  same tumor types
ax = C.panel("e", 22, 199, 40, 31, head="Same tumor types")
for i, t_ in enumerate(TUMS):
    for jf, f in enumerate(["Inflamed/IFN", "Confounder"]):
        y = i + (jf - 0.5) * 0.34
        a_ = MT[(MT.tumor == t_) & (MT.family == f) & (MT.setting == "response")].rho.iloc[0]
        b_ = MT[(MT.tumor == t_) & (MT.family == f) & (MT.setting == "tcga_OS")].rho.iloc[0]
        ax.plot([a_, b_], [y, y], color=K.tint(K.FAM[f], 0.45), lw=1.0, zorder=3)
        K.dots(ax, [a_], [y], K.FAM[f], s=12, z=5)
        K.dots(ax, [b_], [y], K.FAM[f], s=12, filled=False, z=5)
K.ref(ax, 0.5, axis="v")
ax.set_yticks(range(len(TUMS))); ax.set_yticklabels(TUMS, fontsize=K.T_SMALL); ax.set_ylim(len(TUMS) - 0.4, -0.6)
ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1]); K.frame(ax, left=False); ax.set_xlabel("mean ρ of the family")
K.hlegend(C, 47, 194.6, [(K.INK, "response", "o", True), (K.INK, "survival", "o", False)], size=K.T_SMALL, s=7,
          gap=1.5)

# ================================================================== f  independent cohorts against TCGA (round 11)
X = K.D + "/icb/v17_external/"
BT = pd.read_csv(X + "by_type.csv"); PCX = pd.read_csv(X + "per_cohort.csv"); PRX = pd.read_csv(X + "predictions.csv")
PRI = pd.read_csv(X + "primary.csv").iloc[0]
ax = C.panel("f", 82, 199, 36, 31, head="Independent cohorts")
K.ref(ax, 0.5); K.ref(ax, 0.5, axis="v"); K.identity(ax, 0, 1)
for _, q in BT.iterrows():
    K.dots(ax, [q.tcga], [q.ext], K.BLUE, s=13, z=5)
    dy_ = {"CESC": 0.022, "SARC": -0.022}.get(q.type, 0)
    ax.text(q.tcga + 0.025, q.ext + dy_, q.type, fontsize=4.3, color=K.INK2, va="center", ha="left")
ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05); ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
ax.set_xlabel("inflamed mean ρ, TCGA"); ax.set_ylabel("inflamed mean ρ, independent")
ax.text(0.97, 0.03, f"Spearman {PRI.spearman:.2f}".replace("-", "\u2212") + f", P = {K.sci(PRI.p_one_sided)}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=K.T_SMALL - 0.2, color=K.INK)

# ================================================================== g  each independent cohort, and the response cohort
ax = C.panel("g", 140, 199, 38, 31, head="Per independent cohort")
K.legend(ax, [(K.FAM["Inflamed/IFN"], "inflamed"), (K.FAM["Confounder"], "confounder")], x=0.70, y=1.10, dy=0.065, s=6,
         size=K.T_SMALL - 0.3, xpad=0.03)
PCX = PCX.sort_values("tcga_rho_inflamed")
ys = np.arange(len(PCX))
for y, (_, q) in zip(ys, PCX.iterrows()):
    for f_, lab_, off in ((K.FAM["Inflamed/IFN"], "inflamed", 0.17), (K.FAM["Confounder"], "confounder", -0.17)):
        K.hbar_interval(ax, y + off, q[f"rho_{lab_}_lo"], q[f"rho_{lab_}_hi"], K.tint(f_, 0.45), lw=1.0)
        K.dots(ax, [q[f"rho_{lab_}"]], [y + off], f_, s=8, z=5)
yr = len(PCX) + 0.6
for key, f_, off in (("GSE176307 inflamed", K.FAM["Inflamed/IFN"], 0.17), ("GSE176307 confounder", K.FAM["Confounder"], -0.17)):
    rr = PRX[PRX.test.str.startswith(key)].iloc[0]
    K.hbar_interval(ax, yr + off, rr.lo, rr.hi, K.tint(f_, 0.45), lw=1.0)
    K.dots(ax, [rr.estimate], [yr + off], f_, s=9, marker="D", z=5)
K.ref(ax, 0.5, axis="v")
ax.set_yticks(list(ys) + [yr])
ax.set_yticklabels([f"{t} {g}" for t, g in zip(PCX.type, PCX.cohort.str.split("_").str[1])] + ["BLCA GSE176307, response"], fontsize=3.9)
ax.set_ylim(-0.6, yr + 0.7); ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1]); K.frame(ax, left=False)
ax.set_xlabel("mean ρ (95% cluster interval)")
C.save("Figure8_GB14")
print("Figure8_GB14 | rows", nr, "| contrast", round(p0.contrast_inflamed_minus_confounder, 3))
