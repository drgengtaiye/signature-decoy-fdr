"""Figure 7: what a signature's genes add beyond its structure, by endpoint and cancer type (style: pixl_nc.py): gene content and structure
(validation/PLAN_round10_content_and_structure.md; icb/v15_content from decoys/v15_content_structure.py) and the independent
cohorts (validation/PLAN_round11_external_replication.md; icb/v17_external from decoys/v17_external_replication.py).
Panels outside the figure plan are descriptive reads of those tables (validation/NOTE_figures_v17_descriptive_only.md).

rho(s, c) = the share of a signature's 320 co-expressed decoys whose gain over the null in cohort c lies below its own."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt, fm, TextPath

V = D + "/icb/v15_content/"
X_ = D + "/icb/v17_external/"
FAMD = pd.read_csv(D + "/icb/v5_signature_clusters.csv").set_index("signature").family
FAMS = ["Inflamed/IFN", "Other", "Confounder"]
FAMLAB = {"Inflamed/IFN": "inflamed/IFN", "Other": "other", "Confounder": "confounder"}
FTXT = {"Inflamed/IFN": INFL, "Other": MUTED, "Confounder": CONF}       # the grey of the family is too pale for 5 pt words
CELL = {k: pd.read_csv(V + f"cells_{k}.csv.gz") for k in ("response", "tcga_OS")}
SG = {k: pd.read_csv(V + f"signature_{k}.csv", index_col=0) for k in ("response", "tcga_OS")}
FBS = pd.read_csv(V + "family_by_setting.csv")
PC = pd.read_csv(V + "primary_contrast.csv")
MT = pd.read_csv(V + "matched_tumor_types.csv")
AT = pd.read_csv(V + "atlas_by_type.csv")
BT = pd.read_csv(X_ + "by_type.csv")
PCX = pd.read_csv(X_ + "per_cohort.csv")
PRX = pd.read_csv(X_ + "predictions.csv")
PRI = pd.read_csv(X_ + "primary.csv").iloc[0]
SX = pd.read_csv(X_ + "signature_external_survival.csv", index_col=0)
LX = pd.read_csv(X_ + "library_external_survival.csv")
CXS = pd.read_csv(X_ + "cells_external_survival.csv.gz")
CXR = pd.read_csv(X_ + "cells_external_response.csv.gz")
ICB_TUM = {"ICB_Gide": "melanoma", "ICB_Hugo": "melanoma", "ICB_Liu": "melanoma", "ICB_Riaz": "melanoma",
           "ICB_Van_Allen": "melanoma", "ICB_Puch": "melanoma", "ICB_Braun": "kidney", "ICB_Miao1": "kidney",
           "ICB_Jung": "lung", "ICB_Limagne1": "lung", "ICB_Mariathasan": "bladder", "ICB_Kim": "gastric"}
COHORT = {"ICB_Braun": "Braun", "ICB_Gide": "Gide", "ICB_Hugo": "Hugo", "ICB_Jung": "Jung", "ICB_Kim": "Kim",
          "ICB_Limagne1": "Limagne", "ICB_Liu": "Liu", "ICB_Mariathasan": "Mariathasan", "ICB_Miao1": "Miao",
          "ICB_Puch": "Puch", "ICB_Riaz": "Riaz", "ICB_Van_Allen": "Van Allen"}
TAB = {"melanoma": "mel", "lung": "lung", "bladder": "blad", "gastric": "gast", "kidney": "kid"}
TUMS = ["melanoma", "lung", "bladder", "gastric", "kidney"]
SHORT = {"TUMOR-VASCULATURE-UP_Lu": "Tumor vasculature (Lu)", "CD39-CD8Tcell_Chow": "CD39 CD8 (Chow)",
         "Tcell-inflamed_Ayers": "T-cell-inflamed (Ayers)", "WNT-beta-catenin_Takeuchi": "WNT/β-catenin (Takeuchi)"}
RESP, SURV = "o", "^"       # endpoint: circle response, triangle survival; family by colour; filled or open only for a call


def lab(s):
    if s in SHORT:
        return SHORT[s]
    if "_" in s:
        a, b = s.rsplit("_", 1)
        return f"{a.replace('_', ' ')} ({b})"
    return s


def pfmt(p):
    """(sign, value): two significant digits down to 0.001; below, no superscript (p_small: 0.0004, or < 0.0001)"""
    if p >= 0.001:
        return "=", (f"{p:.2g}" if p < 0.01 else f"{p:.2f}")
    return p_small(p)


def width(s, size=5.0):
    return TextPath((0, 0), s, prop=fm.FontProperties(family="Helvetica", size=size)).get_extents().width * 25.4 / 72


pg = Page(240)
WB = {"bbox": {"fc": "white", "ec": "none", "pad": 0.35}}

# ================================================================== a, b  rho over every test, by family
for lt, xa, k, tl in (("a", 8, "response", "All tests, response"), ("b", 29, "tcga_OS", "All tests, OS")):
    pg.panel(1 if lt == "a" else 26.6, 3, lt, tl)
    cc = CELL[k].assign(family=lambda d: d.target.map(FAMD))
    rows = []
    for i, f in enumerate(FAMS):
        ax = pg.ax(xa, 9 + i * 10.4, 15.5, 9.6)
        v = cc[cc.family == f].rho
        hist(ax, v, np.linspace(0, 1, 21), density=True)
        zero(ax, 1, axis="y")
        mf = FBS[(FBS.setting == k) & (FBS.family == f)].rho_bar.iloc[0]
        ax.plot([mf, mf], [0, 5.6], color=FAM[f] if f != "Other" else MUTED, lw=0.6, zorder=3, solid_capstyle="butt")
        ax.set_xlim(0, 1); ax.set_ylim(0, 7.5); ax.set_yticks([0, 5], ["0", "5"] if lt == "a" else ["", ""])
        ax.text(0.02, 7.6, FAMLAB[f], color=FTXT[f], fontsize=5, va="top")
        ax.text(0.02, 5.5, f"n = {len(v)}", fontsize=5, va="top")
        if i == 2:
            ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"]); ax.set_xlabel("ρ")
        else:
            ax.set_xticks([])
        if i == 1 and lt == "a":
            ax.set_ylabel("density")
        h, e = np.histogram(v, np.linspace(0, 1, 21), density=True)
        rows += [{"family": f, "bin_lo": a, "bin_hi": b, "density": d, "tests": len(v), "family_mean_rho": mf} for a, b, d in zip(e[:-1], e[1:], h)]
    pg.src(lt, pd.DataFrame(rows))

# ================================================================== c  response against survival, per signature
pg.panel(47, 3, "c", "Response against survival")
ax = pg.ax(54.5, 9, 27, 31)
r, o = SG["response"], SG["tcga_OS"].reindex(SG["response"].index)
zero(ax, 0.5); zero(ax, 0.5, axis="y")
ax.plot([0.2, 1], [0.2, 1], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
out = []
for f in ["Other", "Confounder", "Inflamed/IFN"]:
    sel = r.index[r.family == f]
    for s in sel:
        cc_ = bool(r.loc[s, "content_call"]) or bool(o.loc[s, "content_call"])
        point(ax, r.loc[s, "rho_bar"], o.loc[s, "rho_bar"], FAM[f], ms=2.8, filled=cc_, zorder=5)
        out.append((s, f, r.loc[s, "rho_bar"], o.loc[s, "rho_bar"], cc_))
    ax.plot(r.loc[sel, "rho_bar"].mean(), o.loc[sel, "rho_bar"].mean(), "+", ms=6.5, mew=1.0, color=FAM[f] if f != "Other" else MUTED, zorder=6)
ax.set_xlim(0.2, 1); ax.set_ylim(0.2, 1); ax.set_xticks([.2, .4, .6, .8, 1], ["0.2", "0.4", "0.6", "0.8", "1.0"]); ax.set_yticks([.2, .4, .6, .8, 1], ["0.2", "0.4", "0.6", "0.8", "1.0"])
ax.set_xlabel("mean ρ, response (12 cohorts)"); ax.set_ylabel("mean ρ, OS (20 types)")
ax.text(0.21, 0.99, "more for survival", ha="left", va="top", fontsize=5, color=MUTED, zorder=2, **WB)
ax.text(0.99, 0.345, "more for response", ha="right", va="center", fontsize=5, color=MUTED, zorder=2, **WB)
p0 = PC[PC.b == "tcga_OS"].iloc[0]
op_, pv_ = pfmt(p0.p_perm)
xr_ = 81.2
line = f"{p0.contrast_inflamed_minus_confounder:.2f} ({p0.lo:.2f}–{p0.hi:.2f}), P {op_} {pv_}"
pg.txt(xr_, 36.7, "inflamed − confounder", ha="right", **WB)
pg.txt(xr_, 38.7, line, ha="right", **WB)
pg.c.plot(56.2, 14.3, "o", ms=2.6, mfc=INK, mec="white", mew=0.3, zorder=7); pg.txt(57.4, 14.3, "content call", **WB)
pg.c.plot(56.2, 16.7, "o", ms=2.4, mfc="white", mec=INK, mew=0.6, zorder=7); pg.txt(57.4, 16.7, "none")
pg.c.plot(56.2, 19.1, "+", ms=4.5, mew=0.8, color=INK, zorder=7); pg.txt(57.4, 19.1, "family mean", **WB)
pg.src("c", pd.DataFrame(out, columns=["signature", "family", "mean_rho_response", "mean_rho_OS", "content_call_either"])
       .assign(inflamed_minus_confounder=p0.contrast_inflamed_minus_confounder, lo=p0.lo, hi=p0.hi, P_permutation=p0.p_perm))

# ================================================================== d  family means, per setting
pg.panel(84.5, 3, "d", "Family means")
AY, AH = 9.6, 30.4
ax = pg.ax(100.5, AY, 22, AH)
ax.set_ylim(AY + AH, AY)
SETS = [("response", "response"), ("tcga_OS", "OS"), ("tcga_DSS", "DSS"), ("tcga_PFI", "PFI"), ("tcga_OS_matched", "OS, same\ntumor types")]
rp = AH / (15 + 4 * 0.45)
out = []
yy = AY + rp / 2
for k, lb in SETS:
    y0_ = yy
    for f in FAMS:
        q = FBS[(FBS.setting == k) & (FBS.family == f)].iloc[0]
        ci(ax, q.lo, q.hi, yy, FAM[f], lw=0.7)
        point(ax, q.rho_bar, yy, FAM[f], ms=2.6, zorder=5)
        pg.txt(123.5, yy, f"{q.rho_bar:.2f} ({q.lo:.2f}–{q.hi:.2f})")
        out.append((k, f, q.n, q.rho_bar, q.lo, q.hi))
        yy += rp
    pg.txt(99.3, y0_ + rp, lb, ha="right", linespacing=1.0)
    yy += 0.45 * rp
zero(ax, 0.5)
ax.set_xlim(0.15, 1.0); ax.set_xticks([.25, .5, .75, 1], ["0.25", "0.50", "0.75", "1.00"]); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("mean ρ")
pg.txt(130.3, 7.6, "ρ (95% interval)", ha="center", weight="bold")
pg.src("d", pd.DataFrame(out, columns=["setting", "family", "signatures", "mean_rho", "lo", "hi"]))
key(pg, 86.5, 6.9, [(INFL, "inflamed", "o"), (OTH, "other", "o")], gap=1.2)
key(pg, 86.5, 9.3, [(CONF, "confounder", "o")])

# ================================================================== e  rho among random-gene decoys
pg.panel(140, 3, "e", "ρ, two decoy classes")
ax = pg.ax(148.5, 9, 28.5, 31)
ax.plot([0.2, 1], [0.2, 1], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
out = []
for k, m in (("tcga_OS", SURV), ("response", RESP)):
    s_ = SG[k]
    for f in FAMS:
        x = s_[s_.family == f]
        point(ax, x.rho_bar, x.rho_bar_rand, FAM[f], m=m, ms=2.4 if m == SURV else 2.6)
        out.append(x.assign(setting=k).reset_index()[["target", "family", "setting", "rho_bar", "rho_bar_rand"]])
ax.set_xlim(0.2, 1); ax.set_ylim(0.2, 1.02); ax.set_xticks([.2, .6, 1], ["0.2", "0.6", "1.0"]); ax.set_yticks([.2, .6, 1], ["0.2", "0.6", "1.0"])
ax.set_xlabel("mean ρ, co-expressed decoys"); ax.set_ylabel("mean ρ, random-gene decoys")
rr = {k: np.corrcoef(SG[k].rho_bar, SG[k].rho_bar_rand)[0, 1] for k in SG}
ax.text(0.99, 0.31, f"response, r = {rr['response']:.2f}", ha="right", va="center", fontsize=5)
ax.text(0.99, 0.25, f"OS, r = {rr['tcga_OS']:.2f}; n = 52", ha="right", va="center", fontsize=5)
key(pg, 149.5, 11, [(INK, "response", RESP)]); key(pg, 149.5, 13.4, [(INK, "OS", SURV)])
pg.src("e", pd.concat(out))

# ================================================================== the aligned block: columns and rows shared by f, g, h and j, k, l
r_cols = AT[AT.setting == "response"].assign(tum=lambda d: d.cohort.map(ICB_TUM))
r_cols = r_cols.sort_values(["tum", "rho_inflamed"], key=lambda s: s.map(TUMS.index) if s.name == "tum" else -s)
o_cols = AT[AT.setting == "tcga_OS"].sort_values("rho_inflamed", ascending=False)
i_cols = PCX.sort_values("tcga_rho_inflamed", ascending=False)
Cs = SG["response"].rho_bar - SG["tcga_OS"].rho_bar.reindex(SG["response"].index)
rows = []
for f in FAMS:
    rows += sorted([s for s in Cs.index if FAMD[s] == f], key=lambda s: -Cs[s])
MR = CELL["response"].pivot_table(index="target", columns="cohort", values="rho").reindex(index=rows)[list(r_cols.cohort)]
MO = CELL["tcga_OS"].pivot_table(index="target", columns="cohort", values="rho").reindex(index=rows)[list(o_cols.cohort)]
MI = CXS.pivot_table(index="target", columns="cohort", values="rho").reindex(index=rows)[list(i_cols.cohort)]
MI["GSE176307"] = CXR.set_index("target").rho.reindex(rows)
CW, RH, FGAP = 2.45, 1.88, 0.6
XR = 24.8; XO = XR + 12 * CW + 1.6; XI = XO + 20 * CW + 1.6
XI_ = {c: XI + j * CW for j, c in enumerate(i_cols.cohort)}
XI_["GSE176307"] = XI + 11 * CW + 0.6
YH = 75.5
ry = {}
y_ = YH
for i, s in enumerate(rows):
    if i and FAMD[s] != FAMD[rows[i - 1]]:
        y_ += FGAP
    ry[s] = y_
    y_ += RH
YE = y_

# ================================================================== f, g, h  calls and expected false per cohort
YT, YB = 56.5, 68.5
YMAX = 36
FDRL = lambda q: "–" if not q.calls else (f"{q.fdr:.0f}" if q.fdr < 1000 else ">999")  # noqa: E731
pg.panel(1, 47.5, "f", "Calls, response")
pg.panel(XO, 47.5, "g", "Calls, overall survival")
pg.panel(XI, 47.5, "h", "Calls, independent")
pg.txt(XR - 1.2, 53.8, "FDR (%)", ha="right", weight="bold")
for lt, x0, tab, cols, nc in (("f", XR, r_cols, "cohort", 12), ("g", XO, o_cols, "cohort", 20), ("h", XI, i_cols, "cohort", 11)):
    wd = (nc * CW) if lt != "h" else (XI_["GSE176307"] + CW - XI)
    ax = pg.ax(x0, YT, wd, YB - YT)
    ax.set_xlim(x0, x0 + wd); ax.set_ylim(-1.3, YMAX); ax.set_xticks([])
    ax.set_yticks([0, 20], ["0", "20"] if lt == "f" else ["", ""])
    if lt == "f":
        ax.set_ylabel("calls")
    out = []
    for j, (_, q) in enumerate(tab.iterrows()):
        xc = x0 + (j + 0.5) * CW
        ax.plot([xc, xc], [q.exp_false, q.calls], color=OTHER, lw=0.5, zorder=2, solid_capstyle="butt")
        point(ax, xc, q.exp_false, COEX, "s", ms=1.9, zorder=4)
        point(ax, xc, q.calls, SIG, "o", ms=2.1, zorder=5)
        pg.txt(xc, 53.8, FDRL(q), rotation=90, ha="center")
        out.append((q[cols], q.n, q.calls, q.exp_false, q.fdr))
    if lt == "h":
        xc = XI_["GSE176307"] + CW / 2
        c_, e_ = int(CXR.clear.sum()), CXR.dclear.sum()
        ax.plot([xc, xc], [e_, c_], color=OTHER, lw=0.5, zorder=2, solid_capstyle="butt")
        point(ax, xc, e_, COEX, "s", ms=1.9, zorder=4)
        point(ax, xc, c_, SIG, "o", ms=2.1, zorder=5)
        out.append(("GSE176307 (response)", int(CXR.n.iloc[0]), c_, e_, np.nan))
    pg.src(lt, pd.DataFrame(out, columns=["cohort", "patients", "calls", "expected_false", "fdr_pct"]))
key(pg, 1, 58, [(SIG, "calls", "o")])
key(pg, 1, 60.6, [(COEX, "expected false", "s")])

# ================================================================== i  content calls by family and setting
pg.panel(138, 47.5, "i", "Content calls")
SI = [("response", "resp."), ("tcga_OS", "OS"), ("tcga_DSS", "DSS"), ("tcga_PFI", "PFI"), ("tcga_OS_matched", "OS, 6"), ("indep", "ind.")]
x0, y0, cw, ch = 148.2, 54.8, 5.0, 2.85
for j, (_, h) in enumerate(SI):
    pg.txt(x0 + cw * j + cw / 2, 53.1, h, ha="center", weight="bold")
out = []
for i, f in enumerate(FAMS + ["all"]):
    pg.txt(x0 - 0.6, y0 + ch * i + ch / 2, FAMLAB.get(f, "all").replace("/IFN", ""), ha="right")
    for j, (k, _) in enumerate(SI):
        q = LX[LX.family == f].iloc[0] if k == "indep" else FBS[(FBS.setting == k) & (FBS.family == f)].iloc[0]
        tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, BLUE(q.content_calls / q.n), f"{int(q.content_calls)}/{int(q.n)}")
        out.append((k, f, int(q.content_calls), int(q.n), int(q.calls)))
pg.src("i", pd.DataFrame(out, columns=["setting", "family", "content_calls", "signatures", "signature_calls"]))

# ================================================================== j, k, l  rho for every signature in every cohort
pg.panel(1, 71.5, "j", "ρ, response (12 cohorts)")
pg.panel(XO, 71.5, "k", "ρ, overall survival (20 cancer types)")
pg.panel(XI, 71.5, "l", "ρ, independent cohorts")
for s in rows:
    pg.txt(23, ry[s] + (RH - 0.15) / 2, lab(s), ha="right")
    pg.box(23.25, ry[s], 1.05, RH, FAM[FAMD[s]], z=2)
for M_, xpos in ((MR, {c: XR + j * CW for j, c in enumerate(MR.columns)}), (MO, {c: XO + j * CW for j, c in enumerate(MO.columns)}), (MI, XI_)):
    for c in M_.columns:
        for s in rows:
            v = M_.loc[s, c]
            if np.isfinite(v):
                tile(pg, xpos[c], ry[s], CW, RH, PUOR(float(v)))
for j, (_, q) in enumerate(r_cols.iterrows()):
    pg.txt(XR + (j + 0.5) * CW, YE + 0.4, f"{COHORT[q.cohort]} ({TAB[ICB_TUM[q.cohort]]})", rotation=90, ha="center", va="top")
for j, (_, q) in enumerate(o_cols.iterrows()):
    pg.txt(XO + (j + 0.5) * CW, YE + 0.4, q.cohort, rotation=90, ha="center", va="top")
for c in i_cols.cohort:
    pg.txt(XI_[c] + CW / 2, YE + 0.4, c.replace("_", " "), rotation=90, ha="center", va="top")
pg.txt(XI_["GSE176307"] + CW / 2, YE + 0.4, "GSE176307, response", rotation=90, ha="center", va="top")
for lt, M_ in (("j", MR), ("k", MO), ("l", MI)):
    pg.src(lt, M_.rename_axis("signature").reset_index().assign(family=lambda d: d.signature.map(FAMD)))

# ================================================================== m  mean rho and content calls
XM = 138.4; MW = 5.2
pg.panel(XM, 71.5, "m", "Mean ρ")
ME = pd.DataFrame({"response": SG["response"].rho_bar, "OS": SG["tcga_OS"].rho_bar, "independent": SX.rho_bar}).reindex(rows)
CC = pd.DataFrame({"response": SG["response"].content_call, "OS": SG["tcga_OS"].content_call, "independent": SX.content_call}).reindex(rows)
for j, h in enumerate(("resp.", "OS", "ind.")):
    pg.txt(XM + MW * j + MW / 2, 74.4, h, ha="center", weight="bold")
for s in rows:
    for j, c in enumerate(ME.columns):
        v = ME.loc[s, c]
        if np.isfinite(v):
            fc = PUOR(float(v))
            tile(pg, XM + MW * j, ry[s], MW, RH, fc)
            pg.txt(XM + MW * j + 0.32, ry[s] + (RH - 0.15) / 2, f"{v:.2f}", color=ink_on(fc), zorder=3)
            if CC.loc[s, c] is True or CC.loc[s, c] == True:  # noqa: E712
                pg.c.plot(XM + MW * j + MW - 0.68, ry[s] + (RH - 0.15) / 2, "o", ms=1.7, mfc=SIG, mec="white", mew=0.25, zorder=4)
pg.src("m", ME.add_prefix("mean_rho_").join(CC.add_prefix("content_call_")).rename_axis("signature").reset_index())

# ================================================================== n  content and structure, aligned with the rows
XN = 155.6
pg.panel(XN, 71.5, "n", "Content, structure")
out = []
for k, xa, xl, tk, h in (("response", 156.4, (-0.02, 0.11), [0, 0.1], "response"), ("tcga_OS", 168.2, (-0.02, 0.04), [0, 0.03], "OS")):
    ax = pg.ax(xa, YH, 9.8, YE - YH)
    ax.set_ylim(YE, YH); ax.set_xlim(*xl); ax.set_yticks([]); ax.set_xticks(tk, ["0", f"{tk[1]:g}"])
    ax.spines["left"].set_visible(False)
    zero(ax, 0)
    q = SG[k]
    for s in rows:
        yc = ry[s] + (RH - 0.15) / 2
        point(ax, q.loc[s, "med_S"], yc, COEX, "s", ms=1.6, zorder=3)
        point(ax, q.loc[s, "med_G"], yc, SIG, "o", ms=1.9, filled=bool(q.loc[s, "content_call"]), zorder=4, mew=0.45 if not q.loc[s, "content_call"] else None)
    pg.txt(xa + 4.9, 74.4, h, ha="center", weight="bold")
    out.append(q.loc[rows, ["family", "med_S", "med_G", "content_call"]].assign(setting=k).rename_axis("signature").reset_index())
pg.src("n", pd.concat(out))
pg.txt(167.3, YE + 5.0, "median gain", size=6, ha="center")
pg.cbar(XM, YE + 3.0, 14, 1.4, PUOR, 0, 1, [0, 0.5, 1], "ρ", labels=["0", "0.5", "1"])
pg.c.plot(XM + 0.6, YE + 11.0, "o", ms=2.0, mfc=SIG, mec="white", mew=0.25, zorder=4)
pg.txt(XM + 1.6, YE + 11.0, "content call")
key(pg, XN + 0.8, YE + 8.5, [(SIG, "content", "o"), (COEX, "structure", "s")], gap=1.2)
for yk, (c_, lb_) in zip((YE + 3, YE + 5.6, YE + 8.2), ((INFL, "inflamed/IFN"), (OTH, "other"), (CONF, "confounder"))):
    key(pg, 8.5, yk, [(c_, lb_)])

# ================================================================== o  same tumor types
YR = 193.0; YA = 199.5; HA = 34.0
pg.panel(1, YR, "o", "Same tumor types")
ax = pg.ax(12.5, YA, 18.5, HA)
out = []
for i, t_ in enumerate(TUMS):
    for jf, f in enumerate(["Inflamed/IFN", "Confounder"]):
        y = i + (jf - 0.5) * 0.34
        a_ = MT[(MT.tumor == t_) & (MT.family == f) & (MT.setting == "response")].rho.iloc[0]
        b_ = MT[(MT.tumor == t_) & (MT.family == f) & (MT.setting == "tcga_OS")].rho.iloc[0]
        ax.plot([a_, b_], [y, y], color=FAM[f], lw=0.6, zorder=3, solid_capstyle="butt")
        point(ax, a_, y, FAM[f], RESP, ms=2.8, zorder=5)
        point(ax, b_, y, FAM[f], SURV, ms=2.8, zorder=5)
        out.append((t_, f, a_, b_))
zero(ax, 0.5)
ax.set_yticks(range(len(TUMS)), TUMS, fontsize=5); ax.set_ylim(len(TUMS) - 0.4, -0.6); ax.tick_params(axis="y", length=0)
ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1], ["0.0", "0.5", "1.0"]); ax.spines["left"].set_visible(False); ax.set_xlabel("mean ρ of the family")
key(pg, 9, YA - 2.8, [(INK, "response", RESP), (INK, "OS", SURV)], gap=1.4)
pg.src("o", pd.DataFrame(out, columns=["tumor", "family", "mean_rho_response", "mean_rho_OS"]))


# ================================================================== p, q  independent cohorts against TCGA
def types_scatter(ax, x, y, names, c, off):
    ax.plot([0, 1.05], [0, 1.05], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
    zero(ax, 0.5); zero(ax, 0.5, axis="y")
    for a, b, n in zip(x, y, names):
        point(ax, a, b, c, SURV, ms=3.0, zorder=5)
        dx, dy, ha = off.get(n, (0.035, 0, "left"))
        ax.text(a + dx, b + dy, n, fontsize=5, va="center", ha=ha, zorder=2, **WB)
    ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05); ax.set_xticks([0, .5, 1], ["0.0", "0.5", "1.0"]); ax.set_yticks([0, .5, 1], ["0.0", "0.5", "1.0"])


pg.panel(33.5, YR, "p", "Independent cohorts")
ax = pg.ax(40.5, YA - 2.5, 22.5, HA + 2.5)
types_scatter(ax, BT.tcga, BT.ext, BT.type, INFL, {"CESC": (0.035, 0.03, "left"), "SARC": (0.035, -0.035, "left"), "HNSC": (-0.035, 0.03, "right"), "LUAD": (-0.035, 0, "right"),
                                                   "SKCM": (-0.035, 0.0, "right"), "BRCA": (0.035, 0.035, "left")})
ax.set_xlabel("inflamed mean ρ, TCGA"); ax.set_ylabel("inflamed mean ρ, independent")
op_, pb = pfmt(PRI.p_one_sided)
ax.text(1.04, 0.015, f"Spearman {PRI.spearman:.2f},".replace("-", "−") + f"\nP {op_} {pb}", ha="right", va="bottom", fontsize=5, linespacing=1.05, zorder=2, **WB)
pg.src("p", BT.assign(spearman=PRI.spearman, p_one_sided=PRI.p_one_sided))

pg.panel(65.5, YR, "q", "Confounder family")
ax = pg.ax(72.5, YA - 2.5, 22.5, HA + 2.5)
sc = PCX.rho_confounder.corr(PCX.tcga_rho_confounder, method="spearman")
types_scatter(ax, PCX.tcga_rho_confounder, PCX.rho_confounder, PCX.type, CONF,
              {"SKCM": (0.04, 0.02, "left"), "STAD": (0.035, -0.035, "left"), "LUAD": (0.03, 0.02, "left"), "HNSC": (-0.035, 0, "right"), "BRCA": (-0.035, 0, "right")})
ax.set_xlabel("confounder mean ρ, TCGA"); ax.set_ylabel("confounder mean ρ, independent")
ax.text(1.04, 0.015, f"Spearman {sc:.2f}, n = {len(PCX)}", ha="right", va="bottom", fontsize=5, zorder=2, **WB)
pg.src("q", PCX[["cohort", "type", "tcga_rho_confounder", "rho_confounder"]].assign(spearman=sc))

# ================================================================== r  each independent cohort, and the response cohort
pg.panel(97.5, YR, "r", "Per independent cohort")
ax = pg.ax(119.5, YA, 22.5, HA)
PCS = PCX.sort_values("tcga_rho_inflamed")
ys = np.arange(len(PCS))
out = []
XNC = 146.8
pg.txt(XNC, YA - 1.2, "n", ha="right", weight="bold")
ymap = lambda v: YA + HA * (1 - (v + 0.6) / (len(PCS) + 0.6 + 0.7))  # noqa: E731
for y, (_, q) in zip(ys, PCS.iterrows()):
    for f_, lb, off in ((INFL, "inflamed", 0.17), (CONF, "confounder", -0.17)):
        ci(ax, q[f"rho_{lb}_lo"], q[f"rho_{lb}_hi"], y + off, f_, lw=0.6)
        point(ax, q[f"rho_{lb}"], y + off, f_, SURV, ms=2.2, zorder=5)
        out.append((q.cohort, "survival", lb, q[f"rho_{lb}"], q[f"rho_{lb}_lo"], q[f"rho_{lb}_hi"], q.n))
    pg.txt(XNC, ymap(y), f"{int(q.n):,}", ha="right")
yr = len(PCS) + 0.25
for kk, f_, off, lb in (("GSE176307 inflamed", INFL, 0.17, "inflamed"), ("GSE176307 confounder", CONF, -0.17, "confounder")):
    q = PRX[PRX.test.str.startswith(kk)].iloc[0]
    ci(ax, q.lo, q.hi, yr + off, f_, lw=0.6)
    point(ax, q.estimate, yr + off, f_, RESP, ms=2.4, zorder=5)
    out.append(("GSE176307", "response", lb, q.estimate, q.lo, q.hi, int(CXR.n.iloc[0])))
pg.txt(XNC, ymap(yr), f"{int(CXR.n.iloc[0])}", ha="right")
zero(ax, 0.5)
ax.set_yticks(list(ys) + [yr], [f"{t} {g}" for t, g in zip(PCS.type, PCS.cohort.str.split("_").str[1])] + ["BLCA GSE176307, response"], fontsize=5)
ax.tick_params(axis="y", length=0)
ax.set_ylim(-0.6, len(PCS) + 0.7); ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1], ["0.0", "0.5", "1.0"]); ax.spines["left"].set_visible(False)
ax.set_xlabel("mean ρ (95% cluster interval)")
key(pg, 119, YA - 2.8, [(INFL, "inflamed", "o"), (CONF, "confounder", "o")], gap=1.0, sw=1.4)
pg.src("r", pd.DataFrame(out, columns=["cohort", "endpoint", "family", "mean_rho", "lo", "hi", "patients"]))

# ================================================================== s  family means per cohort
pg.panel(148.6, YR, "s", "ρ per cohort")
ax = pg.ax(158, YA, 19.2, HA)
GR = [("response", AT[AT.setting == "response"], RESP), ("OS", AT[AT.setting == "tcga_OS"], SURV), ("ind.", PCX, SURV)]
yy = 0
out, tk = [], []
for (lb, t_, m) in GR:
    tk.append(yy + 1)
    for f, col in zip(FAMS, ("rho_inflamed", "rho_other", "rho_confounder")):
        strip(ax, yy, t_[col], FAM[f], m=m, ms=1.8, spread=0.3, horizontal=True)
        out += [(lb, f, c_, v_) for c_, v_ in zip(t_["cohort"], t_[col])]
        yy += 1
    yy += 0.6
zero(ax, 0.5)
ax.set_ylim(yy - 0.6 - 0.4, -0.6); ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1], ["0.0", "0.5", "1.0"]); ax.spines["left"].set_visible(False)
ax.set_yticks(tk, [g[0] for g in GR], fontsize=5); ax.tick_params(axis="y", length=0)
ax.set_xlabel("family mean ρ")
pg.src("s", pd.DataFrame(out, columns=["setting", "family", "cohort", "family_mean_rho"]))

pg.save("Figure7")
