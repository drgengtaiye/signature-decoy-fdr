"""Figure 9 redrawn at the density and in the look of the AML multi-omic figures (pixl_nc.py).
Published prognostic risk models and the pipeline that builds them, against matched decoys
(validation/PLAN_round13_published_risk_models.md with amendments 1-2; icb/v19_audit from applications/v19_*.R and decoys/v19_audit.py).
New panels are descriptive reads of those tables (validation/NOTE_figures_v17_descriptive_only.md)."""
import glob
import json
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
from pixl_nc import *  # noqa: F401,F403
import pixl_nc as X
from pixl_nc import np, pd, plt

A = D + "/icb/v19_audit/"
V = D + "/signature-decoy-fdr/validation/"
U = pd.read_csv(A + "audit_models_units.csv")
U = U[U.n_coex.notna()]
R = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A + "models_[A-Z]*.csv.gz"))])
DEC = R[R.kind == "decoy"]
AM = pd.read_csv(A + "audit_models.csv")
PU = pd.read_csv(A + "audit_pipeline_units.csv")
PA = pd.read_csv(A + "audit_pipeline.csv")
PD = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(A + "pipeline_decoys_*.csv.gz"))])
BYT = pd.read_csv(A + "audit_models_by_type.csv")
BYF = pd.read_csv(A + "audit_models_by_family.csv")
SEARCH = json.load(open(V + "round13_search/search_summary.json"))
PM = pd.read_csv(V + "published_models.csv")
SCR = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(V + "round13_search/screen_*.csv"))])
IND = U[~U.uses_external.astype(bool)]
USED = U[U.uses_external.astype(bool)]
THEME = {"emt": "EMT", "m6a": "m6A"}


def am(pop, q, dec="coex"):
    return AM[AM.population.str.startswith(pop) & (AM.quantity == q) & (AM.decoys == dec)].iloc[0]


pg = Page(240)

# ---------------------------------------------------------------- a  the screen
pg.panel(1, 3, "a", "Search to scored models")
steps = [("found by the search", sum(v["count"] for v in SEARCH.values())), ("screened", len(SCR)),
         ("eligible", int(SCR.eligible.astype(str).str.upper().eq("TRUE").sum())), ("scored", int(PM.scored.sum())),
         ("with decoys (≥ 3 genes)", len(U)), ("independent of the cohort", len(IND))]
ax = pg.ax(25, 8, 16, 31)
y = np.arange(len(steps))[::-1]
ax.plot([s[1] for s in steps], y, color=OTHER, lw=0.5, zorder=2)
point(ax, [s[1] for s in steps], y, INK, ms=3.2)
for yy, s in zip(y, steps):
    ax.text(s[1] * 1.45, yy, f"{s[1]:,}", va="center", ha="left", fontsize=5)
ax.set_xscale("log"); ax.set_xlim(30, 8000); ax.set_xticks([100, 1000], ["100", "1,000"]); ax.minorticks_off()
ax.set_yticks(y, [s[0] for s in steps], fontsize=5); ax.set_ylim(-0.6, len(steps) - 0.4)
ax.set_xlabel("records or models")
pg.src("a", pd.DataFrame(steps, columns=["step", "n"]))

# ---------------------------------------------------------------- b  models by type
pg.panel(44, 3, "b", "Models by type")
g = U.groupby("type")
T = pd.DataFrame({"all": g.size(), "indep.": IND.groupby("type").size(), "valid.": IND.groupby("type").validated.sum()}).fillna(0).astype(int)
x0, y0, cw, ch = 51, 10.4, 6.0, 2.66
for j, c in enumerate(T.columns):
    pg.txt(x0 + cw * j + cw / 2, 8.3, c, ha="center", weight="bold")
for i, (tp, r) in enumerate(T.iterrows()):
    pg.txt(x0 - 0.8, y0 + ch * i + ch / 2, tp, ha="right")
    for j, c in enumerate(T.columns):
        tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, BLUE(r[c] / 12), str(r[c]))
pg.src("b", T.reset_index())

# ---------------------------------------------------------------- c  genes per model
pg.panel(72, 3, "c", "Genes per model")
ax = pg.ax(79, 8, 21, 31)
hist(ax, U.n_genes, np.arange(2.5, 15.5))
zero(ax, U.n_genes.median())
ax.text(U.n_genes.median() + 0.4, ax.get_ylim()[1], f"median {U.n_genes.median():.0f}", va="top", fontsize=5)
ax.text(14.5, ax.get_ylim()[1] * 0.86, f"n = {len(U)}", ha="right", va="top", fontsize=5)
ax.set_xlabel("genes"); ax.set_ylabel("models"); ax.set_xticks([3, 6, 9, 12])
pg.src("c", U[["model", "type", "n_genes"]])

# ---------------------------------------------------------------- d  validated share
pg.panel(103, 3, "d", "Validated")
ax = pg.ax(111, 8, 22, 31)
rT = DEC[(DEC["class"] == "rand") & DEC.model.isin(IND.model)]
vals = {"TCGA": (IND.tcga_validated.mean(), IND.d_tcga.mean(), ((rT.tcga_p < 0.05) & (rT.tcga_worse == 1)).groupby(rT.model).mean().mean()),
        "independent": (IND.validated.mean(), IND.d_coex.mean(), IND.d_rand.mean())}
rows = []
for i, (k, v) in enumerate(vals.items()):
    for j, (val, u) in enumerate(zip(v, ("sig", "coex", "rand"))):
        xx = i + (j - 1) * 0.27
        point(ax, xx, 100 * val, UNIT[u], MARK[u], ms=3.6)
        ax.text(xx, 100 * val + 4, f"{100 * val:.0f}", ha="center", va="bottom", fontsize=5)
        rows.append((k, u, 100 * val))
ax.set_xticks([0, 1], list(vals), fontsize=5.5); ax.set_xlim(-0.55, 1.55); ax.set_ylim(0, 90); ax.set_yticks([0, 25, 50, 75])
ax.set_ylabel("validated (%)")
pg.src("d", pd.DataFrame(rows, columns=["cohort", "unit", "validated_pct"]))

# ---------------------------------------------------------------- e  FDR
pg.panel(136, 3, "e", "False discovery rate")
ax = pg.ax(158, 8, 15.5, 31)
rows = [("independent, co-expressed", am("primary", "fdr"), "coex"), ("independent, random genes", am("primary", "fdr", "rand"), "rand"),
        ("all models", am("all models", "fdr"), "coex"), ("paper used the cohort", am("models whose", "fdr"), "coex")]
out = []
for i, (lab, r, u) in enumerate(rows):
    yy = len(rows) - 1 - i
    ci(ax, 100 * r.lo, 100 * r.hi, yy, UNIT[u], lw=0.7)
    point(ax, 100 * r.estimate, yy, UNIT[u], MARK[u], ms=3.6)
    ax.text(103, yy, f"{100 * r.estimate:.0f}%", va="center", fontsize=5)
    out.append((lab, r.estimate, r.lo, r.hi))
ax.set_yticks(range(len(rows))[::-1], [r[0] for r in rows], fontsize=5)
ax.set_xlim(0, 100); ax.set_xticks([0, 50, 100]); ax.set_xlabel("estimated FDR (%)"); ax.set_ylim(-0.6, len(rows) - 0.4)
pg.src("e", pd.DataFrame(out, columns=["population", "estimate", "lo", "hi"]))

# ---------------------------------------------------------------- f  each model among its decoys
pg.panel(1, 49, "f", "Independent models among their decoys")
ax = pg.ax(9, 54, 49, 110)
o = IND.sort_values("ext_cdir").reset_index(drop=True)
QD = {}
out = []
for yy, r in o.iterrows():
    dc = DEC[(DEC.model == r.model) & (DEC["class"] == "coex")].ext_cdir
    q = np.quantile(dc, [0.05, 0.25, 0.5, 0.75, 0.95])
    QD[r.model] = q[2]
    ax.plot([q[0], q[4]], [yy, yy], color=tint(COEX, 0.45), lw=0.6, solid_capstyle="butt", zorder=2)
    ax.plot([q[1], q[3]], [yy, yy], color=COEX, lw=1.4, solid_capstyle="butt", zorder=3)
    point(ax, r.ext_cdir, yy, SIG, ms=3.0, filled=bool(r.validated), zorder=5)
    pg.txt(8.6, 54 + 110 * (1 - (yy + 1) / (len(o) + 1)), r.type, ha="right")
    out.append((r.model, r.type, r.ext_cdir, r.validated, *q))
zero(ax, 0.5)
ax.set_xlim(0.38, 0.70); ax.set_ylim(-1, len(o)); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xticks([0.4, 0.5, 0.6, 0.7]); ax.set_xlabel("directed C-index, independent cohort")
for yy, (lab, kind) in zip((10, 7, 4), (("validated", "o"), ("not validated", "open-o"), ("decoys, 25–75%, 5–95%", "line"))):
    key(pg, 40, 54 + 110 * (1 - (yy + 1) / (len(o) + 1)), [(COEX if kind == "line" else SIG, lab, kind)])
pg.src("f", pd.DataFrame(out, columns=["model", "type", "cdir", "validated", "d05", "d25", "d50", "d75", "d95"]))

# ---------------------------------------------------------------- g, h  rho
for (lt, xa, S, name) in (("g", 69, IND, "independent"), ("h", 105, USED, "paper used the cohort")):
    pg.panel(xa - 6, 49, lt, f"ρ, {name}")
    ax = pg.ax(xa, 54, 27, 24)
    hist(ax, S.rho, np.linspace(0, 1, 11), weights=100 * np.ones(len(S)) / len(S))
    zero(ax, 10, axis="y")
    ax.set_ylim(0, 52); ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"]); ax.set_yticks([0, 25, 50])
    ax.text(0.03, 50, f"n = {len(S)}\nmean ρ {S.rho.mean():.2f}", va="top", fontsize=5, linespacing=1.15)
    ax.set_xlabel("ρ"); ax.set_ylabel("models (%)")
    pg.src(lt, S[["model", "type", "rho"]])

# ---------------------------------------------------------------- i  model against decoy median
pg.panel(137, 49, "i", "Model and decoy median")
ax = pg.ax(144, 54, 30, 24)
med = IND.model.map(QD)
lo, hi = 0.40, 0.68
ax.plot([lo, hi], [lo, hi], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
ax.set_aspect("auto")
for v, f in ((1, True), (0, False)):
    s = IND.validated == v
    point(ax, med[s], IND.ext_cdir[s], SIG, ms=2.8, filled=f)
ax.set_xlim(0.47, 0.60); ax.set_ylim(lo, hi); ax.set_xticks([0.48, 0.52, 0.56, 0.60]); ax.set_yticks([0.4, 0.5, 0.6])
ax.set_xlabel("decoy median C-index"); ax.set_ylabel("model C-index")
ax.text(0.60, lo + 0.005, f"r = {np.corrcoef(med, IND.ext_cdir)[0, 1]:.2f}, n = {len(IND)}", ha="right", va="bottom", fontsize=5)
pg.src("i", pd.DataFrame({"model": IND.model, "decoy_median": med, "model_cdir": IND.ext_cdir, "validated": IND.validated}))

# ---------------------------------------------------------------- j  TCGA and independent
pg.panel(63, 88, "j", "TCGA and independent")
ax = pg.ax(69, 93, 27, 24)
for v, f in ((1, True), (0, False)):
    s = IND.validated == v
    point(ax, IND.tcga_cdir[s], IND.ext_cdir[s], SIG, ms=2.8, filled=f)
zero(ax, 0.5, axis="y")
ax.set_xlabel("C-index, TCGA"); ax.set_ylabel("C-index, independent"); ax.set_ylim(0.38, 0.70); ax.set_yticks([0.4, 0.5, 0.6, 0.7])
ax.set_xlim(0.5, 0.82); ax.set_xticks([0.5, 0.6, 0.7, 0.8])
ax.text(0.51, 0.695, f"r = {np.corrcoef(IND.tcga_cdir, IND.ext_cdir)[0, 1]:.2f}, n = {len(IND)}", ha="left", va="top", fontsize=5)
pg.src("j", IND[["model", "type", "tcga_cdir", "ext_cdir", "validated"]])

# ---------------------------------------------------------------- k  decoys validated, per model
pg.panel(99, 88, "k", "Decoys validated, by type")
ax = pg.ax(107, 93, 25, 24)
od = IND.groupby("type").d_coex.median().sort_values().index
for i, tp in enumerate(od):
    strip(ax, i, 100 * IND[IND.type == tp].d_coex, COEX, m="s", ms=1.9, spread=0.28, horizontal=True)
ax.set_yticks(range(len(od)), od, fontsize=5); ax.set_ylim(-0.7, len(od) - 0.3)
ax.set_xlim(-3, 100); ax.set_xticks([0, 50, 100]); ax.set_xlabel("co-expressed decoys (%)")
pg.src("k", IND[["model", "type", "d_coex"]])

# ---------------------------------------------------------------- l  audit by type
pg.panel(137, 88, "l", "All models, by type")
cols = [("models", "models", "{:.0f}"), ("validated", "valid.", "{:.0f}"), ("decoy_exp", "exp.", "{:.1f}"), ("rho", "ρ", "{:.2f}"), ("fdr", "FDR", "{:.2f}")]
x0, y0, cw, ch = 145, 95.2, 6.9, 2.66
for j, (c, h, _) in enumerate(cols):
    pg.txt(x0 + cw * j + cw / 2, 93.2, h, ha="center", weight="bold")
for i, r in BYT.iterrows():
    pg.txt(x0 - 0.8, y0 + ch * i + ch / 2, r.type, ha="right")
    for j, (c, _, f) in enumerate(cols):
        fc = {"models": BLUE(r[c] / 20), "validated": BLUE(r[c] / 10), "decoy_exp": PURPLE(r[c] / 4), "rho": BLUE((r[c] - 0.4) / 0.5),
              "fdr": PURPLE(r[c])}[c]
        tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, fc, f.format(r[c]))
pg.src("l", BYT)

# ---------------------------------------------------------------- m  audit by theme family
pg.panel(63, 127, "m", "All models, by theme")
ax = pg.ax(85, 132, 23, 33)
b = BYF.sort_values("models").reset_index(drop=True)
for i, r in b.iterrows():
    ax.plot([r.decoy_exp, r.validated], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, r.decoy_exp, i, COEX, "s", ms=2.6)
    point(ax, r.validated, i, SIG, ms=2.8)
    ax.text(9.8, i, f"{int(r.models)}", va="center", ha="right", fontsize=5)
ax.text(9.8, len(b) - 0.2, "n", va="bottom", ha="right", fontsize=5, weight="bold")
ax.set_yticks(range(len(b)), [f.replace("immune-related", "immune") for f in b.family], fontsize=5)
ax.set_ylim(-0.7, len(b) - 0.2); ax.set_xlim(-0.3, 10); ax.set_xticks([0, 2, 4, 6, 8]); ax.set_xlabel("models")
pg.src("m", BYF)

# ---------------------------------------------------------------- n  pipeline by type
pg.panel(112, 127, "n", "Pipeline, by cancer type")
ax = pg.ax(122, 132, 37, 33)
by = PU.groupby("type").agg(t=("validated", "mean"), d=("d_coex", "mean"), r=("d_rand", "mean")).sort_values("d")
for i, (tp, r) in enumerate(by.iterrows()):
    ax.plot([100 * r.d, 100 * r.t], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, 100 * r.d, i, COEX, "s", ms=2.8)
    point(ax, 100 * r.t, i, SIG, ms=3.0)
    point(ax, 100 * r.r, i, RAND, "D", ms=2.2, zorder=3)
ax.set_yticks(range(len(by)), by.index, fontsize=5); ax.set_ylim(-0.7, len(by) - 0.3)
ax.set_xlim(-2, 100); ax.set_xticks([0, 50, 100]); ax.set_xlabel("pools yielding a validated model (%)")
fb = PA[(PA.population == "theme pools") & (PA.quantity == "fdr") & (PA.decoys == "coex")].iloc[0]
ax.text(100, 3.0, f"FDR {100 * fb.estimate:.0f}%\n({100 * fb.lo:.0f}–{100 * fb.hi:.0f}%)", ha="right", va="center", fontsize=5, linespacing=1.15)
for yy, (c, lab, k) in zip((134, 137, 140), ((SIG, "15 theme", "o"), (COEX, "20 co-expr.", "s"), (RAND, "5 random", "D"))):
    key(pg, 162, yy, [(c, lab, k)])
pg.src("n", by.reset_index().rename(columns={"t": "theme_pools", "d": "coex_decoy_pools", "r": "random_decoy_pools"}))

# ---------------------------------------------------------------- o  pipeline by theme pool
pg.panel(1, 175, "o", "Pipeline, by theme")
ax = pg.ax(19, 180, 33, 49)
bp = PU.groupby("pool").agg(v=("validated", "sum"), e=("d_coex", "sum")).sort_values("e")
for i, (p, r) in enumerate(bp.iterrows()):
    ax.plot([r.e, r.v], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, r.e, i, COEX, "s", ms=3.0)
    point(ax, r.v, i, SIG, ms=3.2)
ax.set_yticks(range(len(bp)), [THEME.get(p, p) for p in bp.index], fontsize=5); ax.set_ylim(-0.7, len(bp) - 0.3)
ax.set_xlim(0, 7); ax.set_xticks([0, 2, 4, 6]); ax.set_xlabel("of 11 cancer types")
pg.src("o", bp.reset_index().rename(columns={"v": "validated_units", "e": "expected_from_decoys"}))

# ---------------------------------------------------------------- p  type by theme
pg.panel(56, 175, "p", "Decoy pools validated, by type and theme")
tt = by.index[::-1]; pp = bp.index[::-1]
x0, y0, cw, ch = 64, 188.0, 5.3, 3.65
for j, p in enumerate(pp):
    pg.txt(x0 + cw * j + cw / 2, y0 - 1.0, THEME.get(p, p), ha="left", va="bottom", rotation=45, rotation_mode="anchor")
M = PU.set_index(["type", "pool"])
for i, tp in enumerate(tt):
    pg.txt(x0 - 0.8, y0 + ch * i + ch / 2, tp, ha="right")
    for j, p in enumerate(pp):
        r = M.loc[(tp, p)]
        tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, BLUE(float(r.d_coex)), f"{100 * r.d_coex:.0f}")
        if r.validated == 1:
            pg.c.plot(x0 + cw * j + cw - 0.75, y0 + ch * i + 0.65, "o", ms=2.0, mfc=SIG, mec="white", mew=0.3, zorder=4)
pg.cbar(90, 231.2, 22, 1.5, BLUE, 0, 100, [0, 50, 100], "co-expressed decoy pools (%)")
pg.c.plot(64.5, 232, "o", ms=2.6, mfc=SIG, mec="white", mew=0.3, zorder=4)
pg.txt(66, 232, "theme model validated")
pg.src("p", PU[["type", "pool", "d_coex", "validated"]])

# ---------------------------------------------------------------- q  pipeline models, TCGA and independent
C_ = PD[(PD["class"] == "coex") & PD.ext_cdir.notna() & PD.tcga_cdir.notna()]
pg.panel(148, 175, "q", "Pipeline models")
ax = pg.ax(155, 180, 23, 23)
ax.hexbin(C_.tcga_cdir, C_.ext_cdir, gridsize=16, cmap=GREYHEX, mincnt=1, linewidths=0, rasterized=True, extent=(0.5, 0.9, 0.3, 0.75))
T_ = PU[PU.ext_cdir.notna()]
for v, f in ((1, True), (0, False)):
    s = T_.validated == v
    point(ax, T_.tcga_cdir[s], T_.ext_cdir[s], SIG, ms=2.2, filled=f, mew=0.4 if not f else 0.2)
zero(ax, 0.5, axis="y")
ax.set_xlim(0.5, 0.9); ax.set_ylim(0.3, 0.75); ax.set_xticks([0.5, 0.7, 0.9]); ax.set_yticks([0.3, 0.5, 0.7])
ax.set_xlabel("C-index, TCGA"); ax.set_ylabel("C-index, ind.")
pg.src("q", pd.concat([C_[["type", "pool", "kind", "tcga_cdir", "ext_cdir", "validated"]],
                       T_.assign(kind="target")[["type", "pool", "kind", "tcga_cdir", "ext_cdir", "validated"]]]))

# ---------------------------------------------------------------- r  pipeline model size
pg.panel(148, 212, "r", "Pipeline model size")
ax = pg.ax(155, 217, 23, 13.5)
bins = np.arange(-0.5, 38.5, 3)
hist(ax, PD[PD["class"] == "coex"].n_model, bins, weights=100 * np.ones((PD["class"] == "coex").sum()) / (PD["class"] == "coex").sum())
h, _ = np.histogram(PU.n_model, bins)
ax.step(bins[:-1], 100 * h / len(PU), where="post", color=SIG, lw=0.8, zorder=4)
ax.set_xlabel("genes in model"); ax.set_ylabel("pools (%)"); ax.set_xlim(-0.5, 37); ax.set_xticks([0, 12, 24, 36])
pg.src("r", pd.concat([PD[PD["class"] == "coex"][["type", "pool", "n_model"]].assign(kind="decoy"), PU[["type", "pool", "n_model"]].assign(kind="target")]))

pg.save("Fig9")
