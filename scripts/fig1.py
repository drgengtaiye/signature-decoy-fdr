"""Figure 1: a matched null for 52 published signatures, and the target-decoy procedure built on it (style: pixl_nc.py).
Panels a-d are the procedure drawn with real values: Chemokine (Messina) and its first decoy set (a-c), the whole library (d).
Tables: v13_fig/ (figprep_example.py, figprep_roc.py, figprep_r.R); calib_v3/ICB_*_v3.csv, calib_v3_units.csv, calib_v3_fdr.csv and v12_pkg/response_fdr.csv hold the unit-level
values of the first decoy set and the ten-set estimate (calib_score_v3.R, calib_v3_fdr.py, applications/v11_*.R).
Panels outside the figure plan are descriptive reads of those tables (validation/NOTE_figures_v17_descriptive_only.md)."""
import glob
import json
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt
from scipy.stats import spearmanr

I = D + "/icb/"
R = pd.read_csv(I + "native_response_matched_v2.csv"); R = R[R.cohort != "ICB_VanDenEnde"].copy()
R["delta"] = R.AUC - R.rand_med
O = pd.read_csv(I + "native_response.csv"); O = O[O.cohort != "ICB_VanDenEnde"]
FIX = O.groupby("cohort")[["rand_med", "rand_p95"]].first()
famd = pd.read_csv(I + "v3_dissociation.csv").set_index("signature").family
RB = pd.read_csv(I + "v13_fig/r_by_cohort_response.csv")
EX = pd.read_csv(I + "v13_fig/example_units.csv")
ROC = pd.read_csv(I + "v13_fig/roc_example.csv")
UN = pd.read_csv(I + "calib_v3_units.csv")
UN = UN[(UN["null"] == "N0") & (UN.stat == "fold")]
TEN = pd.read_csv(I + "v12_pkg/response_fdr.csv")
TEN = TEN[(TEN.stat == "fold") & (TEN.decoys == "coex") & (TEN.boot == "clusters")].iloc[0]
F3 = pd.read_csv(I + "calib_v3_fdr.csv")
F3 = F3[(F3["null"] == "N0") & (F3.stat == "fold") & (F3.pool == "A") & (F3.decoys == "coex")].iloc[0]
T1, T2, CO = "Chemokine_Messina", "CIN70_Carter", "ICB_Gide"
NAME = {T1: "Chemokine (Messina)", T2: "CIN70 (Carter)"}
COH = {"ICB_Braun": "Braun", "ICB_Gide": "Gide", "ICB_Hugo": "Hugo", "ICB_Jung": "Jung", "ICB_Kim": "Kim", "ICB_Limagne1": "Limagne",
       "ICB_Liu": "Liu", "ICB_Mariathasan": "Mariathasan", "ICB_Miao1": "Miao", "ICB_Puch": "Puch", "ICB_Riaz": "Riaz",
       "ICB_Van_Allen": "Van Allen"}
TUM = {"ICB_Gide": "melanoma", "ICB_Hugo": "melanoma", "ICB_Liu": "melanoma", "ICB_Riaz": "melanoma", "ICB_Van_Allen": "melanoma",
       "ICB_Puch": "melanoma", "ICB_Braun": "kidney", "ICB_Miao1": "kidney", "ICB_Jung": "lung", "ICB_Limagne1": "lung",
       "ICB_Mariathasan": "bladder", "ICB_Kim": "gastric"}
ENG = {"gsva": "GSVA", "wmean": "weighted mean", "ssgsea": "ssGSEA"}
ENGC = {"gsva": INK, "wmean": "#B5367A", "ssgsea": "#7A8B1F"}
FIXC = "#8C8C8C"       # the fixed 48-gene null, grey against the matched null in ink
assert check_palette("engine", {ENG[k]: v for k, v in ENGC.items()})
assert check_palette("null", {"fixed": FIXC, "matched": INK})


def swarm(v, width, step):
    """offsets of a one-dimensional beeswarm: each point, in order, takes the free level nearest the centre; levels are 1/2 step apart"""
    v = np.asarray(v, float)
    off = np.zeros(len(v))
    levels = sorted(np.linspace(-width, width, 7), key=abs)
    placed = []
    for k in np.argsort(v):
        for lv in levels:
            if all(abs(v[k] - a) >= step or abs(lv - b) > 1e-9 for a, b in placed):
                off[k] = lv
                break
        else:
            off[k] = levels[len(placed) % len(levels)]
        placed.append((v[k], off[k]))
    return off


def rows_axis(ax, n, grid=True):
    """aligned cohort rows, row 0 at the top; hairlines carry the eye between panels"""
    ax.set_ylim(n - 0.5, -0.5); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    if grid:
        for i in range(n):
            ax.axhline(i, color="#E4E4E4", lw=0.3, zorder=0)


def flow(x0, x1, y=3.1):
    """the thin line that carries the procedure from one panel to the next"""
    pg.c.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="-|>", lw=0.4, color=MUTED, mutation_scale=4, shrinkA=0, shrinkB=0))


pg = Page(240)

# ================================================================ a-d  the procedure, with real values
# ---------------------------------------------------------------- a  signature and its decoys: internal r per cohort
cols = ["cohort", "unit", "kind", "target", "pool", "structure", "null", "r_int"]
CX = pd.concat([pd.read_csv(f, usecols=cols) for f in sorted(glob.glob(I + "calib_v3/ICB_*_v3.csv"))])
CX = CX[(CX["null"] == "N0") & (CX.target == T1) & (CX.pool.isin(["-", "A"]))]
M = CX.pivot(index="unit", columns="cohort", values="r_int")
co_a = M.loc[T1].sort_values(ascending=False).index
pg.panel(1, 3, "a", "Chemokine (Messina) and its decoys")
x0, y0, cw = 4.4, 15.8, 3.0
yy = y0
for j, co in enumerate(co_a):
    pg.txt(x0 + cw * j + (cw - 0.15) / 2, y0 - 0.7, COH[co], ha="center", va="bottom", rotation=90)
order_a = []
for name, c, lab, rh in (("signature", SIG, ("signature",), 2.2), ("coex", COEX, ("32 co-expressed", "decoys"), 0.62),
                         ("rand", RAND, ("8 random-gene", "decoys"), 0.62)):
    un = [T1] if name == "signature" else list(M.loc[M.index.str.contains(f"_A_{name}_")].mean(axis=1).sort_values(ascending=False).index)
    assert len(un) == {"signature": 1, "coex": 32, "rand": 8}[name]
    order_a += un
    h_ = rh * len(un)
    pg.box(x0 - 1.5, yy, 0.9, h_ - 0.12, c, z=2)
    for i, u in enumerate(un):
        for j, co in enumerate(co_a):
            pg.box(x0 + cw * j, yy + rh * i, cw - 0.15, rh - 0.12, BLUE(float(np.clip(M.loc[u, co] / 0.7, 0, 1))), z=2)
    for k_, s_ in enumerate(lab):
        pg.txt(x0 + cw * 12 + 0.6, yy + h_ / 2 + (k_ - (len(lab) - 1) / 2) * 2.2, s_, color=c)
    yy += h_ + 0.7
pg.cbar(x0, yy + 0.3, 14, 1.3, BLUE, 0, 0.7, [0, 0.35, 0.7], "internal mean pairwise r")
pg.src("a", M.loc[order_a, co_a].reset_index().replace(COH).rename(columns=COH))

# ---------------------------------------------------------------- b  read against the matched null, in one cohort
g_ = EX[(EX.target == T1) & (EX.cohort == CO)]
bank, dec = g_[g_.kind == "bank"].fold.to_numpy(), g_[g_.kind == "decoy"].fold.to_numpy()
tv, med_b = g_[g_.kind == "target"].fold.iloc[0], float(np.median(bank))
pg.panel(65, 3, "b", "Null bank, Gide")
flow(45.8, 64.0)
flow(88.5, 99.0)
flow(131.5, 136.0)
ax = pg.ax(72.5, 9, 25.5, 33)
hb = hist(ax, bank, np.arange(0.5, 0.875, 0.0125))
top = hb[0].max()
ax.set_ylim(0, top * 1.46); ax.set_xlim(0.5, 0.9); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9]); ax.set_yticks([0, 50, 100])
ax.plot([med_b, med_b], [0, top * 1.34], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=3)
ax.plot([tv, tv], [0, top * 1.34], color=SIG, lw=0.8, zorder=4, solid_capstyle="butt")
ax.vlines(dec, top * 1.05, top * 1.17, color=COEX, lw=0.5, zorder=4)
ax.annotate("", xy=(tv, top * 1.3), xytext=(med_b, top * 1.3), arrowprops=dict(arrowstyle="<|-|>", lw=0.4, color=INK, mutation_scale=3.2, shrinkA=0, shrinkB=0), zorder=6)
ax.text((tv + med_b) / 2, top * 1.33, f"Δ = {tv - med_b:+.2f}".replace("-", "−"), ha="center", va="bottom", fontsize=5)
ax.text(0.812, top * 1.11, "32\ndecoys", color=COEX, ha="left", va="center", fontsize=5, linespacing=1.0)
ax.text(0.813, top * 0.5, "signature", color=SIG, ha="left", va="center", fontsize=5, rotation=90)
ax.text(0.646, top * 0.62, "1,000\nrandom sets", color=MUTED, ha="left", va="center", fontsize=5, linespacing=1.05)
ax.set_xlabel("folded AUC"); ax.set_ylabel("random sets")
pg.src("b", g_[["cohort", "target", "kind", "unit", "fold"]].assign(cohort="Gide", bank_median=med_b, delta=tv - med_b))

# ---------------------------------------------------------------- c  called over cohorts: the same test for every unit
pg.panel(100, 3, "c", "Calls over 12 cohorts")
ax = pg.ax(108, 9, 27, 33)
uc = UN[(UN.target == T1) & (UN.pool.isin(["-", "A"]))]
for st, c, m, ms in (("coex", COEX, "s", 2.4), ("rand", RAND, "D", 2.2)):
    s_ = uc[uc.structure == st]
    for f in (False, True):
        k_ = s_[s_.call == f]
        point(ax, k_.med_delta, -np.log10(k_.p), c, m, ms=ms, filled=f, mew=0.4 if not f else 0.2)
ts = uc[uc.structure == "signature"].iloc[0]
point(ax, ts.med_delta, -np.log10(ts.p), SIG, "o", ms=4.0, zorder=6)
zero(ax, -np.log10(0.05), axis="y")
ax.plot([0, 0], [0, 3.15], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
ax.text(0.125, -np.log10(0.05) + 0.04, "P < 0.05", fontsize=5, va="bottom", ha="right")
ax.text(ts.med_delta + 0.002, -np.log10(ts.p) + 0.1, f"P = {ts.p:.3f}", color=SIG, fontsize=5, ha="center", va="bottom")
ax.set_xlim(-0.06, 0.125); ax.set_ylim(0, 4.0); ax.set_xticks([-0.05, 0, 0.05, 0.1], ["−0.05", "0", "0.05", "0.10"]); ax.set_yticks([0, 1, 2, 3])
ax.set_xlabel("median Δ over cohorts"); ax.set_ylabel("−log10 signed-rank P")
n_c, n_r = int(uc[uc.structure == "coex"].call.sum()), int(uc[uc.structure == "rand"].call.sum())
key(pg, 109.5, 10.2, [(COEX, f"co-expressed, {n_c}/32 called", "s")])
key(pg, 109.5, 12.8, [(RAND, f"random-gene, {n_r}/8 called", "D")])
pg.src("c", uc[["unit", "kind", "structure", "ncoh", "med_delta", "p", "call"]])

# ---------------------------------------------------------------- d  the estimate over the library
pg.panel(137, 3, "d", "Estimated FDR")
ax = pg.ax(145, 9, 34, 33)
tg = UN[UN.kind == "target"].set_index("target")
da = UN[(UN.kind == "decoy") & (UN.pool == "A") & (UN.structure == "coex")]
ps = da.groupby("target").call.mean()
assert da.groupby("target").size().eq(32).all() and len(ps) == 52
RR = int(tg.call.sum()); SP = float(ps.sum())
assert RR == int(F3.called) and abs(SP - F3.exp_false) < 1e-9 and RR == int(TEN.called)
dd = pd.DataFrame({"p_s": ps, "called": tg.call}).sort_values(["p_s", "called"], ascending=False).reset_index().rename(columns={"index": "target"})
dd["rank"] = np.arange(1, 53)
for f in (False, True):
    k_ = dd[dd.called == f]
    point(ax, k_["rank"], 100 * k_.p_s, SIG, ms=2.6, filled=f, mew=0.5 if not f else 0.2)
ax.set_xlim(-1, 53.5); ax.set_ylim(-1.5, 36); ax.set_xticks([1, 13, 26, 39, 52]); ax.set_yticks([0, 10, 20, 30])
r1, r2 = dd[dd.target == T1].iloc[0], dd[dd.target == T2].iloc[0]
ax.annotate(NAME[T1], xy=(r1["rank"], 100 * r1.p_s), xytext=(r1["rank"] + 3.2, 100 * r1.p_s + 0.2), fontsize=5, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.3, color=MUTED, shrinkA=0, shrinkB=1.2))
ax.annotate(NAME[T2], xy=(r2["rank"], 0), xytext=(52.8, 5.6), fontsize=5, ha="right", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.3, color=MUTED, shrinkA=0, shrinkB=1.2))
ax.set_xlabel("signatures, ranked"); ax.set_ylabel("decoys called (%)")
cx = (149.8, 167.6, 172.6, 178.9)
yh = 10.4
for x_, h_ in zip(cx[1:], ("exp. false", "R", "FDR")):
    pg.txt(x_, yh, h_, ha="right", weight="bold")
ROWS_D = (("first set", SP, RR, 100 * SP / RR), ("10 sets", TEN.exp_false, int(TEN.called), TEN.fdr),
          (f"10 sets × π0 {TEN.pi0:.2f}", TEN.exp_false * TEN.pi0, int(TEN.called), TEN.fdr_pi0))
for i, (lab, e_, r_, f_) in enumerate(ROWS_D):
    y_ = yh + 2.5 * (i + 1)
    pg.txt(cx[0], y_, lab)
    pg.txt(cx[1], y_, f"{e_:.1f}", ha="right"); pg.txt(cx[2], y_, f"{r_}", ha="right"); pg.txt(cx[3], y_, f"{f_:.0f}%", ha="right")
key(pg, 149.8, 20.9, [(SIG, f"called, R = {RR}", "o"), (SIG, "not called", "open-o")], gap=1.8)
pg.src("d", dd.assign(exp_false_first_set=SP, R=RR, fdr_first_set=100 * SP / RR, exp_false_ten_sets=TEN.exp_false, fdr_ten_sets=TEN.fdr,
                      pi0_ten_sets=TEN.pi0, fdr_pi0_ten_sets=TEN.fdr_pi0))

# ================================================================ e-h  the library and its null
SIZ = R.drop_duplicates("signature")[["signature", "n_genes", "engine"]].assign(family=lambda d: d.signature.map(famd))
MEDSZ = float(SIZ.n_genes.median())
LADDER = sorted(R.null_size.unique())
pg.panel(1, 53, "e", "52 signatures, three engines")
ye, he = 58.0, 25.0
ax = pg.ax(19, ye, 40, he)
mm_e = lambda v: ye + he * (v + 1.05) / 3.95  # noqa: E731
for i, e in enumerate(("gsva", "wmean", "ssgsea")):
    g = SIZ[SIZ.engine == e].sort_values("n_genes")
    off = swarm(np.log10(g.n_genes.to_numpy(float)), 0.33, 0.042)
    for f in ("Other", "Confounder", "Inflamed/IFN"):
        s_ = (g.family == f).to_numpy()
        point(ax, g.n_genes[s_], i + off[s_], FAM[f], ms=2.7)
    pg.txt(18.4, mm_e(i), f"{ENG[e]} ({len(g)})", ha="right")
zero(ax, MEDSZ)
ax.text(MEDSZ * 1.08, 2.66, f"median {MEDSZ:g}", fontsize=5, va="center")
ax.vlines(LADDER, -1.0, -0.66, color=INK, lw=0.5)
pg.txt(18.4, mm_e(-0.83), "null sizes", ha="right")
ax.set_xscale("log"); ax.set_xlim(2.7, 560); ax.set_ylim(2.9, -1.05); ax.minorticks_off()
ax.set_xticks([3, 10, 30, 100, 300], ["3", "10", "30", "100", "300"]); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("genes in signature")
for j, (f, lab) in enumerate((("Inflamed/IFN", "inflamed/IFN"), ("Other", "other"), ("Confounder", "confounder"))):
    key(pg, 45.4, 74.4 + 2.5 * j, [(FAM[f], lab, "o")])
pg.src("e", SIZ.assign(null_sizes=", ".join(map(str, LADDER))))

# ---------------------------------------------------------------- f  family and engine
pg.panel(61, 53, "f", "Family and engine")
tb = pd.crosstab(SIZ.engine, SIZ.family).loc[["gsva", "wmean", "ssgsea"], ["Inflamed/IFN", "Confounder", "Other"]]
tb["all"] = tb.sum(axis=1); tb.loc["all"] = tb.sum()
x0, y0, cw, ch = 74.6, 63.0, 5.4, 4.45
for j, (c, h_) in enumerate(zip(tb.columns, ("infl.", "conf.", "other", "all"))):
    pg.txt(x0 + cw * j + cw / 2, y0 - 2.6, h_, ha="center", weight="bold")
    if c in FAM:
        pg.box(x0 + cw * j, y0 - 1.3, cw - 0.15, 0.7, FAM[c], z=2)
for i, (e, r_) in enumerate(tb.iterrows()):
    pg.txt(x0 - 0.8, y0 + ch * i + ch / 2, ENG.get(e, "all"), ha="right")
    for j, c in enumerate(tb.columns):
        tile(pg, x0 + cw * j, y0 + ch * i, cw, ch, BLUE(min(r_[c] / 34, 1)), str(int(r_[c])))
pg.src("f", tb.reset_index().replace(ENG))

# ---------------------------------------------------------------- g  null by size and engine
pg.panel(100, 53, "g", "Null 95th percentile, by size")
ax = pg.ax(108, 58, 33, 24)
Q = (R.groupby(["engine", "null_size"]).agg(p95=("ladder_p95", "median"), n=("ladder_p95", "size")).reset_index())
Qk = Q[Q.n >= 10]
for e in ("gsva", "ssgsea", "wmean"):
    g = Qk[Qk.engine == e].sort_values("null_size")
    ax.plot(g.null_size, g.p95, color=ENGC[e], lw=0.7, zorder=3)
    point(ax, g.null_size, g.p95, ENGC[e], ms=2.6)
old95 = R.cohort.map(FIX.rand_p95).median()
zero(ax, old95, axis="y", color=MUTED)
ax.text(520, old95 + 0.003, "fixed 48-gene null", fontsize=5, va="bottom", ha="right")
ax.set_xscale("log"); ax.set_xlim(3.2, 560); ax.minorticks_off(); ax.set_xticks([4, 10, 40, 100, 400], ["4", "10", "40", "100", "400"])
ax.set_ylim(.58, .77); ax.set_yticks([.6, .65, .7, .75], ["0.60", "0.65", "0.70", "0.75"])
ax.set_xlabel("random gene-set size"); ax.set_ylabel("null 95th pct. (AUC)")
for e, yy_ in (("gsva", 0.767), ("ssgsea", 0.752), ("wmean", 0.737)):
    ax.text(520, yy_, ENG[e], color=ENGC[e], fontsize=5, ha="right", va="center")
pg.src("g", Q.assign(drawn=Q.n >= 10, fixed_48_gene_median_p95=old95))

# ---------------------------------------------------------------- h  matched against fixed
pg.panel(143, 53, "h", "Matched and fixed null")
ax = pg.ax(151, 58, 28, 24)
m = R.assign(delta_fixed=R.AUC - R.cohort.map(FIX.rand_med))
s = m.groupby("signature")[["delta", "delta_fixed"]].median(); s["family"] = s.index.map(famd)
lo_, hi_ = -.025, .13
ax.plot([lo_, hi_], [lo_, hi_], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
for f in ("Other", "Confounder", "Inflamed/IFN"):
    g = s[s.family == f]
    point(ax, g.delta_fixed, g.delta, FAM[f], ms=2.5)
ax.set_xlim(lo_, hi_); ax.set_ylim(lo_, hi_); ax.set_xticks([0, .05, .1], ["0", "0.05", "0.10"]); ax.set_yticks([0, .05, .1], ["0", "0.05", "0.10"])
ax.set_xlabel("median ΔAUC, fixed null"); ax.set_ylabel("median ΔAUC, matched")
rr = spearmanr(s.delta_fixed, s.delta)[0]
ax.text(lo_ + 0.006, hi_, f"Spearman r = {rr:.2f}, n = {len(s)}", fontsize=5, va="top")
pg.src("h", s.reset_index())

# ================================================================ i, j
# ---------------------------------------------------------------- i  co-expression per cohort
RH = 2.75
pg.panel(1, 93, "i", "Internal correlation, by cohort")
tt = RB.drop_duplicates(["target", "cohort"])
oi = tt.groupby("cohort").r_target.median().sort_values(ascending=False).index
ax = pg.ax(14, 98, 42, RH * 12)
rows_axis(ax, 12, grid=False)
outi = []
for i, co in enumerate(oi):
    pg.txt(13.4, 98 + RH * (i + 0.5), COH[co], ha="right")
    for dy, v, c, mk, ms_, nm in ((-0.28, tt[tt.cohort == co].r_target, SIG, "o", 2.2, "signature"),
                                  (0.0, RB[(RB.cohort == co) & (RB["class"] == "coex")].r_decoy, COEX, "s", 1.9, "coex"),
                                  (0.28, RB[(RB.cohort == co) & (RB["class"] == "rand")].r_decoy, RAND, "D", 1.7, "rand")):
        v = v.dropna()      # CRMA (Shukla) is unscored in Miao
        q = np.percentile(v, [5, 25, 50, 75, 95])
        ci(ax, q[0], q[4], i + dy, tint(c, 0.55), lw=0.4)
        ci(ax, q[1], q[3], i + dy, c, lw=0.9)
        point(ax, q[2], i + dy, c, mk, ms=ms_, mew=0.25, zorder=5)
        outi.append((COH[co], nm, len(v), *q))
zero(ax, 0)
ax.set_xlim(-0.12, 0.9); ax.set_xticks([0, .25, .5, .75], ["0", "0.25", "0.50", "0.75"]); ax.set_xlabel("mean pairwise r")
for j, (c, lab, k_) in enumerate(((SIG, "signatures", "o"), (COEX, "co-expressed", "s"), (RAND, "random-gene", "D"))):
    key(pg, 46.0, 122.7 + 2.5 * j, [(c, lab, k_)])
pg.src("i", pd.DataFrame(outi, columns=["cohort", "unit", "n", "p05", "p25", "median", "p75", "p95"]))

# ---------------------------------------------------------------- j  null 95th percentile, every rung in every cohort
pg.panel(59, 93, "j", "Null 95th percentile, by cohort and size")
coh_n = R.groupby("cohort").n.first()
oh = coh_n.sort_values(ascending=False).index
LP = R.groupby(["engine", "null_size", "cohort"]).ladder_p95.first()
x0, y0, cw, ch = 72.2, 104.4, 4.46, 2.72
scale = lambda v: BLUE(float(np.clip((v - 0.55) / 0.23, 0, 1)))  # noqa: E731
for i, co in enumerate(oh):
    pg.txt(x0 - 0.8, y0 + ch * i + ch / 2, COH[co], ha="right")
xx = x0
outj = []
for e in ("gsva", "wmean", "ssgsea"):
    rungs = sorted(R[R.engine == e].null_size.unique())
    pg.txt(xx + cw * len(rungs) / 2, y0 - 4.6, ENG[e], ha="center", weight="bold", color=ENGC[e])
    for j, k in enumerate(rungs):
        pg.txt(xx + cw * j + cw / 2, y0 - 1.5, str(k), ha="center")
        for i, co in enumerate(oh):
            if (e, k, co) in LP.index:
                v = LP.loc[(e, k, co)]
                tile(pg, xx + cw * j, y0 + ch * i, cw, ch, scale(v), f"{v:.2f}")
                outj.append((ENG[e], k, COH[co], v))
            else:
                pg.box(xx + cw * j, y0 + ch * i, cw - 0.15, ch - 0.15, EMPTY, z=2)
    xx += cw * len(rungs) + 1.1
pg.txt(xx + 2.6, y0 - 4.6, "fixed", ha="center", weight="bold"); pg.txt(xx + 2.6, y0 - 1.5, "48", ha="center")
for i, co in enumerate(oh):
    v = FIX.rand_p95[co]
    tile(pg, xx, y0 + ch * i, 5.2, ch, scale(v), f"{v:.2f}")
    outj.append(("fixed", 48, COH[co], v))
pg.cbar(153, 92.2, 16, 1.3, BLUE, 0.55, 0.78, [0.55, 0.65, 0.75], "AUC")
pg.src("j", pd.DataFrame(outj, columns=["engine", "null_size", "cohort", "null_p95"]))

# ================================================================ k-o  every cohort, aligned rows
EE = {}
for T in (T1, T2):
    e = EX[EX.target == T]
    rows = {}
    for c_, g in e.groupby("cohort"):
        b = g[g.kind == "bank"].fold.to_numpy()
        rows[c_] = dict(fold=g[g.kind == "target"].fold.iloc[0], bank=np.percentile(b, [5, 25, 50, 75, 95]),
                        dec=g[g.kind == "decoy"].set_index("unit").fold)
    EE[T] = rows
ROWS = sorted(EE[T1], key=lambda c_: -(EE[T1][c_]["fold"] - EE[T1][c_]["bank"][2]))
yt, rh = 148.6, 3.2
Hh = rh * 12

# ---------------------------------------------------------------- k  the cohorts
pg.panel(1, 140, "k", "Cohorts")
ST_ = I + "v10_store"
n1 = {co: int(sum(json.load(open(f"{ST_}/{co}_meta.json"))["y"])) for co in R.cohort.unique()}
nn = R.groupby("cohort").n.first()
assert int(nn.sum()) == 827
x0 = 13.6
pg.txt(x0 + 1.6, yt - 1.4, "tumor", weight="bold")
for x_, h_ in ((x0 + 15.0, "n"), (x0 + 20.8, "resp.")):
    pg.txt(x_, yt - 1.4, h_, ha="right", weight="bold")
for i, c_ in enumerate(ROWS):
    y_ = yt + rh * i
    pg.txt(x0 - 0.8, y_ + rh / 2, COH[c_], ha="right", weight="bold" if c_ == CO else None)
    pg.box(x0, y_ + 0.08, 1.0, rh - 0.16, TYPE[TUM[c_]], z=2)
    pg.txt(x0 + 1.6, y_ + rh / 2, TUM[c_])
    pg.txt(x0 + 15.0, y_ + rh / 2, str(int(nn[c_])), ha="right")
    pg.txt(x0 + 20.8, y_ + rh / 2, str(n1[c_]), ha="right")
pg.src("k", pd.DataFrame({"cohort": [COH[c_] for c_ in ROWS], "tumor": [TUM[c_] for c_ in ROWS], "patients": [int(nn[c_]) for c_ in ROWS],
                          "responders": [n1[c_] for c_ in ROWS]}))


# ---------------------------------------------------------------- l, m  folded AUC against bank and decoys
def bank_rows(letter, T, x_, w, title, called):
    pg.panel(x_ - 1.5, 140, letter, title)
    ax = pg.ax(x_, yt, w, Hh)
    rows_axis(ax, 12)
    out = []
    for i, c_ in enumerate(ROWS):
        r_ = EE[T][c_]
        q = r_["bank"]
        ci(ax, q[0], q[4], i, NULLC, lw=0.6, zorder=2)
        ax.add_patch(plt.Rectangle((q[1], i - 0.3), q[3] - q[1], 0.6, fc=NULLC, ec="none", zorder=2))
        ax.plot([q[2], q[2]], [i - 0.3, i + 0.3], color=INK, lw=0.7, zorder=3, solid_capstyle="butt")
        ax.vlines(r_["dec"].to_numpy(), i - 0.3, i + 0.3, color=COEX, lw=0.3, zorder=4)
        point(ax, r_["fold"], i, SIG, ms=3.0, filled=called, zorder=6)
        out.append((COH[c_], r_["fold"], *q))
    ax.set_xlim(0.5, 1.0); ax.set_xticks([0.5, 0.75, 1.0], ["0.5", "0.75", "1"]); ax.set_xlabel("folded AUC")
    pg.src(letter, pd.DataFrame(out, columns=["cohort", "signature_fold_auc", "bank_p05", "bank_p25", "bank_p50", "bank_p75", "bank_p95"]).merge(
        EX[(EX.target == T) & (EX.kind == "decoy")].assign(cohort=lambda d: d.cohort.map(COH)).pivot(index="cohort", columns="unit", values="fold").reset_index(), on="cohort"))


call_t = {T: bool(UN[(UN.kind == "target") & (UN.target == T)].call.iloc[0]) for T in (T1, T2)}
bank_rows("l", T1, 37, 33.5, "Chemokine (Messina), AUC", call_t[T1])
bank_rows("m", T2, 74, 32, "CIN70 (Carter), AUC", call_t[T2])
yk = yt - 2.2
xk = key(pg, 37.2, yk, [(NULLC, "null bank", None)], gap=1.6)
pg.line(xk + 0.9, yk - 0.8, xk + 0.9, yk + 0.8, color=COEX, lw=0.5); pg.txt(xk + 2.4, yk, "decoys")
key(pg, xk + 9.0, yk, [(SIG, "called", "o")], gap=1.2)
pg.line(74.9, yk - 0.8, 74.9, yk + 0.8, color=INK, lw=0.7); pg.txt(76.4, yk, "null median")
key(pg, 88.0, yk, [(SIG, "not called", "open-o")])


# ---------------------------------------------------------------- n, o  delta over cohorts
def delta_rows(letter, T, x_, w, title, called):
    pg.panel(x_ - 1.5, 140, letter, title)
    ax = pg.ax(x_, yt, w, Hh)
    rows_axis(ax, 12)
    u_ = UN[(UN.target == T) & (UN.pool == "A") & (UN.structure == "coex")].set_index("unit").call
    out = []
    for i, c_ in enumerate(ROWS):
        r_ = EE[T][c_]
        md = r_["bank"][2]
        d_ = r_["dec"] - md
        off = swarm(d_.to_numpy(), 0.40, 0.0072)
        f_ = u_.reindex(d_.index).to_numpy(bool)
        point(ax, d_.to_numpy()[~f_], i + off[~f_], COEX, "s", ms=1.45, filled=False, mew=0.3)
        point(ax, d_.to_numpy()[f_], i + off[f_], COEX, "s", ms=1.7, zorder=4.5, mew=0.1)
        tv_ = r_["fold"] - md
        point(ax, tv_, i, SIG, ms=3.4, filled=called, zorder=6)
        out.append(dict(cohort=COH[c_], signature_delta=tv_, **d_.to_dict()))
    zero(ax, 0)
    ax.set_xlim(-0.13, 0.31); ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3], ["−0.1", "0", "0.1", "0.2", "0.3"]); ax.set_xlabel("Δ folded AUC")
    ts_ = UN[(UN.kind == "target") & (UN.target == T)].iloc[0]
    pg.txt(x_, yt - 2.2, f"P = {ts_.p:.3f}, decoys called {int(u_.sum())}/{len(u_)}")
    pg.src(letter, pd.DataFrame(out))


delta_rows("n", T1, 109, 33, "Chemokine (Messina), Δ", call_t[T1])
delta_rows("o", T2, 146, 33, "CIN70 (Carter), Δ", call_t[T2])

# ================================================================ p-s
# ---------------------------------------------------------------- p  ROC
pg.panel(1, 197, "p", "ROC in Gide")
ax = pg.ax(9, 203, 28, 29)
dcols = [c for c in ROC.columns if c.startswith("decoy")]
dq = np.percentile(ROC[dcols].to_numpy(float), [5, 50, 95], axis=1)
for c_, lw_ in (("bank_p05", 0.4), ("bank_p95", 0.4), ("bank_p50", 0.7)):
    ax.step(ROC.fpr, ROC[c_], where="post", color=FIXC, lw=lw_, zorder=2)
for k_, lw_ in ((0, 0.4), (2, 0.4), (1, 0.7)):
    ax.step(ROC.fpr, dq[k_], where="post", color=COEX, lw=lw_, zorder=3)
ax.step(ROC.fpr, ROC.target, where="post", color=SIG, lw=1.0, zorder=4)
ax.plot([0, 1], [0, 1], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_xticks([0, .5, 1], ["0", "0.5", "1"]); ax.set_yticks([0, .5, 1], ["0", "0.5", "1"])
ax.set_xlabel("false-positive rate"); ax.set_ylabel("true-positive rate")
auc_t = np.trapezoid(ROC.target, ROC.fpr)
ax.text(1.0, 0.03, f"signature, AUC {auc_t:.2f}", color=SIG, fontsize=5, ha="right", va="bottom")
ax.text(1.0, 0.105, "decoys 5–95%", color=COEX, fontsize=5, ha="right", va="bottom")
ax.text(1.0, 0.18, "null 5–95%", color=MUTED, fontsize=5, ha="right", va="bottom")
pg.src("p", ROC.assign(decoy_p05=dq[0], decoy_p50=dq[1], decoy_p95=dq[2])[["fpr", "target", "decoy_p05", "decoy_p50", "decoy_p95", "bank_p05", "bank_p50", "bank_p95"]])

# ---------------------------------------------------------------- q  null width against cohort size
pg.panel(40, 197, "q", "Null SD, by cohort size")
ax = pg.ax(48, 203, 38, 29)
for co in nn.index:
    v = R[R.cohort == co].rand_sd.to_numpy()
    ax.plot(np.full(len(v), nn[co]) * np.exp(swarm(v, 0.03, 0.0011)), v, "o", ms=1.3, mfc=NULLC, mec="white", mew=0.1, ls="none", zorder=2)
    ax.plot([nn[co] * 0.955, nn[co] * 1.047], [np.median(v)] * 2, color=INK, lw=0.8, zorder=3, solid_capstyle="butt")
fr = np.median([n1[c_] / nn[c_] for c_ in nn.index])
xx_ = np.linspace(24, 250, 100); hn = np.sqrt(1 - 2 / np.pi)
ax.plot(xx_, hn * np.sqrt((xx_ + 1) / (12 * xx_ * fr * xx_ * (1 - fr))), color=INK, lw=0.7, zorder=4)
ax.text(255, 0.06, "Mann–Whitney,\nfolded", fontsize=5, ha="right", va="center", linespacing=1.05)
ax.set_xscale("log"); ax.set_xlim(22, 270); ax.minorticks_off(); ax.set_xticks([30, 60, 120, 240], ["30", "60", "120", "240"])
ax.set_ylim(0, 0.09); ax.set_yticks([0, 0.02, 0.04, 0.06, 0.08], ["0", "0.02", "0.04", "0.06", "0.08"])
ax.set_xlabel("patients in cohort"); ax.set_ylabel("null SD (AUC)")
pg.src("q", R[["cohort", "n", "signature", "engine", "rand_sd"]].assign(cohort=R.cohort.map(COH), median_responder_share=fr))

# ---------------------------------------------------------------- r  tests clearing each null, per cohort
pg.panel(89, 197, "r", "Tests clearing the null")
mm = R.assign(fix95=R.cohort.map(FIX.rand_p95))
mm["cf"], mm["cm"] = mm.AUC > mm.fix95, mm.AUC > mm.rand_p95
pc = mm.groupby("cohort").agg(tests=("cf", "size"), fixed=("cf", "mean"), matched=("cm", "mean"), n_fixed=("cf", "sum"), n_matched=("cm", "sum"))
pc[["fixed", "matched"]] *= 100
pc = pc.sort_values("matched", ascending=False)
allr = pd.Series(dict(tests=len(mm), fixed=100 * mm.cf.mean(), matched=100 * mm.cm.mean(), n_fixed=mm.cf.sum(), n_matched=mm.cm.sum()), name="all tests")
assert f"{allr.fixed:.1f}" == "38.2" and f"{allr.matched:.1f}" == "23.0"
T_ = pd.concat([allr.to_frame().T, pc])
ax = pg.ax(102, 203, 33, 29)
rows_axis(ax, 13, grid=False)
xa, xb = 141.2, 151.2
pg.txt(xa, 201.2, "fixed", ha="right", weight="bold"); pg.txt(xb, 201.2, "matched", ha="right", weight="bold")
rr_h = 29 / 13
for i, (c_, r_) in enumerate(T_.iterrows()):
    a_, b_ = r_.matched, r_.fixed
    ax.plot([a_, b_], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, b_, i, FIXC, "o", ms=2.7, zorder=3)
    point(ax, a_, i, INK, "o", ms=2.7, zorder=4)
    y_ = 203 + rr_h * (i + 0.5)
    pg.txt(101.4, y_, COH.get(c_, c_), ha="right", weight="bold" if i == 0 else None)
    pg.txt(xa, y_, f"{b_:.1f}", ha="right", weight="bold" if i == 0 else None); pg.txt(xb, y_, f"{a_:.1f}", ha="right", weight="bold" if i == 0 else None)
ax.set_xlim(-1.5, 75); ax.set_xticks([0, 25, 50, 75]); ax.set_xlabel("tests above null 95th pct. (%)")
key(pg, 120.6, 203 + rr_h * 11.5, [(FIXC, "fixed 48-gene", "o")])
key(pg, 120.6, 203 + rr_h * 12.5, [(INK, "matched", "o")])
pg.src("r", T_.reset_index().rename(columns={"index": "cohort"}).replace(COH))

# ---------------------------------------------------------------- s  the same tests by rung of the ladder
pg.panel(153.5, 197, "s", "Tests, by set size")
bs = mm.groupby("null_size").agg(tests=("cf", "size"), fixed=("cf", "mean"), matched=("cm", "mean")).reset_index()
x0, y0, cw, ch = 161.3, 203.0, 5.9, 2.9
pg.txt(x0 - 0.8, y0 - 1.8, "size", ha="right", weight="bold")
for j, h_ in enumerate(("tests", "fixed", "match.")):
    pg.txt(x0 + cw * j + cw / 2, y0 - 1.8, h_, ha="center", weight="bold")
for i, r_ in bs.iterrows():
    pg.txt(x0 - 0.8, y0 + ch * i + ch / 2, str(int(r_.null_size)), ha="right")
    pg.txt(x0 + cw / 2, y0 + ch * i + ch / 2, str(int(r_.tests)), ha="center")
    tile(pg, x0 + cw, y0 + ch * i, cw, ch, BLUE(float(r_.fixed) / 0.7), f"{100 * r_.fixed:.0f}")
    tile(pg, x0 + 2 * cw, y0 + ch * i, cw, ch, BLUE(float(r_.matched) / 0.7), f"{100 * r_.matched:.0f}")
pg.src("s", bs.assign(fixed=100 * bs.fixed, matched=100 * bs.matched))

pg.save("Figure1")
print("Fig1 | fixed", round(allr.fixed, 1), "matched", round(allr.matched, 1), "fixed higher in", int((pc.fixed > pc.matched).sum()),
      "| ROC AUC", round(auc_t, 3), "| rho", round(rr, 2), "| median size", MEDSZ, "| Gide target", round(tv, 3), "median", round(med_b, 3),
      "| first set", SP, RR, round(100 * SP / RR, 2), "| ten sets", TEN.exp_false, TEN.called, round(TEN.fdr, 2), "pi0", round(TEN.pi0, 3), round(TEN.fdr_pi0, 2))
