"""Figure 2: immune signatures carry information about response to checkpoint blockade (style: pixl_nc.py).
Panels outside the figure plan are descriptive reads of the tables below (validation/NOTE_figures_v17_descriptive_only.md).

Tables: icb/native_response_matched_v2.csv (v10_apply_signed.py), icb/v3_dissociation.csv (p_resp, family, n_genes),
icb/v4_permutation_inference.csv and icb/v4_permutation_null_draws.csv.gz (perm_inference.py, signed nulls of validation/PLAN_calibration_round5.md),
icb/cf_cohort_stats.csv (cf1_cohort_stats.R; patients and responders)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pixl_nc import *  # noqa: F401,F403
from pixl_nc import np, pd, plt
from matplotlib.colors import Normalize
from scipy.stats import binom
from statsmodels.stats.proportion import proportion_confint

I = D + "/icb/"
R = pd.read_csv(I + "native_response_matched_v2.csv")
R = R[R.cohort != "ICB_VanDenEnde"].copy()
R["delta"] = R.AUC - R.rand_med
R["clear"] = R.AUC > R.rand_p95
X = pd.read_csv(I + "v3_dissociation.csv").set_index("signature")
R["family"] = R.signature.map(X.family)
TUMS = ["melanoma", "kidney", "lung", "bladder", "gastric"]
ICB_TUM = {"ICB_Gide": "melanoma", "ICB_Hugo": "melanoma", "ICB_Liu": "melanoma", "ICB_Riaz": "melanoma", "ICB_Van_Allen": "melanoma",
           "ICB_Puch": "melanoma", "ICB_Braun": "kidney", "ICB_Miao1": "kidney", "ICB_Jung": "lung", "ICB_Limagne1": "lung",
           "ICB_Mariathasan": "bladder", "ICB_Kim": "gastric"}  # figstyle.ICB_TUM
R["tum"] = R.cohort.map(ICB_TUM)
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]
FAMLAB = {"Inflamed/IFN": "inflamed/IFN", "Other": "other", "Confounder": "confounder"}
famd = X.family
sigp = X.p_resp
med = R.groupby("signature").delta.median()
nsig = int((sigp < 0.05).sum())
assert nsig == 20 and len(R) == 623 and (med[sigp < 0.05] > 0).all()
assert np.allclose(med, X.d_resp.reindex(med.index))
P = pd.read_csv(I + "v4_permutation_inference.csv")
DR = pd.read_csv(I + "v4_permutation_null_draws.csv.gz")
HEAD = "12 cohorts (headline)"
ph = P[(P.analysis == HEAD) & (P.statistic == "signatures with P<0.05")].iloc[0]
pc_ = P[(P.analysis == HEAD) & (P.statistic == "inflamed-confounder clearance (pp)")].iloc[0]
dv = DR[(DR.analysis == HEAD) & (DR.statistic == "signatures with P<0.05")].value.to_numpy()
dc = DR[(DR.analysis == HEAD) & (DR.statistic == "inflamed-confounder clearance (pp)")].value.to_numpy()
assert len(dv) == 2000 and int(ph.observed) == nsig
CS = pd.read_csv(I + "cf_cohort_stats.csv").set_index("cohort")

SHORT = {"CD39-CD8Tcell_Chow": "CD39 CD8 (Chow)", "TUMOR-VASCULATURE-UP_Lu": "Tumor vasculature (Lu)", "Tcell-inflamed_Ayers": "T-cell-inflamed (Ayers)",
         "WNT-beta-catenin_Takeuchi": "WNT/β-catenin (Takeuchi)"}  # cns.SHORT
GENES = {}   # signature names are labels, roman on every sheet


def lab(s, italic=True):
    if italic and s in GENES:
        return GENES[s]
    if s in SHORT:
        return SHORT[s]
    a, b = s.rsplit("_", 1)
    return f"{a.replace('_', ' ')} ({b})"


def coh(c):
    return {"ICB_Limagne1": "Limagne", "ICB_Miao1": "Miao"}.get(c, c.replace("ICB_", "").replace("_", " "))


def swarm(vals, step):
    """tied values side by side, centred (for shares that take few values)"""
    vals = np.asarray(vals, float)
    off = np.zeros(len(vals))
    for u in np.unique(np.round(vals, 6)):
        k = np.where(np.round(vals, 6) == u)[0]
        off[k] = (np.arange(len(k)) - (len(k) - 1) / 2) * step
    return off


def qjit(vals, width, nbins=14):
    """deterministic spread for a strip of many points: low-discrepancy offsets, scaled by the local density"""
    vals = np.asarray(vals, float)
    if len(vals) < 2:
        return np.zeros(len(vals))
    h, e = np.histogram(vals, bins=nbins)
    dens = h[np.clip(np.digitize(vals, e) - 1, 0, nbins - 1)] / h.max()

    def vdc(i):
        x, f = 0.0, 0.5
        while i:
            x += f * (i & 1); i >>= 1; f /= 2
        return x
    off = np.zeros(len(vals))
    off[np.argsort(vals)] = [2 * vdc(i + 1) - 1 for i in range(len(vals))]
    return off * width * np.sqrt(dens)


def no_x(ax):
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)


pg = Page(240)

# ================================================================= a  the matrix
rows = []
for f in FAMORD:
    rows += sorted([s for s in med.index if famd[s] == f], key=lambda s: -med[s])
ncoh = R.groupby("cohort").n.first()
cols = sorted(ncoh.index, key=lambda c: (TUMS.index(ICB_TUM[c]), -ncoh[c]))
M = R.pivot_table(index="signature", columns="cohort", values="delta").loc[rows][cols]
CL = R.pivot_table(index="signature", columns="cohort", values="clear").loc[rows][cols]
LIM = 0.15
norm = Normalize(-LIM, LIM)
M0, MT, cw, rh, gap = 27.5, 15.2, 4.5, 2.04, 0.9
MW = cw * len(cols)
ypos, y_ = {}, MT
span = {}
for f in FAMORD:
    fs_ = [s for s in rows if famd[s] == f]
    span[f] = (y_, y_ + rh * len(fs_))
    for s in fs_:
        ypos[s] = y_
        y_ += rh
    y_ += gap
MB = y_ - gap
pg.panel(1, 3, "a", "Gain over the matched null")
for f in FAMORD:
    y0, y1 = span[f]
    pg.box(2.75, y0, 0.9, y1 - y0 - 0.15, FAM[f], z=2)
    pg.txt(1.25, (y0 + y1) / 2, FAMLAB[f], ha="center", rotation=90)
for s in rows:
    pg.txt(M0 - 0.5, ypos[s] + rh / 2, lab(s), ha="right")
    for j, c in enumerate(cols):
        v = M.loc[s, c]
        if np.isnan(v):
            tile(pg, M0 + cw * j, ypos[s], cw, rh, EMPTY)
            continue
        tile(pg, M0 + cw * j, ypos[s], cw, rh, PUOR(norm(v)))
        if CL.loc[s, c] == 1:
            pg.c.plot(M0 + cw * j + (cw - 0.15) / 2, ypos[s] + (rh - 0.15) / 2, "o", ms=1.35, mfc=INK, mec="white", mew=0.25, zorder=4)
for j, c in enumerate(cols):
    pg.txt(M0 + cw * j + cw / 2 - 0.6, MT - 0.7, coh(c), ha="left", va="bottom", rotation=45, rotation_mode="anchor")
start = 0
for t in TUMS:
    n_ = sum(ICB_TUM[c] == t for c in cols)
    pg.box(M0 + cw * start, MB + 0.25, cw * n_ - 0.3, 0.8, TYPE[t], z=2)
    pg.txt(M0 + cw * (start + n_ / 2) - 0.15, MB + 2.5, t if n_ > 1 else t[:3] + ".", ha="center")
    start += n_
pg.txt(6, 6.5, "ΔAUC vs matched null")
pg.cbar(6, 7.8, 18.5, 1.25, PUOR, -LIM, LIM, [-0.15, 0, 0.15])
pg.c.plot(6.7, 13.3, "o", ms=1.6, mfc=INK, mec="white", mew=0.25, zorder=4)
pg.txt(8, 13.3, "clears null 95th pct")
pg.src("a", R[["signature", "family", "cohort", "tum", "n", "AUC", "rand_med", "rand_p95", "delta", "clear", "null_signed"]])

# ================================================================= b  right marginal: median gain per signature, aligned rows
pg.panel(82, 3, "b", "Median gain")
ax = pg.ax(82.5, MT, 12, MB - MT)
zero(ax, 0)
for s in rows:
    point(ax, med[s], ypos[s] + rh / 2, FAM[famd[s]], ms=2.5, filled=bool(sigp[s] < 0.05))
ax.set_ylim(MB, MT); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlim(-0.015, 0.13); ax.set_xticks([0, 0.1], ["0", "0.1"]); ax.set_xlabel("median ΔAUC", labelpad=1)
clearn = R.groupby("signature").clear.sum()
for xx, h_ in ((98.2, "genes"), (103.4, "clear"), (109.3, "P")):
    pg.txt(xx, MT - 1.2, h_, ha="right", weight="bold")
for s in rows:
    yy = ypos[s] + rh / 2
    pg.txt(98.2, yy, f"{X.n_genes[s]:.0f}", ha="right")
    pg.txt(103.4, yy, f"{clearn[s]:.0f}", ha="right")
    pg.txt(109.3, yy, f"{sigp[s]:.3f}", ha="right")
pg.c.plot(96.4, 6.8, "o", ms=2.5, mfc=INK, mec="white", mew=0.3, zorder=4)
pg.txt(97.8, 6.8, "P < 0.05")
pg.c.plot(96.4, 9.6, "o", ms=2.3, mfc="white", mec=INK, mew=0.5, zorder=4)
pg.txt(97.8, 9.6, "P ≥ 0.05")
pg.src("b", pd.DataFrame({"signature": rows, "family": [famd[s] for s in rows], "median_gain": med[rows].values, "genes": X.n_genes[rows].values,
                          "cohorts_clearing": clearn[rows].values, "signed_rank_P": sigp[rows].values}))

# ================================================================= c  clearance by family (old b)
CAX, CBX, CW = 120.5, 155, 23.5
RY2, RY3, RH = 46.0, 89.0, 28.5
pg.panel(112.5, 3, "c", "Tests clearing, by family")
ax = pg.ax(CAX, 8, CW, RH)
PSIG = R.groupby(["family", "signature"]).clear.mean().mul(100).reset_index()
zero(ax, 5, axis="y")
out = []
for i, f in enumerate(FAMORD):
    g = R[R.family == f]; k, n = int(g.clear.sum()), len(g)
    lo, hi = proportion_confint(k, n, method="wilson")
    v = PSIG.loc[PSIG.family == f, "clear"].to_numpy()
    ax.plot(i - 0.13 + swarm(v, 0.095), v, "o", ms=1.9, mfc=FAM[f], mec="white", mew=0.2, ls="none", zorder=3)
    ci(ax, 100 * lo, 100 * hi, i + 0.34, INK, lw=0.7, horizontal=False)
    ax.plot([i + 0.34 - 0.12, i + 0.34 + 0.12], [100 * k / n] * 2, color=INK, lw=1.0, solid_capstyle="butt", zorder=4)
    ax.text(i + 0.02, 67.5, f"{100 * k / n:.0f}%\n{100 * lo:.0f}–{100 * hi:.0f}", ha="center", va="top", fontsize=5, linespacing=1.1)
    out.append((f, k, n, 100 * k / n, 100 * lo, 100 * hi))
ax.text(2.56, 3.6, "chance (5%)", ha="right", va="top", fontsize=5, color=MUTED)
ax.set_xticks(range(3), ["inflamed/\nIFN", "other", "confounder\ncontrol"], fontsize=5)
ax.tick_params(axis="x", length=0)
ax.set_ylim(-3.5, 68); ax.set_xlim(-0.55, 2.58); ax.set_yticks([0, 10, 20, 30, 40, 50]); ax.set_ylabel("tests clearing the null (%)")
pg.src("c", pd.concat([pd.DataFrame(out, columns=["family", "tests_clearing", "tests", "pooled_pct", "wilson_lo", "wilson_hi"]),
                       PSIG.rename(columns={"clear": "signature_pct"})]))

# ================================================================= d  count against the permutation null (old c)
pg.panel(147, 3, "d", "Signatures at P < 0.05")
ax = pg.ax(CBX, 8, CW, RH)
bins = np.arange(-0.5, 31.5)
hist(ax, dv[dv < nsig], bins, weights=np.full((dv < nsig).sum(), 100 / len(dv)))
hist(ax, dv[dv >= nsig], bins, color=SIG, weights=np.full((dv >= nsig).sum(), 100 / len(dv)))
kk = np.arange(0, 31)
pmf = binom.pmf(kk, 52, 0.05) * 100
ax.plot(kk[pmf > 1e-2], pmf[pmf > 1e-2], color=INK, lw=0.6, zorder=3)
ax.plot([nsig, nsig], [0.01, 6.5], color=SIG, lw=0.8, zorder=4)
ax.text(nsig - 0.7, 2.6, f"observed {nsig}\nP = {ph.perm_P:.3f}", ha="right", va="center", fontsize=5, linespacing=1.15)
ax.set_yscale("log"); ax.set_ylim(0.01, 60); ax.set_yticks([0.1, 1, 10], ["0.1", "1", "10"]); ax.minorticks_off()
ax.set_xlim(-1, 30); ax.set_xticks([0, 10, 20, 30]); ax.set_xlabel("P < 0.05, of 52 signatures"); ax.set_ylabel("permutations (%)")
key(pg, CBX + 8.0, 8.8, [(NULLC, "permutation null")], gap=0)
key(pg, CBX + 8.0, 11.2, [(INK, "binomial", "line")], gap=0)
pg.src("d", pd.DataFrame({"signatures_at_P05": kk, "permutations_pct": [(dv == k_).mean() * 100 for k_ in kk], "binomial_pct": pmf}))

# ================================================================= e  family contrast against the same permutations (new)
pg.panel(112.5, RY2, "e", "Inflamed − confounder gap")
ax = pg.ax(CAX, RY2 + 5, CW, RH)
ob = float(pc_.observed)
bins = np.arange(-15, 27, 2)
hist(ax, dc[dc < ob], bins, weights=np.full((dc < ob).sum(), 100 / len(dc)))
hist(ax, dc[dc >= ob], bins, color=SIG, weights=np.full((dc >= ob).sum(), 100 / len(dc)))
h_, e_ = np.histogram(dc, bins)
top_ = 100 * h_.max() / len(dc)
ax.plot([ob, ob], [0, top_ * 1.02], color=SIG, lw=0.8, zorder=4)
ax.text(25, top_ * 1.42, f"observed {ob:.1f}", ha="right", va="top", fontsize=5, linespacing=1.15)
ax.set_xlim(-15, 25); ax.set_xticks([-10, 0, 10, 20], ["−10", "0", "10", "20"]); ax.set_ylim(0, top_ * 1.42); ax.set_yticks([0, 5, 10, 15])
ax.set_xlabel("inflamed − confounder (pct pts)"); ax.set_ylabel("permutations (%)")
pg.src("e", pd.DataFrame({"bin_left_pp": e_[:-1], "permutations_pct": 100 * h_ / len(dc), "observed_pp": ob, "null_median_pp": pc_.perm_null_median,
                          "null_p95_pp": pc_.perm_null_p95, "null_max_pp": pc_.perm_null_max}))

# ================================================================= f  published against random (old d)
pg.panel(147, RY2, "f", "AUC, signatures and null")
ax = pg.ax(CBX, RY2 + 5, CW, RH)
bins = np.arange(0.30, 0.96, 0.025)
hist(ax, R.AUC, bins, color=tint(SIG, 0.85))
hist(ax, R.rand_med, bins, weights=-np.ones(len(R)))
ax.axhline(0, color=INK, lw=0.4, zorder=3)
h1, _ = np.histogram(R.AUC, bins); h0, _ = np.histogram(R.rand_med, bins)
ma, mr = R.AUC.median(), R.rand_med.median()
ax.plot([ma, ma], [0, h1.max() * 1.04], color=INK, lw=0.6, zorder=4); ax.plot([mr, mr], [0, -h0.max() * 1.02], color=INK, lw=0.6, zorder=4)
ax.text(0.96, h1.max() * 1.0, f"published\nmedian {ma:.3f}", ha="right", va="top", fontsize=5, linespacing=1.1)
ax.text(0.96, -h0.max() * 0.5, f"matched random\nmedian {mr:.3f}", ha="right", va="top", fontsize=5, linespacing=1.1)
ax.set_xlim(0.33, 0.96); ax.set_ylim(-h0.max() * 1.08, h1.max() * 1.12); ax.set_xticks([0.4, 0.6, 0.8], ["0.4", "0.6", "0.8"])
ax.set_yticks([-200, -100, 0, 100], ["200", "100", "0", "100"]); ax.set_xlabel("AUC for response"); ax.set_ylabel("signature × cohort tests")
pg.src("f", pd.DataFrame({"bin_left": bins[:-1], "signature_tests": h1, "matched_random_tests": h0}))

# ================================================================= g  every test against the null 95th percentile (new)
pg.panel(112.5, RY3, "g", "AUC and null 95th pct")
ax = pg.ax(CAX, RY3 + 5, CW, RH)
no = R[~R.clear]
ax.hexbin(no.rand_p95, no.AUC, gridsize=(13, 16), cmap=GREYHEX, mincnt=1, linewidths=0, rasterized=True, extent=(0.57, 0.80, 0.5, 0.9))
ax.plot([0.5, 1], [0.5, 1], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=2)
for f in ("Other", "Confounder", "Inflamed/IFN"):
    s_ = R[R.clear & (R.family == f)]
    point(ax, s_.rand_p95, s_.AUC, FAM[f], ms=1.9, mew=0.2)
ax.set_xlim(0.58, 0.80); ax.set_ylim(0.5, 0.96); ax.set_xticks([0.6, 0.7, 0.8], ["0.6", "0.7", "0.8"]); ax.set_yticks([0.5, 0.7, 0.9], ["0.5", "0.7", "0.9"])
ax.text(0.584, 0.955, f"{int(R.clear.sum())} of {len(R)}\nclearing", ha="left", va="top", fontsize=5, linespacing=1.1)
ax.set_xlabel("null 95th percentile"); ax.set_ylabel("AUC")
pg.src("g", R[["signature", "family", "cohort", "rand_p95", "AUC", "clear"]])

# ================================================================= h  median AUC per signature (new)
pg.panel(147, RY3, "h", "Signature median AUC")
ax = pg.ax(CBX, RY3 + 5, CW, RH)
ms = R.groupby("signature").agg(auc=("AUC", "median"), null=("rand_med", "median")).sort_values("auc", ascending=False)
ms["rank"] = np.arange(1, len(ms) + 1)
point(ax, ms["rank"], ms["null"], NULLC, ms=1.7, mew=0.2)
for s, r in ms.iterrows():
    point(ax, r["rank"], r["auc"], FAM[famd[s]], ms=2.0, mew=0.25 if sigp[s] < 0.05 else 0.45, filled=bool(sigp[s] < 0.05))
ax.text(4, ms.auc.max() + 0.002, f"{lab(ms.index[0])} {ms.auc.max():.3f}", ha="left", va="center", fontsize=5)
ax.set_xlim(-1, 54); ax.set_xticks([1, 26, 52]); ax.set_ylim(0.53, 0.70); ax.set_yticks([0.55, 0.6, 0.65, 0.7], ["0.55", "0.60", "0.65", "0.70"])
ax.set_xlabel("signature rank"); ax.set_ylabel("median AUC")
key(pg, CBX + 1.0, RY3 + 5 + RH - 3.4, [(NULLC, "matched random", "o")], gap=0)
pg.src("h", ms.reset_index().assign(family=lambda d: d.signature.map(famd), signed_rank_P=lambda d: d.signature.map(sigp)).rename(columns={"auc": "median_AUC", "null": "median_random_set_median"}))

# ================================================================= i  top marginal: share clearing per cohort, aligned columns
Y_I = MB + 7.6
pg.panel(1, Y_I, "i", "Tests clearing, per cohort")
ax = pg.ax(M0, Y_I + 2.4, MW, 12.5)
pc = {c: 100 * R[R.cohort == c].clear.mean() for c in cols}
zero(ax, 5, axis="y")
for j, c in enumerate(cols):
    point(ax, j + 0.5, pc[c], TYPE[ICB_TUM[c]], ms=3.0)
    ax.text(j + 0.5, pc[c] + 7, f"{pc[c]:.0f}", ha="center", va="bottom", fontsize=5)
ax.set_xlim(0, len(cols)); ax.set_ylim(-3, 84); ax.set_yticks([0, 25, 50]); ax.set_ylabel("% clear", labelpad=1)
no_x(ax)
for k_, t in enumerate(TUMS):
    key(pg, 4.5, Y_I + 4.2 + 2.25 * k_, [(TYPE[t], t, "o")], gap=0)
pg.src("i", pd.DataFrame({"cohort": cols, "tumor": [ICB_TUM[c] for c in cols], "tests": [int((R.cohort == c).sum()) for c in cols],
                          "clearing": [int(R[R.cohort == c].clear.sum()) for c in cols], "pct_clearing": [pc[c] for c in cols]}))

# ================================================================= j  family clearing per cohort (old e)
Cc = []
for c in cols:
    g = R[R.cohort == c]
    Cc.append(dict(cohort=c, tum=ICB_TUM[c], infl=100 * g[g.family == "Inflamed/IFN"].clear.mean(), conf=100 * g[g.family == "Confounder"].clear.mean()))
Cc = pd.DataFrame(Cc).sort_values("infl", kind="stable").reset_index(drop=True)
TAB = {"melanoma": "mel", "kidney": "kidn", "lung": "lung", "bladder": "blad", "gastric": "gast"}
NAME = {"ICB_Limagne1": "Limagne", "ICB_Miao1": "Miao"}
Y_J = Y_I
pg.panel(84.5, Y_J, "j", "Family clearing, per cohort")
ax = pg.ax(101, Y_J + 6.2, 17, 26)
zero(ax, 5)
for i, r in Cc.iterrows():
    ax.plot([r.conf, r.infl], [i, i], color=OTHER, lw=0.5, zorder=2)
    point(ax, r.conf, i, CONF, ms=2.7)
    point(ax, r.infl, i, INFL, ms=2.7)
    yy = Y_J + 6.2 + 26 * (1 - (i + 0.7) / (len(Cc) + 0.4))
    pg.txt(121.8, yy, f"{r.infl:.0f}", ha="right")
    pg.txt(126.4, yy, f"{r.conf:g}", ha="right")
pg.txt(121.8, Y_J + 4.8, "infl.", ha="right", weight="bold"); pg.txt(126.4, Y_J + 4.8, "conf.", ha="right", weight="bold")
ax.set_yticks(range(len(Cc)), [f"{NAME.get(c, coh(c))} ({TAB[t]})" for c, t in zip(Cc.cohort, Cc.tum)], fontsize=5)
ax.tick_params(axis="y", length=0); ax.set_ylim(-0.7, len(Cc) - 0.3); ax.spines["left"].set_visible(False)
ax.set_xlim(-4, 104); ax.set_xticks([0, 50, 100]); ax.set_xlabel("clearing the null (%)")
pg.src("j", Cc.rename(columns={"infl": "inflamed_pct", "conf": "confounder_pct"}))

# ================================================================= k  agreement between cohorts (new)
pg.panel(128, Y_J, "k", "Gains between cohorts")
RHO = M.corr(method="spearman")
x0, y0, cx, cy = 138.8, Y_J + 9.6, 3.33, 2.44
for j, c in enumerate(cols):
    pg.txt(x0 + cx * j + cx / 2 - 0.5, y0 - 0.5, NAME.get(c, coh(c)), ha="left", va="bottom", rotation=62, rotation_mode="anchor")
for i, a_ in enumerate(cols):
    pg.txt(x0 - 0.5, y0 + cy * i + cy / 2, NAME.get(a_, coh(a_)), ha="right")
    for j, b_ in enumerate(cols):
        if i == j:
            tile(pg, x0 + cx * j, y0 + cy * i, cx, cy, EMPTY)
            continue
        v = RHO.loc[a_, b_]
        tile(pg, x0 + cx * j, y0 + cy * i, cx, cy, TEALOR((v + 1) / 2), f"{v + 0.0:.1f}".replace("-0.0", "0.0").replace("0.", ".").replace("-", "−"))   # no leading zero: "−0.6" filled the cell
pg.src("k", RHO.rename_axis("cohort").reset_index())

# ================================================================= l  gains per cohort: the data behind the marginal (new)
Y_L = Y_I + 18.0
pg.panel(1, Y_L, "l", "Gain, per cohort")
ax = pg.ax(M0, Y_L + 2.2, MW, 18.5)
zero(ax, 0, axis="y")
for j, c in enumerate(cols):
    g = R[R.cohort == c]
    off = qjit(g.delta, 0.34)
    for f in ("Other", "Confounder", "Inflamed/IFN"):
        s = (g.family == f).to_numpy()
        ax.plot(j + 0.5 + off[s], g.delta[s], "o", ms=1.45, mfc=FAM[f], mec="white", mew=0.12, ls="none", zorder=3)
    q = g.delta.median()
    ax.plot([j + 0.16, j + 0.84], [q, q], color=INK, lw=0.8, solid_capstyle="butt", zorder=5)
ax.set_xlim(0, len(cols)); ax.set_ylim(-0.12, 0.34); ax.set_yticks([-0.1, 0, 0.1, 0.2, 0.3], ["−0.1", "0", "0.1", "0.2", "0.3"]); ax.set_ylabel("ΔAUC", labelpad=1)
no_x(ax)
for k_, f in enumerate(FAMORD):
    key(pg, 4.5, Y_L + 4.6 + 2.25 * k_, [(FAM[f], FAMLAB[f], "o")], gap=0)
pg.line(4.5, Y_L + 4.6 + 2.25 * 3, 6.3, Y_L + 4.6 + 2.25 * 3, lw=0.8)
pg.txt(6.9, Y_L + 4.6 + 2.25 * 3, "median")
pg.src("l", R[["cohort", "tum", "signature", "family", "delta"]])

# ================================================================= n  share clearing against cohort size (old i)
Y_R3 = Y_J + 40.8
Y_M = Y_R3
pg.panel(84.5, Y_M, "n", "Clearing and cohort size")
ax = pg.ax(92, Y_M + 5.0, 33, 18.5)
nco = ncoh; shr = R.groupby("cohort").clear.mean() * 100
OFF = {"ICB_Gide": (1.09, 0, "left"), "ICB_Puch": (1.09, 0, "left"), "ICB_Jung": (0.91, 0, "right"), "ICB_Kim": (1.09, 0, "left"), "ICB_Limagne1": (1.09, 0, "left"),
       "ICB_Mariathasan": (1.0, 7, "center"), "ICB_Liu": (1.09, 0, "left"), "ICB_Braun": (1.09, 1.8, "left"), "ICB_Hugo": (0.91, 0, "right"), "ICB_Miao1": (0.91, -3.0, "right"),
       "ICB_Riaz": (1.09, 3.6, "left"), "ICB_Van_Allen": (1.09, -3.2, "left")}
zero(ax, 5, axis="y")
for c in cols:
    point(ax, nco[c], shr[c], TYPE[ICB_TUM[c]], ms=2.9)
    fx, dy, ha = OFF[c]
    ax.text(nco[c] * fx, shr[c] + dy, NAME.get(c, coh(c)), fontsize=5, va="center", ha=ha)
rs = pd.DataFrame({"n": nco, "s": shr}).corr(method="spearman").iloc[0, 1]
ax.text(320, 71, f"Spearman r = {rs:.2f}\nn = {len(cols)}".replace("-", "−"), ha="right", va="top", fontsize=5, linespacing=1.1)
ax.set_xscale("log"); ax.set_xticks([30, 60, 120, 240], ["30", "60", "120", "240"]); ax.minorticks_off()
ax.set_xlim(15.5, 320); ax.set_ylim(-3, 72); ax.set_yticks([0, 25, 50]); ax.set_xlabel("patients", labelpad=1); ax.set_ylabel("% of 52 clearing")
pg.src("n", pd.DataFrame({"cohort": cols, "tumor": [ICB_TUM[c] for c in cols], "patients": [nco[c] for c in cols], "pct_clearing": [shr[c] for c in cols]}).assign(spearman_r=rs))

# ================================================================= m  the cohorts, aligned columns (new)
Y_N = Y_R3
pg.panel(1, Y_N, "m", "Cohorts")
T = pd.DataFrame({"patients": ncoh[cols], "responders": CS.n_resp[cols]})
assert (CS.n[cols] == ncoh[cols]).all() and T.patients.sum() == 827
T["response (%)"] = 100 * T.responders / T.patients
g = R.groupby("cohort")
T["median AUC"] = g.AUC.median()[cols]; T["random median"] = g.rand_med.median()[cols]; T["median gain"] = g.delta.median()[cols]
ch = 2.5
for i, (c_, fmt, fc) in enumerate((("patients", "{:.0f}", lambda v: BLUE(np.log(v / 20) / np.log(14))), ("responders", "{:.0f}", lambda v: BLUE(np.log(v / 7) / np.log(11))),
                                   ("response (%)", "{:.0f}", lambda v: PURPLE(v / 60)), ("median AUC", "{:.2f}", lambda v: BLUE((v - 0.5) / 0.25)),
                                   ("random median", "{:.2f}", lambda v: BLUE((v - 0.5) / 0.25)), ("median gain", "{:.2f}", lambda v: PUOR(norm(v))))):
    yy = Y_N + 2.2 + ch * i
    pg.txt(M0 - 0.5, yy + ch / 2, c_, ha="right")
    for j, c in enumerate(cols):
        tile(pg, M0 + cw * j, yy, cw, ch, fc(T.loc[c, c_]), fmt.format(T.loc[c, c_]).replace("-", "−"))
for j, c in enumerate(cols):
    pg.txt(M0 + cw * j + cw / 2 + 0.6, Y_N + 2.2 + ch * 6 + 0.5, coh(c), ha="right", va="top", rotation=45, rotation_mode="anchor")
pg.src("m", T.rename_axis("cohort").reset_index())

# ================================================================= o  tumor type by family (new)
Y_O = Y_R3
pg.panel(128, Y_O, "o", "Tumor type and family")
x0, y0, tw, tc = 141, Y_O + 8.6, 6.2, 2.85
GR = R.groupby(["tum", "family"])
out = []
for k_, (hd, xs) in enumerate((("median gain", x0), ("clearing (%)", x0 + 3 * tw + 0.9))):
    pg.txt(xs + 1.5 * tw, y0 - 5.0, hd, ha="center", weight="bold")
    pg.line(xs + 0.2, y0 - 3.9, xs + 3 * tw - 0.3, y0 - 3.9, lw=0.4)
    for j, f in enumerate(FAMORD):
        pg.txt(xs + tw * j + tw / 2, y0 - 2.0, {"Inflamed/IFN": "infl.", "Other": "other", "Confounder": "conf."}[f], ha="center")
for i, t in enumerate(TUMS):
    pg.txt(x0 - 0.5, y0 + tc * i + tc / 2, t, ha="right")
    for j, f in enumerate(FAMORD):
        g_ = R[(R.tum == t) & (R.family == f)]
        mm = g_.groupby("signature").delta.median().median()
        sh = 100 * g_.clear.mean()
        tile(pg, x0 + tw * j, y0 + tc * i, tw, tc, PUOR(norm(mm)), f"{mm:.3f}".replace("-", "−"))
        tile(pg, x0 + 3 * tw + 0.9 + tw * j, y0 + tc * i, tw, tc, BLUE(sh / 60), f"{sh:.0f}")
        out.append((t, f, mm, int(g_.clear.sum()), len(g_), sh))
pg.src("o", pd.DataFrame(out, columns=["tumor", "family", "family_median_of_signature_medians", "tests_clearing", "tests", "pct_clearing"]))

# ================================================================= p  by tumor type (old f)
Y_P = 202.0
pg.panel(1, Y_P, "p", "Gain, by tumor type")
ax = pg.ax(9.5, Y_P + 5, 58, 26.5)
zero(ax, 0, axis="y")
qs = []
for i, t in enumerate(TUMS):
    for j, f in enumerate(FAMORD):
        v = R[(R.tum == t) & (R.family == f)].groupby("signature").delta.median()
        xx = i + (j - 1) * 0.29
        ax.plot(xx - 0.03 + qjit(v, 0.065, nbins=8), v, "o", ms=1.45, mfc=FAM[f], mec="white", mew=0.12, ls="none", zorder=3)
        q = np.percentile(v, [5, 25, 50, 75, 95])
        ci(ax, q[0], q[4], xx + 0.085, INK, lw=0.35, horizontal=False)
        ci(ax, q[1], q[3], xx + 0.085, INK, lw=1.0, horizontal=False)
        ax.plot(xx + 0.085, q[2], "o", ms=2.1, mfc="white", mec=INK, mew=0.5, ls="none", zorder=5)
        qs.append((t, f, *q))
ax.set_xticks(range(5), [f"{t}\n({R[R.tum == t].cohort.nunique()})" for t in TUMS], fontsize=5); ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.5, 4.5); ax.set_ylim(-0.16, 0.30); ax.set_yticks([-0.1, 0, 0.1, 0.2, 0.3], ["−0.1", "0", "0.1", "0.2", "0.3"]); ax.set_ylabel("median ΔAUC per signature")
key(pg, 31, Y_P + 0.1, [(FAM[f], {"Inflamed/IFN": "inflamed/IFN", "Other": "other", "Confounder": "confounder control"}[f], "o") for f in FAMORD], gap=1.6)
pg.src("p", pd.DataFrame(qs, columns=["tumor", "family", "q05", "q25", "q50", "q75", "q95"]))

# ================================================================= q  volcano (old g)
pg.panel(71.5, Y_P, "q", "Effect and evidence")
ax = pg.ax(80, Y_P + 5, 27, 26.5)
lp = -np.log10(sigp)
ax.plot([0, 0], [0, 3.05], color=INK, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)
ax.axhline(-np.log10(0.05), color=INK, lw=0.5, ls=(0, (2, 1.5)), zorder=1)
ax.text(0.128, -np.log10(0.05) - 0.06, "P = 0.05", ha="right", va="top", fontsize=5)
for f in ("Other", "Confounder", "Inflamed/IFN"):
    s_ = [s for s in sigp.index if famd[s] == f]
    for s in s_:
        ms_ = np.sqrt(4 + 2.2 * np.log2(X.n_genes[s])) * 0.72
        point(ax, med[s], lp[s], FAM[f], ms=ms_, filled=bool(sigp[s] < 0.05), mew=0.25 if sigp[s] < 0.05 else 0.45)
top4 = lp.sort_values(ascending=False).index[:4]
stack = sorted(top4, key=lambda s: (s != "APM_Thompson", -med[s]))
for k, s in enumerate(stack):
    ax.annotate(lab(s), (med[s], lp[s]), xytext=(0.072, 4.42 - k * 0.37), fontsize=5, ha="right", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.25, color=MUTED, shrinkA=0.5, shrinkB=1.2, relpos=(1, 0.5)))
ax.set_xlim(-0.03, 0.13); ax.set_xticks([0, 0.05, 0.1], ["0", "0.05", "0.10"]); ax.set_ylim(0, 4.6); ax.set_yticks([0, 1, 2, 3])
ax.set_xlabel("median ΔAUC"); ax.set_ylabel("−log10 P (signed rank)")
pg.src("q", pd.DataFrame({"signature": sigp.index, "family": famd[sigp.index].values, "median_gain": med[sigp.index].values, "signed_rank_P": sigp.values, "genes": X.n_genes[sigp.index].values}))

# ================================================================= r  every test by family (old h)
pg.panel(110.5, Y_P, "r", "Every test, by family")
ax = pg.ax(118.5, Y_P + 5, 26, 26.5)
zero(ax, 0, axis="y")
out = []
for i, f in enumerate(FAMORD):
    v = R[R.family == f].delta.to_numpy()
    ax.plot(i + qjit(v, 0.3, nbins=16), v, "o", ms=1.3, mfc=FAM[f], mec="white", mew=0.1, ls="none", zorder=3)
    q = np.percentile(v, [25, 50, 75])
    ci(ax, q[0], q[2], i + 0.4, INK, lw=0.8, horizontal=False)
    ax.plot([i - 0.32, i + 0.46], [q[1], q[1]], color=INK, lw=0.8, solid_capstyle="butt", zorder=5)
    out.append((f, len(v), *q))
ax.set_xticks(range(3), [f"{FAMLAB[f]}\n({(R.family == f).sum()})" for f in FAMORD], fontsize=5); ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.5, 2.6); ax.set_ylim(-0.14, 0.36); ax.set_yticks([-0.1, 0, 0.1, 0.2, 0.3], ["−0.1", "0", "0.1", "0.2", "0.3"]); ax.set_ylabel("ΔAUC")
pg.src("r", pd.concat([pd.DataFrame(out, columns=["family", "tests", "q25", "median", "q75"]), R[["signature", "family", "cohort", "delta"]]]))

# ================================================================= s  gain against size (old j)
pg.panel(147.5, Y_P, "s", "Gain and signature size")
ax = pg.ax(155.5, Y_P + 5, 23, 26.5)
zero(ax, 0, axis="y")
for f in ("Other", "Confounder", "Inflamed/IFN"):
    s_ = [s for s in med.index if famd[s] == f]
    for s in s_:
        point(ax, X.n_genes[s], med[s], FAM[f], ms=2.3, filled=bool(sigp[s] < 0.05), mew=0.25 if sigp[s] < 0.05 else 0.45)
rg = pd.DataFrame({"g": X.n_genes[med.index], "m": med}).corr(method="spearman").iloc[0, 1]
ax.text(560, 0.128, f"Spearman r = {rg:.2f}\nn = {len(med)}".replace("-", "−"), ha="right", va="top", fontsize=5, linespacing=1.1)
ax.set_xscale("log"); ax.set_xticks([3, 10, 30, 100, 300], ["3", "10", "30", "100", "300"]); ax.minorticks_off(); ax.set_xlim(2.4, 560)
ax.set_xlabel("genes"); ax.set_ylabel("median ΔAUC"); ax.set_ylim(-0.03, 0.13); ax.set_yticks([0, 0.05, 0.1], ["0", "0.05", "0.10"])
pg.src("s", pd.DataFrame({"signature": med.index, "family": famd[med.index].values, "genes": X.n_genes[med.index].values, "median_gain": med.values, "signed_rank_P": sigp[med.index].values}).assign(spearman_r=rg))

check_palette("tumor type", TYPE)
check_palette("family", FAM)
pg.sources = dict(sorted(pg.sources.items())); pg.letters = sorted(pg.letters)
pg.save("Figure2")
