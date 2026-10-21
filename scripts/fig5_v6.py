"""Figure 5 (v6) - response discrimination is attenuated and heterogeneous in renal cancer.

a  AUC of the inflamed axis for response in each of 13 cohorts, pooled by renal status
b  the contrast estimated where it varies: log-odds of response per s.d., cohort by cohort
c  the same contrast under five specifications
d  ROC curves   e  patient-level logistic fits   f  response rate by quintile
g  the axis in responders and non-responders   h, i  the axis against tumor mutational burden
The computations are those of the previous version (same estimator, same seeds).
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from sklearn.metrics import roc_auc_score, roc_curve
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
from statsmodels.stats.proportion import proportion_confint
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D
P = pd.read_csv(f"{D}/icb/positive_patientlevel.csv").rename(columns={"inf": "infl"})
rng = np.random.default_rng(0)
NR, RN = fs.BLUE, fs.VERM
TUM = {"Kidney": "kidney", "Melanoma": "melanoma", "Lung": "lung", "Bladder": "bladder",
       "Gastric": "gastric"}
NICE = lambda c: c.replace("ICB_", "").replace("_", " ")


def auc_ci(y, x, B=2000):
    a = roc_auc_score(y, x); bs = []
    for _ in range(B):
        i = rng.integers(0, len(y), len(y))
        if 0 < y[i].sum() < len(y): bs.append(roc_auc_score(y[i], x[i]))
    return a, np.percentile(bs, 2.5), np.percentile(bs, 97.5)


rows = []
for (co, tu), g in P.groupby(["cohort", "tumour"]):
    a, lo, hi = auc_ci(g.y.to_numpy(), g.infl.to_numpy())
    rows.append(dict(cohort=co, tumour=tu, n=len(g), renal=int(tu == "Kidney"), auc=a, lo=lo, hi=hi))
F = pd.DataFrame(rows)
pooled = {}
for lab_, sub in [("non-renal", P[P.is_renal == 0]), ("renal", P[P.is_renal == 1])]:
    pooled[lab_] = auc_ci(sub.y.to_numpy(), sub.infl.to_numpy())
PC = pd.read_csv(f"{D}/icb/v4_renal_percohort.csv")
CL = pd.read_csv(f"{D}/icb/v4_renal_cluster.csv")


def cl_row(analysis, term=None):
    q = CL[CL.analysis == analysis]
    if term is not None: q = q[q.term == term]
    return q.iloc[0]


C = K.Canvas(178)

# ============================================================ a  per-cohort AUC
ax = C.ax(4, 12, 76, 60)
nr = F[F.renal == 0].sort_values("auc"); rr = F[F.renal == 1].sort_values("auc")
y = 0
for blk, name, col in [(rr, "renal", RN), (nr, "non-renal", NR)]:
    a_, lo_, hi_ = pooled[name]
    ax.plot([lo_, hi_], [y, y], color=col, lw=1.1, zorder=3)
    ax.fill([lo_, a_, hi_, a_], [y, y + .40, y, y - .40], color=col, lw=0, zorder=4)
    ax.text(.03, y, f"pooled {name}", ha="left", va="center", fontsize=K.T_NOTE,
            fontweight="bold", color=col, transform=ax.get_yaxis_transform())
    ax.text(1.0, y, f"{a_:.2f}", ha="right", va="center", fontsize=K.T_NOTE, color=col,
            fontweight="bold", transform=ax.get_yaxis_transform())
    y0 = y; y += 1.3
    for _, r in blk.iterrows():
        c = fs.TUMCOL[TUM[r.tumour]]
        ax.plot([r.lo, r.hi], [y, y], color=c, lw=.9, zorder=3, solid_capstyle="butt")
        ax.scatter(r.auc, y, s=8 + r.n / 6, c=c, lw=.4, ec="white", zorder=4)
        ax.text(.03, y, NICE(r.cohort), ha="left", va="center", fontsize=K.T_SMALL,
                transform=ax.get_yaxis_transform())
        ax.text(1.0, y, f"{r.auc:.2f}", ha="right", va="center", fontsize=K.T_SMALL,
                color=fs.SLATE, transform=ax.get_yaxis_transform())
        y += 1
    ax.axhspan(y0 - .65, y - .45, color=fs.PAPER, zorder=0, lw=0)
    y += .9
ax.axvline(.5, color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=1)
ax.set_ylim(-.9, y - .5); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlim(-.28, 1.15); ax.set_xticks([.2, .4, .6, .8, 1.0]); ax.spines["bottom"].set_bounds(.1, 1.0)
ax.set_xlabel("AUC for response")
K.frame(ax, left=False)
xl = .03
for t in ["melanoma", "lung", "bladder", "gastric", "kidney"]:
    ax.scatter(xl, 1.045, s=12, c=fs.TUMCOL[t], lw=0, transform=ax.transAxes, clip_on=False)
    ax.text(xl + .022, 1.045, t, fontsize=K.T_SMALL, va="center", transform=ax.transAxes)
    xl += .075 + .0135 * len(t)
C.letter(1.5, 4, "a")

# ============================================================ b  cohort-level contrast
ax = C.ax(118, 12, 44, 60)
y = 0; ypos, ylab, ycol = [], [], []
for name in ["renal", "non-renal"]:
    blk = PC[PC.renal == int(name == "renal")].sort_values("logOR")
    c = RN if name == "renal" else NR
    p_ = cl_row("meta-analysis within group", name)
    ax.plot([p_.ci_lo, p_.ci_hi], [y, y], color=c, lw=1.1, zorder=3)
    ax.fill([p_.ci_lo, p_.estimate, p_.ci_hi, p_.estimate], [y, y + .42, y, y - .42], color=c, lw=0, zorder=4)
    ypos.append(y); ylab.append(f"pooled {name}"); ycol.append(c)
    ax.text(1.02, y, f"OR {np.exp(p_.estimate):.2f}", ha="left", va="center", fontsize=K.T_NOTE,
            color=c, fontweight="bold", transform=ax.get_yaxis_transform())
    y0 = y; y += 1.3
    for _, r in blk.iterrows():
        ax.plot([r.logOR - 1.96 * r.se, r.logOR + 1.96 * r.se], [y, y], color=c, lw=.8, zorder=3)
        ax.scatter(r.logOR, y, s=7 + r.n / 8, c=c, lw=.35, ec="white", zorder=4)
        ypos.append(y); ylab.append(NICE(r.cohort)); ycol.append(fs.INK)
        y += 1
    ax.axhspan(y0 - .65, y - .45, color=fs.PAPER, zorder=0, lw=0)
    y += .9
ax.axvline(0, color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=1)
ax.set_yticks(ypos); ax.set_yticklabels(ylab, fontsize=K.T_SMALL)
for t, c, lb in zip(ax.get_yticklabels(), ycol, ylab):
    t.set_color(c)
    if lb.startswith("pooled"): t.set_fontweight("bold"); t.set_fontsize(K.T_NOTE)
ax.tick_params(axis="y", length=0)
ax.set_ylim(-.9, y - .5); ax.set_xlim(-1.45, 2.55); ax.set_xticks([-1, 0, 1, 2])
ax.set_xlabel("log odds of response per s.d.")
K.frame(ax, left=False)
C.letter(92, 4, "b")


# ============================================================ c  five specifications
ax = C.ax(42, 90, 30, 30)
SPECS = [("meta-regression (primary)", "meta-regression (DL, cohort level)", fs.BLUE, "D"),
         ("exact permutation", "exact permutation of the renal label across cohorts", fs.SLATE, "o"),
         ("GEE, 13 clusters", "GEE, cluster-robust by cohort", fs.GREY, "D"),
         ("random slope by cohort", "patient level, random slope by cohort (glmer)", fs.GREY, "D"),
         ("patient level, cohort FE", "patient level, cohort fixed effects (published)", fs.VERM, "D")]
yy = np.arange(len(SPECS))[::-1]
for k_, (lab_, key_, c, mk) in enumerate(SPECS):
    r = cl_row(key_); e_ = r.estimate
    lo = r.ci_lo if np.isfinite(r.ci_lo) else (e_ - 1.96 * r.se if np.isfinite(r.se) else np.nan)
    hi = r.ci_hi if np.isfinite(r.ci_hi) else (e_ + 1.96 * r.se if np.isfinite(r.se) else np.nan)
    if np.isfinite(lo): ax.plot([lo, hi], [yy[k_]] * 2, color=c, lw=1.0)
    if mk == "o": ax.scatter(e_, yy[k_], s=18, facecolor="white", edgecolor=c, lw=.9, zorder=4)
    else: ax.scatter(e_, yy[k_], s=16, c=c, marker="D", lw=.3, ec="white", zorder=4)
    ax.text(1.03, yy[k_], fs.fmt_p(r.p, "$P$ = "), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=K.T_SMALL)
ax.axvline(0, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_yticks(yy); ax.set_yticklabels([s_[0] for s_ in SPECS], fontsize=K.T_SMALL); ax.tick_params(axis="y", length=0)
ax.set_xlim(-1.08, .3); ax.set_xticks([-1, -.5, 0]); ax.set_ylim(-.6, len(SPECS) - .4)
ax.set_xlabel("renal − non-renal\n(log odds per s.d.)", fontsize=K.T_SMALL)
K.frame(ax, left=False)
C.letter(1.5, 84, "c")

# ============================================================ d  ROC
ax = C.ax(102, 90, 30, 30)
for (co, tu), g in P.groupby(["cohort", "tumour"]):
    fpr, tpr, _ = roc_curve(g.y, g.infl); ax.plot(fpr, tpr, color=RN if tu == "Kidney" else NR, lw=.5, alpha=.35)
for lab_, sub, c in [("non-renal", P[P.is_renal == 0], NR), ("renal", P[P.is_renal == 1], RN)]:
    fpr, tpr, _ = roc_curve(sub.y, sub.infl); ax.plot(fpr, tpr, color=c, lw=1.6)
ax.plot([0, 1], [0, 1], color=fs.INK, lw=.5, ls=(0, (3, 2)))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate"); K.frame(ax)
K.key(ax, [(NR, "non-renal"), (RN, "renal")], x=.52, y=.20, dy=.12, marker="s", s=10)
C.letter(90, 84, "d")

# ============================================================ e  patient-level fits
ax = C.ax(146, 90, 32, 30)
lims = (P.infl.quantile(.005), P.infl.quantile(.995))
for k, (lab_, sub, c) in enumerate([("non-renal", P[P.is_renal == 0], NR), ("renal", P[P.is_renal == 1], RN)]):
    m = smf.logit("y ~ infl", data=sub).fit(disp=0)
    xs = np.linspace(sub.infl.quantile(.02), sub.infl.quantile(.98), 120)
    ax.plot(xs, 100 * m.predict(pd.DataFrame({"infl": xs})), color=c, lw=1.5)
    q = pd.qcut(sub.infl, 5, labels=False); rq = sub.groupby(q).agg(x=("infl", "median"), p=("y", "mean"))
    ax.scatter(rq.x, 100 * rq.p, s=16, facecolor="white", edgecolor=c, lw=.9, zorder=5)
ax.set_xlim(*lims); ax.set_ylim(0, 70); ax.set_yticks([0, 20, 40, 60])
ax.set_xlabel("inflamed axis (s.d.)"); ax.set_ylabel("response probability (%)"); K.frame(ax)
C.letter(134, 84, "e")

# ============================================================ f  response by quintile
ax = C.ax(12, 142, 34, 30)
w = .2; qtab = {}
for k, (lab_, sub, c) in enumerate([("non-renal", P[P.is_renal == 0], NR), ("renal", P[P.is_renal == 1], RN)]):
    q = pd.qcut(sub.infl, 5, labels=False); t = sub.groupby(q).y.agg(["mean", "size", "sum"]); qtab[lab_] = t
    ci = np.array([proportion_confint(int(k_), int(n_), method="wilson") for k_, n_ in zip(t["sum"], t["size"])]) * 100
    xq = np.arange(5) + (k - .5) * w
    ax.vlines(xq, ci[:, 0], ci[:, 1], color=c, lw=.9, zorder=2)
    ax.plot(xq, 100 * t["mean"], color=c, lw=.7, zorder=3)
    ax.scatter(xq, 100 * t["mean"], s=16, c="white", ec=c, lw=1.0, zorder=4)
ax.set_xticks(range(5)); ax.set_xticklabels(["Q1", "Q2", "Q3", "Q4", "Q5"]); ax.tick_params(axis="x", length=0)
ax.set_xlabel("quintile of the axis"); ax.set_ylabel("responders (%)"); ax.set_ylim(0, 70); ax.set_yticks([0, 20, 40, 60]); ax.set_xlim(-.5, 4.5)
K.key(ax, [(NR, f"non-renal ({int(qtab['non-renal']['size'].sum())})"), (RN, f"renal ({int(qtab['renal']['size'].sum())})")],
      x=.05, y=.95, dy=.13, marker="s", s=10)
K.frame(ax)
C.letter(1.5, 136, "f")

# ============================================================ g  score by response
ax = C.ax(60, 142, 44, 30)
rng2 = np.random.default_rng(0)
pos = {("non-renal", 0): 0, ("non-renal", 1): 1, ("renal", 0): 2.5, ("renal", 1): 3.5}
for (lab_, yv), x0 in pos.items():
    sub = P[(P.is_renal == (lab_ == "renal")) & (P.y == yv)]; c = NR if lab_ == "non-renal" else RN; v = sub.infl.to_numpy()
    parts = ax.violinplot([v], positions=[x0], widths=.80, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(c if yv else "white"); b.set_edgecolor(c); b.set_alpha(.35 if yv else 1); b.set_linewidth(.7)
    ax.scatter(x0 + rng2.uniform(-.12, .12, len(v)), v, s=1.4, c=c, lw=0, alpha=.45, zorder=3, rasterized=True)
    q1, md, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([x0, x0], [q1, q3], color=fs.INK, lw=1.3, zorder=5, solid_capstyle="butt")
    ax.scatter(x0, md, s=12, c="white", ec=fs.INK, lw=.7, zorder=6)
ymax = P.infl.quantile(.998); PV = {}
for lab_, x0 in [("non-renal", .5), ("renal", 3.0)]:
    sub = P[P.is_renal == (lab_ == "renal")]
    p = mannwhitneyu(sub.infl[sub.y == 1], sub.infl[sub.y == 0]).pvalue; PV[lab_] = p
    ax.plot([x0 - .5, x0 - .5, x0 + .5, x0 + .5], [ymax + .15, ymax + .32, ymax + .32, ymax + .15], color=fs.INK, lw=.55)
    ax.text(x0, ymax + .45, fs.fmt_p(p, "$P$ = "), ha="center", va="bottom", fontsize=K.T_SMALL)
ax.set_xticks(list(pos.values())); ax.set_xticklabels(["non-resp.", "resp.", "non-resp.", "resp."], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0)
for lab_, x0, c in [("non-renal", .5, NR), ("renal", 3.0, RN)]:
    ax.text(x0, -.2, lab_, ha="center", transform=ax.get_xaxis_transform(), fontsize=K.T_SMALL, color=c, fontweight="bold")
ax.set_ylabel("inflamed axis (s.d.)"); ax.set_ylim(P.infl.min() - .4, ymax + 1.3); ax.set_xlim(-.6, 4.1); K.frame(ax)
C.letter(50, 136, "g")

# ============================================================ h, i  against TMB
TMB = pd.read_csv(f"{D}/icb/positive_tmb.csv")
T = TMB.copy(); T["tum"] = T.cohort.map(P.drop_duplicates("cohort").set_index("cohort").tumour)
ax = C.ax(116, 142, 24, 30)
cols3 = ["AUC_inflamed", "AUC_TMB", "AUC_combined"]
for _, r in T.iterrows():
    c = RN if r.tum == "Kidney" else NR
    ax.plot([0, 1, 2], r[cols3].astype(float), color=c, lw=.6, alpha=.6); ax.scatter([0, 1, 2], r[cols3].astype(float), s=7, c=c, lw=0, zorder=3)
wm = [np.average(T[c], weights=T.n) for c in cols3]
for k, mu in enumerate(wm): ax.plot([k - .25, k + .25], [mu] * 2, color=fs.INK, lw=1.4, zorder=4)
ax.axhline(.5, color=fs.LIGHT, lw=.6)
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["axis", "TMB", "both"], fontsize=K.T_SMALL); ax.tick_params(axis="x", length=0)
ax.set_xlim(-.4, 2.4); ax.set_ylim(.3, 1.0); ax.set_ylabel("AUC for response"); K.frame(ax)
C.letter(106, 136, "h")
ax = C.ax(162, 142, 16, 30)
T2 = T.sort_values("r_inf_tmb").reset_index(drop=True)
zr = np.arctanh(T2.r_inf_tmb.clip(-.999, .999)); se_ = 1 / np.sqrt((T2.n - 3).clip(lower=1))   # Fisher-z 95% interval
ax.hlines(np.arange(len(T2)), np.tanh(zr - 1.96 * se_), np.tanh(zr + 1.96 * se_), color=[fs.TUMCOL[TUM[t]] for t in T2.tum], lw=.9)
ax.scatter(T2.r_inf_tmb, np.arange(len(T2)), s=8 + T2.n / 8, c=[fs.TUMCOL[TUM[t]] for t in T2.tum], lw=.4, ec="white", zorder=3)
ax.axvline(0, color=fs.INK, lw=.5); wr = np.average(T2.r_inf_tmb, weights=T2.n)
ax.axvline(wr, color=fs.SLATE, lw=.7, ls=(0, (3, 2)))
ax.set_yticks(np.arange(len(T2))); ax.set_yticklabels([f"{NICE(c)} ({n})" for c, n in zip(T2.cohort, T2.n)], fontsize=5.0)
ax.tick_params(axis="y", length=0); ax.set_xlim(-.65, .85); ax.set_xticks([-.5, 0, .5]); ax.set_xticklabels(["−0.5", "0", "0.5"])
ax.set_xlabel("r with TMB", fontsize=K.T_SMALL); K.frame(ax, left=False)
C.letter(146, 136, "i")

C.save("Figure5")
MR = cl_row("meta-regression (DL, cohort level)"); PM = cl_row("exact permutation of the renal label across cohorts")
print(f"Figure5 v6 | pooled AUC {pooled['non-renal'][0]:.3f} / {pooled['renal'][0]:.3f} | MR {MR.estimate:+.3f}±{MR.se:.3f} P {MR.p:.4f}; perm {PM.p:.4f}")
print("  quintiles NR", np.round(100 * qtab['non-renal']['mean'].values), "R", np.round(100 * qtab['renal']['mean'].values),
      f"| MW {PV['non-renal']:.2g} {PV['renal']:.3f} | TMB wm {[round(x, 3) for x in wm]} r {wr:+.3f}")
print('  f Wilson', {k: [(int(a_), int(b_)) for a_, b_ in zip(v['sum'], v['size'])] for k, v in qtab.items()}, '| i r CI', [(c_, round(r_, 2), round(float(np.tanh(z_ - 1.96 * e_)), 2), round(float(np.tanh(z_ + 1.96 * e_)), 2)) for c_, r_, z_, e_ in zip(T2.cohort, T2.r_inf_tmb, zr, se_)])
