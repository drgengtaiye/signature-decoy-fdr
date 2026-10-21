"""Figure 4 (v6) - in randomized trials the benefit tracks the treatment received.

a, b  every signature's score x arm interaction in IMmotion150 and PRINCE (shared row order)
c     signatures fitted inside each arm: hazard ratio per s.d., checkpoint against comparator
d     the checkpoint-versus-comparator hazard ratio across the inflamed score
e     the immune hit count against the no-interaction simulation (primary inference)
f-i   progression-free survival within each arm, split at that arm's median score
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator
from scipy.stats import norm
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D

SHORT = {"GEP18_TIS (Ayers 2017)": "GEP-18/TIS (Ayers)", "IFNG_6gene (Ayers 2017)": "IFNG 6-gene (Ayers)",
         "Chemokine12 (Messina 2012)": "Chemokine-12 (Messina)", "TLS (Cabrita 2020)": "TLS (Cabrita)",
         "TLS_9gene (Prabhakaran)": "TLS 9-gene (Prabhakaran)",
         "EMT_stroma (Chakravarthy 2018)": "EMT/stroma (Chakravarthy)",
         "F_TBRS (Mariathasan 2018)": "F-TBRS (Mariathasan)", "CTRL_Hypoxia (Buffa)": "Hypoxia (Buffa)",
         "CTRL_Proliferation": "Proliferation", "CTRL_Ribosome": "Ribosome",
         "barriernet_Immune": "NetImmune", "barriernet_Barrier": "BarrierIndex",
         "CD8_Tcells": "CD8 T cells", "T_cells": "T cells", "NK_cells": "NK cells", "B_cells": "B cells",
         "Plasma_cells": "Plasma cells", "IFN_typeI": "IFN type I", "M2_TAM": "M2 TAM",
         "MDSC_like": "MDSC-like"}
short = lambda s: SHORT.get(s, s.replace("_", " "))
FC = {"Immune": fs.BLUE, "Control": fs.VERM, "Other": fs.GREY}

A4 = pd.read_csv(f"{D}/pancan/v4_arm_split.csv")
A4S = A4[A4.row_type == "signature"].copy()
A4S["se"] = np.abs(A4S.coef) / np.maximum(norm.ppf(1 - A4S.p / 2), 1e-6)
POOL = A4S[A4S.contrast == "pooled_atezo_vs_sunitinib"].set_index("signature")
PRIN = A4S[A4S.contrast == "nivolumab_vs_non_checkpoint"].set_index("signature")
AS = pd.read_csv(f"{D}/pancan/arm_specific.csv")
KM = {"IMmotion150": pd.read_csv(f"{D}/pancan/immotion150_km.csv", index_col=0),
      "PRINCE": pd.read_csv(f"{D}/pancan/v2_prince_km.csv", index_col=0)}

C = K.Canvas(228)

# ============================================================ a, b  forests
order = POOL.sort_values("coef").index.tolist(); ny = len(order)
XL = (-1.55, 1.25); FY, FH = 14, 80
for k, (S, name, sub, x0, let, lx) in enumerate([
        (POOL, "IMmotion150", "renal · atezolizumab ± bev. vs sunitinib", 38, "a", 1.5),
        (PRIN, "PRINCE", "pancreatic · nivolumab vs non-checkpoint", 79, "b", 75)]):
    ax = C.ax(x0, FY, 36, FH)
    for i, s in enumerate(order):
        r = S.loc[s]; c = FC[r.family]
        lo_, hi_ = r.coef - 1.96 * r.se, r.coef + 1.96 * r.se
        ax.plot([max(lo_, XL[0]), min(hi_, XL[1])], [i, i], color=c, lw=.75,
                alpha=.55 if r.family == "Other" else .95, solid_capstyle="butt")
        ax.scatter(r.coef, i, s=14 if r.family != "Other" else 9, c=c, lw=.35, ec="white",
                   zorder=4, marker="D" if r.p < .05 else "o")
    ax.axvline(0, color=fs.INK, lw=.55, ls=(0, (3, 2)))
    ax.set_ylim(ny - .4, -.6); ax.set_xlim(*XL); ax.set_xticks([-1, 0, 1])
    ax.set_xlabel("interaction coefficient")
    ax.set_yticks(range(ny))
    if k == 0:
        ax.set_yticklabels([short(s) for s in order], fontsize=K.T_SMALL)
        ax.tick_params(axis="y", length=0, pad=9.5)
        for i, s in enumerate(order):
            ax.add_patch(Rectangle((-.075, i - .32), .045, .64, transform=ax.get_yaxis_transform(),
                                   fc=FC[POOL.loc[s, "family"]], ec="none", clip_on=False))
    else:
        ax.set_yticklabels([]); ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False); K.frame(ax, left=False)
    C.title(x0, FY - 7.5, name, fontweight="bold")
    C.text(x0, FY - 4.2, sub, size=K.T_SMALL, color=fs.SLATE)
    C.letter(lx, FY - 8, let)
# marker key under b
kax = C.ax(79, FY + FH + 10, 36, 3); kax.axis("off"); kax.set_xlim(0, 1); kax.set_ylim(0, 1)
for j, (mk, lab_) in enumerate([("D", "$P$ < 0.05"), ("o", "$P$ ≥ 0.05")]):
    kax.scatter([.04 + j * .5], [.5], marker=mk, s=10, c=fs.SLATE, lw=0)
    kax.text(.10 + j * .5, .5, lab_, fontsize=K.T_SMALL, va="center")

# ============================================================ c  inside each arm
ax = C.ax(133, 10, 44, 30)
MK = {"IMmotion150": "o", "PRINCE": "^"}
for tr in ["IMmotion150", "PRINCE"]:
    for f in ["Immune", "Control"]:
        g = AS[(AS.trial == tr) & (AS.family == f)]
        off = -.07 if tr == "IMmotion150" else .07
        for _, r in g.iterrows():
            ax.plot([off, 1 + off], np.log2([r.hr_ici, r.hr_ctl]), color=FC[f], lw=.4,
                    alpha=.45 if f == "Immune" else .85, zorder=2)
        ax.scatter(np.full(len(g), off), np.log2(g.hr_ici), s=8, c=FC[f], marker=MK[tr], lw=0, zorder=3)
        ax.scatter(np.full(len(g), 1 + off), np.log2(g.hr_ctl), s=8, c=FC[f], marker=MK[tr], lw=0, zorder=3)
imm = AS[AS.family == "Immune"]
for i, col in enumerate(["hr_ici", "hr_ctl"]):
    ax.plot([i - .28, i + .28], [np.log2(imm[col].median())] * 2, color=fs.INK, lw=1.4, zorder=5)
ax.axhline(0, color=fs.LIGHT, lw=.6, zorder=0)
ax.set_xticks([0, 1]); ax.set_xticklabels(["checkpoint\narm", "comparator\narm"], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0)
ax.set_xlim(-.45, 1.45); ax.set_ylim(-1.1, 1.1); ax.set_yticks([-1, -.5, 0, .5, 1])
ax.set_ylabel("log$_2$ HR per s.d. (PFS)")
K.frame(ax)
K.key(ax, [(fs.SLATE, "IMmotion150")], x=.05, y=.94, s=10, marker="o")
K.key(ax, [(fs.SLATE, "PRINCE")], x=.05, y=.84, s=10, marker="^")
C.letter(121, 4, "c")

# ============================================================ d  effect across the score
ax = C.ax(133, 54, 44, 28)
TE = []
for tr, col, ls, band in [("PRINCE", fs.SLATE, (0, (3, 1.6)), .12), ("IMmotion150", fs.BLUE, "-", .18)]:
    Kd = KM[tr].copy()
    Kd["z"] = (Kd.s - Kd.s.mean()) / Kd.s.std(); Kd["zx"] = Kd.z * Kd.arm
    cph = CoxPHFitter().fit(Kd[["t", "e", "z", "arm", "zx"]], "t", "e")
    b = cph.params_; V = cph.variance_matrix_
    zs = np.linspace(Kd.z.quantile(.03), Kd.z.quantile(.97), 80)
    est = b["arm"] + b["zx"] * zs
    se = np.sqrt(V.loc["arm", "arm"] + zs ** 2 * V.loc["zx", "zx"] + 2 * zs * V.loc["arm", "zx"])
    ax.fill_between(zs, np.exp(est - 1.96 * se), np.exp(est + 1.96 * se), color=col, alpha=band, lw=0)
    ax.plot(zs, np.exp(est), color=col, lw=1.2, ls=ls)
    TE.append(dict(trial=tr, n=len(Kd), p_zx=cph.summary.loc["zx", "p"],
                   hr_m1=np.exp(b["arm"] - b["zx"]), hr_p1=np.exp(b["arm"] + b["zx"])))
    ax.text(zs[-1] + .08, np.exp(est[-1]) * (1.25 if tr == "PRINCE" else .80), tr, color=col,
            fontsize=K.T_SMALL, va="center", ha="left")
ax.axhline(1, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_yscale("log"); ax.yaxis.set_minor_locator(NullLocator())
ax.set_yticks([.25, .5, 1, 2, 4]); ax.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
ax.set_ylim(.18, 5.5); ax.set_xlim(-2.2, 3.1); ax.set_xticks([-2, -1, 0, 1, 2])
ax.spines["bottom"].set_bounds(-2.2, 2.2)
ax.set_xlabel("inflamed score (s.d. within trial)")
ax.set_ylabel("HR, checkpoint vs comparator")
ax.text(-2.1, 4.6, "favors comparator", fontsize=K.T_SMALL, color=fs.SLATE, va="top")
ax.text(-2.1, .22, "favors checkpoint", fontsize=K.T_SMALL, color=fs.SLATE, va="bottom")
K.frame(ax)
C.letter(121, 48, "d")


# ============================================================ e  no-interaction simulation
V7 = pd.read_csv(f"{D}/pancan/v7_interaction_null_valid.csv")
V7D = pd.read_csv(f"{D}/pancan/v7_null_draws_hit_immune.csv"); KHEAD = 9
KS = sorted(V7D.K.unique())
for j, (tr, col) in enumerate([("IMmotion150", fs.BLUE), ("PRINCE", fs.SLATE)]):
    ax = C.ax(133, 95 + j * 12.5, 44, 9.5)
    o_ = float(V7[(V7.trial == tr) & (V7.statistic == "hit_immune") & (V7.method == f"joint_parametric_pc{KHEAD}")].observed.iloc[0])
    for K_ in KS:                                         # every score subspace, K = 3-20; bars at the drawn K
        v = V7D.loc[(V7D.trial == tr) & (V7D.K == K_), "hit_immune"].to_numpy()
        pc = np.bincount(v, minlength=18)[:18] / len(v) * 100
        if K_ == KHEAD:
            ax.bar(np.arange(18), np.where(pc > 0, pc, np.nan), 1.0, color=[col if k_ < o_ else fs.VERM for k_ in range(18)], alpha=.8, lw=0, zorder=2)
        else:
            ax.step(np.arange(18), np.where(pc > 0, pc, np.nan), where="mid", color=fs.INK, lw=.4, alpha=.55, zorder=3)
    pk = V7[(V7.trial == tr) & (V7.statistic == "hit_immune") & V7.method.str.startswith("joint_parametric")].perm_p
    ax.plot([o_, o_], [.03, 45], color=fs.VERM, lw=.9, zorder=4)
    ax.text(17.5, 400, f"{tr}: {o_:.0f} hits, $P$ = {pk.min():.3f}–{pk.max():.3f}", ha="right", va="top", fontsize=5.0, color=fs.INK)
    ax.set_yscale("log"); ax.set_ylim(.03, 600); ax.set_yticks([.1, 10]); ax.set_yticklabels(["0.1", "10"]); ax.minorticks_off()
    ax.set_xlim(-.6, 17.6); ax.set_xticks([0, 5, 10, 15]); K.frame(ax)
    if j == 0: ax.set_xticklabels([])
ax.set_xlabel("immune hits (of 17)")
C.text(125.5, 106, "replicates (%)", size=K.T_LAB, rotation=90, ha="center", va="center")
C.letter(121, 90, "e")

# ============================================================ f-i  survival within each arm
KMS = [("f", "IMmotion150", 1, "Atezolizumab ± bev."), ("g", "IMmotion150", 0, "Sunitinib"),
       ("h", "PRINCE", 1, "Nivolumab-containing"), ("i", "PRINCE", 0, "Non-checkpoint")]
XS = [14, 59, 104, 149]; KY, KH = 136, 31
KMOUT = []
for (let, tr, arm, name), x0 in zip(KMS, XS):
    ax = C.ax(x0, KY, 30, KH)
    g = KM[tr][KM[tr].arm == arm].copy(); hi = g.s > g.s.median(); tmax = 24
    meds = {}
    for msk, col, lb in [(hi, fs.BLUE, "high"), (~hi, fs.VERM, "low")]:
        kf = KaplanMeierFitter().fit(g.t[msk], g.e[msk]); sf = kf.survival_function_
        ax.step(sf.index, sf.iloc[:, 0], where="post", color=col, lw=1.1)
        cens = g[msk & (g.e == 0) & (g.t <= tmax)].t
        ax.scatter(cens, kf.survival_function_at_times(cens).values, marker="|", s=8, lw=.45,
                   c=col, zorder=4)
        meds[lb] = kf.median_survival_time_
    lr = logrank_test(g.t[hi], g.t[~hi], g.e[hi], g.e[~hi]).p_value
    KMOUT.append(dict(trial=tr, arm=name, n_high=int(hi.sum()), n_low=int((~hi).sum()),
                      med_high=meds["high"], med_low=meds["low"], logrank_p=lr))
    ax.set_xlim(0, tmax); ax.set_ylim(0, 1.03); ax.set_xticks([0, 6, 12, 18, 24]); ax.set_yticks([0, .5, 1])
    ax.set_xlabel("months", labelpad=1)
    if let == "f": ax.set_ylabel("progression-free survival")
    ax.text(.97, .97, fs.fmt_p(lr, "$P$ = "), transform=ax.transAxes, ha="right", va="top",
            fontsize=K.T_NOTE, color=fs.INK if lr < .05 else fs.SLATE)
    K.frame(ax)
    C.title(x0 - 3, KY - 8.5, name, fontweight="bold")
    C.text(x0 - 3, KY - 4.8, tr, size=K.T_SMALL, color=fs.SLATE)
    C.letter(x0 - 12, KY - 9, let)
    for j, (msk, col, lb) in enumerate([(hi, fs.BLUE, "high"), (~hi, fs.VERM, "low")]):
        yy = KY + KH + 10.5 + j * 3.4
        for t_ in [0, 6, 12, 18, 24]:
            C.text(x0 + 30 * t_ / tmax, yy, f"{int((g.t[msk] >= t_).sum())}", size=K.T_SMALL,
                   ha="center", va="center")
        if let == "f":
            C.fig.add_artist(__import__("matplotlib.lines", fromlist=["Line2D"]).Line2D(
                [(x0 - 11.5) / C.w, (x0 - 8.5) / C.w], [1 - yy / C.h] * 2, color=col, lw=1.3))
            C.text(x0 - 8, yy, lb, size=K.T_SMALL, va="center")
    if let == "f":
        C.text(x0 - 11.5, KY + KH + 7.2, "at risk", size=K.T_SMALL, color=fs.SLATE, va="center",
               style="italic")

# ============================================================ j  each arm against its own matched null
# The interaction above is a contrast between arms. This asks a simpler question inside each
# arm: do the signatures beat random gene sets of their own size, scored the same way, in
# these patients? By construction 5% of random sets clear their own 95th percentile.
W8 = pd.read_csv(f"{D}/pancan/v8_within_arm_null.csv")
W9 = pd.read_csv(f"{D}/validation/ispy2_within_arm_null.csv")
W9 = W9.rename(columns={"trial": "trial"})
WA = pd.concat([W8[["trial", "arm", "family", "clears"]], W9[["trial", "arm", "family", "clears"]]])
LAB = {"IMmotion150": "IMmotion150\natezo ± bev\nvs sunitinib", "PRINCE": "PRINCE\nnivolumab\nvs sotigalimab",
       "ISPY2_pembro": "I-SPY2\npaclitaxel + pembro\nvs paclitaxel", "ISPY2_durva": "I-SPY2\ndurva + olaparib\nvs control"}
ORDER = ["IMmotion150", "PRINCE", "ISPY2_pembro", "ISPY2_durva"]
ax = C.ax(20, 196, 150, 24)
bw = .34
for i, tr in enumerate(ORDER):
    for j_, (arm_, col) in enumerate([("checkpoint", fs.BLUE), ("comparator", fs.VERM)]):
        g = WA[(WA.trial == tr) & (WA.arm == arm_) & (WA.family == "Immune")]
        if not len(g):
            continue
        v = 100 * g.clears.mean()
        ax.bar(i + (j_ - .5) * bw, v, bw * .92, color=col, lw=0, zorder=3)
        ax.text(i + (j_ - .5) * bw, max(v + 1.6, 8.5), f"{int(g.clears.sum())}/{len(g)}", ha="center",
                va="bottom", fontsize=K.T_SMALL, color=col)
        c_ = WA[(WA.trial == tr) & (WA.arm == arm_) & (WA.family == "Control")]
        if len(c_):
            ax.scatter([i + (j_ - .5) * bw], [100 * c_.clears.mean()], s=13, marker="_",
                       c=fs.INK, lw=1.1, zorder=4)
# the same comparison one gene at a time, with the arms matched for size
TC = pd.read_csv(f"{D}/icb/tc_gene_arms.csv")
for i, tr in enumerate(ORDER):
    for j_, (arm_, col) in enumerate([("checkpoint", fs.BLUE), ("comparator", fs.VERM)]):
        g = TC[(TC.trial == tr) & (TC.arm == arm_) & (TC.gene_set == "immune")]
        if len(g):
            ax.scatter([i + (j_ - .5) * bw], [100 * float(g.frac_matched.iloc[0])], s=20, facecolors="white",
                       edgecolors=col, lw=.9, zorder=5)
ax.axhline(5, color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=2)
ax.text(.99, .97, "dashed, the 5% a random set clears by construction\ndashes, confounder controls",
        transform=ax.transAxes, fontsize=K.T_SMALL, color=fs.SLATE, ha="right", va="top")
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([LAB[t] for t in ORDER], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0); ax.set_xlim(-.6, 3.6); ax.set_ylim(0, 100)
ax.set_ylabel("immune signatures clearing\ntheir null in that arm (%)")
K.frame(ax)
K.key(ax, [(fs.BLUE, "checkpoint arm"), (fs.VERM, "comparator arm")], x=.02, y=.97, dy=.12, s=16)
ax.scatter([.028], [.73], s=16, facecolors="white", edgecolors=fs.SLATE, lw=.9,
           transform=ax.transAxes, clip_on=False, zorder=6)
ax.text(.055, .73, "the same test one gene at a time (arms matched for size)", transform=ax.transAxes,
        fontsize=K.T_SMALL, color=fs.SLATE, va="center")
C.letter(8, 190, "j")

C.save("Figure4")
print("Figure4 v6 |", pd.DataFrame(TE).round(3).to_string(index=False))
print(pd.DataFrame(KMOUT).round(3).to_string(index=False))
print(f"  arm medians HR {imm.hr_ici.median():.2f} vs {imm.hr_ctl.median():.2f}; immune tests {len(imm)}")
for S, nm in [(POOL, "IMmotion150"), (PRIN, "PRINCE")]:
    im = S[S.family == "Immune"]; ct = S[S.family == "Control"]
    print(f"  {nm}: immune<0 {(im.coef<0).sum()}/17, P<.05 & <0 {((im.p<.05)&(im.coef<0)).sum()}/17, "
          f"controls P<.05 {(ct.p<.05).sum()}/3; n sigs {len(S)}")
