"""Supplementary figure - the randomized trials in detail (panels moved from the former Fig. 2).

a  every patient's inflamed score, by trial and arm
b  the interaction against size-matched random gene sets
c  the IMmotion150 checkpoint arm refitted as its two regimens against the same comparator
d  the two trials against each other, signature by signature
(the no-interaction simulation is Fig. 4e)
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy.stats import norm, spearmanr, binomtest
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D
FC = {"Immune": fs.BLUE, "Control": fs.VERM, "Other": fs.GREY}
KM = {"IMmotion150": pd.read_csv(f"{D}/pancan/immotion150_km.csv", index_col=0),
      "PRINCE": pd.read_csv(f"{D}/pancan/v2_prince_km.csv", index_col=0)}
N4 = pd.read_csv(f"{D}/pancan/v4_matched_interaction_null.csv")
RUNG = N4[N4.row_type == "rung"]; NSIG = N4[N4.row_type == "signature"]; NSUM = N4[N4.row_type == "summary"]
A4 = pd.read_csv(f"{D}/pancan/v4_arm_split.csv")
A4S = A4[A4.row_type == "signature"].copy()
POOL = A4S[A4S.contrast == "pooled_atezo_vs_sunitinib"].set_index("signature")
PRIN = A4S[A4S.contrast == "nivolumab_vs_non_checkpoint"].set_index("signature")
V7 = pd.read_csv(f"{D}/pancan/v7_interaction_null_valid.csv")
V7D = pd.read_csv(f"{D}/pancan/v7_null_draws_hit_immune.csv")
KS = sorted(V7D.K.unique()); KHEAD = 9


def v7(trial, stat, K_=KHEAD):
    return V7[(V7.trial == trial) & (V7.statistic == stat) & (V7.method == f"joint_parametric_pc{K_}")].iloc[0]


def v7range(trial, stat):
    g = V7[(V7.trial == trial) & (V7.statistic == stat) & V7.method.str.startswith("joint_parametric_pc")].perm_p
    return g.min(), g.max()


C = K.Canvas(112)

# ---------------------------------------------------------------- a  score by arm
ax = C.ax(40, 8, 44, 38)
ROWS = [("IMmotion150", 1, "atezolizumab ± bev.", fs.BLUE), ("IMmotion150", 0, "sunitinib", fs.GREY),
        ("PRINCE", 1, "nivolumab-containing", fs.BLUE), ("PRINCE", 0, "non-checkpoint", fs.GREY)]
rng = np.random.default_rng(11)
for i, (tr, arm, lab_, col) in enumerate(ROWS):
    g = KM[tr][KM[tr].arm == arm]
    ax.scatter(g.s, i + rng.uniform(-.23, .23, len(g)), s=4, c=col, lw=0, alpha=.75, zorder=3)
    q1, md, q3 = np.percentile(g.s, [25, 50, 75])
    ax.plot([q1, q3], [i, i], color=fs.INK, lw=1.0, solid_capstyle="butt", zorder=4)
    ax.scatter([md], [i], s=14, marker="|", c=fs.INK, lw=1.2, zorder=5)
ax.set_yticks(range(4))
ax.set_yticklabels([f"{r[2]} (n = {len(KM[r[0]][KM[r[0]].arm == r[1]])})" for r in ROWS], fontsize=K.T_SMALL)
ax.tick_params(axis="y", length=0)
ax.axhline(1.5, color=fs.LIGHT, lw=.6)
for yy, tr in ((.5, "IMmotion150"), (2.5, "PRINCE")):
    ax.text(1.02, yy, tr, transform=ax.get_yaxis_transform(), rotation=90, ha="left", va="center",
            fontsize=K.T_SMALL, fontweight="bold")
ax.axvline(0, color=fs.LIGHT, lw=.6, zorder=1)
ax.set_ylim(3.6, -.6); ax.set_xlabel("inflamed score (s.d. within trial)")
K.frame(ax, left=False)
C.letter(1.5, 3, "a")

# ---------------------------------------------------------------- b  size-matched null
ax = C.ax(110, 18, 66, 28)
for tr, col, ls in [("IMmotion150", fs.BLUE, "-"), ("PRINCE", fs.SLATE, (0, (3, 1.6)))]:
    g = RUNG[(RUNG.trial == tr) & (RUNG.n_draws == 200)].sort_values("size")
    ax.plot(g["size"], 100 * g.frac_p05_checkpoint, color=col, lw=1.1, ls=ls, zorder=4)
    ax.scatter(g["size"], 100 * g.frac_p05_checkpoint, s=8, c=col, lw=0, zorder=5)
    if tr == "PRINCE":
        ax.plot(g["size"], 100 * g.frac_p05, color=col, lw=.8, ls=(0, (1, 1.4)), alpha=.8, zorder=3)
ax.axhline(5, color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=1)
sz = NSIG[NSIG.trial == "IMmotion150"]["size"]
ax.scatter(sz, np.full(len(sz), -.7), marker="|", s=18, lw=.8, c=fs.SLATE, clip_on=False, zorder=6)
ax.set_xscale("log"); ax.set_xlim(3.3, 520); ax.minorticks_off()
ax.set_xticks([4, 10, 40, 100, 400]); ax.set_xticklabels(["4", "10", "40", "100", "400"])
ax.set_ylim(-1.4, 15.5)
ax.set_xlabel("random gene-set size"); ax.set_ylabel("random sets with\ninteraction $P$ < 0.05 (%)")
KEY = [(fs.BLUE, "-", "IMmotion150, checkpoint-favoring"), (fs.SLATE, (0, (3, 1.6)), "PRINCE, checkpoint-favoring"),
       (fs.SLATE, (0, (1, 1.4)), "PRINCE, either direction")]
for j, (col, ls, lab_) in enumerate(KEY):
    yk = 1.34 - j * .115
    ax.plot([.0, .09], [yk, yk], transform=ax.transAxes, color=col, lw=1.0, ls=ls, clip_on=False)
    ax.text(.11, yk, lab_, transform=ax.transAxes, va="center", fontsize=K.T_SMALL)
K.frame(ax)
C.letter(93, 3, "b")
# ---------------------------------------------------------------- c  regimens separately (was d)
CONTRASTS = [("atezo_mono_vs_sunitinib", "atezolizumab vs sunitinib"),
             ("atezo_bev_vs_sunitinib", "atezolizumab + bev. vs sunitinib"),
             ("pooled_atezo_vs_sunitinib", "pooled checkpoint arm vs sunitinib")]
ax = C.ax(40, 64, 60, 34)
rng_e = np.random.default_rng(1)
for i, (ck, lab_) in enumerate(CONTRASTS):
    g = A4S[A4S.contrast == ck]
    for f in ["Other", "Control", "Immune"]:
        h = g[g.family == f]
        ax.scatter(h.coef, np.full(len(h), i) + rng_e.uniform(-.17, .17, len(h)), s=10 if f != "Other" else 7,
                   c=FC[f], lw=.3, ec="white", alpha=1 if f != "Other" else .75, zorder=4 if f == "Immune" else 3)
    ax.plot([g[g.family == "Immune"].coef.median()] * 2, [i - .32, i + .32], color=fs.INK, lw=1.3, zorder=6)
ax.axvline(0, color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=1)
ax.set_ylim(2.6, -.6); ax.set_xlim(-1.15, .85)
ax.set_yticks(range(3)); ax.set_yticklabels([l for _, l in CONTRASTS], fontsize=K.T_SMALL)
ax.tick_params(axis="y", length=0)
ax.set_xlabel("interaction coefficient")
K.frame(ax, left=False)
C.letter(1.5, 58, "c")

# ---------------------------------------------------------------- d  the trials agree
ax = C.ax(128, 62, 36, 36)
m = POOL[["coef", "family"]].join(PRIN[["coef"]], rsuffix="_p").rename(columns={"coef": "coef_r"})
for f in ["Other", "Control", "Immune"]:
    g = m[m.family == f]
    ax.scatter(g.coef_r, g.coef_p, s=12 if f != "Other" else 8, c=FC[f], lw=.35, ec="white",
               zorder=4 if f == "Immune" else 3, alpha=1 if f != "Other" else .75)
lim = .78
ax.plot([-lim, lim], [-lim, lim], color=fs.INK, lw=.5, ls=(0, (3, 2)), zorder=2)
ax.axhline(0, color=fs.LIGHT, lw=.6, zorder=1); ax.axvline(0, color=fs.LIGHT, lw=.6, zorder=1)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xticks([-.5, 0, .5]); ax.set_yticks([-.5, 0, .5])
ax.set_xlabel("coefficient, IMmotion150"); ax.set_ylabel("coefficient, PRINCE")
K.key(ax, [(fs.BLUE, "immune"), (fs.VERM, "confounder control"), (fs.GREY, "other")], x=.05, y=.95, dy=.09, s=10)
K.frame(ax)
C.letter(112, 58, "d")

C.save("SF_trials")
ph150, phprn = v7("IMmotion150", "hit_immune"), v7("PRINCE", "hit_immune")
r150, rprn = v7range("IMmotion150", "hit_immune"), v7range("PRINCE", "hit_immune")
ra, pa = spearmanr(m.coef_r, m.coef_p); mi = m[m.family == "Immune"]; ri, pi = spearmanr(mi.coef_r, mi.coef_p)
sg = int(((mi.coef_r < 0) & (mi.coef_p < 0)).sum() + ((mi.coef_r > 0) & (mi.coef_p > 0)).sum())
ims = {ck: A4S[(A4S.contrast == ck) & (A4S.family == "Immune")] for ck, _ in CONTRASTS}
r_mc = np.corrcoef(ims["atezo_mono_vs_sunitinib"].sort_values("signature").coef,
                   ims["atezo_bev_vs_sunitinib"].sort_values("signature").coef)[0, 1]
print(f"S-trials | hits IMmotion150 {int(ph150.observed)}/17 P {r150[1]:.3f} ({r150[0]:.3f}-{r150[1]:.3f}); "
      f"PRINCE {int(phprn.observed)}/17 P {rprn[1]:.3f} | trials rho all {ra:+.2f} P {pa:.1e}; immune {ri:+.2f} P {pi:.2f}; same sign {sg}/17")
for ck, lb in CONTRASTS:
    g = A4S[A4S.contrast == ck]; im = g[g.family == "Immune"]
    print(f"  {lb}: n {int(g.n_treated.iloc[0])}/{int(g.n_comparator.iloc[0])} events {int(g.events.iloc[0])} "
          f"<0 {(im.coef<0).sum()}/17 P<.05 {((im.p<.05)&(im.coef<0)).sum()}/17 median {im.coef.median():+.2f}")
print(f"  regimen correlation r {r_mc:.2f}; mono sign test {binomtest(16, 17, .5, alternative='greater').pvalue:.1e}")
