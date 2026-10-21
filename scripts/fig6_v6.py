"""Figure 6 (v6) - why this was not seen: the conventional validation measures something else.

a-c  the immune axis in three gastric cohorts, split at the median, numbers at risk
d    hazard ratio per s.d. in each cohort and pooled
e    the same log-rank test on 1,000 random axes per cohort
f    how many cohorts a random axis is significant in, against independence
g    52 signatures x 20 TCGA cancer types against their matched null (hero; row order as Fig. 3a)
h    how often each family clears its null, by endpoint
i    median gain per signature, by endpoint
j    stage plus age against the best signature, per cancer type
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy.stats import t as tdist
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.stats import binom, norm
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import v4kit as K
import figstyle as fs
K.apply()
D = K.D; RR = f"{D}/revision/results"
COH = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254 (ACRG)")]
RND = pd.read_csv(f"{RR}/v2_random_logrank.csv")


def load(tag):
    d = pd.read_csv(f"{RR}/B2_scores_{tag}.csv")
    d = d[d.time.notna() & (d.time > 0) & d.event.notna() & d.NetImmune.notna()].copy()
    if tag == "TCGA": d["time"] = d.time / 30.44
    return d


PR = pd.read_csv(f"{D}/icb/native_prognosis_matched.csv")
PR["delta"] = PR.C - PR.rand_med; PR["beats"] = PR.C > PR.rand_p95
V3 = pd.read_csv(f"{D}/icb/v3_dissociation.csv").set_index("signature")
PR["family"] = PR.signature.map(V3.family); famd = V3.family
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]
os_ = PR[PR.endpoint == "OS"]
C = K.Canvas(292)

# ============================================================ a-c  survival
XS = [13, 53, 93]; KMP = {}
for j, ((tag, name), x0) in enumerate(zip(COH, XS)):
    ax = C.ax(x0, 12, 30, 24)
    d = load(tag); hi = d.NetImmune > d.NetImmune.median()
    tmax = {"TCGA": 120, "GSE84437": 150, "GSE62254": 100}[tag]
    ticks = np.arange(0, tmax + 1, 30 if tmax > 100 else 25)
    for k, (m, c) in enumerate([(hi, fs.BLUE), (~hi, fs.VERM)]):
        km = KaplanMeierFitter().fit(d.time[m], d.event[m]); sf = km.survival_function_
        ax.step(sf.index, sf.iloc[:, 0], where="post", color=c, lw=1.0)
        cens = d.time[m & (d.event == 0)]
        ax.scatter(cens, km.survival_function_at_times(cens), marker="|", s=5, c=c, lw=.4, zorder=4)
        for t, n in zip(ticks, [int((d.time[m] >= t).sum()) for t in ticks]):
            C.text(x0 + 30 * (t + tmax * .03) / (tmax * 1.05), 43.5 + k * 3.0, f"{n}", size=5.0, ha="center", va="center", color=c)
    p = logrank_test(d.time[hi], d.time[~hi], d.event[hi], d.event[~hi]).p_value; KMP[tag] = p
    ax.set_xlim(-tmax * .03, tmax * 1.02); ax.set_ylim(0, 1.03); ax.set_xticks(ticks); ax.set_yticks([0, .5, 1])
    ax.set_xlabel("months", labelpad=1)
    if j == 0: ax.set_ylabel("overall survival")
    ax.text(.97, .97, fs.fmt_p(p, "$P$ = "), transform=ax.transAxes, ha="right", va="top", fontsize=K.T_SMALL)
    K.frame(ax); ax.set_title(name, fontsize=K.T_SMALL, loc="left", pad=2, fontweight="bold")
    C.letter(x0 - 11.5, 4, "abc"[j])
C.fig.add_artist(Line2D([2 / C.w, 5.5 / C.w], [1 - 43.5 / C.h] * 2, color=fs.BLUE, lw=1.3))
C.fig.add_artist(Line2D([2 / C.w, 5.5 / C.w], [1 - 46.5 / C.h] * 2, color=fs.VERM, lw=1.3))
C.text(6, 43.5, "high", size=5.0, va="center"); C.text(6, 46.5, "low", size=5.0, va="center")

# ============================================================ d  pooled HR
rows = []
for tag, name in COH:
    d = load(tag); d["z"] = (d.NetImmune - d.NetImmune.mean()) / d.NetImmune.std()
    s = CoxPHFitter().fit(d[["time", "event", "z"]], "time", "event").summary.loc["z"]
    rows.append(dict(cohort=name, n=len(d), ev=int(d.event.sum()), b=s["coef"], se=s["se(coef)"]))
H = pd.DataFrame(rows)
w = 1 / H.se**2; bf = (w * H.b).sum() / w.sum(); Qs = (w * (H.b - bf)**2).sum(); k = len(H)
tau2 = max(0, (Qs - (k - 1)) / (w.sum() - (w**2).sum() / w.sum()))
ws = 1 / (H.se**2 + tau2); bd = (ws * H.b).sum() / ws.sum(); sd = np.sqrt(1 / ws.sum())
I2 = max(0, (Qs - (k - 1)) / Qs); p_meta = 2 * norm.sf(abs(bd / sd))
ax = C.ax(150, 12, 28, 24)
for yy, (_, r) in zip([3.2, 2.2, 1.2], H.iterrows()):
    lo, hi_ = np.exp(r.b - 1.96 * r.se), np.exp(r.b + 1.96 * r.se)
    ax.plot([lo, hi_], [yy] * 2, color=fs.SLATE, lw=1.0)
    ax.scatter(np.exp(r.b), yy, s=10 + r.ev / 7, marker="s", c=fs.SLATE, lw=0, zorder=4)
E, L_, U_ = np.exp(bd), np.exp(bd - 1.96 * sd), np.exp(bd + 1.96 * sd)
ax.fill([L_, E, U_, E], [0, .30, 0, -.30], color=fs.BLUE, lw=0, zorder=4)
ax.axvline(1, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_xlim(.58, 1.1); ax.set_ylim(-.6, 3.7); ax.set_yticks([3.2, 2.2, 1.2, 0])
ax.set_yticklabels([n for _, n in COH] + ["pooled"], fontsize=5.0)
ax.get_yticklabels()[-1].set_color(fs.BLUE); ax.get_yticklabels()[-1].set_fontweight("bold")
ax.tick_params(axis="y", length=0); ax.set_xticks([.6, .8, 1.0])
ax.set_xlabel("HR per s.d. (OS)"); K.frame(ax, left=False)
C.letter(128, 4, "d")

# ============================================================ e  random axes
EST = {}
for j, ((tag, name), x0) in enumerate(zip(COH, XS)):
    ax = C.ax(x0, 66, 30, 22)
    u = RND[(RND.cohort == tag) & (RND.type == "random")]; o = RND[(RND.cohort == tag) & (RND.type == "immune axis")].iloc[0]
    x = -np.log10(u.p); bins = np.arange(0, 6.01, .25); cut = -np.log10(.05)
    n_, b_, patches = ax.hist(x.clip(upper=5.9), bins=bins, color=fs.LIGHT, ec="white", lw=.25)
    for pt, left in zip(patches, b_[:-1]):
        if left >= cut - 1e-9: pt.set_facecolor(fs.GREY)
    ax.axvline(cut, color=fs.INK, lw=.55, ls=(0, (3, 2)))
    ox = -np.log10(o.p)
    ax.annotate("", xy=(ox, 0), xytext=(ox, max(n_) * .55), arrowprops=dict(arrowstyle="-|>", color=fs.BLUE, lw=1.0, mutation_scale=7))
    EST[tag] = (100 * (u.p < .05).mean(), int((u.p <= o.p).sum()) + 1)
    ax.set_xlim(0, 6.0); ax.set_ylim(0, max(n_) * 1.12); ax.set_xticks([0, 2, 4, 6]); ax.set_xticklabels(["0", "2", "4", "≥6"])
    ax.set_xlabel("log-rank −log$_{10}$ $P$")
    if j == 0: ax.set_ylabel("random axes")
    ax.set_title(name, fontsize=K.T_SMALL, loc="left", pad=2); K.frame(ax)
C.letter(1.5, 58, "e")

# ============================================================ f  joint significance
ax = C.ax(146, 66, 32, 22)
Wm = RND[RND.type == "random"].pivot(index="set", columns="cohort", values="p")
cnt = (Wm < .05).sum(1); obs = np.array([(cnt == k_).mean() for k_ in range(4)]) * 100
exp_ = binom.pmf(range(4), 3, .05) * 100
ax.bar(np.arange(4), obs, .62, color=fs.GREY, lw=0)
ax.scatter(np.arange(4), exp_, s=16, facecolor="white", edgecolor=fs.INK, lw=.7, zorder=4)
ax.annotate("", xy=(3, obs[3] + 5), xytext=(3, 40), arrowprops=dict(arrowstyle="-|>", color=fs.BLUE, lw=.9, mutation_scale=7))
ax.text(3, 42, "immune\naxis", ha="center", va="bottom", fontsize=5.0, color=fs.BLUE)
ax.set_xticks(range(4)); ax.set_xlim(-.55, 3.55); ax.set_ylim(0, 92); ax.set_yticks([0, 30, 60, 90]); ax.tick_params(axis="x", length=0)
ax.set_xlabel("cohorts with $P$ < 0.05"); ax.set_ylabel("random axes (%)"); K.frame(ax)
C.letter(128, 58, "f")

# ============================================================ g  the pan-cancer matrix
R = pd.read_csv(f"{D}/icb/native_response_matched_v2.csv"); R = R[R.cohort != "ICB_VanDenEnde"]
rmed = (R.AUC - R.rand_med).groupby(R.signature).median()
rows_ = []
for f in FAMORD: rows_ += sorted([s for s in rmed.index if famd[s] == f], key=lambda s: -rmed[s])
M = os_.pivot_table(index="signature", columns="cancer", values="delta").loc[rows_]
Bt = os_.pivot_table(index="signature", columns="cancer", values="beats").loc[M.index, M.columns]
pmed = M.median(1); nr, nc = M.shape
MX, MY, MW, MH = 38, 116, 68, 114; LIM = .1
ax = C.ax(MX, MY, MW, MH)
ax.imshow(M.values, aspect="auto", cmap=fs.DIV, norm=Normalize(-LIM, LIM), interpolation="nearest")
for i in range(nr + 1): ax.axhline(i - .5, color="white", lw=.3)
for j in range(nc + 1): ax.axvline(j - .5, color="white", lw=.3)
yy_, xx_ = np.where(Bt.values == 1); ax.scatter(xx_, yy_, s=2.2, c=fs.INK, lw=0, zorder=4)
ax.set_xticks(range(nc)); ax.set_xticklabels(M.columns, rotation=90, fontsize=5.0); ax.xaxis.tick_top()
ax.set_yticks(range(nr)); ax.set_yticklabels([K.lab(s) for s in M.index], fontsize=5.0)
for t, s in zip(ax.get_yticklabels(), M.index): t.set_color(fs.INK if famd[s] != "Other" else fs.SLATE)
ax.tick_params(length=0, pad=1.5)
for sp in ax.spines.values(): sp.set_visible(False)
st = C.ax(4.5, MY, 1.6, MH); st.set_xlim(0, 1); st.set_ylim(nr - .5, -.5); st.axis("off")
start = 0
for f in FAMORD:
    n = int(sum(famd[s] == f for s in M.index))
    st.add_patch(Rectangle((0, start - .5), 1, n, fc=fs.FAM[f], ec="white", lw=.5))
    st.text(-.9, start + n / 2 - .5, fs.FAMLAB[f].replace(" control", ""), rotation=90, ha="center", va="center",
            fontsize=K.T_SMALL, color=fs.FAM[f], fontweight="bold")
    start += n
rb = C.ax(MX + MW + 2, MY, 14, MH)
rb.barh(range(nr), pmed.values, height=.74, color=[fs.FAM[famd[s]] for s in M.index], lw=0)
rb.set_ylim(nr - .5, -.5); rb.set_yticks([]); rb.axvline(0, color=fs.INK, lw=.45)
rb.set_xlim(-.012, .035); rb.set_xticks([0, .03]); rb.set_xticklabels(["0", "0.03"])
rb.xaxis.tick_top(); rb.xaxis.set_label_position("top"); rb.set_xlabel("median ΔC", fontsize=K.T_SMALL, labelpad=2)
rb.tick_params(axis="x", labelsize=5.0)
for sp in ("left", "bottom", "right"): rb.spines[sp].set_visible(False)
rb.spines["top"].set_visible(True)
cax = C.ax(8, 103, 22, 1.8)
cb = C.fig.colorbar(plt.cm.ScalarMappable(Normalize(-LIM, LIM), fs.DIV), cax=cax, orientation="horizontal", ticks=[-.1, 0, .1])
cb.outline.set_linewidth(.3); cax.tick_params(labelsize=5.0, length=1.5, width=.4, pad=1)
cax.set_title("ΔC-index vs matched null", fontsize=K.T_SMALL, pad=2)
C.text(8, 108.5, "•  clears its null", size=5.0)
C.letter(1.5, 97, "g")

# ============================================================ h  clearance by family
ax = C.ax(140, 118, 38, 28)
clr = {}; rng6 = np.random.default_rng(3)
for i, ep in enumerate(["OS", "DSS", "PFI"]):
    s_ = PR[PR.endpoint == ep]
    for j_, f_ in enumerate(FAMORD):
        t = s_[s_.family == f_].groupby("signature").beats.mean() * 100; xc = i + (j_ - 1) * .27
        h_ = tdist.ppf(.975, len(t) - 1) * t.sem()                       # mean over signatures, t-based 95% interval
        ax.scatter(xc + rng6.uniform(-.07, .07, len(t)), t, s=3, c=fs.FAM[f_], lw=0, alpha=.35, zorder=2)
        ax.vlines(xc, t.mean() - h_, t.mean() + h_, color=fs.FAM[f_], lw=1.2, zorder=3)
        ax.scatter([xc], [t.mean()], s=18, c="white", ec=fs.FAM[f_], lw=1.1, zorder=4)
        clr.setdefault(ep, {})[f_] = (round(t.mean(), 1), round(t.mean() - h_, 1), round(t.mean() + h_, 1))
ax.axhline(5, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_xticks(range(3)); ax.set_xticklabels(["OS", "DSS", "PFI"]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("cancer types\ncleared (%)"); ax.set_ylim(-1, 42); ax.set_yticks([0, 10, 20, 30, 40]); ax.set_xlim(-.5, 2.5)
for j_, f_ in enumerate(FAMORD):
    xk = [-.02, .36, .64][j_]
    ax.scatter([xk], [1.10], s=12, c=fs.FAM[f_], marker="s", lw=0, transform=ax.transAxes, clip_on=False)
    ax.text(xk + .04, 1.10, fs.FAMLAB[f_].replace(" control", ""), transform=ax.transAxes, fontsize=5.0, va="center")
K.frame(ax)
C.letter(128, 110, "h")

# ============================================================ i  median gain by endpoint
ax = C.ax(140, 158, 38, 26)
rng = np.random.default_rng(0)
for i, ep in enumerate(["OS", "DSS", "PFI"]):
    v = PR[PR.endpoint == ep].groupby(["signature", "family"]).delta.median().reset_index()
    for f_ in FAMORD:
        g = v[v.family == f_]
        ax.scatter(np.full(len(g), i) + rng.uniform(-.2, .2, len(g)), g.delta, s=6, c=fs.FAM[f_], lw=0, alpha=.85)
    ax.plot([i - .32, i + .32], [v.delta.median()] * 2, color=fs.INK, lw=1.1)
ax.axhline(0, color=fs.LIGHT, lw=.6)
ax.set_xticks(range(3)); ax.set_xticklabels(["OS", "DSS", "PFI"]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("median ΔC-index"); ax.set_xlim(-.5, 2.5); ax.set_ylim(-.03, .075); K.frame(ax)
C.letter(128, 152, "i")

# ============================================================ j  stage + age
ax = C.ax(144, 196, 34, 34)
C2 = pd.read_csv(f"{D}/pancan/ceiling2_results.csv")
clin = C2[C2.endpoint == "OS"].groupby("cancer").C_clinical.first()
b = pd.concat([clin.rename("clin"), os_.groupby("cancer").C.max().rename("best"), os_.groupby("cancer").n.first().rename("n")], axis=1).dropna()
m0, m1 = .5, .9
ax.plot([m0, m1], [m0, m1], color=fs.INK, lw=.55, ls=(0, (3, 2)), zorder=2)
win = b.clin > b.best
ax.scatter(b.clin[win], b.best[win], s=np.clip(b.n[win] / 20, 6, 30), c=fs.TEAL, lw=.35, ec="white", zorder=3)
ax.scatter(b.clin[~win], b.best[~win], s=np.clip(b.n[~win] / 20, 6, 30), c=fs.BLUE, lw=.35, ec="white", zorder=3)
ax.set_xlim(m0, m1); ax.set_ylim(m0, m1); ax.set_aspect("equal")
ax.set_xticks([.5, .7, .9]); ax.set_yticks([.5, .7, .9])
ax.set_xlabel("C-index, stage + age"); ax.set_ylabel("C-index, best signature")
ax.text(.52, .88, "signature better", fontsize=5.0, color=fs.BLUE, va="top")
ax.text(.895, .505, "stage + age better", fontsize=5.0, color=fs.TEAL, ha="right", va="bottom")
K.frame(ax)
C.letter(128, 190, "j")

# ================================================= k, l  survival under checkpoint blockade
# The benchmark above compares response under treatment with survival without it. These two
# panels put survival on both sides: the same signatures, the same matched null, in patients
# who did receive checkpoint blockade (10 cohorts with overall and 9 with progression-free
# survival) and in the 20 TCGA cancer types who did not.
ST = pd.read_csv(f"{D}/icb/v4_standardised_effects.csv").set_index("signature")
famd = V3.family; ORD3 = ["Other", "Confounder", "Inflamed/IFN"]
IS = pd.read_csv(f"{D}/icb/native_icbsurv_matched.csv")
IS["family"] = IS.signature.map(famd)
IS["z"] = (IS.C - IS.rand_med) / IS.rand_sd
ZI = {ep: IS[IS.endpoint == ep].groupby("signature").z.median() for ep in ("OS", "PFS")}
SET = [("response\nunder ICB", ST.z_resp),
       ("progression-free\nsurvival under ICB", ZI["PFS"]),
       ("overall survival\nunder ICB", ZI["OS"]),
       ("survival without\nICB (TCGA)", ST.z_prog)]
ax = C.ax(22, 248, 74, 34)
rng2 = np.random.default_rng(4)
for f in ORD3:
    med = []
    for i, (nm, v) in enumerate(SET):
        w = v.reindex(famd[famd == f].index.values).dropna()
        ax.scatter(i + rng2.uniform(-.13, .13, len(w)), w, s=5, c=fs.FAM[f], lw=0, alpha=.35, zorder=2)
        med.append(w.median())
    ax.plot(range(len(SET)), med, color=fs.FAM[f], lw=1.5, zorder=4)
    ax.scatter(range(len(SET)), med, s=22, c=fs.FAM[f], lw=.4, ec="white", zorder=5)

ax.axhline(0, color=fs.INK, lw=.55, ls=(0, (3, 2)))
ax.set_xticks(range(len(SET))); ax.set_xticklabels([n for n, _ in SET], fontsize=K.T_SMALL)
ax.tick_params(axis="x", length=0); ax.set_xlim(-.4, 3.05); ax.set_ylim(-2.6, 5.2)
ax.set_ylabel("effect (null s.d.)"); K.frame(ax)
K.key(ax, [(fs.FAM[f], fs.FAMLAB[f]) for f in ["Inflamed/IFN", "Other", "Confounder"]], x=.02, y=.97, dy=.10, s=16)
C.letter(8, 242, "k")

ax = C.ax(126, 248, 34, 34)
zo, zp = ZI["OS"], ST.z_prog
both = [s_ for s_ in famd.index if s_ in zo.index and s_ in zp.index]
for f in ORD3:
    g = [s_ for s_ in both if famd.get(s_) == f]
    ax.scatter(zp.reindex(g), zo.reindex(g), s=12, c=fs.FAM[f], lw=.3, ec="white", zorder=3)
lim = (-1.6, 2.6)
ax.plot(lim, lim, color=fs.INK, lw=.5, ls=(0, (3, 2))); ax.axhline(0, color=fs.LIGHT, lw=.6); ax.axvline(0, color=fs.LIGHT, lw=.6)
ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_xticks([-1, 0, 1, 2]); ax.set_yticks([-1, 0, 1, 2])
ax.set_xlabel("without ICB (null s.d.)"); ax.set_ylabel("overall survival\nunder ICB (null s.d.)")
ax.text(.03, .97, "36 of 52 signatures\nabove the line", transform=ax.transAxes, fontsize=K.T_SMALL,
        color=fs.SLATE, va="top")
K.frame(ax)
C.letter(112, 242, "l")
print(f"k,l | ICB survival cohorts OS {IS[IS.endpoint=='OS'].cohort.nunique()} PFS {IS[IS.endpoint=='PFS'].cohort.nunique()}"
      f" | inflamed median z: PFS {ZI['PFS'].reindex(famd[famd == 'Inflamed/IFN'].index.values).median():+.2f}"
      f" OS {ZI['OS'].reindex(famd[famd == 'Inflamed/IFN'].index.values).median():+.2f}"
      f" | above identity: {int((zo.reindex(both).values > zp.reindex(both).values).sum())}/{len(both)}")

C.save("Figure6")
print("Figure6 v6 | KM", {k: f"{v:.2g}" for k, v in KMP.items()}, f"| HR {np.round(np.exp(H.b), 2).tolist()} pooled {E:.3f} ({L_:.3f}-{U_:.3f}) P {p_meta:.1e} I2 {100*I2:.0f}%")
print("  random", {k: (round(v[0], 1), v[1]) for k, v in EST.items()}, f"| joint {np.round(obs, 1).tolist()} exp {np.round(exp_, 3).tolist()}")
print(f"  matrix mean {os_.delta.mean():+.4f} above0 {100*(os_.delta>0).mean():.1f}% clear {100*os_.beats.mean():.1f}% | clearance {clr} | stage+age {int(win.sum())}/{len(b)}")
