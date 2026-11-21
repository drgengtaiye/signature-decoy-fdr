"""Shared panels of Supplementary Figures 6-8: an application of decoyfdr read from icb/v11_pkg."""
import numpy as np, pandas as pd
import v4kit as K
import figstyle as fs
C_FOLD, C_DIR, C_SIG, C_DEC = fs.GREY, fs.BLUE, fs.TEAL, fs.SLATE


def bars(ax, F, combos, ymax):
    """Expected true (solid) and false (hatched) calls with the estimated FDR above; combos: (label, row)."""
    for x, (lab, r, cc) in enumerate(combos):
        et, ef = r.called - r.exp_false, r.exp_false
        ax.bar(x, et, .7, color=cc, lw=0, zorder=3)
        ax.bar(x, ef, .7, bottom=et, color="white", edgecolor=cc, hatch="/////", lw=.6, zorder=3)
        ax.text(x, r.called + ymax * .02, f"{r.fdr:.0f}%" if r.called else "–", ha="center", va="bottom",
                fontsize=K.T_SMALL, color=cc, fontweight="bold")
    ax.set_xticks(range(len(combos))); ax.set_xticklabels([c[0] for c in combos], fontsize=5.0)
    ax.tick_params(axis="x", length=0, pad=1.5); ax.set_ylim(0, ymax); ax.set_xlim(-.6, len(combos) - .4)
    K.frame(ax)


def matching(ax, units, info):
    """Each target's internal r against its co-expressed and random decoys' median, per cohort."""
    pass


def per_signature(ax, calls, stat, fam=None, cls="coex", label_called=True, maxn=60):
    """Each target's median delta (dot) against the 5th-95th percentile and median of its decoys' (bar)."""
    C = calls[calls.stat == stat]
    t = C[C.kind == "target"].set_index("target")
    d = C[(C.kind == "decoy") & (C["class"] == cls)]
    q = d.groupby("target").med_delta.quantile([.05, .5, .95]).unstack()
    order = t.med_delta.sort_values().index
    for i, s in enumerate(order):
        if s in q.index:
            ax.plot([q.loc[s, .05], q.loc[s, .95]], [i, i], color=fs.LIGHT, lw=1.6, solid_capstyle="butt", zorder=2)
            ax.plot([q.loc[s, .5]] * 2, [i - .35, i + .35], color=fs.SLATE, lw=.6, zorder=3)
        c = fs.FAM.get(fam.get(s), fs.SLATE) if fam is not None else fs.SLATE
        ax.scatter(t.loc[s, "med_delta"], i, s=10 if t.loc[s, "call"] else 6, color=c,
                   marker="o" if t.loc[s, "call"] else "o", facecolors=c if t.loc[s, "call"] else "white",
                   edgecolors=c, lw=.6, zorder=4)
    ax.axvline(0, color=fs.INK, lw=.5)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([K.lab(s) if "_" in s else s for s in order], fontsize=3.6)
    ax.tick_params(axis="y", length=0, pad=1); ax.set_ylim(-.8, len(order) - .2)
    K.frame(ax, left=False)
    return order


def globalnull(ax, summ, stats):
    top = max(summ.realized.max(), summ.estimated.max())
    for gi, st in enumerate(stats):
        r = summ[summ.stat == st].iloc[0]
        ax.scatter(gi - .15, r.realized, s=14, color=fs.INK, zorder=4)
        ax.scatter(gi + .15, r.estimated, s=14, marker="D", color=C_SIG, zorder=4)
        ax.text(gi, max(r.realized, r.estimated) + top * .08, f"ratio {r.ratio:.2f}\n({r.ratio_lo:.2f}–{r.ratio_hi:.2f})",
                ha="center", va="bottom", fontsize=K.T_SMALL)
    ax.set_xticks(range(len(stats))); ax.set_xticklabels(["folded" if s == "fold" else "directed" for s in stats],
                                                         fontsize=K.T_SMALL)
    ax.tick_params(axis="x", length=0); ax.set_xlim(-.6, len(stats) - .4); ax.set_ylim(0, top * 1.55)
    ax.set_ylabel("false calls per library\n(outcomes permuted)")
    K.key(ax, [(fs.INK, "signatures called"), (C_SIG, "decoy estimate")], x=.04, y=.14, s=9, size=K.T_SMALL)
    K.frame(ax)


def bars2(ax, combos, ymax):
    """Expected true (solid) and false (hatched) calls by the plain estimate; above each bar the plain estimate
    and, beneath it, the pi0-corrected one. combos: (label, row, color), row from decoy_fdr()."""
    for x, (lab, r, cc) in enumerate(combos):
        et, ef = r.called - r.exp_false, r.exp_false
        ax.bar(x, max(et, 0), .7, color=cc, lw=0, zorder=3)
        ax.bar(x, ef, .7, bottom=max(et, 0), color="white", edgecolor=cc, hatch="/////", lw=.6, zorder=3)
        txt = f"{r.fdr:.0f}%\n{r.fdr_pi0:.0f}%" if r.called else "–"
        ax.text(x, max(r.called, r.exp_false) + ymax * .02, txt, ha="center", va="bottom", fontsize=K.T_SMALL, color=cc,
                fontweight="bold", linespacing=1.0)
    ax.set_xticks(range(len(combos))); ax.set_xticklabels([c[0] for c in combos], fontsize=5.0)
    ax.tick_params(axis="x", length=0, pad=1.5); ax.set_ylim(0, ymax); ax.set_xlim(-.6, len(combos) - .4)
    K.frame(ax)
