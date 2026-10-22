"""Shared style for the v10 main figures: the v8/v9 module (cns.py) with a deeper palette and the forms the redesign
adds (nested interval bands, filled densities, two-dimensional density contours, alluvial ribbons, Kaplan-Meier
curves with numbers at risk, circular tracks).

Palette, validated all-pairs with the dataviz checker (worst CVD deltaE 10.7 deutan, 8.7 tritan; normal-vision 18.9;
every hue >= 3:1 on white): BLUE #1f62b0, VERM #d24b2a, TEAL #169073. The roles are those of cns.py:
family (inflamed/IFN BLUE, confounder VERM, other GRAY), unit (signature INK, co-expressed decoy TEAL, random-gene
decoy LGRAY), statistic (folded SLATE, directed BLUE).
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.colors import LinearSegmentedColormap, to_rgb
import cns

cns.BLUE, cns.ORANGE, cns.AQUA = "#1f62b0", "#d24b2a", "#169073"
cns.GRAY, cns.LGRAY, cns.SLATE = "#a3abb3", "#cdd2d7", "#4c5864"
cns.FAM = {"Inflamed/IFN": cns.BLUE, "Confounder": cns.ORANGE, "Other": cns.GRAY}
cns.STAT = {"fold": cns.SLATE, "dir": cns.BLUE}
cns.DIV = LinearSegmentedColormap.from_list("div", [cns.BLUE, "#8fb3dc", "#f4f4f2", "#eda48f", "#b8371c"])
cns.T_TICK, cns.T_LAB, cns.T_NOTE, cns.T_SMALL, cns.T_HEAD, cns.T_LET = 5.5, 6.0, 5.5, 5.0, 6.4, 8.5
from cns import *                                           # noqa: E402,F401,F403  (after the palette is set)
from cns import tint, frame, dots, ref, lab, Canvas          # noqa: E402

TEAL, VERM = cns.AQUA, cns.ORANGE
SEQ_BLUE = LinearSegmentedColormap.from_list("seqb", ["#f4f6f8", "#c4d6ec", "#6f9bd1", cns.BLUE, "#123d70"])
SEQ_TEAL = LinearSegmentedColormap.from_list("seqt", ["#f4f6f5", "#bfe1d6", "#5fb69c", cns.AQUA, "#0b5442"])


def apply():
    cns.apply()
    plt.rcParams.update({"font.size": cns.T_NOTE, "axes.labelsize": cns.T_LAB, "xtick.labelsize": cns.T_TICK,
                         "ytick.labelsize": cns.T_TICK, "legend.fontsize": cns.T_NOTE, "axes.linewidth": 0.5,
                         "xtick.major.width": 0.5, "ytick.major.width": 0.5})


# ---------------------------------------------------------------- intervals and densities
def band(ax, pos, v=None, q=None, color=None, width=1.0, orient="v", dot=True, z=3, dot_s=9, dot_c=None):
    """Nested interval: 5th-95th percentile (thin, light), 25th-75th (thick, mid) and the median (ringed dot).
    Give the values v or the quantiles q = (q05, q25, q50, q75, q95)."""
    color = color or cns.INK
    if q is None:
        q = np.nanpercentile(np.asarray(v, float), [5, 25, 50, 75, 95])
    q05, q25, q50, q75, q95 = q
    lw1, lw2 = 1.1 * width, 3.0 * width
    seg = (lambda a, b, lw, c: ax.plot([pos, pos], [a, b], color=c, lw=lw, solid_capstyle="butt", zorder=z)) if orient == "v" \
        else (lambda a, b, lw, c: ax.plot([a, b], [pos, pos], color=c, lw=lw, solid_capstyle="butt", zorder=z))
    seg(q05, q95, lw1, tint(color, 0.55))
    seg(q25, q75, lw2, tint(color, 0.18))
    if dot:
        xy = ([pos], [q50]) if orient == "v" else ([q50], [pos])
        ax.scatter(*xy, s=dot_s, color=dot_c or "white", edgecolors=color, linewidths=0.7, zorder=z + 1)
    return q


def density(ax, v, color, lo=None, hi=None, bw=None, scale=None, base=0.0, orient="h", alpha=0.85, lw=0.7, n=256,
            outline=True, z=3, fill=True, weights=None):
    """Filled kernel density with a darker outline. orient 'h': density along y over x; 'v': along x over y.
    scale=None keeps the density's own height; otherwise the peak is scaled to that height above base."""
    from scipy.stats import gaussian_kde
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3 or np.ptp(v) == 0:
        return None
    lo = np.min(v) if lo is None else lo; hi = np.max(v) if hi is None else hi
    xs = np.linspace(lo, hi, n); d = gaussian_kde(v, bw_method=bw, weights=weights)(xs)
    if scale is not None:
        d = d / d.max() * scale
    if orient == "h":
        if fill:
            ax.fill_between(xs, base, base + d, color=tint(color, 1 - alpha) if alpha < 1 else color, lw=0, zorder=z)
        if outline:
            ax.plot(xs, base + d, color=color, lw=lw, zorder=z + 0.1)
    else:
        if fill:
            ax.fill_betweenx(xs, base, base + d, color=tint(color, 1 - alpha) if alpha < 1 else color, lw=0, zorder=z)
        if outline:
            ax.plot(base + d, xs, color=color, lw=lw, zorder=z + 0.1)
    return xs, d


def contour2d(ax, x, y, color, xlim, ylim, levels=(0.25, 0.5, 0.75, 0.9), n=120, bw=None, z=2, lines=True):
    """Highest-density regions of (x, y): filled tints of one hue, the densest darkest; a hairline at each level."""
    from scipy.stats import gaussian_kde
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    k = gaussian_kde(np.vstack([x[m], y[m]]), bw_method=bw)
    gx, gy = np.meshgrid(np.linspace(*xlim, n), np.linspace(*ylim, n))
    dz = k(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
    s = np.sort(dz.ravel())[::-1]; cs = np.cumsum(s) / s.sum()
    th = sorted(s[np.searchsorted(cs, lv)] for lv in levels)          # density thresholds, outermost first
    fills = [tint(color, f) for f in np.linspace(0.82, 0.25, len(th))]
    ax.contourf(gx, gy, dz, levels=list(th) + [dz.max() * 1.01], colors=fills, zorder=z, antialiased=True)
    if lines:
        ax.contour(gx, gy, dz, levels=th, colors=[color], linewidths=0.35, zorder=z + 0.1)


def hist_steps(ax, v, bins, color, alpha=0.9, density_=True, z=3, lw=0.7):
    """Histogram drawn as a filled step outline (no gaps between bins)."""
    h, e = np.histogram(v, bins=bins, density=density_)
    xs = np.repeat(e, 2)[1:-1]; ys = np.repeat(h, 2)
    ax.fill_between(xs, 0, ys, color=tint(color, 1 - alpha), lw=0, zorder=z)
    ax.plot(np.r_[e[0], xs, e[-1]], np.r_[0, ys, 0], color=color, lw=lw, zorder=z + 0.1)
    return h, e


# ---------------------------------------------------------------- alluvial
def ribbon(ax, x0, y0a, y0b, x1, y1a, y1b, color, alpha=0.55, z=2):
    """A flow from [y0a, y0b] at x0 to [y1a, y1b] at x1, with cubic edges."""
    xm = (x0 + x1) / 2
    verts = [(x0, y0a), (xm, y0a), (xm, y1a), (x1, y1a), (x1, y1b), (xm, y1b), (xm, y0b), (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), fc=tint(color, 1 - alpha), ec="none", zorder=z))


def alluvial(ax, rows, axes_cats, colors_of, xs=None, width=0.06, gap=0.012, label_fn=None, label_size=None,
             label_min=1, flow_color=None):
    """rows: list of tuples (one category per stratum) - one per unit. axes_cats: the category order of each stratum.
    colors_of(stratum_index, category) -> colour of the stratum block; flow_color(row) -> colour of the unit's ribbon.
    Draws blocks and the ribbons between neighbouring strata; returns block positions {(k, cat): (y0, y1)}."""
    from collections import Counter
    n = len(rows); ns = len(axes_cats)
    xs = xs or list(np.linspace(0, 1, ns))
    pos = {}
    for k, cats in enumerate(axes_cats):
        cnt = Counter(r[k] for r in rows); used = [c for c in cats if cnt.get(c, 0)]
        tot_gap = gap * (len(used) - 1); y = 1.0
        for c in used:
            h = cnt[c] / n * (1 - tot_gap)
            pos[(k, c)] = (y - h, y); y -= h + gap
    for k in range(ns - 1):
        top_l = {c: pos[(k, c)][1] for (kk, c) in pos if kk == k}
        top_r = {c: pos[(k + 1, c)][1] for (kk, c) in pos if kk == k + 1}
        order = sorted(range(n), key=lambda i: (axes_cats[k].index(rows[i][k]), axes_cats[k + 1].index(rows[i][k + 1])))
        # right side ordered by the left category so ribbons do not cross inside a block
        order_r = sorted(range(n), key=lambda i: (axes_cats[k + 1].index(rows[i][k + 1]), axes_cats[k].index(rows[i][k])))
        unit_h = (1 - gap * (len([c for c in axes_cats[k] if (k, c) in pos]) - 1)) / n
        unit_hr = (1 - gap * (len([c for c in axes_cats[k + 1] if (k + 1, c) in pos]) - 1)) / n
        yl, yr = {}, {}
        for i in order:
            c = rows[i][k]; yl[i] = (top_l[c] - unit_h, top_l[c]); top_l[c] -= unit_h
        for i in order_r:
            c = rows[i][k + 1]; yr[i] = (top_r[c] - unit_hr, top_r[c]); top_r[c] -= unit_hr
        for i in range(n):
            col = flow_color(rows[i]) if flow_color else colors_of(k, rows[i][k])
            ribbon(ax, xs[k] + width / 2, yl[i][1], yl[i][0], xs[k + 1] - width / 2, yr[i][1], yr[i][0], col, alpha=0.45)
    for (k, c), (y0, y1) in pos.items():
        ax.add_patch(Rectangle((xs[k] - width / 2, y0), width, y1 - y0, fc=colors_of(k, c), ec="white", lw=0.4, zorder=3))
    ax.set_xlim(xs[0] - width, xs[-1] + width); ax.set_ylim(-0.01, 1.01); ax.axis("off")
    return pos


# ---------------------------------------------------------------- survival
def km_risk(canvas, x, y, w, h, df, col, tlab, color, xmax, xticks, xlabel, risk_h=7.0, show_y=True, unit_lab=None):
    """Kaplan-Meier curves for patients above (dark) and below (light) the median of df[col], the log-rank P, and the
    numbers at risk beneath."""
    ax = canvas.ax(x, y, w, h)
    hi = df[col].to_numpy() > np.median(df[col].to_numpy())
    lo_c = tint(color, 0.5)
    km(ax, df.time[hi], df.event[hi], color, lw=1.0)
    km(ax, df.time[~hi], df.event[~hi], lo_c, lw=1.0)
    p = logrank_p(df.time, df.event, hi)
    ax.set_xlim(0, xmax); ax.set_xticks(xticks); ax.set_ylim(0, 1.04); ax.set_yticks([0, 0.5, 1])
    ax.set_xticklabels([f"{t_:g}" for t_ in xticks])
    ax.set_xlabel(xlabel)
    if show_y:
        ax.set_ylabel("event-free")
    else:
        ax.set_yticklabels([])
    ax.text(0.04, 0.05, f"log-rank P = {sci(p)}", transform=ax.transAxes, fontsize=cns.T_SMALL, color=cns.INK)
    ax.text(0.97, 0.98, tlab, transform=ax.transAxes, ha="right", va="top", fontsize=cns.T_NOTE, color=color,
            fontweight="bold")
    for j_, (cc, ll) in enumerate(((color, "above median"), (lo_c, "below"))):
        yy = 0.30 - j_ * 0.11
        ax.plot([0.04, 0.12], [yy, yy], color=cc, lw=1.1, transform=ax.transAxes)
        ax.text(0.15, yy, ll, transform=ax.transAxes, va="center", fontsize=cns.T_SMALL)
    frame(ax)
    # numbers at risk
    rk = canvas.ax(x, y + h + 9.0, w, risk_h); rk.set_xlim(0, xmax); rk.set_ylim(-0.5, 1.5); rk.axis("off")
    t = df.time.to_numpy()
    for j_, (m, cc) in enumerate(((hi, color), (~hi, lo_c))):
        for tt in xticks:
            rk.text(tt, 1 - j_, str(int((t[m] >= tt).sum())), ha="center", va="center", fontsize=cns.T_SMALL - 0.3,
                    color=cc if j_ == 0 else cns.INK2)
    if show_y:
        rk.text(-0.10 * xmax, 1.9, "at risk", ha="right", va="center", fontsize=cns.T_SMALL, color=cns.INK2)
        for j_, ll in enumerate(("above", "below")):
            rk.text(-0.10 * xmax, 1 - j_, ll, ha="right", va="center", fontsize=cns.T_SMALL - 0.3, color=cns.INK2)
    return ax, p


# ---------------------------------------------------------------- misc
def cbar(canvas, x, y, w, h, cmap, vmin, vmax, ticks, label=None, horizontal=True, fmt=None):
    from matplotlib.colors import Normalize
    cax = canvas.ax(x, y, w, h)
    cb = canvas.fig.colorbar(plt.cm.ScalarMappable(Normalize(vmin, vmax), cmap), cax=cax, ticks=ticks,
                             orientation="horizontal" if horizontal else "vertical")
    cb.outline.set_linewidth(0.3)
    cax.tick_params(labelsize=cns.T_SMALL, length=1.2, pad=1, width=0.4)
    if fmt:
        cb.set_ticklabels([fmt(t) for t in ticks])
    if label:
        if horizontal:
            cax.set_title(label, fontsize=cns.T_SMALL, pad=2)
        else:
            cax.set_title(label, fontsize=cns.T_SMALL, pad=2, loc="left")
    return cax


def clean(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def sci(p):
    """P value as text: two significant digits, or a power of ten below 0.001."""
    if p >= 0.001:
        return f"{p:.2g}" if p < 0.01 else f"{p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return f"{m}×10$^{{{int(e)}}}$"


def forest_bands(ax, calls, stat, fam=None, cls="coex", maxlab=4.2, color_default=None, stripes=True):
    """Each target's median delta (dot; filled if called) against its decoys' nested 5-95% and 25-75% bands and median
    (tick), one column per target sorted by the target's median delta."""
    c = calls[calls.stat == stat]
    t = c[c.kind == "target"].set_index("target")
    d = c[(c.kind == "decoy") & (c["class"] == cls)]
    qd = d.groupby("target").med_delta.quantile([.05, .25, .5, .75, .95]).unstack()
    order = t.med_delta.sort_values(ascending=False).index
    for i, s in enumerate(order):
        if stripes and i % 2 == 0:
            ax.axvspan(i - .5, i + .5, color=cns.GRID, lw=0, zorder=0)
        if s in qd.index:
            band(ax, i, q=qd.loc[s].to_numpy(), color=cns.AQUA, width=0.9, dot=False)
            ax.plot([i - .28, i + .28], [qd.loc[s, .5]] * 2, color=cns.AQUA, lw=0.8, zorder=4)
        col = cns.FAM.get(fam.get(s), cns.GRAY) if fam is not None else (color_default or cns.INK)
        dots(ax, [i], [t.loc[s, "med_delta"]], col, s=12 if t.loc[s, "call"] else 9, filled=bool(t.loc[s, "call"]), z=5)
    ref(ax, 0)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([lab(s) for s in order], rotation=90, fontsize=maxlab)
    ax.tick_params(axis="x", length=0, pad=1); ax.set_xlim(-0.6, len(order) - 0.4)
    return order
