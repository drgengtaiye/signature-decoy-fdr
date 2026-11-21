"""Shared style of the earlier figure set (cns10.py builds on it; cohort labels and the log-rank P used by the figures).

180 mm canvas, axes placed in millimetres from the top-left; sans-serif type at 5-8 pt; hairline axes; marks separated
by a white gap or ring instead of outlines; no hatching. Colour does one job per figure:

  family       inflamed/IFN BLUE, confounder ORANGE, other GRAY (de-emphasis)
  unit         signature INK, co-expressed decoy AQUA, random-gene decoy LGRAY
  statistic    folded SLATE, directed BLUE
  estimate     filled marker = plain estimate, open marker = pi0-corrected, black = realized
  planted      signal in the library BLUE, elsewhere in signatures / mixed ORANGE, outside the library AQUA

The categorical hues are slots 1-3 of the reference palette (blue, orange, aqua), validated all-pairs (worst CVD
deltaE 9.2, normal-vision 24.0); aqua sits below 3:1 on white, so it always carries a legend or a direct label.
"""
import sys
sys.path.insert(0, "/path/to/analysis_root")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import figstyle as fs

D = "/path/to/analysis_root"
OUTDIR = f"{D}/figs_pub"
MM = 1 / 25.4
W = 180.0

INK, INK2, MUTED = "#1a1f24", "#52514e", "#8a9097"
RULE, GRID = "#c4c9ce", "#eceff1"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GRAY, LGRAY, SLATE = "#9aa3ac", "#c9cfd5", "#5f6b76"
FAM = {"Inflamed/IFN": BLUE, "Confounder": ORANGE, "Other": GRAY}
FAMLAB = {"Inflamed/IFN": "inflamed/IFN", "Confounder": "confounder", "Other": "other"}
STAT = {"fold": SLATE, "dir": BLUE}
STATLAB = {"fold": "folded", "dir": "directed"}

T_TICK, T_LAB, T_NOTE, T_SMALL, T_HEAD, T_LET = 5.5, 6.0, 5.5, 4.9, 6.3, 8.0


def tint(c, f=0.5):
    """Lighter step of the same hue (mixed with white)."""
    r, g, b = to_rgb(c)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def apply():
    fs.apply()
    plt.rcParams.update({
        "font.size": T_NOTE, "axes.labelsize": T_LAB, "axes.titlesize": T_HEAD,
        "xtick.labelsize": T_TICK, "ytick.labelsize": T_TICK, "legend.fontsize": T_NOTE,
        "axes.linewidth": 0.45, "xtick.major.width": 0.45, "ytick.major.width": 0.45,
        "xtick.minor.width": 0.35, "ytick.minor.width": 0.35,
        "xtick.major.size": 2.0, "ytick.major.size": 2.0, "xtick.minor.size": 1.2, "ytick.minor.size": 1.2,
        "axes.labelpad": 1.8, "xtick.major.pad": 1.4, "ytick.major.pad": 1.4,
        "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
        "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
    })


class Canvas:
    """A 180 mm wide figure with axes placed in millimetres from the top-left."""

    def __init__(self, h_mm, w_mm=W):
        self.w, self.h = w_mm, h_mm
        self.fig = plt.figure(figsize=(w_mm * MM, h_mm * MM))

    def ax(self, x, y, w, h, **kw):
        a = self.fig.add_axes([x / self.w, 1 - (y + h) / self.h, w / self.w, h / self.h], **kw)
        frame(a)
        return a

    def letter(self, x, y, s):
        self.fig.text(x / self.w, 1 - y / self.h, s, fontsize=T_LET, fontweight="bold", va="top", ha="left", color=INK)

    def head(self, x, y, s, size=T_HEAD):
        """Short panel header (a few words), set beside the panel letter."""
        self.fig.text(x / self.w, 1 - y / self.h, s, fontsize=size, fontweight="bold", va="top", ha="left", color=INK)

    def text(self, x, y, s, size=T_NOTE, color=None, **kw):
        kw.setdefault("va", "top"); kw.setdefault("ha", "left")
        return self.fig.text(x / self.w, 1 - y / self.h, s, fontsize=size, color=color or INK, **kw)

    def panel(self, letter, x, y, w, h, head=None, **kw):
        """Axes at (x, y, w, h) mm with its letter 6 mm left of and 4 mm above the axes, and an optional header."""
        a = self.ax(x, y, w, h, **kw)
        self.letter(x - 7.5, y - 5.2, letter)
        if head:
            self.head(x - 2.5, y - 5.0, head)
        return a

    def save(self, stem):
        for ext in ("pdf", "png"):
            self.fig.savefig(f"{OUTDIR}/{stem}.{ext}", format=ext, dpi=450)
        plt.close(self.fig)


def frame(ax, left=True, bottom=True):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left); ax.spines["bottom"].set_visible(bottom)
    if not left:
        ax.tick_params(axis="y", length=0)
    if not bottom:
        ax.tick_params(axis="x", length=0)


def ref(ax, v, axis="h", label=None, color=RULE, lw=0.6, pos=0.98, ha="right", va="bottom", z=1):
    """Solid hairline reference (identity, nominal rate, zero), optionally labelled at its end."""
    if axis == "h":
        ax.axhline(v, color=color, lw=lw, zorder=z)
        if label:
            ax.text(pos, v, label, transform=ax.get_yaxis_transform(), ha=ha, va=va, fontsize=T_SMALL, color=INK2)
    else:
        ax.axvline(v, color=color, lw=lw, zorder=z)
        if label:
            ax.text(v, pos, label, transform=ax.get_xaxis_transform(), ha="left", va="top", fontsize=T_SMALL, color=INK2)


def identity(ax, lo, hi, color=RULE, lw=0.6):
    ax.plot([lo, hi], [lo, hi], color=color, lw=lw, zorder=1)


def dots(ax, x, y, c, s=12, filled=True, marker="o", ring=True, z=4, alpha=1.0, lw=0.7, **kw):
    """Markers with a white surface ring; open markers keep the hue on the edge."""
    if filled:
        return ax.scatter(x, y, s=s, c=[c] if isinstance(c, str) or (isinstance(c, tuple) and len(c) in (3, 4)) else c,
                          marker=marker, edgecolors="white" if ring else "none", linewidths=0.35 if ring else 0,
                          zorder=z, alpha=alpha, **kw)
    return ax.scatter(x, y, s=s, facecolors="white", edgecolors=c, linewidths=lw, marker=marker, zorder=z,
                      alpha=alpha, **kw)


def hbar_interval(ax, y, lo, hi, c, lw=0.9, z=3):
    ax.plot([lo, hi], [y, y], color=c, lw=lw, solid_capstyle="butt", zorder=z)


def vbar_interval(ax, x, lo, hi, c, lw=0.9, z=3):
    ax.plot([x, x], [lo, hi], color=c, lw=lw, solid_capstyle="butt", zorder=z)


def legend(ax, items, x, y, dy=0.075, size=T_NOTE, s=10, xpad=0.035, transform=None):
    """Compact legend: items are (colour, label) or (colour, label, marker, filled)."""
    tr = transform or ax.transAxes
    for i, it in enumerate(items):
        c, lab = it[0], it[1]
        mk = it[2] if len(it) > 2 else "o"
        filled = it[3] if len(it) > 3 else True
        yy = y - i * dy
        if mk == "line":
            ax.plot([x - 0.012, x + 0.018], [yy, yy], color=c, lw=0.9, transform=tr, clip_on=False, zorder=10)
        elif mk == "box":
            ax.scatter([x], [yy], s=s * 1.3, marker="s", c=[c], linewidths=0, transform=tr, clip_on=False, zorder=10)
        else:
            if filled:
                ax.scatter([x], [yy], s=s, marker=mk, c=[c], edgecolors="white", linewidths=0.3, transform=tr,
                           clip_on=False, zorder=10)
            else:
                ax.scatter([x], [yy], s=s, marker=mk, facecolors="white", edgecolors=c, linewidths=0.7, transform=tr,
                           clip_on=False, zorder=10)
        ax.text(x + xpad, yy, lab, transform=tr, va="center", ha="left", fontsize=size, color=INK, zorder=10)


def stacked_calls(ax, x, called, exp_false, c, width=0.62, text=None):
    """A column of calls split into expected true (solid) and expected false (light tint), white gap between."""
    et = max(called - exp_false, 0)
    ax.bar(x, et, width, color=c, edgecolor="white", linewidth=0.6, zorder=3)
    ax.bar(x, exp_false, width, bottom=et, color=tint(c, 0.62), edgecolor="white", linewidth=0.6, zorder=3)
    if text:
        ax.text(x, called + 0.6, text, ha="center", va="bottom", fontsize=T_SMALL, color=INK, linespacing=1.0)


def swarm(v, width=0.32, nb=16):
    v = np.asarray(v, float); off = np.zeros(len(v))
    if len(v) < 2:
        return off
    bins = np.linspace(np.nanmin(v), np.nanmax(v) + 1e-9, nb)
    idx = np.digitize(v, bins)
    for b in np.unique(idx):
        k = np.where(idx == b)[0]
        k = k[np.argsort(v[k])]
        steps = np.array([0] + [s for i in range(1, len(k)) for s in (i, -i)])[:len(k)]
        off[k] = steps * width / max(4, nb / 3)
    return np.clip(off, -width, width)


SHORT = {"TUMOR-VASCULATURE-UP_Lu": "Tumor vasculature (Lu)", "CD39-CD8Tcell_Chow": "CD39 CD8 (Chow)",
         "Tcell-inflamed_Ayers": "T-cell-inflamed (Ayers)", "WNT-beta-catenin_Takeuchi": "WNT/β-catenin (Takeuchi)"}


def lab(s):
    if s in SHORT:
        return SHORT[s]
    if "_" in s:
        a, b = s.rsplit("_", 1)
        return f"{a.replace('_', ' ')} ({b})"
    return s


def pct(v, nd=0):
    return f"{v:.{nd}f}%"


def minus(s):
    return s.replace("-", "−")


COHORT = {"ICB_Braun": "Braun", "ICB_Gide": "Gide", "ICB_Hugo": "Hugo", "ICB_Jung": "Jung", "ICB_Kim": "Kim",
          "ICB_Limagne1": "Limagne", "ICB_Liu": "Liu", "ICB_Mariathasan": "Mariathasan", "ICB_Miao1": "Miao",
          "ICB_Puch": "Puch", "ICB_Riaz": "Riaz", "ICB_Van_Allen": "Van Allen"}


def hlegend(canvas, x, y, items, size=T_NOTE, s=9, gap=3.0):
    """A one-line legend in figure millimetres: items are (colour, label) or (colour, label, marker, filled)."""
    fig = canvas.fig
    ax = fig.add_axes([0, 0, 1, 1], facecolor="none"); ax.set_xlim(0, canvas.w); ax.set_ylim(canvas.h, 0); ax.axis("off")
    ax.set_zorder(20)
    r = fig.canvas.get_renderer() if hasattr(fig.canvas, "get_renderer") else None
    xx = x
    for it in items:
        c, lab = it[0], it[1]
        mk = it[2] if len(it) > 2 else "o"
        filled = it[3] if len(it) > 3 else True
        if mk == "line":
            ax.plot([xx - 1.2, xx + 1.4], [y, y], color=c, lw=0.9)
        elif filled:
            ax.scatter([xx], [y], s=s, marker=mk, c=[c], edgecolors="white", linewidths=0.3)
        else:
            ax.scatter([xx], [y], s=s, marker=mk, facecolors="white", edgecolors=c, linewidths=0.7)
        t = ax.text(xx + 1.8, y, lab, fontsize=size, va="center", ha="left", color=INK)
        w_mm = len(lab) * size * 0.19 + 1.0          # approximate width of the label in mm at this size
        xx += 1.8 + w_mm + gap
    return ax

from matplotlib.colors import LinearSegmentedColormap as _LSC
DIV = _LSC.from_list("div", [BLUE, "#9fc2ea", "#f3f3f1", "#f0a58a", "#d6452c"])


def key(ax, items, x, y, dy=.085, size=T_NOTE, marker="o", s=10, xpad=.035):
    """Compatibility with v4kit.key: (colour, label) items with one marker shape."""
    legend(ax, [(c, lab, "box" if marker == "s" else marker) for c, lab in items], x, y, dy=dy, size=size, s=s, xpad=xpad)


def forest_columns(ax, calls, stat, fam=None, cls="coex", maxlab=4.2, color_default=None):
    """Each target's median delta (dot; filled if called) against its decoys' 5-95% range (bar) and median (tick),
    one column per target sorted by the target's median delta."""
    c = calls[calls.stat == stat]
    t = c[c.kind == "target"].set_index("target")
    d = c[(c.kind == "decoy") & (c["class"] == cls)]
    qd = d.groupby("target").med_delta.quantile([.05, .5, .95]).unstack()
    order = t.med_delta.sort_values(ascending=False).index
    for i, s in enumerate(order):
        if s in qd.index:
            vbar_interval(ax, i, qd.loc[s, .05], qd.loc[s, .95], tint(AQUA, .45), lw=2.0)
            ax.plot([i - .3, i + .3], [qd.loc[s, .5]] * 2, color=AQUA, lw=0.7, zorder=3)
        col = FAM.get(fam.get(s), GRAY) if fam is not None else (color_default or INK)
        dots(ax, [i], [t.loc[s, "med_delta"]], col, s=11 if t.loc[s, "call"] else 8, filled=bool(t.loc[s, "call"]), z=5)
    ref(ax, 0)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([lab(s) for s in order], rotation=90, fontsize=maxlab)
    ax.tick_params(axis="x", length=0, pad=1); ax.set_xlim(-0.8, len(order) - 0.2)
    return order


def globalnull_hist(ax, per_perm, stat="fold", summary=None):
    """Realized false calls per permutation (bars) with the realized and estimated means (lines)."""
    x = per_perm[per_perm.stat == stat]
    mx = int(min(max(np.ceil(np.percentile(x.R, 99)), 4), 14))
    bins = np.arange(-0.5, mx + 1.5, 1)
    ax.hist(x.R, bins=bins, color=LGRAY, edgecolor="white", linewidth=0.5, density=True, zorder=2)
    ax.axvline(x.R.mean(), color=INK, lw=1.0, zorder=3)
    ax.axvline(x.E.mean(), color=AQUA, lw=1.0, zorder=4)
    ax.set_xlim(-0.6, mx + 0.6)
    ax.set_xlabel("false calls per library"); ax.set_ylabel("permutations (share)")
    if summary is not None:
        r = summary[summary.stat == stat].iloc[0]
        ax.text(0.97, 0.97, f"ratio {r.ratio:.2f}\n({r.ratio_lo:.2f}–{r.ratio_hi:.2f})", transform=ax.transAxes,
                ha="right", va="top", fontsize=T_SMALL, linespacing=1.1)


# ---------------------------------------------------------------- richer forms
def km_curve(time, event):
    """Kaplan-Meier estimate: step times and survival."""
    t = np.asarray(time, float); e = np.asarray(event, int)
    o = np.argsort(t); t, e = t[o], e[o]
    ut = np.unique(t[e == 1]); s = 1.0; xs, ys = [0.0], [1.0]
    n_at = len(t)
    for u in ut:
        d = int(((t == u) & (e == 1)).sum()); n_risk = int((t >= u).sum())
        s *= 1 - d / n_risk; xs.append(u); ys.append(s)
    return np.array(xs), np.array(ys)


def km(ax, time, event, color, ls="-", lw=1.0, label=None, censor=True):
    xs, ys = km_curve(time, event)
    ax.step(xs, ys, where="post", color=color, lw=lw, ls=ls, label=label, zorder=3)
    if censor:
        t = np.asarray(time, float); e = np.asarray(event, int)
        ct = t[e == 0]
        if len(ct):
            sy = np.interp(ct, xs, ys, left=1.0)
            idx = np.searchsorted(xs, ct, side="right") - 1
            ax.scatter(ct, ys[np.clip(idx, 0, len(ys) - 1)], marker="|", s=6, color=color, lw=0.5, zorder=4)


def logrank_p(time, event, group):
    """Two-group log-rank P (chi-square, 1 df)."""
    from scipy.stats import chi2
    t = np.asarray(time, float); e = np.asarray(event, int); g = np.asarray(group, bool)
    o_e = 0.0; v = 0.0
    for u in np.unique(t[e == 1]):
        at = t >= u; n = at.sum(); n1 = (at & g).sum()
        d = ((t == u) & (e == 1)).sum(); d1 = ((t == u) & (e == 1) & g).sum()
        if n > 1:
            o_e += d1 - d * n1 / n
            v += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    return float(chi2.sf(o_e ** 2 / v, 1)) if v > 0 else 1.0


def violin(ax, v, pos, color, width=0.36, side="both", alpha=0.85, bw=None, points=False, z=3):
    """KDE violin (or half-violin) with median tick and interquartile bar."""
    from scipy.stats import gaussian_kde
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3 or np.ptp(v) == 0:
        return
    k = gaussian_kde(v, bw_method=bw)
    ys = np.linspace(v.min(), v.max(), 120); dens = k(ys); dens = dens / dens.max() * width
    if side in ("both", "left"):
        ax.fill_betweenx(ys, pos - dens, pos, color=color, alpha=alpha, lw=0, zorder=z)
    if side in ("both", "right"):
        ax.fill_betweenx(ys, pos, pos + dens, color=color, alpha=alpha, lw=0, zorder=z)
    q1, md, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([pos, pos], [q1, q3], color=INK, lw=1.4, solid_capstyle="butt", zorder=z + 1)
    ax.scatter([pos], [md], s=6, color="white", edgecolors=INK, linewidths=0.4, zorder=z + 2)
    if points:
        ax.scatter(pos + np.random.default_rng(1).uniform(-width * 0.25, width * 0.25, len(v)), v, s=1.2, color=INK, alpha=0.4,
                   lw=0, zorder=z + 1)


def ridge(ax, y0, v, color, height=0.8, alpha=0.75, lo=None, hi=None, bw=None, lw=0.5, fill=True, z=3, line=None):
    """One ridgeline at baseline y0 (density scaled to height); fill=False draws the outline only."""
    from scipy.stats import gaussian_kde
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3:
        return
    lo = np.min(v) if lo is None else lo; hi = np.max(v) if hi is None else hi
    xs = np.linspace(lo, hi, 200); d = gaussian_kde(v, bw_method=bw)(xs); d = d / d.max() * height
    if fill:
        ax.fill_between(xs, y0, y0 + d, color=color, alpha=alpha, lw=0, zorder=z)
        ax.plot(xs, y0 + d, color=line or "white", lw=lw, zorder=z + 0.01)
    else:
        ax.fill_between(xs, y0, y0 + d, color="white", alpha=0.0, lw=0, zorder=z)
        ax.plot(xs, y0 + d, color=color, lw=lw, zorder=z + 0.01)


def marginal(canvas, x, y, w, h, xdata_list, ydata_list, colors, size=4.0, xlim=None, ylim=None):
    """Marginal densities above (x) and to the right (y) of a main axes at (x, y, w, h) mm."""
    from scipy.stats import gaussian_kde
    top = canvas.fig.add_axes([x / canvas.w, 1 - y / canvas.h, w / canvas.w, size / canvas.h]); top.axis("off")
    rgt = canvas.fig.add_axes([(x + w) / canvas.w, 1 - (y + h) / canvas.h, size / canvas.w, h / canvas.h]); rgt.axis("off")
    for xd, yd, c in zip(xdata_list, ydata_list, colors):
        xd = np.asarray(xd, float); yd = np.asarray(yd, float)
        gx = np.linspace(*(xlim or (np.nanmin(xd), np.nanmax(xd))), 150)
        gy = np.linspace(*(ylim or (np.nanmin(yd), np.nanmax(yd))), 150)
        if len(xd) > 3:
            top.fill_between(gx, 0, gaussian_kde(xd[np.isfinite(xd)])(gx), color=c, alpha=0.45, lw=0)
        if len(yd) > 3:
            rgt.fill_betweenx(gy, 0, gaussian_kde(yd[np.isfinite(yd)])(gy), color=c, alpha=0.45, lw=0)
    if xlim: top.set_xlim(*xlim)
    if ylim: rgt.set_ylim(*ylim)
    top.set_ylim(bottom=0); rgt.set_xlim(left=0)
    return top, rgt


def km_pair(canvas, x, y, w, h, df, gap=6.0, xmax=None, xticks=None, xlabel="time", tlab="signature", dlab="decoy",
            tcol=None, dcol=None, leg_y=0.30):
    """Two Kaplan-Meier axes side by side: patients split at the median of the signature's score (left) and of one of
    its co-expressed decoys' (right); the log-rank P is printed for orientation."""
    tcol = tcol or BLUE; dcol = dcol or AQUA
    w1 = (w - gap) / 2; axes = []
    for k, (col, c, lab_) in enumerate((("target", tcol, tlab), ("decoy", dcol, dlab))):
        ax = canvas.ax(x + k * (w1 + gap), y, w1, h)
        hi = df[col].to_numpy() > np.median(df[col].to_numpy())
        km(ax, df.time[hi], df.event[hi], c, lw=1.0)
        km(ax, df.time[~hi], df.event[~hi], tint(c, 0.45), lw=1.0)
        p = logrank_p(df.time, df.event, hi)
        ax.text(0.97, 0.97, lab_, transform=ax.transAxes, ha="right", va="top", fontsize=T_SMALL, color=INK)
        for j_, (cc_, ll_) in enumerate(((c, "above median"), (tint(c, 0.45), "below"))):
            yy_ = leg_y - j_ * 0.11
            ax.plot([0.55, 0.63], [yy_, yy_], color=cc_, lw=1.0, transform=ax.transAxes)
            ax.text(0.66, yy_, ll_, transform=ax.transAxes, va="center", fontsize=T_SMALL, color=INK)
        ax.text(0.03, 0.04, f"log-rank P = {p:.1g}" if p < 0.001 else f"log-rank P = {p:.2f}", transform=ax.transAxes,
                fontsize=T_SMALL, color=INK2)
        ax.set_ylim(0, 1.03); ax.set_yticks([0, 0.5, 1])
        if xmax: ax.set_xlim(0, xmax)
        if xticks is not None: ax.set_xticks(xticks)
        ax.set_xlabel(xlabel)
        if k == 0:
            ax.set_ylabel("event-free")
        else:
            ax.set_yticklabels([])
        frame(ax)
        axes.append((ax, p))
    return axes
