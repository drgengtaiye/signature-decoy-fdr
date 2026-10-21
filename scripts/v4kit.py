"""Shared layout kit for the v4 figure set (journal-standard redesign).

One message per figure, four to six panels, statistics in the legends rather than on
the panels, type never below 5 pt at print size. Axes are placed in millimetres from the
top-left of a 180 mm canvas so rows align exactly across panels.
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric")
import numpy as np
import matplotlib.pyplot as plt
import figstyle as fs

D = "/groups/adv2105_gp/yichen/Yi/com/gastric"
MM = fs.MM
W = 180.0

# type sizes at print scale (pt)
T_TICK, T_LAB, T_NOTE, T_SMALL, T_TITLE, T_LET = 6.0, 6.6, 5.8, 5.3, 6.9, 9.0


def apply():
    fs.apply()
    plt.rcParams.update({
        "font.size": T_NOTE, "axes.labelsize": T_LAB, "axes.titlesize": T_TITLE,
        "xtick.labelsize": T_TICK, "ytick.labelsize": T_TICK, "legend.fontsize": T_NOTE,
        "axes.linewidth": 0.55, "xtick.major.width": 0.55, "ytick.major.width": 0.55,
        "xtick.major.size": 2.4, "ytick.major.size": 2.4,
        "axes.labelpad": 2.0, "xtick.major.pad": 1.6, "ytick.major.pad": 1.6,
    })


class Canvas:
    """A 180 mm wide figure with axes placed in millimetres from the top-left."""

    def __init__(self, h_mm, w_mm=W):
        self.w, self.h = w_mm, h_mm
        self.fig = plt.figure(figsize=(w_mm * MM, h_mm * MM))

    def ax(self, x, y, w, h, **kw):
        return self.fig.add_axes([x / self.w, 1 - (y + h) / self.h, w / self.w, h / self.h], **kw)

    def letter(self, x, y, s):
        """Bold panel letter with its top-left corner at (x, y) mm."""
        self.fig.text(x / self.w, 1 - y / self.h, s, fontsize=T_LET, fontweight="bold",
                      va="top", ha="left", color=fs.INK)

    def text(self, x, y, s, size=T_NOTE, color=None, **kw):
        kw.setdefault("va", "top"); kw.setdefault("ha", "left")
        return self.fig.text(x / self.w, 1 - y / self.h, s, fontsize=size,
                             color=color or fs.INK, **kw)

    def title(self, x, y, s, size=T_TITLE, color=None, **kw):
        return self.text(x, y, s, size=size, color=color or fs.INK, **kw)

    def save(self, stem):
        fs.savefig(self.fig, stem)
        plt.close(self.fig)


def frame(ax, left=True, bottom=True):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left); ax.spines["bottom"].set_visible(bottom)


def minus(s):
    """Typographic minus for numbers placed as text."""
    return s.replace("-", "−")


def swarm(v, width=.30, nb=18):
    """Deterministic beeswarm offsets for a one-dimensional sample."""
    v = np.asarray(v, float); off = np.zeros(len(v))
    if len(v) < 2:
        return off
    bins = np.linspace(v.min(), v.max() + 1e-9, nb)
    idx = np.digitize(v, bins)
    for b in np.unique(idx):
        k = np.where(idx == b)[0]
        k = k[np.argsort(v[k])]
        steps = np.array([0] + [s for i in range(1, len(k)) for s in (i, -i)])[:len(k)]
        off[k] = steps * width / max(4, nb / 3)
    return np.clip(off, -width, width)


def key(ax, items, x, y, dy=.085, size=T_NOTE, marker="o", s=12, xpad=.035):
    """Compact key of (colour, label) in axes coordinates."""
    for i, (c, lab) in enumerate(items):
        yy = y - i * dy
        ax.scatter([x], [yy], s=s, c=[c], marker=marker, lw=0, transform=ax.transAxes,
                   clip_on=False, zorder=10)
        ax.text(x + xpad, yy, lab, transform=ax.transAxes, va="center", ha="left",
                fontsize=size, color=fs.INK, zorder=10)


SHORT = {"TUMOR-VASCULATURE-UP_Lu": "Tumor vasculature Lu",
         "CD39-CD8Tcell_Chow": "CD39-CD8 Chow",
         "Teff-IFNG_Fehrenbacher": "Teff-IFNG Fehrenbacher",
         "Tcell-inflamed_Ayers": "T-cell-inflamed Ayers",
         "WNT-beta-catenin_Takeuchi": "WNT/β-catenin Takeuchi"}


def lab(s):
    return SHORT.get(s, s.replace("_", " "))
