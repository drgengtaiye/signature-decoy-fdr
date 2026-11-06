"""The look of the AML multi-omic figures (Yi/pixlgen/integration/cns_figures), for the main figures of the decoy paper.

  page     180 mm wide, at most 240 mm tall, everything placed in mm from the top-left; panels small and many
  type     real Helvetica; letter 8 bold, title 6.5 bold, axis 6, tick 5.5, everything else 5; nothing under 5 pt
  colour   the two sides of the paper in the two project hues: a signature or published model is a red circle, a co-expressed
           decoy a navy square; a random-gene decoy an orange diamond; a null bank is grey. Filled = called or validated, open =
           not. Family: inflamed/IFN teal, confounder purple, other grey. Ramps in one hue (BLUE, PURPLE) or two opposing
           hues centred on the reference (PUOR, TEALOR); never a rainbow, never grey tiles
  matrix   tiles parted by white hairlines with the value printed inside in white or ink
  values   a point, not a bar: on a dotted reference, with n and values as plain columns at the right
  bars     only for histograms, in grey
"""
from __future__ import annotations
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath

D = "/groups/adv2105_gp/yichen/Yi/com/gastric"
OUT = os.environ.get("PIXL_OUT", f"{D}/figs_pub/pixl")
os.makedirs(OUT, exist_ok=True)

# ---- Font: the real faces only. ~/.fonts holds a different cut under the same family name, so it is dropped in memory.
FACES = "/groups/adv2105_gp/yichen/font/helvetica_faces/Helvetica-%s.ttf"
fm.fontManager.ttflist = [f for f in fm.fontManager.ttflist if f.name != "Helvetica"]
for face in ("Regular", "Bold", "Oblique", "BoldOblique"):
    fm.fontManager.addfont(FACES % face)
for props, face in (({}, "Regular"), ({"weight": "bold"}, "Bold"), ({"style": "italic"}, "Oblique")):
    found = fm.findfont(fm.FontProperties(family="Helvetica", **props), fallback_to_default=False)
    assert found.endswith(f"helvetica_faces/Helvetica-{face}.ttf"), (props, found)

INK, MUTED, OTHER, EMPTY = "#1a1a1a", "#595959", "#a6a6a6", "#f0f0f0"
plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 5, "axes.linewidth": 0.5, "axes.labelsize": 6, "axes.labelpad": 2,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 2,
    "xtick.minor.size": 1.2, "ytick.minor.size": 1.2, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "axes.spines.top": False, "axes.spines.right": False, "text.hinting": "no_hinting", "pdf.fonttype": 42,
    "ps.fonttype": 42, "svg.fonttype": "none", "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK,
    "ytick.color": INK, "axes.edgecolor": INK, "legend.frameon": False, "legend.fontsize": 5, "axes.unicode_minus": True,
    "mathtext.fontset": "custom", "mathtext.rm": "Helvetica", "mathtext.it": "Helvetica:italic",
    "mathtext.bf": "Helvetica:bold", "mathtext.bfit": "Helvetica:bold:italic", "mathtext.default": "regular",
})

# ---- Roles. Every set checked with check_palette below: OKLab delta E x100 >= 15 in normal vision, >= 8 under CVD simulation.
SIG, COEX, RAND, NULLC = "#E7213A", "#3C5488", "#E18727", "#BDBDBD"     # signature/model, co-expressed decoy, random-gene decoy, null bank
UNIT = {"sig": SIG, "coex": COEX, "rand": RAND, "null": NULLC}
MARK = {"sig": "o", "coex": "s", "rand": "D", "null": "o"}
INFL, CONF, OTH = "#00A087", "#9B6FD6", "#BDBDBD"                       # family
FAM = {"Inflamed/IFN": INFL, "Confounder": CONF, "Other": OTH}
EST, NOM, PI0 = "#00A087", "#BDBDBD", "#8A6100"                         # decoy estimate, nominal, pi0-corrected (where drawn together)
SCEN = {"P1": "#00A087", "P2": "#E3B341", "P3": "#7A0F3C"}              # planted scenarios in Fig. 5
TYPE = {"melanoma": "#4A3B8F", "kidney": "#C2549C", "lung": "#3B8ED0", "bladder": "#E3B341", "gastric": "#8A6100"}
# Ramps
BLUE = LinearSegmentedColormap.from_list("agree", ["#EEF4FA", "#9CC3E0", "#1F6FAE", "#08306B"])
PURPLE = LinearSegmentedColormap.from_list("seq", ["#F4F1F8", "#B9B0DC", "#6A51A3", "#3F007D"])
PUOR = LinearSegmentedColormap.from_list("puor0", ["#3F1F87", "#8C74BE", "#F7F7F7", "#E89A3C", "#9C4A00"])
TEALOR = LinearSegmentedColormap.from_list("div", ["#00596B", "#4FA8B8", "#F5F5F5", "#D88A3D", "#8C3D00"])
GREYHEX = LinearSegmentedColormap.from_list("grey_seq", ["#E4E4E4", "#8C8C8C", "#111111"])


def ink_on(fc):
    """white or ink, whichever has the better contrast on this fill"""
    lin = lambda c: c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4  # noqa: E731
    lum = sum(w * lin(c) for w, c in zip((0.2126, 0.7152, 0.0722), to_rgb(fc)))
    return "white" if 1.05 / (lum + 0.05) > (lum + 0.05) / (0.0103 + 0.05) else INK


def tint(c, f):
    r = np.array(to_rgb(c))
    return tuple(1 - (1 - r) * f)


class Page:
    """A sheet whose content is placed in mm from the top-left. pg.src(letter, frame) collects the source data of a panel."""

    def __init__(self, H, W=180.0):
        assert W <= 180.01 and H <= 240.01, "Nature Communications page: at most 180 x 240 mm"
        self.W, self.H = W, H
        self.fig = plt.figure(figsize=(W / 25.4, H / 25.4))
        self.fig.patch.set_facecolor("white")
        self.c = self.fig.add_axes([0, 0, 1, 1])
        self.c.set_xlim(0, W); self.c.set_ylim(H, 0); self.c.axis("off"); self.c.set_zorder(5)
        self.letters, self.sources = [], {}

    def ax(self, x, y, w, h, **kw):
        a = self.fig.add_axes([x / self.W, 1 - (y + h) / self.H, w / self.W, h / self.H], **kw)
        a.patch.set_visible(False)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        return a

    def box(self, x, y, w, h, fc, ec="none", lw=0, z=1):
        self.c.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=z))

    def txt(self, x, y, s, size=5.0, ha="left", va="center", color=INK, **kw):
        return self.c.text(x, y, s, fontsize=size, ha=ha, va=va, color=color, **kw)

    def panel(self, x, y, letter, title=None):
        self.letters.append(letter)
        self.txt(x, y, letter, size=8, weight="bold")
        if title:
            self.txt(x + 3.6, y + 0.1, title, size=6.5, weight="bold")

    def line(self, x0, y0, x1, y1, color=INK, lw=0.4, **kw):
        self.c.plot([x0, x1], [y0, y1], color=color, lw=lw, solid_capstyle="butt", **kw)

    def cbar(self, x, y, w, h, cmap, lo, hi, ticks, label=None, labels=None):
        a = self.ax(x, y, w, h)
        a.imshow(np.linspace(0, 1, 200)[None, :], aspect="auto", cmap=cmap, extent=(lo, hi, 0, 1), rasterized=True)
        a.set_yticks([])
        a.set_xticks(ticks, labels or [minus(t) for t in ticks], fontsize=5)
        a.tick_params(length=1.5, pad=1, width=0.4)
        for s in a.spines.values():
            s.set_visible(False)
        if label:
            self.txt(x + w + 1.5, y + h / 2, label, size=5)
        return a

    def src(self, letter, df):
        self.sources[letter] = df if isinstance(df, pd.DataFrame) else pd.DataFrame(df)

    def save(self, name):
        self.fig.savefig(f"{OUT}/{name}.pdf", dpi=600)
        self.fig.savefig(f"{OUT}/{name}.png", dpi=300)
        plt.close(self.fig)
        if self.sources:
            os.makedirs(f"{OUT}/source", exist_ok=True)
            with pd.ExcelWriter(f"{OUT}/source/{name}_source.xlsx") as xw:
                for k, df in self.sources.items():
                    df.to_excel(xw, sheet_name=str(k)[:31], index=False)
        print("wrote", f"{OUT}/{name}.pdf", self.W, "x", self.H, "panels", len(self.letters), "".join(self.letters))


def minus(v):
    return str(v).replace("-", "−")


def tile(pg, x, y, w, h, fc, s=None, size=5.0, weight=None):
    """one matrix cell, parted from its neighbours by paper; the value printed in the contrasting ink"""
    pg.box(x, y, w - 0.15, h - 0.15, fc, z=2)
    if s is not None:
        pg.txt(x + (w - 0.15) / 2, y + (h - 0.15) / 2, s, size=size, ha="center", color=ink_on(fc), zorder=3,
               **({"weight": weight} if weight else {}))


def zero(ax, at=0.0, axis="x", color=INK):
    """the dotted reference every value is read against"""
    (ax.axvline if axis == "x" else ax.axhline)(at, color=color, lw=0.5, ls=(0, (1.0, 1.5)), zorder=1)


def point(ax, x, y, c, m="o", ms=3.4, zorder=4, filled=True, mew=None):
    """one marker: filled with a white rim, or open with a coloured rim (not called, not validated)"""
    if filled:
        ax.plot(x, y, marker=m, ms=ms, mfc=c, mec="white", mew=0.3 if mew is None else mew, ls="none", zorder=zorder)
    else:
        ax.plot(x, y, marker=m, ms=ms, mfc="white", mec=c, mew=0.6 if mew is None else mew, ls="none", zorder=zorder)


def ci(ax, x0, x1, y, c=INK, lw=0.6, horizontal=True, zorder=3):
    a, b = ([x0, x1], [y, y]) if horizontal else ([y, y], [x0, x1])
    ax.plot(a, b, color=c, lw=lw, solid_capstyle="butt", zorder=zorder)


def strip(ax, k, vals, c, m="o", ms=3.0, spread=0.16, median=True, horizontal=False, filled=True):
    """points of one group with a black median line"""
    vals = np.asarray(vals, float); vals = vals[np.isfinite(vals)]
    off = np.zeros(len(vals))
    if len(vals) > 1:
        off[np.argsort(vals)] = np.linspace(-spread, spread, len(vals)) * np.array([1, -1] * (len(vals) // 2 + 1))[: len(vals)]
    a, b = (vals, k + off) if horizontal else (k + off, vals)
    point(ax, a, b, c, m=m, ms=ms, filled=filled)
    if median and len(vals):
        q = float(np.median(vals))
        if horizontal:
            ax.plot([q, q], [k - 0.28, k + 0.28], color=INK, lw=0.8, zorder=5, solid_capstyle="butt")
        else:
            ax.plot([k - 0.28, k + 0.28], [q, q], color=INK, lw=0.8, zorder=5, solid_capstyle="butt")


def key(pg, x, y, items, size=5.0, sw=1.8, gap=2.4):
    """a horizontal key. items: (colour, label) or (colour, label, kind), kind 'o','s','D','^','open-o','open-s',
    'line','dot'. Returns the end x."""
    fp = fm.FontProperties(family="Helvetica", size=size)
    for it in items:
        c, lab = it[0], it[1]
        kind = it[2] if len(it) > 2 else None
        if kind and kind.startswith("open-"):
            pg.c.plot(x + sw / 2, y, kind[-1], ms=3.0, mfc="white", mec=c, mew=0.6, zorder=3)
        elif kind in ("o", "s", "D", "^"):
            pg.c.plot(x + sw / 2, y, kind, ms=3.2, mfc=c, mec="white", mew=0.3, zorder=3)
        elif kind == "line":
            pg.line(x, y, x + sw, y, color=c, lw=0.8)
        elif kind == "dot":
            pg.line(x, y, x + sw, y, color=c, lw=0.6, ls=(0, (1.0, 1.5)))
        else:
            pg.box(x, y - sw / 2 * 0.8, sw, sw * 0.8, c, z=2)
        pg.txt(x + sw + 0.6, y, lab, size=size)
        w = TextPath((0, 0), lab.replace("$\\it{", "").replace("}$", ""), prop=fp).get_extents().width * 25.4 / 72
        x += sw + 0.6 + w + gap
    return x


def hist(ax, v, bins, color="#BDBDBD", **kw):
    """grey histogram bars with thin white partitions; the only bars on a sheet"""
    return ax.hist(v, bins=bins, color=color, edgecolor="white", linewidth=0.3, **kw)


# ---- Palette check, as in cns_common.py
_M = {"protan": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
      "deutan": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
      "tritan": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]]}


def _oklab(lin):
    lin = np.clip(lin, 0, 1)
    a = np.array([[.4122214708, .5363325363, .0514459929], [.2119034982, .6806995451, .1073969566], [.0883024619, .2817188376, .6299787005]]) @ lin
    return np.array([[.2104542553, .7936177850, -.0040720468], [1.9779984951, -2.4285922050, .4505937099], [.0259040371, .7827717662, -.8086757660]]) @ np.cbrt(a)


def check_palette(name, colours, normal_floor=15, cvd_floor=8):
    """every pair differs by OKLab delta E x100 >= 15 in normal vision and >= 8 under simulated CVD"""
    import itertools

    def lin(h):
        c = np.array(to_rgb(h))
        return np.where(c <= .04045, c / 12.92, ((c + .055) / 1.055) ** 2.4)
    bad = []
    for (an, a), (bn, b) in itertools.combinations(colours.items(), 2):
        n = np.linalg.norm(_oklab(lin(a)) - _oklab(lin(b))) * 100
        c = min(np.linalg.norm(_oklab(np.array(m) @ lin(a)) - _oklab(np.array(m) @ lin(b))) * 100 for m in _M.values())
        if n < normal_floor or c < cvd_floor:
            bad.append((an, bn, round(float(n), 1), round(float(c), 1)))
    print(f"palette {name}: {'pass' if not bad else 'FAIL ' + str(bad)}")
    return not bad


if __name__ == "__main__":
    check_palette("unit", UNIT)
    check_palette("family", FAM)
    check_palette("unit+family", {**UNIT, **{k: v for k, v in FAM.items() if k != "Other"}})
    check_palette("scenario", SCEN)
    check_palette("tumour type", TYPE)
