"""Publication style for the manuscript figures.

Nature Communications specifications: 180 mm double column, 88 mm single column,
sans-serif face, 5-7 pt type, 0.25-1 pt rules. Everything here exists so that the five
main figures read as one set rather than five separate scripts.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

MM = 1 / 25.4
W2 = 180 * MM          # double column
W1 = 88 * MM           # single column

# colour-blind-safe; blue carries the immune/positive story, vermilion the confounders
BLUE   = "#0F6FC5"     # immune / inflamed
VERM   = "#C1401A"     # confounders and controls
TEAL   = "#1B7F79"     # stroma
GREY   = "#8E9AA6"     # other / neutral marks
INK    = "#1A1F24"     # text and axes
LIGHT  = "#D8DEE4"     # null bands, guide rules
SLATE  = "#3F4C5A"     # full-transcriptome / learned models

_fam = [f.name for f in font_manager.fontManager.ttflist]
SANS = next((f for f in ("Helvetica", "Arial", "Nimbus Sans", "Liberation Sans",
                         "DejaVu Sans") if f in _fam), "DejaVu Sans")


def apply():
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,
        # mathtext must match the body face, or every rho and subscript arrives in a
        # different typeface than the label beside it
        "mathtext.fontset": "custom", "mathtext.rm": SANS, "mathtext.it": f"{SANS}:italic",
        "mathtext.bf": f"{SANS}:bold", "mathtext.default": "it",
        "font.family": SANS, "font.size": 6.2,
        "axes.labelsize": 6.4, "axes.titlesize": 7.0,
        "xtick.labelsize": 5.9, "ytick.labelsize": 5.9,
        "legend.fontsize": 5.6, "figure.dpi": 400, "savefig.dpi": 400,
        "axes.linewidth": 0.5, "axes.edgecolor": INK, "axes.labelcolor": INK,
        "text.color": INK, "xtick.color": INK, "ytick.color": INK,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "legend.handlelength": 1.1,
        "legend.handletextpad": 0.5, "legend.labelspacing": 0.35,
        "lines.linewidth": 0.9, "lines.markersize": 3.2,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.01,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def panel(ax, letter, dx=-0.085, dy=1.055):
    """Bold panel letter, placed identically on every axes in the set."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8.2,
            fontweight="bold", va="top", ha="left", color=INK)


def title(ax, text, pad=4):
    ax.set_title(text, fontsize=6.8, loc="left", pad=pad, color=INK)


def despine(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save(fig, stem, outdir="."):
    for ext in ("pdf", "png"):
        fig.savefig(f"{outdir}/{stem}.{ext}", format=ext)
    plt.close(fig)


# ---------------------------------------------------------------------------------
# Shared encodings (v2 figure set). Every main figure draws its colours from here so
# that a family, a cell type or a tumour type is the same colour on every page.
# ---------------------------------------------------------------------------------
from matplotlib.colors import LinearSegmentedColormap as _LSC

FIGDIR = "/groups/adv2105_gp/yichen/Yi/com/gastric/figs_pub"
AMBER  = "#E3A33B"     # 'missed' / intermediate calls only
PAPER  = "#F4F6F8"     # neutral midpoint, panel wash

DIV = _LSC.from_list("div", ["#0F6FC5", "#89B6E0", "#F4F6F8", "#E2A08A", "#C1401A"])
SEQ = _LSC.from_list("seq", ["#F4F6F8", "#9CC3E6", "#0F6FC5", "#0A3D6E"])

# signature families (manuscript naming -> colour)
FAM = {"Inflamed/IFN": BLUE, "Immune": BLUE, "Other": GREY,
       "Confounder": VERM, "Control": VERM}
FAMLAB = {"Inflamed/IFN": "inflamed/IFN", "Immune": "immune", "Other": "other",
          "Confounder": "confounder control", "Control": "confounder control"}

# cell types of the gastric atlas, grouped lymphoid / myeloid / stroma / epithelium
LYM = ["Tcell", "NK", "Bcell", "Plasma"]; MYE = ["Macro", "DC"]
STR = ["Fibro", "Endo", "Peri"];          EPI = ["Epi"]
CTORDER = LYM + MYE + STR + EPI
CT = {"Tcell": "#1B6FA8", "NK": "#2E9BC9", "Bcell": "#6FB8DC", "Plasma": "#9FD1E8",
      "Macro": "#E08A2E", "DC": "#F2B85C",
      "Fibro": "#7A8B99", "Endo": "#4E6472", "Peri": "#A9B6C1", "Epi": "#C1401A"}
CTGROUP = {t: g for g, ts in [("lymphoid", LYM), ("myeloid", MYE), ("stroma", STR),
                              ("epithelium", EPI)] for t in ts}
GROUPCOL = {"lymphoid": CT["Tcell"], "myeloid": CT["Macro"], "stroma": CT["Fibro"],
            "epithelium": CT["Epi"]}

# tumour types of the checkpoint cohorts
TUMCOL = {"melanoma": "#2E6FA8", "kidney": "#C1401A", "lung": "#1B7F79",
          "bladder": "#4E6472", "gastric": "#E08A2E", "pancreas": "#8A6BB0"}
ICB_TUM = {"ICB_Gide": "melanoma", "ICB_Hugo": "melanoma", "ICB_Liu": "melanoma",
           "ICB_Riaz": "melanoma", "ICB_Van_Allen": "melanoma", "ICB_Puch": "melanoma",
           "ICB_Braun": "kidney", "ICB_Miao1": "kidney", "ICB_Jung": "lung",
           "ICB_Limagne1": "lung", "ICB_Mariathasan": "bladder", "ICB_Kim": "gastric"}

_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def fmt_p(p, lead="P = "):
    """Nature style: two significant figures, 'a × 10⁻ᵇ' below 0.001."""
    if p is None or p != p:
        return f"{lead}n.a."
    if p >= 0.001:
        return f"{lead}{p:.2g}" if p < 0.1 else f"{lead}{p:.2f}"
    m, e = f"{p:.1e}".split("e")
    return f"{lead}{m} × 10{str(int(e)).translate(_SUP)}"


def note(ax, x, y, text, ha="left", va="top", size=5.2, color=None, **kw):
    """Secondary annotation in axes coordinates, always in the muted ink."""
    return ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=size,
                   color=color or SLATE, **kw)


def clean(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def panel_at(fig, ax, letter, dx_mm=5.5, dy_mm=1.2):
    """Panel letter placed a fixed physical distance left of / above the axes box, so
    letters in one row line up regardless of each axes' tick-label width."""
    bb = ax.get_position()
    W, H = fig.get_size_inches()
    x = bb.x0 - dx_mm * MM / W
    y = bb.y1 + dy_mm * MM / H
    fig.text(x, y, letter, fontsize=8.2, fontweight="bold", va="bottom", ha="left",
             color=INK)


def swatch_legend(ax, items, x=0.02, y=0.98, dy=0.075, size=5.0, box=0.028,
                  marker="s"):
    """Compact legend of (colour, label) drawn in axes coordinates - avoids the default
    legend box colliding with data when panels are small."""
    for i, (c, lab) in enumerate(items):
        yy = y - i * dy
        ax.scatter([x + box / 2], [yy], s=9, c=[c], marker=marker, lw=0,
                   transform=ax.transAxes, clip_on=False, zorder=10)
        ax.text(x + box + 0.012, yy, lab, transform=ax.transAxes, va="center",
                ha="left", fontsize=size, color=SLATE, zorder=10)


def savefig(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIGDIR}/{stem}.{ext}", dpi=400)
