"""Supplementary Figure 9: gene content against structure, single-cohort estimates and libraries.   python figS9.py
(validation/PLAN_round10_content_and_structure.md; icb/v15_content from decoys/v15_content_structure.py).

rho(s, c) = the share of a signature's 320 co-expressed decoys whose gain over the null in cohort c lies below its own
(0.5 when its genes add nothing to a module of its structure); the other displays of rho are in Figure 7 (fig7.py).
a  single-cohort estimates against the inflamed family's rho
b  libraries: structure (decoys' per-test clearance) against content (mean rho of called targets)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import cns10 as K
from matplotlib.colors import LinearSegmentedColormap, Normalize
K.apply()
V = K.D + "/icb/v15_content/"
FAMD = pd.read_csv(K.D + "/icb/v5_signature_clusters.csv").set_index("signature").family
FAMS = ["Inflamed/IFN", "Other", "Confounder"]
CELL = {k: pd.read_csv(V + f"cells_{k}.csv.gz") for k in ("response", "tcga_OS")}
SIG = {k: pd.read_csv(V + f"signature_{k}.csv", index_col=0) for k in ("response", "tcga_OS")}
FBS = pd.read_csv(V + "family_by_setting.csv")
PC = pd.read_csv(V + "primary_contrast.csv")
MT = pd.read_csv(V + "matched_tumor_types.csv")
AT = pd.read_csv(V + "atlas_by_type.csv")
LIB = pd.read_csv(V + "libraries.csv")
ICB_TUM = {"ICB_Gide": "melanoma", "ICB_Hugo": "melanoma", "ICB_Liu": "melanoma", "ICB_Riaz": "melanoma",
           "ICB_Van_Allen": "melanoma", "ICB_Puch": "melanoma", "ICB_Braun": "kidney", "ICB_Miao1": "kidney",
           "ICB_Jung": "lung", "ICB_Limagne1": "lung", "ICB_Mariathasan": "bladder", "ICB_Kim": "gastric"}
RHO = LinearSegmentedColormap.from_list("rho", [K.AQUA, "#9fd3c3", "#f4f4f2", "#9aa3ad", "#2b3440"])
SETL = {"response": "response", "tcga_OS": "overall survival"}



C = K.Canvas(70)
# ================================================================== f  single-cohort estimates
ax = C.panel("a", 16, 12, 62, 46, head="Single-cohort estimates")
CAP = 150
FLAB = {"SKCM": (0.70, 22), "CESC": (0.64, 9), "Gide": (1.06, 27), "SARC": (0.97, 44), "KIRC": (0.157, 158),
        "LUSC": (0.29, 158), "STAD": (0.23, 127)}
for k_, mk, col, fill in (("tcga_OS", "o", K.SLATE, False), ("response", "D", K.INK, True)):
    q = AT[(AT.setting == k_) & (AT.calls > 0)]
    for _, rr in q.iterrows():
        y = min(rr.fdr, CAP)
        K.dots(ax, [rr.rho_inflamed], [y], col, s=4 + 0.45 * rr.calls, filled=fill, marker=mk, z=5)
        name = K.COHORT.get(rr.cohort, rr.cohort)
        if name in FLAB:
            tx, ty = FLAB[name]
            lead = abs(tx - rr.rho_inflamed) > 0.05 or abs(ty - y) > 10
            ax.annotate(name, (rr.rho_inflamed, y), (tx, ty), ha="center" if name in ("KIRC", "LUSC", "Gide", "SARC")
                        else ("left" if tx > rr.rho_inflamed else "right"), va="center", fontsize=4.3, color=K.INK2,
                        arrowprops=dict(arrowstyle="-", lw=0.35, color=K.INK2, shrinkA=1, shrinkB=2) if lead else None)
ax.axhline(CAP, color=K.RULE, lw=0.4, ls=(0, (1, 1.5)))
ax.set_xlim(0.1, 1.12); ax.set_ylim(0, CAP + 14); ax.set_xticks([.25, .5, .75, 1]); ax.set_yticks([0, 50, 100, 150])
ax.set_yticklabels(["0", "50", "100", "≥150"])
ax.set_xlabel("mean ρ, inflamed family"); ax.set_ylabel("estimated FDR (%)")
K.legend(ax, [(K.SLATE, "TCGA type", "o", False), (K.INK, "response cohort", "D", True)], x=0.55, y=0.80, dy=0.085,
         s=7, size=K.T_SMALL)

# ================================================================== g  libraries
ax = C.panel("b", 106, 12, 70, 46, head="Libraries")
LL = {"response": "immune, response", "tcga_OS": "immune, OS", "tcga_DSS": "immune, DSS", "tcga_PFI": "immune, PFI",
      "breast_pc1_origin_out": "breast, PC1", "breast_mean": "breast, mean z", "hallmark_OS": "Hallmark, OS"}
OFF = {"response": (4, 0.0), "tcga_OS": (4, 0.008), "tcga_DSS": (4, 0.004), "tcga_PFI": (4, -0.012),
       "breast_pc1_origin_out": (-4, 0.0), "breast_mean": (-4, 0.0), "hallmark_OS": (4, 0.0)}
for _, q in LIB[LIB.library.isin(LL)].iterrows():
    x_, y_ = 100 * q.decoy_clear, q.rho_bar_called
    K.dots(ax, [x_], [y_], K.INK, s=14, z=5)
    dx, dy = OFF[q.library]
    ax.text(x_ + dx * 0.35, y_ + dy, f"{LL[q.library]} {q.fdr_reported:.0f}%", ha="left" if dx > 0 else "right",
            va="center", fontsize=4.3, color=K.INK2)
K.ref(ax, 5, axis="v")
ax.set_xlim(3, 17); ax.set_ylim(0.6, 0.92); ax.set_xticks([5, 10, 15])
ax.set_xlabel("decoys above the null's 95th pct (%)"); ax.set_ylabel("mean ρ of the called")

C.save("SupplementaryFigure9")
print("SupplementaryFigure9")
