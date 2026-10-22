"""Run a supplementary-figure script of the GB v6/v7 set in the palette of the GB v10 main figures (cns10.py), saving under
the GB10 suffix:  python run_restyled_v10.py <script.py>
The script's own layout and data are unchanged; figstyle's colour constants and hatching are replaced before it runs."""
import sys, runpy
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/figscripts")
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric")
import figstyle as fs
import cns10 as K
import matplotlib.axes
fs.BLUE, fs.VERM, fs.TEAL, fs.GREY, fs.INK, fs.LIGHT, fs.SLATE, fs.AMBER = (K.BLUE, K.ORANGE, K.AQUA, K.GRAY, K.INK, K.LGRAY,
                                                                            K.SLATE, K.ORANGE)
fs.FAM.update({"Inflamed/IFN": K.BLUE, "Immune": K.BLUE, "Other": K.GRAY, "Confounder": K.ORANGE})
fs.DIV = K.DIV
_save = fs.savefig
fs.savefig = lambda fig, stem: _save(fig, stem.replace("_GB7", "_GB10").replace("_GB6", "_GB10") if stem.startswith("FigS") else stem + "_restyled_unused")
_bar = matplotlib.axes.Axes.bar


def bar(self, *a, **kw):                     # no hatching: an expected-false segment becomes a light tint of its hue
    if kw.get("hatch"):
        kw.pop("hatch"); ec = kw.get("edgecolor", K.GRAY)
        kw["color"] = K.tint(ec, 0.62); kw["edgecolor"] = "white"
    return _bar(self, *a, **kw)


matplotlib.axes.Axes.bar = bar
runpy.run_path(sys.argv[1], run_name="__main__")
