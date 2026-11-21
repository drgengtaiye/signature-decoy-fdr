"""Run a supplementary-figure script drawn with figstyle.py (figS2.py, figS4.py to figS8.py) in the palette of cns10.py:
    python run_restyled.py <script.py>
The script's own layout and data are unchanged; figstyle's colour constants and hatching are replaced before it runs."""
import os, sys, runpy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/path/to/analysis_root")
import figstyle as fs
import cns10 as K
import matplotlib.axes
fs.BLUE, fs.VERM, fs.TEAL, fs.GREY, fs.INK, fs.LIGHT, fs.SLATE, fs.AMBER = (K.BLUE, K.ORANGE, K.AQUA, K.GRAY, K.INK, K.LGRAY,
                                                                            K.SLATE, K.ORANGE)
fs.FAM.update({"Inflamed/IFN": K.BLUE, "Immune": K.BLUE, "Other": K.GRAY, "Confounder": K.ORANGE})
fs.DIV = K.DIV
_bar = matplotlib.axes.Axes.bar


def bar(self, *a, **kw):                     # no hatching: an expected-false segment becomes a light tint of its hue
    if kw.get("hatch"):
        kw.pop("hatch"); ec = kw.get("edgecolor", K.GRAY)
        kw["color"] = K.tint(ec, 0.62); kw["edgecolor"] = "white"
    return _bar(self, *a, **kw)


matplotlib.axes.Axes.bar = bar
runpy.run_path(sys.argv[1], run_name="__main__")
