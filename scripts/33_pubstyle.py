"""33_pubstyle.py -- publication restyle hook for the Genes & Diseases figure rebuild.
Applies Times-metric serif type (TeX Gyre Termes), 600-dpi raster output and removal of
descriptive panel headings, then executes the target figure script unchanged.
Usage: python3 33_pubstyle.py <mode> <script.py> [args...]
  mode = all   -> every set_title/suptitle suppressed (composite/standalone figure scripts)
  mode = long  -> only headings longer than KEEP chars suppressed; short identity labels
                  (sample names, gene names, programme names) are kept
  mode = none  -> nothing suppressed (script's titles are all data labels)
  mode = keep:A,B -> everything suppressed except the exact titles listed
"""
import sys, os, matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

TG = "/usr/share/texmf/fonts/opentype/public/tex-gyre"
for f in sorted(os.listdir(TG)):
    if f.startswith("texgyretermes"):
        fm.fontManager.addfont(os.path.join(TG, f))
SERIF = "TeX Gyre Termes"

MODE = sys.argv[1]
KEEP = 26  # characters; identity labels are short, descriptive headings are not

plt.rcParams.update({
    "font.family": "serif", "font.serif": [SERIF, "Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "custom", "mathtext.rm": SERIF, "mathtext.it": f"{SERIF}:italic",
    "mathtext.bf": f"{SERIF}:bold",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 600, "savefig.dpi": 600,
    "lines.antialiased": True, "patch.antialiased": True,
    "axes.titlesize": 7.5, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "font.size": 7.5,
})

# the target scripts set their own fonts and dpi; lock the keys this restyle owns
LOCK = {"font.family", "font.serif", "font.sans-serif", "mathtext.fontset", "mathtext.rm",
        "mathtext.it", "mathtext.bf", "figure.dpi", "savefig.dpi", "pdf.fonttype",
        "ps.fonttype", "svg.fonttype"}
_si = matplotlib.RcParams.__setitem__
def __setitem__(self, k, v):
    if k in LOCK and self is plt.rcParams: return
    return _si(self, k, v)
matplotlib.RcParams.__setitem__ = __setitem__

_t, _s = Axes.set_title, Figure.suptitle
KEEPLIST = MODE[5:].split(",") if MODE.startswith("keep:") else []
def drop(x):
    x = str(x)
    if MODE == "none": return False
    if MODE == "all": return True
    if KEEPLIST: return x not in KEEPLIST
    return len(x) > KEEP
def set_title(self, label="", *a, **k):
    return _t(self, "" if drop(label) else label, *a, **k)
def suptitle(self, t="", *a, **k):
    return _s(self, "" if drop(t) else t, *a, **k)
Axes.set_title, Figure.suptitle = set_title, suptitle

_sf = Figure.savefig
def savefig(self, fname, *a, **k):
    if isinstance(fname, str) and fname.lower().endswith((".png", ".jpg", ".tif", ".tiff")):
        k["dpi"] = max(int(k.get("dpi") or 0), 600)
    return _sf(self, fname, *a, **k)
Figure.savefig = savefig

script = sys.argv[2]
sys.argv = sys.argv[2:]
g = {"__name__": "__main__", "__file__": os.path.abspath(script)}
exec(compile(open(script).read(), script, "exec"), g)
