"""28_svg_export_hook.py -- run a figure script so every PDF it saves is also saved as an SVG with live text.
Usage: python3 28_svg_export_hook.py <figure_script.py> [script args...]  (set SVG_ON_PNG=1 to also export PNG-only plots)
Text is kept as <text> elements (svg.fonttype='none') so fonts and sizes can be edited in PowerPoint,
Illustrator or Inkscape; scatter layers with >3,000 points are rasterized in the SVG to keep files small."""
import os, sys, runpy, matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mf
from matplotlib.collections import PathCollection, PolyCollection, QuadMesh
_orig = mf.Figure.savefig
def _save(self, fname, *a, **k):
    out = _orig(self, fname, *a, **k)
    if isinstance(fname, str) and (fname.endswith(".pdf") or (os.environ.get("SVG_ON_PNG") and fname.endswith(".png"))):
        heavy = []
        for ax in self.axes:
            for c in ax.collections:
                n = len(c.get_offsets()) if isinstance(c, PathCollection) else (len(c.get_paths()) if isinstance(c, PolyCollection) else 0)
                if isinstance(c, QuadMesh) or n > 3000:
                    if not c.get_rasterized(): heavy.append(c); c.set_rasterized(True)
            for im in ax.images: pass  # images are embedded as bitmaps already
        with matplotlib.rc_context({"svg.fonttype": "none"}):
            k2 = {kk: v for kk, v in k.items() if kk != "dpi"}
            _orig(self, fname[:-4] + ".svg", dpi=300, **k2)
        for c in heavy: c.set_rasterized(False)
    return out
mf.Figure.savefig = _save
script = sys.argv[1]; sys.argv = sys.argv[1:]
runpy.run_path(script, run_name="__main__")
