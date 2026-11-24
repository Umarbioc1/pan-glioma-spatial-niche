"""31b_fix_png_fallbacks.py -- replace pptxgenjs' SVG-as-PNG fallbacks with real PNGs (figs_v2), so viewers without SVG support show the figure.
Usage: python3 31b_fix_png_fallbacks.py <in.pptx> <out.pptx> <png_dir>"""
import sys, zipfile, re, io
from PIL import Image
src, dst, pngdir = sys.argv[1:4]
order = ["Figure_1","Figure_2","Figure_3","Figure_4","Figure_5","Figure_6","Figure_7","Figure_S1","Figure_S2","Figure_S3","Figure_S4","Figure_S5","Figure_S6","Figure_S7","Figure_S8","Figure_S9"]
zin = zipfile.ZipFile(src); zout = zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED)
for it in zin.infolist():
    data = zin.read(it.filename)
    m = re.match(r"ppt/media/image-(\d+)-1\.png$", it.filename)
    if m:
        key = order[int(m.group(1)) - 2]
        im = Image.open(f"{pngdir}/{key}.png").convert("RGB"); im.thumbnail((2400, 2400))
        b = io.BytesIO(); im.save(b, "PNG", optimize=True); data = b.getvalue()
    if not it.filename.endswith("/"): zout.writestr(it, data)
zout.close(); print("ok", dst)
