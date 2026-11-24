"""30_compose_editable_svgs.py -- build one editable SVG per article figure (v2 numbering).
Composite figures made by 16_article_figures.py embed earlier plots as bitmaps; here each such bitmap is replaced by the
vector SVG of the same plot (re-drawn with identical code by 28_svg_export_hook.py / 29_replot_map_panels.py), so every
label, tick and title is live text. Fonts are set to Arial. Spots in spatial maps stay rasterized inside the SVG.
Usage: python3 30_compose_editable_svgs.py <outdir>
"""
import re, sys, os, json
OUT = sys.argv[1] if len(sys.argv) > 1 else "edit_figs/svg"; os.makedirs(OUT, exist_ok=True)
P = "edit_figs/panels"; RR = "/tmp/rr"
# composite (16_article_figures output name) -> panels in drawing order; optional white "blank" box in source-PNG pixels
COMPOSITE = {
    "Figure_1": [("%s/compartment_maps.svg" % P, None, (1529, 1425))],
    "Figure_2": [("%s/niche_maps.svg" % P, (1330, 1180, 1529, 1425), (1529, 1425))],
    "Figure_3": [("%s/lr_out/lr_niche_dotplot.svg" % RR, None, None)],
    "Figure_5": [("%s/bulk_out/scores_by_type_grade.svg" % RR, None, None), ("%s/cgga_out/cgga_scores_by_type_grade.svg" % RR, None, None),
                 ("%s/cgga_out/meta_forest.svg" % RR, None, None), ("%s/cgga_out/cgga_KM_HypoxicMES_N6.svg" % RR, None, None)],
    "Figure_S1": [("%s/program_maps_greenwald.svg" % P, None, None), ("%s/program_maps_ravi.svg" % P, None, None)],
    "Figure_S2": [("%s/niche_profiles_programscore.svg" % P, None, None), ("%s/niche_profiles_programscore_harmony.svg" % P, None, None)],
    "Figure_S3": [("%s/prog/progression_niche_score_maps.svg" % RR, None, None)],
    "Figure_S4": [("%s/icb_out/gse121810_arm_scores.svg" % RR, None, None)],
}
SRC = {k: f"manuscript_figs/{k}.svg" for k in ["Figure_1", "Figure_2", "Figure_3", "Figure_4", "Figure_5", "Figure_6", "Figure_S1", "Figure_S2", "Figure_S3", "Figure_S4"]}
SRC.update({"Figure_5_new": "edit_figs/raw/Figure_5_new.svg", "Figure_S5": "edit_figs/raw/Figure_S5.svg", "Figure_S6": "edit_figs/raw/Figure_S6.svg",
            "Figure_S7": "edit_figs/raw/Figure_S7.svg", "Figure_S8": "edit_figs/raw/Figure_S8.svg", "Figure_S9": "edit_figs/raw/Figure_S9.svg"})
RENUM = json.load(open("figs_v2/renumbering_map.json"))

def panel_group(path, tag, x, top, w, h, blank, png_size):
    s = open(path).read()
    vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', s).group(1).split()]
    W, H = vb[2], vb[3]
    inner = re.search(r"<svg[^>]*>(.*)</svg>", s, re.S).group(1)
    inner = re.sub(r"<metadata>.*?</metadata>", "", inner, flags=re.S)
    ids = set(re.findall(r'id="([^"]+)"', inner))
    for i in sorted(ids, key=len, reverse=True):  # namespace ids so panels cannot collide
        inner = inner.replace(f'id="{i}"', f'id="{tag}{i}"').replace(f"url(#{i})", f"url(#{tag}{i})").replace(f'href="#{i}"', f'href="#{tag}{i}"')
    sx, sy = w / W, h / H
    if abs(sx - sy) / sx > 0.02: print(f"  warning: aspect mismatch in {path}: {sx:.4f} vs {sy:.4f}")
    extra = ""
    if blank:
        pw, ph = png_size; kx, ky = W / pw, H / ph
        extra = f'<rect x="{blank[0]*kx:.2f}" y="{blank[1]*ky:.2f}" width="{(blank[2]-blank[0])*kx:.2f}" height="{(blank[3]-blank[1])*ky:.2f}" fill="#ffffff"/>'
    return f'<g id="{tag}panel" transform="translate({x:.3f} {top:.3f}) scale({sx:.6f} {sy:.6f})">{inner}{extra}</g>', W / H

def image_boxes(s):
    out = []
    for m in re.finditer(r"<image\b[^>]*/>", s, re.S):
        t = m.group(0)
        g = lambda a: float(re.search(rf'\b{a}="([-0-9.e]+)"', t).group(1))
        x, y, w, h = g("x"), g("y"), g("width"), g("height")
        top = -y if "scale(1 -1)" in t else y
        out.append((m.start(), m.end(), x, top, w, h))
    return out

for old, path in SRC.items():
    s = open(path).read()
    if old in COMPOSITE:
        boxes = image_boxes(s); panels = COMPOSITE[old]; used = set(); repl = []
        for k, (pp, blank, png_size) in enumerate(panels):
            ps = open(pp).read(); vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', ps).group(1).split()]; ar = vb[2] / vb[3]
            cand = [b for b in boxes if b[0] not in used and abs((b[4] / b[5]) / ar - 1) < 0.02]
            if not cand: raise SystemExit(f"no image slot for {pp} in {old}")
            b = cand[0]; used.add(b[0])
            grp, _ = panel_group(pp, f"p{k}_", b[2], b[3], b[4], b[5], blank, png_size); repl.append((b[0], b[1], grp))
        for a, z, grp in sorted(repl, reverse=True): s = s[:a] + grp + s[z:]
    s = re.sub(r"font-family:\s*'?DejaVu Sans'?", "font-family:Arial", s).replace("'DejaVu Sans'", "Arial").replace("DejaVu Sans", "Arial")
    new = RENUM[old]
    open(f"{OUT}/{new}.svg", "w").write(s)
    print(new, "from", old, "| live text:", len(re.findall(r"<text", s)), "| bitmaps left:", len(re.findall(r"<image", s)), f"| {len(s)/1e6:.2f} MB")
