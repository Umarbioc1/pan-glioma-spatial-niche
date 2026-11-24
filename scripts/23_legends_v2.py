"""23_legends_v2.py -- build manuscript_figs/legends.json (v2 numbering) from ms_v2.md legends and figs_v2 images."""
import re, json
from PIL import Image
txt = open("ms_v2.md").read()
leg = txt.split("## Figure legends")[1].split("### Supplementary tables")[0]
out = []
for m in re.finditer(r"\*\*(Figure (S?\d+)\.[^*]*)\*\*\s*(.+)", leg):
    title, num, body = m.group(1).strip(), m.group(2), m.group(3).strip()
    key = f"Figure_{num}"; png = f"figs_v2/{key}.png"
    w, h = Image.open(png).size
    out.append(dict(key=key, title=title, legend=body, png=png, w=w, h=h))
json.dump(out, open("v2_out/legends.json", "w"), ensure_ascii=False, indent=1)
print(len(out), [o["key"] for o in out])
