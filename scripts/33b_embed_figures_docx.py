"""33_embed_figures_docx.py -- reading/editing version of the Neuro-Oncology manuscript with Figures 1-6 placed in the Results
(at the end of the subsection that first discusses each figure); full captions stay in the Figure Captions section.
Usage: python3 33_embed_figures_docx.py <manuscript.md> <figure_png_dir> <out.md>"""
import sys, re
src, figdir, out = sys.argv[1:4]
s = open(src).read()
caps = dict(re.findall(r"\*\*Figure (\d)\. ([^*]+)\*\*", s.split("## Figure Captions")[1]))
anchor = {1: "### Ten spatial niches", 2: "### Niche 6 is a hypoxia", 3: "### Most niche-6 tumor genes", 4: "### ANGPTL4 is produced",
          5: "### The niche signature tracks", 6: "### Niche-resolved prioritization"}
for n, a in anchor.items():
    i = s.index(a)
    block = f"![**Figure {n}.** {caps[str(n)].strip()} Full legend in Figure Captions.]({figdir}/Figure_{n}.png){{width=6.5in}}\n\n"
    s = s[:i] + block + s[i:]
open(out, "w").write(s); print("figures placed:", len(anchor))
