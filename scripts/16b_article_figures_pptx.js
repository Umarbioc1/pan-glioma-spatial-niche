// One slide per article figure (image) with editable title and legend text; source figures in 06_manuscript/figures
const pptxgen = require("pptxgenjs"); const fs = require("fs");
const pres = new pptxgen(); pres.layout = "LAYOUT_WIDE"; pres.title = "Article figures";
pres.theme = { headFontFace: "Cambria", bodyFontFace: "Calibri" };
const L = JSON.parse(fs.readFileSync("manuscript_figs/legends.json"));
for (const f of L) {
  const s = pres.addSlide();
  s.addText(f.title, { x: 0.4, y: 0.2, w: 12.5, h: 0.5, fontSize: 18, bold: true, color: "1B2A41", isTextBox: true, fontFace: "Cambria" });
  const r = Math.min(7.6 / f.w, 6.2 / f.h); const w = f.w * r, h = f.h * r;
  s.addImage({ path: f.png, x: 0.4 + (7.6 - w) / 2, y: 0.85 + (6.2 - h) / 2, w, h, altText: f.title });
  s.addText(f.legend, { x: 8.2, y: 0.85, w: 4.8, h: 6.2, fontSize: 10, color: "333333", valign: "top", isTextBox: true, fontFace: "Calibri" });
}
pres.writeFile({ fileName: "manuscript_figs/Article_Figures_editable.pptx" }).then(() => console.log("ok"));
