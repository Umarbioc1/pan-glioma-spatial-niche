// 31_editable_figures_pptx.js -- PowerPoint with every article figure as an SVG whose text stays editable.
// In PowerPoint (Microsoft 365 / 2019+): select a figure -> right-click -> Convert to Shape -> Ungroup -> edit any text.
// Usage: node 31_editable_figures_pptx.js <svg_dir> <legends.json> <out.pptx>
const pptxgen = require("pptxgenjs"); const fs = require("fs"); const path = require("path");
const [svgDir, legF, outF] = process.argv.slice(2);
const legends = JSON.parse(fs.readFileSync(legF, "utf8"));
const byKey = Object.fromEntries(legends.map(l => [l.key, l]));
const order = ["Figure_1","Figure_2","Figure_3","Figure_4","Figure_5","Figure_6","Figure_7","Figure_S1","Figure_S2","Figure_S3","Figure_S4","Figure_S5","Figure_S6","Figure_S7","Figure_S8","Figure_S9"];
const pres = new pptxgen(); pres.layout = "LAYOUT_WIDE"; pres.title = "Editable article figures";
pres.theme = { headFontFace: "Arial", bodyFontFace: "Arial" };
const W = 13.333, H = 7.5;
// how-to slide
let s = pres.addSlide();
s.addText("Editable figures: how to change fonts and sizes", { x: 0.6, y: 0.5, w: 12, h: 0.7, fontSize: 26, bold: true, color: "1F2A44", isTextBox: true });
s.addText([
  { text: "Each figure is an SVG graphic with live text (not a picture of text). Fonts are Arial.", options: { bullet: true, breakLine: true } },
  { text: "PowerPoint (Microsoft 365 or 2019+): click the figure, then right-click → Convert to Shape (or Graphics Format → Convert to Shape).", options: { bullet: true, breakLine: true } },
  { text: "Right-click again → Group → Ungroup (repeat if needed). Every label, tick, title and legend entry is now a text box you can restyle.", options: { bullet: true, breakLine: true } },
  { text: "To restyle many labels at once, Ctrl-click them (or drag a box around them), then set the font and size on the Home tab.", options: { bullet: true, breakLine: true } },
  { text: "Spatial spot maps and heatmap cells stay as images inside the graphic; their titles, axes and legends are editable.", options: { bullet: true, breakLine: true } },
  { text: "For Illustrator or Inkscape, open the matching .svg in the editable_svg folder: all text is selectable there too.", options: { bullet: true } },
], { x: 0.8, y: 1.5, w: 11.8, h: 4.8, fontSize: 16, color: "333333", paraSpaceAfter: 10, valign: "top", isTextBox: true });
for (const k of order) {
  const f = path.join(svgDir, k + ".svg"); const svg = fs.readFileSync(f, "utf8");
  const vb = svg.match(/viewBox="([^"]+)"/)[1].split(/\s+/).map(Number); const ar = vb[2] / vb[3];
  const sl = pres.addSlide();
  const lab = k.replace("Figure_", "Figure ");
  const title = byKey[k] ? byKey[k].title : lab;
  sl.addText(title, { x: 0.3, y: 0.12, w: W - 0.6, h: 0.45, fontSize: 14, bold: true, color: "1F2A44", isTextBox: true, margin: 0 });
  const boxW = W - 0.6, boxH = H - 0.8; let w = boxW, h = w / ar; if (h > boxH) { h = boxH; w = h * ar; }
  sl.addImage({ path: f, x: (W - w) / 2, y: 0.65, w, h, altText: title, objectName: lab });
  sl.addNotes((byKey[k] ? byKey[k].title + " " + byKey[k].legend : lab) + "\n\nTo edit: right-click the figure → Convert to Shape → Ungroup.");
}
pres.writeFile({ fileName: outF }).then(f => console.log("written", f));
// NOTE: pptxgenjs stores the SVG itself as the "PNG" fallback; 31b_fix_png_fallbacks.py replaces those with real PNG renders.
