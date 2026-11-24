// Editable PPTX record of all project figures. Rebuild after each analysis step:
//   python3 scripts/08_prep_figure_deck.py && node scripts/08_build_figure_deck.js
const pptxgen = require("pptxgenjs");
const fs = require("fs");
const { applyTheme } = require(process.env.SKILL + "/scripts/apply_theme.js");
const D = JSON.parse(fs.readFileSync("deck/deck_data.json"));
const OUT = "deck/Glioma_Spatial_Atlas_Figures.pptx";
const THEME = { name: "Glioma Atlas", headFontFace: "Cambria", bodyFontFace: "Calibri",
  colors: { dk1: "1E1E24", lt1: "FFFFFF", dk2: "1B2A41", lt2: "EEF1F5",
    accent1: "B5172F", accent2: "2A8C82", accent3: "D9A441", accent4: "5B6C8F", accent5: "8E5572", accent6: "7A8B99",
    hlink: "2A8C82", folHlink: "5B6C8F" } };
const pres = new pptxgen(); pres.layout = "LAYOUT_WIDE";
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "Pan-glioma spatial niche atlas: figure record"; pres.author = "Umar Farooq";
const C = pres.SchemeColor;
const W = 13.333, H = 7.5;
pres.defineSlideMaster({ title: "TITLE", background: { color: THEME.colors.dk2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.3, w: 11.7, h: 1.6, fontSize: 40, bold: true, color: C.background1, valign: "bottom" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 4.1, w: 11.7, h: 1.4, fontSize: 18, color: "C9D3E0", valign: "top" }, text: "" } } ] });
pres.defineSlideMaster({ title: "SECTION", background: { color: THEME.colors.dk2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.6, w: 11.7, h: 1.2, fontSize: 36, bold: true, color: C.background1, valign: "bottom" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 3.95, w: 11.7, h: 1.2, fontSize: 16, color: "C9D3E0", valign: "top" }, text: "" } } ],
  slideNumber: { x: 12.4, y: 7.0, w: 0.6, h: 0.3, fontSize: 10, color: "C9D3E0" } });
pres.defineSlideMaster({ title: "FIGURE", background: { color: THEME.colors.lt1 },
  margin: [0.5, 0.5, 0.55, 0.5],
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.3, w: 12.3, h: 0.75, fontSize: 28, bold: true, color: C.text2, valign: "middle" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.5, y: 6.45, w: 12.3, h: 0.6, fontSize: 11, color: "4A5568", valign: "top" }, text: "" } },
    { text: { text: "Pan-glioma spatial niche atlas · figure record", options: { x: 0.5, y: 7.08, w: 8, h: 0.3, fontSize: 9, color: "8A94A6" } } } ],
  slideNumber: { x: 12.4, y: 7.08, w: 0.6, h: 0.3, fontSize: 9, color: "8A94A6" } });

let sec = "";
function section(title, sub) { sec = title; pres.addSection({ title });
  const s = pres.addSlide({ masterName: "SECTION", sectionTitle: sec });
  s.addText(title, { placeholder: "title" }); s.addText(sub, { placeholder: "body" }); }
function fig(title, caption, notes) { const s = pres.addSlide({ masterName: "FIGURE", sectionTitle: sec });
  s.addText(title, { placeholder: "title" }); s.addText(caption, { placeholder: "body" }); s.addNotes(notes); return s; }
function place(s, im, x, y, w, h, name) {
  const r = Math.min(w / im.w, h / im.h); const iw = im.w * r, ih = im.h * r;
  s.addImage({ path: im.file, x: x + (w - iw) / 2, y: y + (h - ih) / 2, w: iw, h: ih, altText: name, objectName: name }); }
const BOX = { x: 0.5, y: 1.2, w: 12.3, h: 5.15 };
function imgSlide(title, caption, notes, ims, name) {
  const s = fig(title, caption, notes); const n = ims.length, gap = 0.3, w = (BOX.w - gap * (n - 1)) / n;
  ims.forEach((im, i) => place(s, im, BOX.x + i * (w + gap), BOX.y, w, BOX.h, name + (n > 1 ? " part " + (i + 1) : "")));
  return s; }
// RdBu_r colour for z-scores in [-1.5, 1.5]
function rdbu(v) { const t = Math.max(-1, Math.min(1, v / 1.5));
  const lerp = (a, b, u) => Math.round(a + (b - a) * u);
  const hex = (r, g, b) => [r, g, b].map(c => c.toString(16).padStart(2, "0")).join("").toUpperCase();
  return t >= 0 ? hex(lerp(247, 178, t), lerp(247, 24, t), lerp(247, 43, t)) : hex(lerp(247, 33, -t), lerp(247, 102, -t), lerp(247, 172, -t)); }

// ---------- Title ----------
pres.addSection({ title: "Cover" }); sec = "Cover";
const t = pres.addSlide({ masterName: "TITLE", sectionTitle: sec });
t.addText("Pan-glioma spatial niche atlas", { placeholder: "title" });
t.addText("Editable record of all project figures · public Visium data (Greenwald 2024, Ravi 2022)\nUmar Farooq · updated " + new Date().toISOString().slice(0, 10), { placeholder: "body" });

// ---------- 1 QC ----------
section("1  Data and quality control", "40 Visium sections passed QC (109,415 spots): Greenwald 2024 (13 GBM + 6 IDH-mutant) and Ravi 2022 (16 GBM + 2 IDH-mutant + 3 cortex controls); 36 tumour sections used for niches");
imgSlide("Spatial UMI depth: Greenwald 2024 cohort", "Total UMI per spot after QC (min 1 count/spot; mitochondrial filter). Source: 01_load_and_qc.py",
  "File: 04_results/01_QC/QC_spatial_UMI_Greenwald2024.png. Script: 03_scripts/01_load_and_qc.py", D.img.qc_g, "QC Greenwald");
imgSlide("Spatial UMI depth: Ravi 2022 cohort", "Total UMI per spot after QC; includes non-tumour cortex sections UKF256_C, UKF259_C, UKF265_C. Source: 01_load_and_qc.py",
  "File: 04_results/01_QC/QC_spatial_UMI.png (Ravi). Script: 03_scripts/01_load_and_qc.py", D.img.qc_r, "QC Ravi");

// ---------- 2 Programs ----------
section("2  Malignant and TME programs", "Gavish 2023 3CA meta-program and cell-state scores per spot");
imgSlide("Program scores per spot: Greenwald 2024", "Left and right panels are the two halves of one figure (sections in order). Source: 02a_program_scores.py",
  "File: 04_results/02_programs/program_maps_all_samples.png (Greenwald). Split in two for slide fit.", D.img.prog_g, "Program maps Greenwald");
imgSlide("Program scores per spot: Ravi 2022", "Left and right panels are the two halves of one figure. Source: 02a_program_scores.py",
  "File: 04_results/02_programs/program_maps_all_samples.png (Ravi). Split in two for slide fit.", D.img.prog_r, "Program maps Ravi");

// ---------- 3 CNV ----------
section("3  CNV-defined compartments", "Malignant / TME-rich / Brain; hallmarks: +7/−10 (GBM), −1p/−19q (oligodendroglioma); reference: Ravi cortex UKF256_C + UKF265_C");
imgSlide("Tissue compartments from inferred CNV", "Spots labelled Malignant, TME-rich or Brain from CNV hallmark signatures (threshold = reference 95th percentile). Source: 03_cnv_compartments.py",
  "File: 04_results/03_CNV_compartments/compartment_maps.png. infercnvpy, hg38 positions from infercna.", D.img.cnv, "CNV compartment maps");

// ---------- 4 Deconvolution ----------
section("4  Cell-type deconvolution", "NNLS against the GBmap core reference (20 level-3 types), per-section platform refit");
{ const s = fig("Deconvolution agrees with independent marker scores", "Pearson r between estimated cell-type fraction and an independent marker-gene score per spot. Median fit r (observed vs predicted expression) = 0.90. Source: 05_deconvolution.py",
    "Values from the validation run of 03_scripts/05_deconvolution.py (per-section refit version). Chart is native: right-click > Edit Data to change.");
  s.addChart(pres.charts.BAR, [{ name: "Pearson r", labels: D.deconv_valid.types, values: D.deconv_valid.r }], {
    x: 0.7, y: 1.35, w: 7.4, h: 4.9, barDir: "bar", chartColors: [THEME.colors.accent2], showValue: true, dataLabelPosition: "outEnd",
    dataLabelFormatCode: "0.00", dataLabelFontSize: 12, dataLabelFontFace: "+mn-lt", dataLabelColor: "1E1E24",
    valAxisMinVal: 0, valAxisMaxVal: 1, valAxisLabelFontSize: 11, catAxisLabelFontSize: 13, catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
    catAxisLabelColor: "1E1E24", valAxisLabelColor: "4A5568", valGridLine: { color: "E2E8F0", size: 0.75 }, catGridLine: { style: "none" },
    catAxisOrientation: "maxMin", showLegend: false, showTitle: true, title: "Agreement with marker score (r)", titleFontSize: 14, titleFontFace: "+mn-lt", titleColor: "1E1E24" });
  s.addText([{ text: "0.90", options: { fontSize: 60, bold: true, color: C.accent2, breakLine: true } },
             { text: "median fit r, observed vs predicted spot expression", options: { fontSize: 14, color: C.text2 } }],
    { x: 8.6, y: 1.6, w: 4.0, h: 2.0, isTextBox: true, valign: "top", objectName: "Fit callout" });
  s.addText("Microglia vs blood-derived macrophages separate only partly at Visium resolution; myeloid fractions are RNA fractions and likely underestimate cell numbers.",
    { x: 8.6, y: 3.9, w: 4.0, h: 2.0, fontSize: 14, color: "4A5568", isTextBox: true, valign: "top", objectName: "Limitation note" }); }

// ---------- 5 Niches ----------
section("5  Lineage composition patterns and niches", "LCPs (k-means, k=14 per compartment) → niches (k=10) from LCP composition within 2 rings; Harmony by cohort");
{ const s = fig("Cell-type composition of the 10 niches", "Mean z-scored deconvolved fraction per niche (red = enriched, blue = depleted). Native table: values and colours are editable. Source: 06_lcp_niches_deconv.py",
    "Data: 04_results/06_niches_deconv/niche_program_profile.csv. Colour scale RdBu, clipped at ±1.5.");
  const hdr = [{ text: "Cell type", options: { bold: true, fill: { color: THEME.colors.dk2 }, color: "FFFFFF", align: "left" } }]
    .concat(D.profile.niches.map(n => ({ text: n, options: { bold: true, fill: { color: THEME.colors.dk2 }, color: "FFFFFF", align: "center" } })));
  const rows = [hdr];
  D.profile.types.forEach((ty, i) => rows.push([{ text: ty, options: { align: "left", bold: true, color: "1E1E24" } }]
    .concat(D.profile.z[i].map(v => ({ text: v.toFixed(1), options: { align: "center", fill: { color: rdbu(v) }, color: Math.abs(v) > 1.0 ? "FFFFFF" : "1E1E24" } })))));
  s.addTable(rows, { x: 0.5, y: 1.15, w: 12.3, colW: [2.3].concat(Array(10).fill(1.0)), rowH: 0.24, fontSize: 10, fontFace: "Calibri",
    border: { type: "solid", pt: 0.5, color: "FFFFFF" }, margin: [0, 0.05, 0, 0.05], valign: "middle", objectName: "Niche composition table" }); }
{ const s = fig("Niche frequency differs between GBM and IDH-mutant glioma", "Percent of spots in each niche, by subtype. Niche 6 (MES-like + TAM-BDM) is almost GBM-only; niche 7 (brain/infiltration) dominates IDH-mutant. Source: 06_lcp_niches_deconv.py",
    "Data: 04_results/06_niches_deconv/niche_frequency_by_subtype.csv. Native chart: right-click > Edit Data.");
  const ser = Object.entries(D.freq.series).map(([k, v]) => ({ name: k, labels: D.freq.niches, values: v }));
  s.addChart(pres.charts.BAR, ser, { x: 0.6, y: 1.2, w: 12.1, h: 5.1, barDir: "col", barGrouping: "clustered", barGapWidthPct: 60,
    chartColors: [THEME.colors.accent1, THEME.colors.accent2], showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: "0.0",
    dataLabelFontSize: 10, dataLabelFontFace: "+mn-lt", dataLabelColor: "4A5568",
    catAxisLabelFontSize: 12, valAxisLabelFontSize: 11, catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
    catAxisLabelColor: "1E1E24", valAxisLabelColor: "4A5568", valGridLine: { color: "E2E8F0", size: 0.75 }, catGridLine: { style: "none" },
    showValAxisTitle: true, valAxisTitle: "% of spots", valAxisTitleFontSize: 12, valAxisTitleFontFace: "+mn-lt", valAxisTitleColor: "4A5568",
    showLegend: true, legendPos: "t", legendFontSize: 12, legendFontFace: "+mn-lt" }); }
imgSlide("Niche maps across all sections", "Each spot coloured by niche (deconvolution-based niches, final version). Source: 06_lcp_niches_deconv.py",
  "File: 04_results/06_niches_deconv/niche_maps.png", D.img.niche_maps, "Niche maps");

// ---------- 6 Ligand–receptor ----------
section("6  Niche ligand–receptor signalling", "CellChatDB (1,912 pairs); spot × first-ring co-expression; niche vs rest Mann–Whitney per section; rank-based meta-analysis");
imgSlide("Niche-enriched ligand–receptor pairs", "Colour = median AUC (niche vs rest of section); size = fraction of sections with AUC > 0.6; ring = significant (FDR < 0.05, median AUC > 0.6, ≥60% of sections). Source: 07b_lr_subtype_and_plot.py",
  "File: 04_results/07_ligand_receptor/lr_niche_dotplot.png. Full tables: lr_niche_meta.csv, lr_niche_significant.csv", D.img.lr_dot, "LR dot plot");
function lrTable(title, caption, notes, data) {
  const s = fig(title, caption, notes);
  const head = ["Niche", "Interaction", "Median AUC", "Sections AUC>0.6", "Sections (G/R)", "AUC GBM", "AUC IDH-mut", "Likely sender", "Likely receiver", "FDR"]
    .map(h => ({ text: h, options: { bold: true, fill: { color: THEME.colors.dk2 }, color: "FFFFFF", align: "center" } }));
  const rows = [head].concat(data.map((r, i) => r.map((v, j) => ({ text: v, options: { align: j === 1 ? "left" : "center", color: "1E1E24",
    fill: { color: i % 2 ? "FFFFFF" : THEME.colors.lt2 } } }))));
  s.addTable(rows, { x: 0.5, y: 1.15, w: 12.3, colW: [0.65, 2.75, 0.95, 1.1, 1.15, 0.9, 1.0, 1.35, 1.35, 1.1], rowH: 0.215, fontSize: 10, fontFace: "Calibri",
    border: { type: "solid", pt: 0.5, color: "FFFFFF" }, margin: [0, 0.05, 0, 0.05], valign: "middle", objectName: "LR table" });
}
lrTable("Niche 6 is a GBM-specific hypoxic/ECM signalling hub", "22 significant pairs; enriched in GBM sections (AUC 0.60–0.77) but not in the few IDH-mutant sections containing niche 6 (AUC ≈ 0.50–0.58). Sender/receiver = top GBmap type expressing the gene (guide only).",
  "Data: lr_niche_significant.csv merged with lr_significant_by_subtype.csv. G/R = Greenwald / Ravi sections.", D.lr6);
lrTable("Significant interactions in the other niches", "SPP1 → CD44 in the myeloid niche 3 is enriched in both GBM (0.67) and IDH-mutant (0.65) tumours. Neuronal and oligodendrocyte niches show normal-brain signalling (internal controls).",
  "Data: lr_niche_significant.csv merged with lr_significant_by_subtype.csv.", D.lr_other);

// ---------- 7 NCDEGs ----------
if (D.ncdeg_tab) {
section("7  Niche-dependent genes (NCDEGs)", "Within malignant-rich and myeloid-rich spots: expression ~ niche + admixture fractions + log UMI; section = replicate (one-sample t-test across sections)");
imgSlide("Malignant cells in niche 6 switch on hypoxia genes", "A: 47 malignant NCDEGs in niche 6 (FDR < 0.05, |β| ≥ 0.2, ≥75% sign consistency). B: after also adjusting for malignant state mix, effects shrink ~3-fold; ANGPTL4 remains up in 14/15 sections. Source: 09_ncdeg.py, 09b_ncdeg_plot.py",
  "Files: 04_results/09_NCDEGs/ncdeg_niche6_malignant.png, niche6_malignant_unadjusted_vs_stateadjusted.csv", D.img.ncdeg, "NCDEG niche 6");
{ const s = fig("Top niche-6 malignant NCDEGs, before and after state adjustment", "β = niche-6 effect (log1p CP10k). State-adjusted model adds MES/NPC/OPC-like share within malignant cells. Myeloid NCDEGs: none significant (3–12 sections per niche; underpowered).",
    "Data: niche6_malignant_unadjusted_vs_stateadjusted.csv. Native table.");
  const head = ["Gene", "Direction", "β (admixture-adj.)", "FDR", "β (+ state-adj.)", "p (state-adj.)", "Sign consistency (state-adj.)"]
    .map(h => ({ text: h, options: { bold: true, fill: { color: THEME.colors.dk2 }, color: "FFFFFF", align: "center" } }));
  const rows = [head].concat(D.ncdeg_tab.map((r, i) => r.map((v, j) => ({ text: v, options: { align: j === 0 ? "left" : "center",
    bold: r[0] === "ANGPTL4", color: j === 1 ? (v === "Up" ? THEME.colors.accent1 : "2A6FB0") : "1E1E24", fill: { color: i % 2 ? "FFFFFF" : THEME.colors.lt2 } } }))));
  s.addTable(rows, { x: 1.2, y: 1.15, w: 10.9, colW: [1.6, 1.1, 1.75, 1.1, 1.75, 1.6, 2.0], rowH: 0.215, fontSize: 10, fontFace: "Calibri",
    border: { type: "solid", pt: 0.5, color: "FFFFFF" }, margin: [0, 0.05, 0, 0.05], valign: "middle", objectName: "NCDEG table" }); }
{ const s = fig("Niche-6 NCDEGs map to hypoxia and MES programs", "Hypergeometric overlap of niche-6 malignant NCDEGs with Gavish 2023 cancer meta-programs (tested-gene universe). Source: 09_ncdeg.py",
    "Data: 04_results/09_NCDEGs/ncdeg_metaprogram_enrichment.csv (all cell-type MPs); chart shows cancer MPs only. Native chart.");
  s.addChart(pres.charts.BAR, [{ name: "-log10 FDR", labels: D.mp_enr.labels, values: D.mp_enr.values }], {
    x: 0.8, y: 1.25, w: 11.7, h: 5.0, barDir: "bar", chartColors: [THEME.colors.accent1], showValue: true, dataLabelPosition: "outEnd",
    dataLabelFormatCode: "0.0", dataLabelFontSize: 11, dataLabelFontFace: "+mn-lt", dataLabelColor: "1E1E24",
    catAxisLabelFontSize: 12, valAxisLabelFontSize: 11, catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
    catAxisLabelColor: "1E1E24", valAxisLabelColor: "4A5568", valGridLine: { color: "E2E8F0", size: 0.75 }, catGridLine: { style: "none" },
    catAxisOrientation: "maxMin", showLegend: false, showValAxisTitle: true, valAxisTitle: "−log10 FDR", valAxisTitleFontSize: 12, valAxisTitleFontFace: "+mn-lt", valAxisTitleColor: "4A5568" }); }
}

// ---------- 8 Bulk validation ----------
if (D.cox_tab) {
section("8  Bulk validation: TCGA LGG + GBM", "674 tumours (PanCancer Atlas, cBioPortal API); niche signatures scored as mean z of log2 RSEM; WHO 2021 groups from IDH/1p19q subtype");
imgSlide("Niche scores track WHO type and grade", "Hypoxic-MES (N6) and myeloid (N3) scores rise from oligodendroglioma → IDH-mut astrocytoma G2–G4 → GBM IDH-wt (Kruskal p = 2e-70 and 1e-50); neuronal and oligodendroglial scores fall. Astro G4 n = 6. Source: 11_bulk_tcga_validation.py",
  "Files: 04_results/11_bulk_TCGA/scores_by_type_grade.png, score_means_by_type.csv, score_type_grade_tests.csv", D.img.bulk_box, "Bulk scores by type");
imgSlide("Hypoxic-MES score separates survival across gliomas, not within type", "Median split. Across all gliomas the difference mostly reflects WHO type/grade; within IDH-mut astrocytoma (p = 0.71) and GBM IDH-wt (p = 0.92) there is no separation. Source: 11_bulk_tcga_validation.py",
  "File: 04_results/11_bulk_TCGA/KM_HypoxicMES_N6.png", D.img.bulk_km, "Bulk KM");
{ const s = fig("Survival models: niche-6 score and ANGPTL4", "HR per SD. 'All gliomas' adjusted for age, grade and WHO type. '+ hypoxia program' additionally adjusts for the Gavish hypoxia meta-program. Unadjusted for multiple testing; needs CGGA replication.",
    "Data: 04_results/11_bulk_TCGA/cox_results.csv (also PFS, Buffa hypoxia and MES-program adjustments). Native table.");
  const head = ["Group", "Variable", "Adjusted for", "n", "Events", "HR per SD (95% CI)", "p"]
    .map(h => ({ text: h, options: { bold: true, fill: { color: THEME.colors.dk2 }, color: "FFFFFF", align: "center" } }));
  const rows = [head].concat(D.cox_tab.map((r, i) => r.map((v, j) => ({ text: v, options: { align: j < 3 ? "left" : "center",
    bold: j === 6 && parseFloat(v) < 0.05, color: j === 6 && parseFloat(v) < 0.05 ? THEME.colors.accent1 : "1E1E24", fill: { color: i % 2 ? "FFFFFF" : THEME.colors.lt2 } } }))));
  s.addTable(rows, { x: 0.6, y: 1.15, w: 12.1, colW: [2.4, 1.9, 2.2, 0.8, 0.9, 2.6, 1.3], rowH: 0.23, fontSize: 10, fontFace: "Calibri",
    border: { type: "solid", pt: 0.5, color: "FFFFFF" }, margin: [0, 0.05, 0, 0.05], valign: "middle", objectName: "Cox table" }); }
}

// ---------- 8b CGGA + meta ----------
if (D.meta_tab) {
imgSlide("CGGA replicates the type and grade pattern", "CGGA-693 + CGGA-325 primary tumours. N6 and N3 rise from oligodendroglioma through IDH-mut astrocytoma G2→G4 to GBM IDH-wt; N1 and N10 fall. Source: 13_bulk_cgga_validation.py",
  "File: 04_results/13_bulk_CGGA/cgga_scores_by_type_grade.png", D.img.cgga_box, "CGGA scores");
imgSlide("In CGGA, the N6 score separates survival within IDH-mut astrocytoma", "Median split within each cohort (primary tumours). Astrocytoma p = 0.0035 (not grade-adjusted here; see Cox meta-analysis); GBM IDH-wt p = 0.81. Source: 13_bulk_cgga_validation.py",
  "File: 04_results/13_bulk_CGGA/cgga_KM_HypoxicMES_N6.png", D.img.cgga_km, "CGGA KM");
imgSlide("Meta-analysis across three cohorts (1,159 patients)", "N6 is prognostic in lower-grade gliomas, independent of age and grade, and beyond the generic hypoxia program, but not beyond the MES program and not within GBM. The neuronal niche is protective. ANGPTL4 alone does not replicate. Source: 13b_meta_forest.py",
  "File: 04_results/13_bulk_CGGA/meta_forest.png; table meta_TCGA_CGGA.csv", D.img.meta_forest, "Meta forest");
{ const s = fig("Meta-analysis table: overall survival", "Fixed-effect inverse-variance meta-analysis of per-SD Cox HRs (adjusted for age; + grade within lower-grade types; + type in 'All gliomas'). Per-cohort HRs: T = TCGA, C693/C325 = CGGA. Oligodendroglioma: TCGA + CGGA-693 only.",
    "Data: 04_results/13_bulk_CGGA/meta_TCGA_CGGA.csv. Native table.");
  const head = ["Variable", "Group", "n", "Events", "HR meta (95% CI)", "p", "I²", "Per-cohort HR"]
    .map(h => ({ text: h, options: { bold: true, fill: { color: THEME.colors.dk2 }, color: "FFFFFF", align: "center" } }));
  const rows = [head].concat(D.meta_tab.map((r, i) => r.map((v, j) => ({ text: v, options: { align: j < 2 || j === 7 ? "left" : "center",
    bold: j === 5 && parseFloat(v) < 0.05, color: j === 5 && parseFloat(v) < 0.05 ? THEME.colors.accent1 : "1E1E24", fill: { color: i % 2 ? "FFFFFF" : THEME.colors.lt2 } } }))));
  s.addTable(rows, { x: 0.5, y: 1.15, w: 12.3, colW: [1.7, 2.1, 0.7, 0.8, 2.0, 0.9, 0.6, 3.5], rowH: 0.24, fontSize: 10, fontFace: "Calibri",
    border: { type: "solid", pt: 0.5, color: "FFFFFF" }, margin: [0, 0.05, 0, 0.05], valign: "middle", objectName: "Meta table" }); }
}

// ---------- 8c Progression case ----------
if (D.img.prog) {
imgSlide("Case study: one IDH-mut astrocytoma sampled in 2015, 2017 and 2019 (GSE318238)", "Niche-signature scores per spot (z across all three sections). A structured N6 + myeloid territory appears only in 2017; the 2019 section is mostly white matter (oligodendrocyte 62%), so it is not a tumour-core comparison. Grade per time point unconfirmed; n = 1 patient. Source: 14b_progression_scores.py",
  "Files: 04_results/14_progression_GSE318238/. Copy-number compartments failed (probe-based CytAssist vs 3' reference) and niches were not projected.", [D.img.prog[0]], "Progression case");
}

// ---------- 8d ICB exploratory ----------
if (D.img.icb) {
imgSlide("Exploratory: niche scores after neoadjuvant anti-PD-1 (GSE121810)", "Recurrent GBM, adjuvant (n = 10) vs neoadjuvant (n = 12) pembrolizumab (22 of 29 tumours in NCBI's processed table). Positive controls trend as reported (IFN-γ/T-cell up, p = 0.22; cell cycle down, p = 0.11); no niche score differs (all p > 0.3). Underpowered; outcome data not in GEO. Source: 15_icb_gse121810.py",
  "Files: 04_results/15_ICB_GSE121810/. Mann-Whitney, BH FDR across 16 scores (all FDR > 0.77).", D.img.icb, "ICB exploratory");
}

// ---------- 8e IvyGAP anatomical validation ----------
if (D.img.ivy) {
imgSlide("IvyGAP: niche 6 maps to pseudopalisading tissue around necrosis", "122 laser-microdissected samples, 10 GBMs. N6 highest in PAN (β +1.39 vs CT, FDR 2e-30; top region in 8/8 tumours); ANGPTL4 β +1.90; 8 of 9 tier-1 genes peak in PAN (VMP1 in MVP); ANGPTL4 integrin receptors peak in MVP. N6 and ANGPTL4 track hypoxia (r 0.90, 0.86), so this localises the niche but cannot separate it from hypoxia. Source: 17_ivygap_validation.py",
  "Files: 04_results/17_IvyGAP/ (scores, region means, within-tumour contrasts). OLS with tumour fixed effect; paired Wilcoxon on tumour means; BH.", D.img.ivy, "IvyGAP regions");
}

// ---------- 8f Hypoxia adjustment + GBmap source ----------
if (D.img.S6) {
imgSlide("ANGPTL4 niche-6 enrichment survives hypoxia adjustment", "Malignant-rich spots, 15 sections. + hypoxia program (leave-one-out): ANGPTL4 β 0.47 → 0.29 (FDR 0.01, 14/15 sections); VEGFA, ADM, PGK1, LDHA lose enrichment. + hypoxia + state: ANGPTL4 is rank 1 of 1,041 up genes (p 7e-4), no gene passes transcriptome-wide FDR. Source: 09c_ncdeg_hypoxia.py",
  "Files: 04_results/09_NCDEGs/hypoxia_adjusted/. r(niche 6, hypoxia score) median 0.43 within sections: conservative test.", D.img.S6, "Hypoxia adjustment");
}
if (D.img.S7) {
imgSlide("GBmap: ANGPTL4 source and receptor-expressing cells", "ANGPTL4: hypoxic MES-like 34% vs hypoxia-independent MES-like 21% of cells (p 0.02; VEGFA 80 vs 46%), also hypoxic monocytes/BDM, astrocytes, mural cells. Receptors: ITGA5 arterial/tip endothelium; ITGB1, SDC2, CDH11 pericytes/SMC; SDC4, CD44 MES-like, astrocytes, myeloid. Source: 18_gbmap_angptl4_source.py",
  "Files: 04_results/18_GBmap_ANGPTL4/. GBmap core reference subsample (≤400 cells per level-3 type, 7,368 cells, 16 studies).", D.img.S7, "GBmap ANGPTL4");
}

if (D.img.S8) {
imgSlide("CPTAC: protein-level support (99 GBMs)", "Protein vs mRNA: niche signatures ρ 0.81–0.90, ANGPTL4 0.80; VMP1 0.25, NRN1 0.20 (weak). ANGPTL4 protein ~ N6 (without ANGPTL4) ρ 0.64, hypoxia 0.76, ITGA5 0.59. Protein model: N6 β 0.16 (ns) with hypoxia; 0.38 (p 0.015) after + MES. mRNA: N6 β 0.35 (p 0.006). Bulk cannot separate niche from hypoxia. Source: 19_cptac_figure.py",
  "Data: cBioPortal gbm_cptac_2021, computed in-browser (19a_fetch_cptac_cbioportal.js); results 04_results/19_CPTAC/cptac_results.json.", D.img.S8, "CPTAC protein");
}

if (D.img.S9) {
imgSlide("GLASS: niche signatures through recurrence (214 pairs)", "Type ordering reproduced (Kruskal p 7e-23). Astrocytoma (36 pairs): N6 up in 61% (p 0.10); ANGPTL4 up in 69% (p 0.03), 80% with grade increase (p 0.01; hypoxia-adjusted p 0.02); none pass FDR. IDH-wt: oligodendroglial gain at recurrence (FDR 0.002), hypoxia-adjusted N6 slightly down (FDR 0.006). Primary N6 not prognostic in GLASS (HR 1.10, p 0.38). Source: 20_glass_longitudinal.py",
  "Data: Synapse syn17038081 (gene TPM, clinical_surgeries, clinical_cases, analysis_rna_silver_set). Results: 04_results/20_GLASS/.", D.img.S9, "GLASS longitudinal");
}

// ---------- 9 Targets ----------
if (D.target_tab) {
section("9  Niche-resolved targets", "Niche-6 NCDEGs and ligand/receptor genes × state-adjusted effect × Open Targets (GBM association, tractability, drugs, DepMap glioma lines)");
{ const s = fig("Prioritised niche-6 targets", "Tier 1 = niche-specific, still up after tumour-state adjustment (p < 0.05) and tractable. OT = Open Targets GBM association score (low = little prior GBM evidence). DepMap = median Chronos effect in glioma lines (< −1 = dependency). Source: 12_target_prioritisation.py",
    "Data: 04_results/12_targets/niche6_target_prioritisation.csv (48 genes, full drug lists and notes). Native table.");
  const head = ["Gene", "Tier", "Niche-6 role", "β state-adj.", "p state-adj.", "OT GBM score", "DepMap glioma", "Most advanced drug"]
    .map(h => ({ text: h, options: { bold: true, fill: { color: THEME.colors.dk2 }, color: "FFFFFF", align: "center" } }));
  const rows = [head].concat(D.target_tab.map((r, i) => r.map((v, j) => ({ text: v, options: { align: j === 0 || j === 7 ? "left" : "center",
    bold: j === 0 && r[1] === "1", color: "1E1E24", fill: { color: r[1] === "1" ? (i % 2 ? "FFFFFF" : THEME.colors.lt2) : "F7F7F9" } } }))));
  s.addTable(rows, { x: 0.5, y: 1.15, w: 12.3, colW: [1.2, 0.7, 1.3, 1.15, 1.15, 1.2, 1.3, 4.3], rowH: 0.24, fontSize: 10, fontFace: "Calibri",
    border: { type: "solid", pt: 0.5, color: "FFFFFF" }, margin: [0, 0.05, 0, 0.05], valign: "middle", objectName: "Target table" }); }
}

// ---------- Article v2 ----------
if (D.img.fig5v2) {
section("Article v2  ANGPTL4 axis and figure renumbering", "Manuscript reorganised around the state- and hypoxia-independent ANGPTL4 signal");
imgSlide("Figure 5 (v2): the ANGPTL4 niche-6 axis", "Hypoxia adjustment (rank 1 of 1,041 in hypoxia + state model), IvyGAP PAN localisation, GBmap source (hypoxic MES-like) and receptors (ITGA5 endothelium; ITGB1/SDC2/CDH11 pericytes), CPTAC protein ρ 0.64 with N6. Source: 21_figure5_angptl4_axis.py",
  "File: 06_manuscript/figures/Figure_5.png (v2 numbering). Panels: A hypoxia adjustment, B IvyGAP boxplots, C IvyGAP heatmap, D GBmap % ANGPTL4, E receptor dot plot, F CPTAC protein correlations.", D.img.fig5v2, "Figure 5 v2");
if (D.renum) { const s = fig("Figure renumbering: v1 to v2 manuscript", "Old file name (v1) and new number in the reorganised manuscript. Supplementary tables re-ordered the same way (see README sheet).",
    "Map stored in 06_manuscript/renumbering_map.json.");
  const head = ["Old (v1)", "New (v2)"].map(h => ({ text: h, options: { bold: true, fill: { color: THEME.colors.dk2 }, color: "FFFFFF", align: "center" } }));
  const rows = [head].concat(D.renum.map((r, i) => r.map(v => ({ text: v, options: { align: "center", color: "1E1E24", fill: { color: i % 2 ? "FFFFFF" : THEME.colors.lt2 } } }))));
  s.addTable(rows, { x: 3.4, y: 1.15, w: 6.5, colW: [3.6, 2.9], rowH: 0.29, fontSize: 11, fontFace: "Calibri",
    border: { type: "solid", pt: 0.5, color: "FFFFFF" }, valign: "middle", objectName: "Renumbering table" }); }
}

// ---------- Appendix ----------
section("Appendix  Original and superseded figures", "Kept for the record; not used for the main results");
imgSlide("Niche profiles: original matplotlib figure (final version)", "Same data as the editable niche table in section 5. Source: 06_lcp_niches_deconv.py",
  "File: 04_results/06_niches_deconv/niche_profiles.png", D.img.niche_prof, "Niche profiles deconvolution PNG");
imgSlide("Superseded: niches from program scores", "First version, built from program scores instead of deconvolution. MES-like/hypoxic niche was already GBM-specific here. Source: 04_lcp_niches.py",
  "File: 04_results/04_LCP_niches/niche_profiles.png", D.img.old_prog, "Superseded niches program score");
imgSlide("Superseded: program-score niches with Harmony", "Program-score niches after Harmony cohort correction. Source: 04_lcp_niches.py (Harmony)",
  "File: 04_results/04_LCP_niches_harmony/niche_profiles.png", D.img.old_harm, "Superseded niches harmony");

(async () => { await pres.writeFile({ fileName: OUT }); await applyTheme(OUT, THEME); console.log("written", OUT); })();
