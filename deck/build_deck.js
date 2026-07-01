/*
 * ICIpredict results deck — built with pptxgenjs.
 * Palette and figures are matched to the analysis figures (teal/orange/red).
 */
const pptxgen = require("pptxgenjs");
const sharp = require("sharp");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const path = require("path");
const {
  FaDna, FaUserInjured, FaChartBar, FaFlask, FaBalanceScale, FaShieldVirus,
  FaLightbulb, FaBookOpen, FaExclamationTriangle, FaArrowRight, FaMicroscope,
  FaCheckCircle, FaDatabase,
} = require("react-icons/fa");

const FIG = "C:/Projects/ICIpredict/deck/figs";

// ---- palette (matches results figures) ----
const TEAL = "0E7C86", ORANGE = "E07B39", RED = "C0392B";
const SLATE = "14252B", PANEL = "EEF4F5", INK = "1F2D33", MUTED = "5C6B72";
const WHITE = "FFFFFF", LINE = "D5E0E2";
const HFONT = "Cambria", BFONT = "Calibri";

const pres = new pptxgen();
pres.defineLayout({ name: "W", width: 13.333, height: 7.5 });
pres.layout = "W";
pres.author = "ICIpredict";
pres.title = "Predicting immunotherapy response from genomic data";
const W = 13.333, H = 7.5;

// ---- icon helper ----
async function icon(Comp, color, size = 256) {
  const svg = ReactDOMServer.renderToStaticMarkup(
    React.createElement(Comp, { color: "#" + color, size: String(size) }));
  const png = await sharp(Buffer.from(svg)).png().toBuffer();
  return "image/png;base64," + png.toString("base64");
}
async function imgDims(p) {
  const m = await sharp(p).metadata();
  return { w: m.width, h: m.height };
}
// place an image into a box (contain), return actual placed rect
async function fitImage(slide, file, bx, by, bw, bh, opts = {}) {
  const { w, h } = await imgDims(path.join(FIG, file));
  const r = Math.min(bw / w, bh / h);
  const iw = w * r, ih = h * r;
  const ix = bx + (bw - iw) / 2, iy = by + (bh - ih) / 2;
  slide.addImage({ path: path.join(FIG, file), x: ix, y: iy, w: iw, h: ih, ...opts });
  return { x: ix, y: iy, w: iw, h: ih };
}

const shadow = () => ({ type: "outer", color: "000000", blur: 7, offset: 3, angle: 90, opacity: 0.12 });

function footer(slide, n) {
  slide.addText([
    { text: "ICIpredict", options: { bold: true, color: TEAL } },
    { text: "   ·   Samstein et al. 2019  (cBioPortal tmb_mskcc_2018, n = 1,630)", options: { color: MUTED } },
  ], { x: 0.55, y: H - 0.42, w: 10, h: 0.3, fontSize: 9, fontFace: BFONT, margin: 0 });
  slide.addText(String(n), { x: W - 1.0, y: H - 0.42, w: 0.5, h: 0.3, fontSize: 9,
    color: MUTED, align: "right", fontFace: BFONT });
}

async function sectionHeader(slide, kicker, title, IconComp) {
  slide.addImage({ data: await icon(IconComp, TEAL), x: 0.55, y: 0.5, w: 0.42, h: 0.42 });
  slide.addText(kicker.toUpperCase(), { x: 1.12, y: 0.45, w: 11, h: 0.25, fontSize: 11,
    color: TEAL, bold: true, charSpacing: 2, fontFace: BFONT, margin: 0 });
  slide.addText(title, { x: 1.12, y: 0.68, w: 11.6, h: 0.62, fontSize: 25, bold: true,
    color: INK, fontFace: HFONT, margin: 0 });
}

// ============================================================ SLIDE 1: TITLE
async function titleSlide() {
  const s = pres.addSlide();
  s.background = { color: SLATE };
  // faint DNA motif
  s.addImage({ data: await icon(FaDna, "1E3A42"), x: 9.7, y: 0.6, w: 3.4, h: 6.3, transparency: 22 });
  s.addText("GENOMICS × IMMUNO-ONCOLOGY", { x: 0.7, y: 1.35, w: 9, h: 0.4, fontSize: 13,
    color: ORANGE, bold: true, charSpacing: 3, fontFace: BFONT, margin: 0 });
  s.addText("Predicting whether a tumor will\nrespond to immunotherapy", { x: 0.66, y: 1.8, w: 9.6,
    h: 1.9, fontSize: 41, bold: true, color: WHITE, fontFace: HFONT, lineSpacingMultiple: 1.02, margin: 0 });
  s.addText("An interpretable, cancer-type-aware survival model that integrates tumor mutational burden with gene-level mutations — and a test of whether genes add value beyond TMB.",
    { x: 0.7, y: 3.95, w: 8.7, h: 1.0, fontSize: 15, color: "CFE0E3", fontFace: BFONT, margin: 0 });
  // meta chips
  const chips = [["1,630", "ICI-treated patients"], ["11", "cancer types"],
    ["69", "gene features"], ["0.55 → 0.64", "C-index: TMB → integrated"]];
  let cx = 0.7;
  for (const [big, lab] of chips) {
    const cw = 0.55 + big.length * 0.16 + 1.3;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: cx, y: 5.35, w: cw, h: 0.95,
      fill: { color: "1B3138" }, rectRadius: 0.08, line: { color: "26454E", width: 1 } });
    s.addText(big, { x: cx + 0.2, y: 5.45, w: cw - 0.3, h: 0.42, fontSize: 20, bold: true,
      color: ORANGE, fontFace: HFONT, margin: 0 });
    s.addText(lab, { x: cx + 0.2, y: 5.86, w: cw - 0.3, h: 0.32, fontSize: 10.5,
      color: "AFC6CB", fontFace: BFONT, margin: 0 });
    cx += cw + 0.25;
  }
  s.addText("Data: Samstein et al., Nature Genetics 2019 · pulled live from the cBioPortal REST API",
    { x: 0.7, y: 6.75, w: 12, h: 0.3, fontSize: 10.5, color: "8FAAB0", italic: true, fontFace: BFONT, margin: 0 });
  s.addNotes("Title. The clinical question: can we predict immunotherapy (immune checkpoint inhibitor) response from a tumor's DNA? We build an integrated, cancer-type-aware survival model on 1,630 ICI-treated patients and rigorously test whether gene-level data adds anything beyond tumor mutational burden (TMB).");
}

// ============================================================ SLIDE 2: QUESTION
async function questionSlide() {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  await sectionHeader(s, "The problem", "TMB is the FDA biomarker — but it is far from enough", FaMicroscope);

  s.addText([
    { text: "Checkpoint inhibitors (anti–PD-1 / CTLA-4) transform outcomes — but only for some patients.", options: { breakLine: true, bullet: true } },
    { text: "Tumor mutational burden (TMB ≥ 10 mut/Mb) is the FDA's tissue-agnostic biomarker, yet many high-TMB tumors still fail to respond (Rizvi 2015).", options: { breakLine: true, bullet: true } },
    { text: "The literature points to a multi-signal, cancer-type-specific answer: TMB + neoantigens + resistance genes (B2M, JAK1/2) + cancer context.", options: { bullet: true } },
  ], { x: 0.7, y: 1.7, w: 6.7, h: 2.6, fontSize: 14.5, color: INK, fontFace: BFONT,
       paraSpaceAfter: 10, lineSpacingMultiple: 1.03 });

  // research question card
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.7, y: 4.55, w: 6.7, h: 2.25,
    fill: { color: PANEL }, rectRadius: 0.1, shadow: shadow() });
  s.addText("RESEARCH QUESTION", { x: 1.0, y: 4.78, w: 6, h: 0.3, fontSize: 11, bold: true,
    color: TEAL, charSpacing: 2, fontFace: BFONT, margin: 0 });
  s.addText([
    { text: "Does integrating TMB with specific gene mutations and cancer-type context predict overall survival after immunotherapy ", options: {} },
    { text: "better than TMB alone", options: { bold: true, color: TEAL } },
    { text: " — and can gene context explain the high-TMB non-responders?", options: {} },
  ], { x: 1.0, y: 5.12, w: 6.1, h: 1.55, fontSize: 15.5, color: INK, fontFace: HFONT,
       italic: true, lineSpacingMultiple: 1.05, margin: 0 });

  // right: stat cards
  const cards = [
    [FaUserInjured, "1,630", "patients treated with\nimmune checkpoint inhibitors"],
    [FaDna, "11", "cancer types — pan-cancer,\nnot one tissue"],
    [FaChartBar, "810", "deaths (49.7%) — mature\noverall-survival follow-up"],
  ];
  let y = 1.75;
  for (const [Ic, big, lab] of cards) {
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 7.75, y, w: 5.0, h: 1.55,
      fill: { color: WHITE }, rectRadius: 0.09, line: { color: LINE, width: 1 }, shadow: shadow() });
    s.addShape(pres.shapes.OVAL, { x: 8.05, y: y + 0.42, w: 0.72, h: 0.72, fill: { color: PANEL } });
    s.addImage({ data: await icon(Ic, TEAL), x: 8.19, y: y + 0.56, w: 0.44, h: 0.44 });
    s.addText(big, { x: 9.0, y: y + 0.2, w: 3.6, h: 0.6, fontSize: 30, bold: true, color: ORANGE,
      fontFace: HFONT, margin: 0, valign: "middle" });
    s.addText(lab, { x: 9.0, y: y + 0.82, w: 3.55, h: 0.6, fontSize: 11.5, color: MUTED,
      fontFace: BFONT, margin: 0 });
    y += 1.72;
  }
  footer(s, 2);
  s.addNotes("The clinical problem and our question. TMB is approved but imperfect — high-TMB non-responders exist. We ask whether a multivariable, cancer-type-aware genomic model beats TMB alone on overall survival, using a large pan-cancer cohort.");
}

// ============================================================ SLIDE 3: DATA & METHODS
async function dataSlide() {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  await sectionHeader(s, "Data & approach", "A reproducible, leakage-controlled survival pipeline", FaDatabase);

  s.addText([
    { text: "Cohort", options: { bold: true, color: TEAL, breakLine: true } },
    { text: "Samstein 2019 (MSK-IMPACT), pulled live from the cBioPortal REST API; OS endpoint.", options: { breakLine: true } },
    { text: "Features", options: { bold: true, color: TEAL, breakLine: true } },
    { text: "log-TMB, age, sex, biopsy site, ICI drug class, cancer type, + 69 gene-mutation flags restricted to the IMPACT341 core (a 0 = wild-type, never “unsequenced”).", options: { breakLine: true } },
    { text: "Models", options: { bold: true, color: TEAL, breakLine: true } },
    { text: "Baselines (FDA cutoff, Samstein type-specific rule, Cox on TMB) vs integrated models (elastic-net Cox, random survival forest).", options: { breakLine: true } },
    { text: "Rigor", options: { bold: true, color: TEAL, breakLine: true } },
    { text: "Repeated stratified 5-fold CV; the elastic-net penalty is tuned inside each fold (nested CV); Harrell C, IPCW-C, time-AUC; paired bootstrap with Holm correction.", options: {} },
  ], { x: 0.7, y: 1.7, w: 5.85, h: 4.9, fontSize: 13, color: INK, fontFace: BFONT,
       paraSpaceAfter: 7, lineSpacingMultiple: 1.0 });

  // right: fig1
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 6.85, y: 1.7, w: 5.95, h: 4.55,
    fill: { color: WHITE }, rectRadius: 0.08, line: { color: LINE, width: 1 }, shadow: shadow() });
  await fitImage(s, "fig1_tmb_by_cancer_type.png", 7.0, 1.85, 5.65, 3.95);
  s.addText("TMB spans orders of magnitude across cancer types — a single global cutoff misfits most. This is why cancer-type context matters.",
    { x: 7.0, y: 5.78, w: 5.65, h: 0.5, fontSize: 10, color: MUTED, italic: true, fontFace: BFONT, margin: 0 });
  footer(s, 3);
  s.addNotes("Data and methods. Everything is pulled programmatically from cBioPortal for reproducibility. We restrict gene flags to the IMPACT341 core so a 0 truly means wild-type. Crucially, the elastic-net penalty is tuned inside each CV fold (nested CV) so the reported numbers carry no hyperparameter-selection optimism.");
}

// ============================================================ SLIDE 4: RESULT 1 — discrimination
async function result1Slide() {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  await sectionHeader(s, "Result 1 — discrimination", "Integrated model beats TMB — but context, not genes, drives the gain", FaChartBar);

  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.6, y: 1.65, w: 7.5, h: 5.1,
    fill: { color: WHITE }, rectRadius: 0.08, line: { color: LINE, width: 1 }, shadow: shadow() });
  await fitImage(s, "fig2_cindex_comparison.png", 0.75, 1.8, 7.2, 4.8);

  // right column: big stat + decomposition
  s.addText("Cross-validated C-index", { x: 8.35, y: 1.75, w: 4.4, h: 0.3, fontSize: 12,
    color: MUTED, bold: true, fontFace: BFONT, margin: 0 });
  s.addText([
    { text: "0.55", options: { color: MUTED } },
    { text: "  →  ", options: { color: INK } },
    { text: "0.64", options: { color: ORANGE } },
  ], { x: 8.3, y: 2.05, w: 4.5, h: 0.9, fontSize: 46, bold: true, fontFace: HFONT, margin: 0 });
  s.addText("TMB alone  →  integrated model", { x: 8.35, y: 2.95, w: 4.4, h: 0.3, fontSize: 11.5,
    color: MUTED, fontFace: BFONT, margin: 0 });

  const rows = [
    ["TMB alone", "0.55", MUTED],
    ["+ clinical & cancer type", "0.62", TEAL],
    ["+ TMB on top", "0.64", TEAL],
    ["+ genes (full model)", "0.64", ORANGE],
  ];
  let y = 3.6;
  for (const [lab, val, col] of rows) {
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 8.3, y, w: 4.5, h: 0.62,
      fill: { color: PANEL }, rectRadius: 0.06 });
    s.addText(lab, { x: 8.5, y, w: 3.2, h: 0.62, fontSize: 12, color: INK, fontFace: BFONT,
      valign: "middle", margin: 0 });
    s.addText(val, { x: 11.5, y, w: 1.1, h: 0.62, fontSize: 17, bold: true, color: col,
      align: "right", fontFace: HFONT, valign: "middle", margin: 0 });
    y += 0.72;
  }
  s.addText("Most of the lift over TMB is clinical + cancer-type context — an honest decomposition, not a black box.",
    { x: 8.3, y: y + 0.05, w: 4.5, h: 0.8, fontSize: 11, color: MUTED, italic: true, fontFace: BFONT, margin: 0 });
  footer(s, 4);
  s.addNotes("Result 1. The integrated model lifts the C-index from 0.55 (TMB alone) to 0.64. But the decomposition is the honest story: clinical + cancer-type covariates already reach 0.62, TMB adds a bit, and genes add only a little more. We don't oversell the genomics.");
}

// ============================================================ SLIDE 5: RESULT 2 — do genes add value
async function result2Slide() {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  await sectionHeader(s, "Result 2 — the honest test", "Do genes add value beyond TMB + clinical? Barely — and we say so.", FaBalanceScale);

  // comparison table (native)
  const head = (t) => ({ text: t, options: { bold: true, color: WHITE, fill: { color: TEAL }, fontFace: BFONT, align: "center", valign: "middle" } });
  const cell = (t, o = {}) => ({ text: t, options: { color: o.color || INK, fontFace: BFONT, valign: "middle", align: o.align || "center", bold: o.bold, fill: o.fill } });
  const data = [
    [head("Comparison (ΔC)"), head("ΔC"), head("95% CI"), head("Holm p"), head("Verdict")],
    [cell("Integrated  vs  TMB alone", { align: "left", bold: true }), cell("+0.099", { color: TEAL, bold: true }), cell("[+0.08, +0.12]"), cell("0.004"), cell("robust ✓", { color: TEAL, bold: true })],
    [cell("Integrated  vs  TMB + clinical", { align: "left", bold: true }), cell("+0.008", { color: ORANGE, bold: true }), cell("[−0.003, +0.017]"), cell("0.26"), cell("not sig.", { color: ORANGE, bold: true })],
    [cell("Forest  vs  TMB + clinical", { align: "left", bold: true }), cell("+0.005", { color: MUTED }), cell("[−0.005, +0.015]"), cell("0.32"), cell("null", { color: MUTED, bold: true })],
  ];
  s.addTable(data, { x: 0.7, y: 1.8, w: 7.5, h: 2.5, colW: [2.75, 1.0, 1.55, 0.9, 1.3],
    rowH: [0.5, 0.6, 0.6, 0.6], fontSize: 11.5, border: { pt: 0.5, color: LINE },
    align: "center", valign: "middle" });

  s.addText([
    { text: "The win over TMB alone is real and robust. But the pure genomic increment over TMB + clinical is ", options: {} },
    { text: "tiny (ΔC ≈ +0.008), its CI crosses zero, and it is not significant after multiplicity correction", options: { bold: true } },
    { text: " (a paired per-fold test detects a consistent but negligible +0.006).", options: {} },
  ], { x: 0.7, y: 4.5, w: 7.5, h: 1.1, fontSize: 13, color: INK, fontFace: BFONT, lineSpacingMultiple: 1.05 });

  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.7, y: 5.65, w: 7.5, h: 1.1,
    fill: { color: PANEL }, rectRadius: 0.09 });
  s.addImage({ data: await icon(FaCheckCircle, TEAL), x: 0.95, y: 5.92, w: 0.5, h: 0.5 });
  s.addText("On this cohort, gene flags add no robust gain beyond TMB + clinical + cancer type. Reporting that honestly — instead of overclaiming — is the scientifically correct call.",
    { x: 1.65, y: 5.75, w: 6.4, h: 0.9, fontSize: 12, color: INK, italic: true, fontFace: BFONT, valign: "middle", margin: 0 });

  // right: AUC over time figure
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 8.45, y: 1.8, w: 4.35, h: 4.95,
    fill: { color: WHITE }, rectRadius: 0.08, line: { color: LINE, width: 1 }, shadow: shadow() });
  await fitImage(s, "fig3_auc_over_time.png", 8.55, 2.4, 4.15, 3.9);
  s.addText("Time-dependent AUC: integrated models lead at every horizon.", { x: 8.55, y: 1.95, w: 4.15, h: 0.45,
    fontSize: 10.5, color: MUTED, italic: true, fontFace: BFONT, margin: 0 });
  footer(s, 5);
  s.addNotes("Result 2 is the honest test. We compare four model pairs with a paired bootstrap and Holm multiplicity correction. The big advantage over TMB alone is robust. The pure genomic increment over TMB+clinical is only ~+0.01, Coxnet-specific, and marginal — so we label it suggestive, not established.");
}

// ============================================================ SLIDE 6: RESULT 3 — masking
async function result3Slide() {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  await sectionHeader(s, "Result 3 — the integration insight", "Why high-TMB non-responders evade the eye: TMB masks resistance genes", FaShieldVirus);

  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.55, y: 1.6, w: 8.05, h: 5.2,
    fill: { color: WHITE }, rectRadius: 0.08, line: { color: LINE, width: 1 }, shadow: shadow() });
  await fitImage(s, "fig7_resistance_tmb_confounding.png", 0.65, 1.72, 7.85, 4.95);

  // right explanation + stats
  s.addText("Resistance mutations hide inside high TMB", { x: 8.8, y: 1.7, w: 4.1, h: 0.7,
    fontSize: 16, bold: true, color: INK, fontFace: HFONT, margin: 0 });
  s.addText([
    { text: "STK11/KEAP1/B2M mutations occur preferentially in ", options: {} },
    { text: "high-TMB", options: { bold: true, color: TEAL } },
    { text: " tumors (9.8 vs 5.2 mut/Mb, p≈10⁻²⁸). Because high TMB is favorable, their harm is ", options: {} },
    { text: "masked", options: { bold: true, color: RED } },
    { text: ".", options: {} },
  ], { x: 8.8, y: 2.45, w: 4.15, h: 1.5, fontSize: 12.5, color: INK, fontFace: BFONT,
       lineSpacingMultiple: 1.04, margin: 0 });

  const hr = [
    ["STK11 (NSCLC)", "1.31", "1.42", "0.041"],
    ["KEAP1 (NSCLC)", "1.35", "1.54", "0.014"],
    ["STK11 (pan-cancer*)", "1.34", "1.50", "0.006"],
  ];
  s.addText([
    { text: "Hazard ratio   ", options: { bold: true } },
    { text: "without → with TMB", options: { italic: true, color: MUTED } },
  ], { x: 8.8, y: 3.95, w: 4.15, h: 0.3, fontSize: 11, color: INK, fontFace: BFONT, margin: 0 });
  let y = 4.3;
  for (const [g, b, a, p] of hr) {
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 8.8, y, w: 4.15, h: 0.6, fill: { color: PANEL }, rectRadius: 0.06 });
    s.addText(g, { x: 8.95, y, w: 2.0, h: 0.6, fontSize: 11, color: INK, bold: true, fontFace: BFONT, valign: "middle", margin: 0 });
    s.addText([{ text: b + "  ", options: { color: MUTED } }, { text: "→ " + a, options: { color: RED, bold: true } }],
      { x: 10.7, y, w: 1.55, h: 0.6, fontSize: 12.5, fontFace: HFONT, valign: "middle", align: "right", margin: 0 });
    s.addText("p=" + p, { x: 12.25, y, w: 0.62, h: 0.6, fontSize: 9.5, color: MUTED, fontFace: BFONT, valign: "middle", align: "right", margin: 0 });
    y += 0.7;
  }
  s.addText("Holding cancer type fixed and adding TMB unmasks STK11/KEAP1 as resistance markers — the mechanism behind “TMB necessary but not sufficient.”   *type-adjusted",
    { x: 8.8, y: y + 0.02, w: 4.15, h: 0.95, fontSize: 10, color: MUTED, italic: true, fontFace: BFONT, margin: 0 });
  footer(s, 6);
  s.addNotes("Result 3 — the novel insight. Rizvi 2015 saw high-TMB tumors that didn't respond but couldn't explain it by eye. We show why: resistance-gene mutations co-occur with high (favorable) TMB, so their harmful effect is masked by positive confounding. Hold cancer type fixed and add TMB to the model, and STK11/KEAP1 emerge as resistance markers (HR ~1.4-1.5). This is the mechanistic reason an integrated model can beat TMB alone.");
}

// ============================================================ SLIDE 7: WHAT IT LEARNED
async function learnedSlide() {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  await sectionHeader(s, "Interpretability", "What the model learned — and it lines up with the biology", FaFlask);

  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.6, y: 1.65, w: 7.0, h: 5.1,
    fill: { color: WHITE }, rectRadius: 0.08, line: { color: LINE, width: 1 }, shadow: shadow() });
  await fitImage(s, "fig4_coxnet_hazard_ratios.png", 0.72, 1.78, 6.75, 4.85);

  // right: two stacked cards (resistance / benefit)
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 7.85, y: 1.72, w: 4.95, h: 2.4,
    fill: { color: "FBEDE9" }, rectRadius: 0.1, shadow: shadow() });
  s.addImage({ data: await icon(FaExclamationTriangle, RED), x: 8.15, y: 2.0, w: 0.45, h: 0.45 });
  s.addText("Worse survival (resistance)", { x: 8.75, y: 1.97, w: 3.9, h: 0.5, fontSize: 15, bold: true, color: RED, fontFace: HFONT, valign: "middle", margin: 0 });
  s.addText("STK11 · KEAP1 · SMARCA4 · TP53 · glioma context — known cold-tumor / ICI-resistance biology.",
    { x: 8.15, y: 2.6, w: 4.5, h: 1.35, fontSize: 13, color: INK, fontFace: BFONT, lineSpacingMultiple: 1.08, margin: 0 });

  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 7.85, y: 4.32, w: 4.95, h: 2.4,
    fill: { color: "E5F2F1" }, rectRadius: 0.1, shadow: shadow() });
  s.addImage({ data: await icon(FaCheckCircle, TEAL), x: 8.15, y: 4.6, w: 0.45, h: 0.45 });
  s.addText("Better survival (benefit)", { x: 8.75, y: 4.57, w: 3.9, h: 0.5, fontSize: 15, bold: true, color: TEAL, fontFace: HFONT, valign: "middle", margin: 0 });
  s.addText("Higher TMB · melanoma & renal-cell context · VHL — the most ICI-responsive settings in the clinic.",
    { x: 8.15, y: 5.2, w: 4.5, h: 1.35, fontSize: 13, color: INK, fontFace: BFONT, lineSpacingMultiple: 1.08, margin: 0 });
  footer(s, 7);
  s.addNotes("Interpretability. The elastic-net Cox coefficients are biologically sensible: STK11, KEAP1, SMARCA4, TP53 and glioma context push toward worse survival (resistance), while higher TMB, melanoma/renal context and VHL push toward benefit. The model recovers known immuno-oncology biology rather than memorizing noise.");
}

// ============================================================ SLIDE 8: BIBLIOGRAPHY SYNTHESIS
async function biblioSlide() {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  await sectionHeader(s, "Synthesis", "One model that ties the group's five sources together", FaBookOpen);

  const cards = [
    [FaChartBar, "Lee & Samstein — TMB", "We reproduce TMB's signal and show the FDA cutoff and Samstein's type-specific rule both under-perform a cancer-type-aware multivariable model."],
    [FaDna, "Rizvi — mutational landscape", "POLE/MMR track benefit; and we give a concrete reason high-TMB non-responders evade the eye — masking by TMB and cancer type."],
    [FaShieldVirus, "Jamieson & Maker — resistance", "STK11/KEAP1 carry the expected resistance direction once TMB and type are modeled — consistent with acquired-resistance biology."],
  ];
  let x = 0.7;
  for (const [Ic, title, body] of cards) {
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 1.95, w: 3.85, h: 4.4,
      fill: { color: WHITE }, rectRadius: 0.1, line: { color: LINE, width: 1 }, shadow: shadow() });
    s.addShape(pres.shapes.OVAL, { x: x + 0.35, y: 2.35, w: 0.95, h: 0.95, fill: { color: PANEL } });
    s.addImage({ data: await icon(Ic, TEAL), x: x + 0.58, y: 2.58, w: 0.5, h: 0.5 });
    s.addText(title, { x: x + 0.32, y: 3.5, w: 3.25, h: 0.75, fontSize: 15.5, bold: true, color: INK, fontFace: HFONT, margin: 0 });
    s.addText(body, { x: x + 0.32, y: 4.3, w: 3.25, h: 1.95, fontSize: 12.5, color: MUTED, fontFace: BFONT, lineSpacingMultiple: 1.06, margin: 0 });
    x += 4.1;
  }
  footer(s, 8);
  s.addNotes("Synthesis. The project isn't a one-off model — it operationalizes the group's five annotated sources into a single tested hypothesis, reproducing TMB's signal, explaining Rizvi's non-responders, and recovering Jamieson & Maker's resistance genes.");
}

// ============================================================ SLIDE 9: LIMITATIONS / NEXT
async function limitationsSlide() {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  await sectionHeader(s, "Honesty & next steps", "What this does not show — and where it goes next", FaLightbulb);

  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.7, y: 1.85, w: 5.85, h: 4.85,
    fill: { color: "FBEDE9" }, rectRadius: 0.1, shadow: shadow() });
  s.addText("Limitations", { x: 1.0, y: 2.1, w: 5, h: 0.4, fontSize: 17, bold: true, color: RED, fontFace: HFONT, margin: 0 });
  s.addText([
    { text: "Targeted panel, not whole-exome; no PD-L1, RNA / IFN-γ signature, or neoantigen calling.", options: { bullet: true, breakLine: true } },
    { text: "Mutation flag = any non-silent variant present — not demonstrated functional loss.", options: { bullet: true, breakLine: true } },
    { text: "OS, single institution (MSK), retrospective, no ICI-vs-control arm — can't separate predictive from prognostic.", options: { bullet: true, breakLine: true } },
    { text: "Genomic increment is small and Coxnet-specific: suggestive, not established.", options: { bullet: true } },
  ], { x: 1.0, y: 2.6, w: 5.3, h: 3.9, fontSize: 13, color: INK, fontFace: BFONT, paraSpaceAfter: 11, lineSpacingMultiple: 1.03 });

  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 6.85, y: 1.85, w: 5.85, h: 4.85,
    fill: { color: "E5F2F1" }, rectRadius: 0.1, shadow: shadow() });
  s.addText("Next steps", { x: 7.15, y: 2.1, w: 5, h: 0.4, fontSize: 17, bold: true, color: TEAL, fontFace: HFONT, margin: 0 });
  s.addText([
    { text: "External validation on independent ICI whole-exome cohorts (Hugo, Van Allen, Liu).", options: { bullet: true, breakLine: true } },
    { text: "Add PD-L1, IFN-γ expression signature, and neoantigen burden as features.", options: { bullet: true, breakLine: true } },
    { text: "Strictly fold-local gene selection + functional (loss-of-function) annotation.", options: { bullet: true, breakLine: true } },
    { text: "Per-cancer-type models where sample size allows (NSCLC, melanoma, bladder).", options: { bullet: true } },
  ], { x: 7.15, y: 2.6, w: 5.3, h: 3.9, fontSize: 13, color: INK, fontFace: BFONT, paraSpaceAfter: 11, lineSpacingMultiple: 1.03 });
  footer(s, 9);
  s.addNotes("Limitations and next steps. We're explicit about what the data can't show: targeted panel, OS-only, single center, no control arm, and a small genomic increment. Next: external WES validation, richer features (PD-L1, RNA, neoantigens), fold-local + functional gene annotation, and per-cancer models.");
}

// ============================================================ SLIDE 10: TAKEAWAYS
async function takeawaySlide() {
  const s = pres.addSlide();
  s.background = { color: SLATE };
  s.addImage({ data: await icon(FaDna, "1E3A42"), x: 10.0, y: 4.3, w: 3.1, h: 3.1, transparency: 28 });
  s.addText("KEY TAKEAWAYS", { x: 0.7, y: 0.85, w: 9, h: 0.4, fontSize: 13, color: ORANGE, bold: true, charSpacing: 3, fontFace: BFONT, margin: 0 });
  s.addText("What we can — and cannot — claim", { x: 0.66, y: 1.25, w: 11, h: 0.7, fontSize: 30, bold: true, color: WHITE, fontFace: HFONT, margin: 0 });

  const items = [
    ["1", "Integration beats TMB alone.", "C-index 0.55 → 0.64 — but most of that gain is clinical + cancer-type context, honestly decomposed."],
    ["2", "Gene flags add almost nothing — and we say so.", "The increment over TMB + clinical is ~+0.008, CI crosses zero, not significant after multiplicity correction. We don't overclaim."],
    ["3", "TMB masks resistance genes.", "STK11/KEAP1 only reveal their harm once TMB and cancer type are modeled jointly — explaining high-TMB non-responders."],
  ];
  let y = 2.35;
  for (const [n, h, b] of items) {
    s.addShape(pres.shapes.OVAL, { x: 0.75, y: y + 0.05, w: 0.78, h: 0.78, fill: { color: TEAL } });
    s.addText(n, { x: 0.75, y: y + 0.05, w: 0.78, h: 0.78, fontSize: 26, bold: true, color: WHITE, align: "center", valign: "middle", fontFace: HFONT, margin: 0 });
    s.addText([
      { text: h + "  ", options: { bold: true, color: WHITE } },
      { text: b, options: { color: "BFD4D9" } },
    ], { x: 1.8, y: y - 0.02, w: 9.4, h: 1.05, fontSize: 15, fontFace: BFONT, lineSpacingMultiple: 1.05, valign: "middle", margin: 0 });
    y += 1.32;
  }
  s.addText("Method, data, and figures are fully reproducible: python main.py  →  results/report.md",
    { x: 0.7, y: 6.7, w: 12, h: 0.35, fontSize: 11, color: "8FAAB0", italic: true, fontFace: BFONT, margin: 0 });
  s.addNotes("Takeaways. Three honest claims: (1) integration beats TMB alone, mostly via clinical/cancer-type context; (2) the pure genomic increment is small and we don't oversell it; (3) the real conceptual contribution is showing TMB masks resistance genes, which explains the high-TMB non-responders. Everything reproduces from one command.");
}

(async () => {
  await titleSlide();
  await questionSlide();
  await dataSlide();
  await result1Slide();
  await result2Slide();
  await result3Slide();
  await learnedSlide();
  await biblioSlide();
  await limitationsSlide();
  await takeawaySlide();
  await pres.writeFile({ fileName: "C:/Projects/ICIpredict/results/ICIpredict_results.pptx" });
  console.log("DECK WRITTEN");
})();
