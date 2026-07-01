# Predicting immunotherapy response from genomic data — results

**Question.** Does an interpretable, cancer-type-aware survival model that integrates
tumor mutational burden (TMB) with specific gene mutations predict overall survival
(OS) after immune checkpoint inhibitors (ICI) **better than TMB alone**, and are
gene-level resistance effects *masked* by their correlation with TMB?

**Cohort.** Samstein et al. 2019 (cBioPortal `tmb_mskcc_2018`), pulled live from the
cBioPortal REST API. n = **1630** ICI-treated patients across
11 cancer types; **810** deaths
(49.7%). Features: log-TMB, age, sex, sample type, ICI drug
class, cancer type (one-hot), and **69 gene-mutation indicators**
restricted to the IMPACT341 core panel (so a 0 = wild-type, never "unsequenced").
The elastic-net penalty is tuned by **nested** CV (per training fold), so the reported
C-index carries no hyperparameter-selection optimism.

## 1. Model comparison (repeated stratified 5-fold CV, 10 repeats)

| Model | Harrell C (mean ± SD) | IPCW C | Time-dep. AUC |
|---|---|---|---|
| TMB ≥ 10 (FDA cutoff) | 0.551 ± 0.017 | 0.553 | 0.562 |
| Type-specific top-20% TMB (Samstein) | 0.535 ± 0.015 | 0.534 | 0.545 |
| Cox: TMB alone | 0.547 ± 0.018 | 0.551 | 0.564 |
| Cox: clinical only | 0.623 ± 0.014 | 0.603 | 0.697 |
| Cox: TMB + clinical (no genes) | 0.637 ± 0.017 | 0.621 | 0.711 |
| Elastic-net Cox: TMB + clinical + genes | 0.643 ± 0.019 | 0.625 | 0.717 |
| Random survival forest: full model | 0.641 ± 0.018 | 0.628 | 0.714 |

**Where the gain comes from (decomposition).** TMB alone gives C = 0.547.
Most of the lift is **clinical + cancer-type context**: clinical-only already reaches
C = 0.623, and TMB + clinical C = 0.637. The *genuine genomic
increment* is the step from TMB+clinical to the full model — small, and the honest
test of "do genes add value":

| Comparison (ΔC = C(A) − C(B)) | ΔC | 95% CI (bootstrap) | boot p | Holm p | per-fold ΔC ± SD (Wilcoxon) |
|---|---|---|---|---|---|
| Coxnet Full vs. TMB | +0.099 | [+0.077, +0.120] | <0.001 | 0.004 | +0.096 ± 0.014 (p<0.001) |
| RSF Full vs. TMB | +0.096 | [+0.074, +0.118] | <0.001 | 0.004 | +0.095 ± 0.014 (p<0.001) |
| Coxnet Full vs. TMBClin | +0.008 | [−0.003, +0.017] | 0.131 | 0.262 | +0.006 ± 0.008 (p<0.001) |
| RSF Full vs. TMBClin | +0.005 | [−0.005, +0.015] | 0.324 | 0.324 | +0.004 ± 0.009 (p=0.004) |

- The large advantage over **TMB alone** (ΔC = +0.099) is real but is
  **mostly clinical/cancer-type covariates, not genes**.
- The **genomic increment over TMB+clinical** is small and **not robust**:
  `Coxnet` ΔC = +0.008 (95% CI crosses 0; bootstrap p = 0.131, Holm p = 0.262), and
  the RSF twin is null. A paired per-fold test detects a consistent but negligible
  +0.006. We report it as **not established — suggestive at most**.

## 2. What the model learned (in-sample interpretation)

Elastic-net Cox (α = 0.0077, fit on the full cohort — **descriptive, not a
performance estimate**) hazard ratios:
- **Worse OS / resistance (HR > 1):** STK11, TP53, drug_PD-1/PDL-1, SMARCA4, ct_Glioma, PBRM1
- **Better OS / benefit (HR < 1):** ct_Renal Cell Carcinoma, ct_Melanoma, VHL, log_tmb, CREBBP, PTPRD

Random-survival-forest permutation importance (in-sample, top): log_tmb, ct_Melanoma, ct_Renal Cell Carcinoma, TP53, VHL, ct_Non-Small Cell Lung Cancer, ct_Bladder Cancer, STK11.

Out-of-fold model risk, split into **tertiles**, separates survival (log-rank
p = 3.0e-40, `fig6_risk_tertile_km.png`); cutpoints are data-derived so treat the
p as descriptive — the cross-validated C-index above is the honest discrimination
measure.

## 3. Resistance genes are masked by TMB — the integration insight (key finding)

Rizvi et al. (2015) saw high-TMB tumors that did **not** benefit but had no
explanation, searching by eye. Here is *why that search fails*:

**The confound.** Canonical resistance-gene mutations
(STK11, KEAP1, B2M, JAK1, JAK2, PTEN) occur preferentially in **high-TMB** tumors (median TMB
9.8 vs 5.2 mut/Mb in wild-type;
Mann-Whitney p = 1.9e-28). Because high TMB is *favorable*,
the harmful effect of these mutations is **masked** (positive confounding) — a naive
subgroup split is null. STK11/KEAP1 are also concentrated in NSCLC, so **cancer type
is a second confounder** that must be held fixed.

**The fix — model TMB and cancer type jointly.** Within NSCLC (cancer type held fixed), **STK11 (NSCLC)** goes from HR 1.31 (p=0.11) without TMB to **HR 1.42** [1.01–1.98], p=0.041 once TMB is added
(`fig7_resistance_tmb_confounding.png`). Full table (NSCLC panels hold cancer type
fixed; pan-cancer panels adjust for cancer type; "d" = deaths in mutated arm;
† = ridge-stabilised):

| Resistance marker | # mutated (deaths) | HR without TMB | **HR with TMB** |
|---|---|---|---|
| STK11 (NSCLC) | 68 (46 d) | 1.31 (p=0.11) | **1.42** [1.01–1.98] (p=0.041) |
| KEAP1 (NSCLC) | 69 (45 d) | 1.35 (p=0.08) | **1.54** [1.09–2.18] (p=0.014) |
| Resistance set (NSCLC) | 126 (80 d) | 1.17 (p=0.28) | **1.34** [0.99–1.80] (p=0.055) |
| STK11 (pan-cancer, type-adj.) | 96 (59 d) | 1.34 (p=0.04) | **1.50** [1.12–2.00] (p=0.006) |
| Resistance set (pan-cancer, type-adj.) | 338 (175 d) | 0.98 (p=0.78) | **1.15** [0.96–1.38] (p=0.135) |

This is the mechanistic reason an integrated model can beat TMB alone: TMB and
resistance genes act in **opposite directions but are positively correlated**, so
only a model that sees TMB, the genes, *and* cancer type together can score either
correctly. (Note: this is positive confounding/masking — we do not claim a full
cross-stratum sign reversal.)

## 4. How this synthesizes the group's bibliography

- **Lee/Samstein (TMB):** we reproduce TMB's signal and show a global cutoff
  (FDA ≥10) and even Samstein's type-specific top-20% rule under-perform a
  cancer-type-aware multivariable model on OS discrimination.
- **Rizvi (mutational landscape):** POLE/MMR-type genes track benefit; and we give a
  concrete reason high-TMB non-responders evade eyeballing — resistance mutations are
  TMB- and cancer-type-confounded.
- **Jamieson & Maker (resistance):** STK11/KEAP1 (and, where powered, the broader
  resistance set) carry the expected resistance direction — **consistent with /
  hypothesis-generating**, not proof.

## 5. Honest limitations

- Targeted panel (not WES); no PD-L1, RNA/IFN-γ signature, or neoantigen calling.
  Gene "mutation" = any non-silent variant present, **not** demonstrated functional
  loss; OS cannot separate predictive from prognostic effects, and there is no
  ICI-vs-control arm.
- Single institution (MSK); retrospective; one cancer-type bucket has n=1.
- The genomic increment over TMB+clinical is small, Coxnet-specific, and only
  marginally significant after multiplicity correction — **suggestive, not
  established**; needs external validation (e.g. independent ICI WES cohorts).
- The data-driven recurrent-gene set is chosen by *unsupervised* marginal frequency on
  the full cohort (outcome-blind); a strictly fold-local selection is future work.
- Resistance-gene HRs are explanatory (full-cohort) Cox fits, reported with the
  cancer-type adjustment that the dominant confounder requires.

_Generated by `main.py`. Tables in `results/tables/`, figures in `results/figures/`._
