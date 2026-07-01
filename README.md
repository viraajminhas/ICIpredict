# ICIpredict — predicting immunotherapy response from genomic data

> **Research question.** Does an *interpretable, cancer-type-aware survival model
> that integrates tumor mutational burden (TMB) with specific gene-level
> mutations* predict overall-survival benefit from immune checkpoint inhibitors
> (ICI) **better than TMB alone** — and can gene context explain the "high-TMB
> non-responders" that Rizvi et al. (2015) could not?

This project synthesizes the group's annotated bibliography into one runnable,
falsifiable experiment on **real clinical data**.

📊 **[Download the results slide deck (.pptx)](https://github.com/viraajminhas/ICIpredict/releases/latest)**
· or view it in the repo at [`results/ICIpredict_results.pptx`](results/ICIpredict_results.pptx)
· written write-up: [`results/report.md`](results/report.md)

## Data

[Samstein et al., *Nature Genetics* 2019](https://www.nature.com/articles/s41588-018-0312-8)
— "Tumor mutational load predicts survival after immunotherapy across multiple
cancer types." cBioPortal study `tmb_mskcc_2018`: **1,661 ICI-treated patients**,
11 cancer types, MSK-IMPACT targeted panel, overall-survival outcomes. Pulled
live (and cached) from the **cBioPortal REST API** — no manual download needed.

After QC (positive follow-up, valid TMB): **n = 1,630**, 810 deaths (49.7%).

## What it does

| Stage | File | Output |
|---|---|---|
| Download + cache clinical & mutation data | `src/cbioportal.py` | `data/raw/*.csv` |
| Build survival labels + feature matrix | `src/build_features.py` | `data/processed/*` |
| Models, cross-validated comparison, bootstrap test | `src/survival_models.py` | — |
| Coxnet HRs, RSF importance, **high-TMB analysis** | `src/analysis.py` | — |
| Figures | `src/figures.py` | `results/figures/*.png` |
| Orchestration + auto-report | `main.py` | `results/report.md`, `results/tables/*` |

### Models compared (endpoint = overall survival, scored as risk)

- **FDA_TMB10** — single tissue-agnostic cutoff TMB ≥ 10 mut/Mb
- **TypeSpecTMB** — Samstein's top-20% TMB *within each cancer type* (CV-correct)
- **Cox_TMB** — Cox on log(TMB) only *(the "TMB alone" baseline)*
- **Cox_Clin / Cox_TMB_Clin** — clinical only / clinical + TMB (no genes)
- **Coxnet_Full** — elastic-net Cox on TMB + clinical + 69 genes *(interpretable)*
- **RSF_Full** — random survival forest on the full matrix *(nonlinear, interactions)*

### Methodological care

- Gene indicators restricted to the **IMPACT341 core** (genes on *all three* panel
  versions) so a `0` means wild-type, never "unsequenced" — a real confounder fix.
- Repeated stratified 5-fold CV; out-of-fold risk averaged over repeats.
- Model comparison by **paired patient bootstrap** of the C-index difference.
- Metrics: Harrell's C, IPCW-C, and time-dependent AUC.

## Run it

```bash
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python main.py            # full run (~few minutes)
.venv\Scripts\python main.py --quick    # fast smoke test
```

Read `results/report.md` for the auto-generated write-up, figures, and the
bibliography synthesis. See that file's "Honest limitations" section — effect
sizes are modest, consistent with the literature that no single genomic axis
fully predicts ICI response.
