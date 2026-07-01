"""
End-to-end pipeline:
  download -> features -> model selection -> cross-validated comparison ->
  bootstrap significance -> interpretation -> high-TMB analysis -> figures + report.

Run:
    python main.py            # full run (~a few minutes)
    python main.py --quick    # fewer CV repeats / bootstraps for a fast smoke test
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings

import numpy as np
import pandas as pd

try:                                    # make Windows consoles UTF-8 safe
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore", category=RuntimeWarning)

from src import config as C
from src import analysis, figures
from src import survival_models as sm
from src.build_features import build


def main(quick: bool = False, force_download: bool = False):
    t0 = time.time()
    if quick:
        C.N_REPEATS, C.N_BOOTSTRAP = 2, 300

    print("=" * 70, "\n[1/7] Building dataset", "\n" + "=" * 70)
    ds = build(force_download=force_download)
    print(f"  n={ds.n}  events={int(ds.event.sum())} ({ds.event.mean()*100:.1f}%)  "
          f"features={ds.X.shape[1]}  genes={len(ds.gene_cols)}")

    print("\n[2/7] Selecting elastic-net Cox regularization (inner CV)")
    alpha = sm.select_coxnet_alpha(ds)

    print("\n[3/7] Cross-validated model comparison "
          f"({C.N_SPLITS}x{C.N_REPEATS} folds; Coxnet penalty tuned per fold)")
    specs = sm.build_specs(ds)
    results = sm.cross_validate(ds, specs)
    summaries = [results[s.name].summary() for s in specs]
    summ_df = pd.DataFrame(summaries)
    summ_df.to_csv(C.TABLES / "model_comparison.csv", index=False)
    print(summ_df[["model", "cindex_mean", "cindex_sd", "ipcw_mean",
                   "auc_mean"]].round(4).to_string(index=False))

    print("\n[4/7] Model comparison: paired bootstrap (Holm-corrected) + fold test")
    tmb_oof = results["Cox_TMB"].oof_risk()
    tmbclin_oof = results["Cox_TMB_Clin"].oof_risk()
    best_name = max(["Coxnet_Full", "RSF_Full"],
                    key=lambda m: results[m].summary()["cindex_mean"])
    comparisons, fold_tests = {}, {}
    pairs = {"Cox_TMB": tmb_oof, "Cox_TMB_Clin": tmbclin_oof}
    for integrated in ["Coxnet_Full", "RSF_Full"]:
        oof = results[integrated].oof_risk()
        for ref_name, ref_oof in pairs.items():
            tag = f"{integrated}_vs_{'TMB' if ref_name=='Cox_TMB' else 'TMBClin'}"
            comparisons[tag] = sm.bootstrap_cindex_diff(ds, oof, ref_oof)
            fold_tests[tag] = sm.paired_fold_test(results[integrated], results[ref_name])
    # Holm-Bonferroni across all four comparisons (addresses multiplicity)
    holm = sm.holm_correct({k: v["p_value"] for k, v in comparisons.items()})
    for k in comparisons:
        comparisons[k]["p_holm"] = holm[k]
    for k, v in comparisons.items():
        ft = fold_tests[k]
        print(f"  {k}: dC={v['delta']:+.4f} 95%CI[{v['ci_low']:+.4f},{v['ci_high']:+.4f}] "
              f"p={v['p_value']:.4f} p_holm={v['p_holm']:.4f} | "
              f"fold dC={ft['delta_mean']:+.4f}+-{ft['delta_sd']:.4f} "
              f"wilcoxon_p={ft['wilcoxon_p']:.4f}")
    with open(C.TABLES / "bootstrap_comparisons.json", "w") as fh:
        json.dump({"bootstrap": comparisons, "fold_paired": fold_tests}, fh, indent=2)

    print("\n[5/7] Interpretation (Coxnet hazard ratios + RSF importance)")
    hr_df = analysis.interpret_coxnet(ds, alpha)
    hr_df.to_csv(C.TABLES / "coxnet_hazard_ratios.csv", index=False)
    imp_df = analysis.rsf_importance(ds, n_repeats=5 if quick else 10)
    imp_df.to_csv(C.TABLES / "rsf_importance.csv", index=False)
    print("  top resistance (HR>1):",
          ", ".join(hr_df[hr_df.coef > 0].head(5)["feature"]))
    print("  top benefit (HR<1):   ",
          ", ".join(hr_df[hr_df.coef < 0].head(5)["feature"]))

    print("\n[6/7] Resistance-gene / TMB confounding analysis")
    conf = analysis.resistance_tmb_confounding(ds)
    conf["hr_table"].to_csv(C.TABLES / "resistance_hr_table.csv", index=False)
    print(f"  resistance-mutated TMB median {conf['tmb_median_mut']:.1f} vs "
          f"WT {conf['tmb_median_wt']:.1f} (Mann-Whitney p={conf['tmb_mannwhitney_p']:.1e})")
    for _, r in conf["hr_table"].iterrows():
        print(f"  {r['marker']:22s} unadj HR={r['unadj_hr']:.2f} (p={r['unadj_p']:.3f})"
              f"  ->  TMB-adj HR={r['adj_hr']:.2f} (p={r['adj_p']:.4f})")

    print("\n[7/7] Figures + report")
    km, _, km_p = analysis.risk_strata_km(ds, results[best_name].oof_risk())
    auc_specs = [s for s in specs
                 if s.name in ("Cox_TMB", "Cox_TMB_Clin", "Coxnet_Full", "RSF_Full")]
    auc_curves = sm.auc_over_time(ds, auc_specs, C.EVAL_TIMES,
                                  n_repeats=2 if quick else 5)
    labels = {s.name: s.label for s in specs}

    figures.fig_tmb_by_type(ds.meta)
    figures.fig_cindex_compare(summaries)
    figures.fig_auc_time(auc_curves, labels)
    figures.fig_coxnet_hr(hr_df)
    figures.fig_rsf_importance(imp_df)
    figures.fig_risk_strata_km(km, km_p)
    figures.fig_resistance_confounding(conf)

    write_report(ds, summ_df, comparisons, fold_tests, hr_df, imp_df, conf, alpha,
                 best_name, km_p)
    print(f"\nDone in {time.time()-t0:.0f}s. See results/report.md")


def _fmt_p(v):
    """Bootstrap p: report as '<1/(B+1)' when it hits the resolution floor."""
    floor = 1.0 / (v.get("n_boot", C.N_BOOTSTRAP) + 1)
    return f"<{floor:.1g}" if v["p_value"] <= 2 * floor + 1e-12 else f"{v['p_value']:.3f}"


def write_report(ds, summ_df, comparisons, fold_tests, hr_df, imp_df, conf, alpha,
                 best_name, km_p):
    s = summ_df.set_index("model")
    c_tmb = s.loc["Cox_TMB", "cindex_mean"]
    c_clin = s.loc["Cox_Clin", "cindex_mean"]
    c_tmbclin = s.loc["Cox_TMB_Clin", "cindex_mean"]
    best = comparisons[f"{best_name}_vs_TMB"]
    best_clin = comparisons[f"{best_name}_vs_TMBClin"]

    def row(m):
        r = s.loc[m]
        return (f"| {r['label']} | {r['cindex_mean']:.3f} ± {r['cindex_sd']:.3f} "
                f"| {r['ipcw_mean']:.3f} | {r['auc_mean']:.3f} |")

    def comp_row(tag):
        v, ft = comparisons[tag], fold_tests[tag]
        return (f"| {tag.replace('_', ' ').replace('vs', 'vs.')} "
                f"| {v['delta']:+.3f} | [{v['ci_low']:+.3f}, {v['ci_high']:+.3f}] "
                f"| {_fmt_p(v)} | {v['p_holm']:.3f} "
                f"| {ft['delta_mean']:+.3f} ± {ft['delta_sd']:.3f} (p={ft['wilcoxon_p']:.3f}) |")

    res_hr = hr_df[hr_df.coef > 0].head(6)["feature"].tolist()
    ben_hr = hr_df[hr_df.coef < 0].head(6)["feature"].tolist()
    top_imp = imp_df.head(8)["feature"].tolist()

    ht = conf["hr_table"]

    def hr_line(r):
        dag = " †" if r.get("penalized") else ""
        return (f"| {r['marker']}{dag} | {r['n_mut']:.0f} ({r['n_mut_events']:.0f} d) "
                f"| {r['unadj_hr']:.2f} (p={r['unadj_p']:.2f}) "
                f"| **{r['adj_hr']:.2f}** [{r['adj_lo']:.2f}–{r['adj_hi']:.2f}] "
                f"(p={r['adj_p']:.3f}) |")

    # pick a clean within-NSCLC headline marker if present
    nsclc_rows = ht[ht["marker"].str.contains("NSCLC")]
    headline = (nsclc_rows.iloc[0] if len(nsclc_rows) else
                (ht.iloc[0] if len(ht) else None))
    if headline is not None:
        head_txt = (f"Within NSCLC (cancer type held fixed), **{headline['marker']}** "
                    f"goes from HR {headline['unadj_hr']:.2f} (p={headline['unadj_p']:.2f}) "
                    f"without TMB to **HR {headline['adj_hr']:.2f}** "
                    f"[{headline['adj_lo']:.2f}–{headline['adj_hi']:.2f}], "
                    f"p={headline['adj_p']:.3f} once TMB is added")
    else:
        head_txt = "No resistance marker met the minimum-events threshold"

    md = f"""# Predicting immunotherapy response from genomic data — results

**Question.** Does an interpretable, cancer-type-aware survival model that integrates
tumor mutational burden (TMB) with specific gene mutations predict overall survival
(OS) after immune checkpoint inhibitors (ICI) **better than TMB alone**, and are
gene-level resistance effects *masked* by their correlation with TMB?

**Cohort.** Samstein et al. 2019 (cBioPortal `tmb_mskcc_2018`), pulled live from the
cBioPortal REST API. n = **{ds.n}** ICI-treated patients across
{ds.meta['CANCER_TYPE'].nunique()} cancer types; **{int(ds.event.sum())}** deaths
({ds.event.mean()*100:.1f}%). Features: log-TMB, age, sex, sample type, ICI drug
class, cancer type (one-hot), and **{len(ds.gene_cols)} gene-mutation indicators**
restricted to the IMPACT341 core panel (so a 0 = wild-type, never "unsequenced").
The elastic-net penalty is tuned by **nested** CV (per training fold), so the reported
C-index carries no hyperparameter-selection optimism.

## 1. Model comparison (repeated stratified {C.N_SPLITS}-fold CV, {C.N_REPEATS} repeats)

| Model | Harrell C (mean ± SD) | IPCW C | Time-dep. AUC |
|---|---|---|---|
{chr(10).join(row(m) for m in ['FDA_TMB10','TypeSpecTMB','Cox_TMB','Cox_Clin','Cox_TMB_Clin','Coxnet_Full','RSF_Full'])}

**Where the gain comes from (decomposition).** TMB alone gives C = {c_tmb:.3f}.
Most of the lift is **clinical + cancer-type context**: clinical-only already reaches
C = {c_clin:.3f}, and TMB + clinical C = {c_tmbclin:.3f}. The *genuine genomic
increment* is the step from TMB+clinical to the full model — small, and the honest
test of "do genes add value":

| Comparison (ΔC = C(A) − C(B)) | ΔC | 95% CI (bootstrap) | boot p | Holm p | per-fold ΔC ± SD (Wilcoxon) |
|---|---|---|---|---|---|
{chr(10).join(comp_row(t) for t in ['Coxnet_Full_vs_TMB','RSF_Full_vs_TMB','Coxnet_Full_vs_TMBClin','RSF_Full_vs_TMBClin'])}

- The large advantage over **TMB alone** (ΔC = {best['delta']:+.3f}) is real but is
  **mostly clinical/cancer-type covariates, not genes**.
- The **genomic increment over TMB+clinical** is small and **model-dependent**:
  `Coxnet` ΔC = {best_clin['delta']:+.3f} (Holm p = {best_clin['p_holm']:.3f}), while
  the RSF twin is essentially null. After multiplicity correction it is marginal —
  we report it as *suggestive, not established*.

## 2. What the model learned (in-sample interpretation)

Elastic-net Cox (α = {alpha:.4f}, fit on the full cohort — **descriptive, not a
performance estimate**) hazard ratios:
- **Worse OS / resistance (HR > 1):** {', '.join(res_hr) or '—'}
- **Better OS / benefit (HR < 1):** {', '.join(ben_hr) or '—'}

Random-survival-forest permutation importance (in-sample, top): {', '.join(top_imp)}.

Out-of-fold model risk, split into **tertiles**, separates survival (log-rank
p = {km_p:.1e}, `fig6_risk_tertile_km.png`); cutpoints are data-derived so treat the
p as descriptive — the cross-validated C-index above is the honest discrimination
measure.

## 3. Resistance genes are masked by TMB — the integration insight (key finding)

Rizvi et al. (2015) saw high-TMB tumors that did **not** benefit but had no
explanation, searching by eye. Here is *why that search fails*:

**The confound.** Canonical resistance-gene mutations
({', '.join(conf['genes'])}) occur preferentially in **high-TMB** tumors (median TMB
{conf['tmb_median_mut']:.1f} vs {conf['tmb_median_wt']:.1f} mut/Mb in wild-type;
Mann-Whitney p = {conf['tmb_mannwhitney_p']:.1e}). Because high TMB is *favorable*,
the harmful effect of these mutations is **masked** (positive confounding) — a naive
subgroup split is null. STK11/KEAP1 are also concentrated in NSCLC, so **cancer type
is a second confounder** that must be held fixed.

**The fix — model TMB and cancer type jointly.** {head_txt}
(`fig7_resistance_tmb_confounding.png`). Full table (NSCLC panels hold cancer type
fixed; pan-cancer panels adjust for cancer type; "d" = deaths in mutated arm;
† = ridge-stabilised):

| Resistance marker | # mutated (deaths) | HR without TMB | **HR with TMB** |
|---|---|---|---|
{chr(10).join(hr_line(r) for _, r in ht.iterrows()) if len(ht) else '| (none met the ≥10-events threshold) | | | |'}

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
"""
    (C.RESULTS / "report.md").write_text(md, encoding="utf-8")
    print(f"  report -> {C.RESULTS / 'report.md'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--force-download", action="store_true")
    args = ap.parse_args()
    main(quick=args.quick, force_download=args.force_download)
