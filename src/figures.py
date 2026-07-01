"""Publication-style figures. All write to results/figures/."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config as C

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 150, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold",
})
TEAL, ORANGE, GREY, RED, BLUE = "#0E7C86", "#E07B39", "#888888", "#C0392B", "#2C5F8A"


def _save(fig, name):
    fp = C.FIGURES / name
    fig.tight_layout()
    fig.savefig(fp, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure -> {fp.name}")


def fig_tmb_by_type(meta: pd.DataFrame):
    # drop cancer types with too few samples to form a real box (e.g. n=1)
    counts = meta["CANCER_TYPE"].value_counts()
    keep = counts[counts >= 5].index
    sub = meta[meta["CANCER_TYPE"].isin(keep)]
    order = (sub.groupby("CANCER_TYPE")["tmb"].median().sort_values(ascending=False)
             .index.tolist())
    data = [sub.loc[sub["CANCER_TYPE"] == t, "tmb"].clip(upper=60) for t in order]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    bp = ax.boxplot(data, vert=True, patch_artist=True, showfliers=False,
                    widths=0.6)
    for patch in bp["boxes"]:
        patch.set(facecolor=TEAL, alpha=0.55)
    for med in bp["medians"]:
        med.set(color="black", linewidth=1.4)
    ax.axhline(C.FDA_TMB_CUTOFF, color=RED, ls="--", lw=1.2,
               label=f"FDA cutoff = {C.FDA_TMB_CUTOFF:.0f} mut/Mb")
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([t.replace(" Cancer", "").replace("Carcinoma", "Ca.")
                        for t in order], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("TMB (mut/Mb, clipped at 60)")
    ax.set_title("TMB varies widely across cancer types — one global cutoff misfits most")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, "fig1_tmb_by_cancer_type.png")


def fig_cindex_compare(summaries: list[dict]):
    df = pd.DataFrame(summaries).sort_values("cindex_mean")
    colors = [ORANGE if m in ("Coxnet_Full", "RSF_Full") else
              (BLUE if m == "Cox_TMB" else GREY) for m in df["model"]]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    y = np.arange(len(df))
    ax.barh(y, df["cindex_mean"], xerr=df["cindex_sd"], color=colors,
            alpha=0.85, capsize=3, error_kw=dict(lw=1, alpha=0.6))
    ax.axvline(0.5, color="black", lw=0.8, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=9)
    ax.set_xlabel("Cross-validated Harrell's C-index (mean ± SD over folds)")
    ax.set_xlim(0.45, max(0.7, df["cindex_mean"].max() + 0.05))
    ax.set_title("Cancer-type + clinical context drives most of the gain over TMB;\n"
                 "genes add a small further increment", fontsize=10.5)
    for yi, (v, s) in enumerate(zip(df["cindex_mean"], df["cindex_sd"])):
        ax.text(v + s + 0.003, yi, f"{v:.3f}", va="center", fontsize=8)
    _save(fig, "fig2_cindex_comparison.png")


def fig_auc_time(auc_curves: dict, labels: dict):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    palette = {"Cox_TMB": BLUE, "Cox_TMB_Clin": GREY,
               "Coxnet_Full": TEAL, "RSF_Full": ORANGE}
    for name, d in auc_curves.items():
        m = ~np.isnan(d["auc_mean"])
        ax.plot(d["times"][m], d["auc_mean"][m], "-o", color=palette.get(name, GREY),
                label=labels.get(name, name), lw=2, ms=5)
    ax.axhline(0.5, color="black", lw=0.8, ls=":")
    ax.set_xlabel("Months since start of immunotherapy")
    ax.set_ylabel("Time-dependent AUC")
    ax.set_title("Discrimination of survival over time")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, "fig3_auc_over_time.png")


def fig_coxnet_hr(hr_df: pd.DataFrame, top: int = 18):
    d = hr_df.head(top).iloc[::-1]
    colors = [RED if c > 0 else TEAL for c in d["coef"]]
    fig, ax = plt.subplots(figsize=(7.5, 6))
    y = np.arange(len(d))
    ax.barh(y, d["hazard_ratio"] - 1.0, left=1.0, color=colors, alpha=0.85)
    ax.axvline(1.0, color="black", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(d["feature"], fontsize=8)
    ax.set_xlabel("Hazard ratio (>1 = worse OS / resistance,  <1 = better OS / benefit)")
    ax.set_title("Elastic-net Cox: gene-level drivers of ICI outcome")
    _save(fig, "fig4_coxnet_hazard_ratios.png")


def fig_rsf_importance(imp_df: pd.DataFrame, top: int = 20):
    d = imp_df.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 6))
    y = np.arange(len(d))
    ax.barh(y, d["importance"], xerr=d["importance_sd"], color=TEAL, alpha=0.85,
            capsize=2, error_kw=dict(lw=0.8, alpha=0.5))
    ax.set_yticks(y)
    ax.set_yticklabels(d["feature"], fontsize=8)
    ax.set_xlabel("Permutation importance (drop in C-index)")
    ax.set_title("Random survival forest: most predictive features")
    _save(fig, "fig5_rsf_importance.png")


def fig_risk_strata_km(km: dict, p_value: float):
    fig, ax = plt.subplots(figsize=(7, 4.8))
    colors = {"Low risk": TEAL, "Intermediate": ORANGE, "High risk": RED}
    for g, kmf in km.items():
        kmf.plot_survival_function(ax=ax, color=colors[g], ci_show=True, lw=2)
    ax.set_xlabel("Months")
    ax.set_ylabel("Overall survival probability")
    ax.set_title(f"Model-predicted risk tertiles stratify survival "
                 f"(log-rank p = {p_value:.1e})")
    ax.legend(frameon=False, fontsize=9)
    _save(fig, "fig6_risk_tertile_km.png")


def fig_resistance_confounding(res: dict):
    """Two panels: (A) resistance mutations co-occur with high TMB (the confound);
    (B) adjusting for TMB unmasks them as resistance markers."""
    df = res["df"]
    tab = res["hr_table"].iloc[::-1].reset_index(drop=True)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5),
                                   gridspec_kw={"width_ratios": [1, 1.4]})

    # Panel A: TMB by resistance status
    mut = df.loc[df["res"] == 1, "tmb"].clip(upper=60)
    wt = df.loc[df["res"] == 0, "tmb"].clip(upper=60)
    bp = axA.boxplot([wt, mut], patch_artist=True, showfliers=False, widths=0.55)
    for patch, col in zip(bp["boxes"], [TEAL, RED]):
        patch.set(facecolor=col, alpha=0.55)
    for med in bp["medians"]:
        med.set(color="black", lw=1.4)
    axA.set_xticks([1, 2]); axA.set_xticklabels(["Resistance\nwild-type",
                                                 "Resistance\nmutated"])
    axA.set_ylabel("TMB (mut/Mb, clipped at 60)")
    axA.set_title(f"A. Resistance mutations co-occur with HIGH TMB\n"
                  f"(median {res['tmb_median_wt']:.1f} vs "
                  f"{res['tmb_median_mut']:.1f}, Mann-Whitney p = "
                  f"{res['tmb_mannwhitney_p']:.1e})", fontsize=9)

    # Panel B: forest plot, HR without vs with TMB in the model
    y = np.arange(len(tab))
    axB.errorbar(tab["unadj_hr"], y - 0.13, fmt="o", color=GREY, ms=6,
                 label="Without TMB in model (effect masked)")
    axB.errorbar(tab["adj_hr"], y + 0.13,
                 xerr=[tab["adj_hr"] - tab["adj_lo"], tab["adj_hi"] - tab["adj_hr"]],
                 fmt="s", color=RED, ms=6, capsize=3,
                 label="With TMB in model (effect unmasked)")
    axB.axvline(1.0, color="black", lw=1)
    axB.set_yticks(y)
    axB.set_yticklabels([f"{m}\n(n={n}, mut={k})" for m, n, k in
                         zip(tab["marker"], tab["n"], tab["n_mut"])], fontsize=7.5)
    for yi, p in zip(y, tab["adj_p"]):
        star = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
        axB.text(tab["adj_hi"].iloc[yi] + 0.05, yi + 0.13, star,
                 va="center", fontsize=9, color=RED)
    axB.set_xlabel("Hazard ratio (>1 = resistance / worse OS)")
    axB.set_title("B. Adding TMB (holding cancer type fixed) unmasks\n"
                  "resistance-gene effects that confounding hides", fontsize=9)
    axB.legend(frameon=False, fontsize=8, loc="lower right")
    _save(fig, "fig7_resistance_tmb_confounding.png")
