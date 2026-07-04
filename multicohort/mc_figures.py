"""Figures for the multi-cohort LOCO benchmark."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from mc_data import PROC, COHORTS
from mc_signatures import SIGNATURES

ROOT = PROC.parent
FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
TEAL, ORANGE, RED, GREY = "#0E7C86", "#E07B39", "#C0392B", "#888888"
plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 150, "font.size": 9,
                     "axes.titleweight": "bold"})

SHORT = {c: c.split("_iatlas_")[-1].split("_")[0].title() if "_iatlas_" in c else c
         for c in COHORTS}


def _auc(y, s):
    return roc_auc_score(y, s) if len(np.unique(y)) == 2 else np.nan


def fig_stability(df, sig_cols):
    cohorts = [c for c in df["cohort"].unique()]
    M = np.full((len(sig_cols), len(cohorts)), np.nan)
    for i, sig in enumerate(sig_cols):
        for j, c in enumerate(cohorts):
            g = df[df["cohort"] == c]
            if g["response"].nunique() == 2 and min(g["response"].value_counts()) >= 5:
                M[i, j] = _auc(g["response"], g[sig])
    order = np.argsort(-np.nanmean(M, axis=1))
    M, sig_cols = M[order], [sig_cols[i] for i in order]

    fig, ax = plt.subplots(figsize=(8.5, 5))
    im = ax.imshow(M, cmap="RdBu_r", vmin=0.30, vmax=0.70, aspect="auto")
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels([f"{SHORT.get(c,c)}\n(n={len(df[df.cohort==c])})" for c in cohorts],
                       fontsize=8, rotation=30, ha="right")
    ax.set_yticks(range(len(sig_cols))); ax.set_yticklabels(sig_cols, fontsize=8)
    for i in range(len(sig_cols)):
        for j in range(len(cohorts)):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(M[i, j]-0.5) > 0.13 else "black")
    ax.set_title("Same biomarker, wildly different AUC across cohorts\n"
                 "(single-signature response AUC per cohort)")
    cb = fig.colorbar(im, ax=ax, shrink=0.8); cb.set_label("AUC", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "mc_fig1_signature_instability.png"); plt.close(fig)
    print("  -> mc_fig1_signature_instability.png")


def fig_gap(bench):
    d = bench.copy()
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    y = np.arange(len(d)); h = 0.38
    ax.barh(y + h/2, d["within_cohort_AUC"], h, color=GREY, alpha=0.9,
            label="Within-cohort CV (optimistic)")
    ax.barh(y - h/2, d["LOCO_AUC"], h, color=TEAL, alpha=0.95,
            label="Leave-one-cohort-out (honest)")
    ax.axvline(0.5, color="black", lw=0.8, ls=":")
    ax.set_yticks(y); ax.set_yticklabels(d["method"], fontsize=8)
    ax.set_xlim(0.45, 0.72); ax.set_xlabel("AUC")
    ax.set_title("The generalization gap: within-cohort performance does not transfer")
    for yi, (w, l) in enumerate(zip(d["within_cohort_AUC"], d["LOCO_AUC"])):
        ax.text(w+0.004, yi+h/2, f"{w:.3f}", va="center", fontsize=7)
        ax.text(l+0.004, yi-h/2, f"{l:.3f}", va="center", fontsize=7, color=TEAL)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "mc_fig2_generalization_gap.png"); plt.close(fig)
    print("  -> mc_fig2_generalization_gap.png")


def fig_per_cohort(per):
    d = per.sort_values("LOCO_AUC")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    colors = [RED if v < 0.5 else (ORANGE if v < 0.6 else TEAL) for v in d["LOCO_AUC"]]
    ax.barh(range(len(d)), d["LOCO_AUC"], color=colors, alpha=0.9)
    ax.axvline(0.5, color="black", lw=0.9, ls=":")
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([SHORT.get(c, c) for c in d.index], fontsize=8)
    ax.set_xlabel("LOCO AUC (held-out cohort)"); ax.set_xlim(0.4, 0.85)
    ax.set_title("Held-out cohort AUC varies from below-chance to strong\n"
                 "(pooled model on signatures, trained on the other cohorts)")
    for i, v in enumerate(d["LOCO_AUC"]):
        ax.text(v+0.005, i, f"{v:.2f}", va="center", fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "mc_fig3_per_cohort_loco.png"); plt.close(fig)
    print("  -> mc_fig3_per_cohort_loco.png")


if __name__ == "__main__":
    df = pd.read_csv(PROC / "multicohort.csv", index_col=0)
    sig_cols = [c for c in SIGNATURES if c in df.columns]
    fig_stability(df, sig_cols)
    fig_gap(pd.read_csv(PROC / "loco_benchmark.csv"))
    per = pd.read_csv(PROC / "loco_per_cohort.csv", index_col=0)
    per.columns = ["LOCO_AUC"]
    fig_per_cohort(per)
