"""
The motivating benchmark: do ICI-response biomarkers / models that look good
WITHIN a cohort survive LEAVE-ONE-COHORT-OUT (LOCO) transfer?

Reports, for each method:
  * within-cohort 5-fold CV AUC (the optimistic number usually published)
  * LOCO AUC  (train on all-but-one cohort, test on the held-out cohort)
and per-signature per-cohort AUC to expose biomarker instability.
"""
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from mc_data import PROC
from mc_signatures import SIGNATURES

warnings.filterwarnings("ignore")
RNG = 20260704
MIN_CLASS = 12          # a held-out cohort needs >= this many of each class to score


def load():
    df = pd.read_csv(PROC / "multicohort.csv", index_col=0)
    sig_cols = [c for c in SIGNATURES if c in df.columns]
    gene_cols = [c for c in df.columns if c.startswith("g_")]
    return df, sig_cols, gene_cols


def _auc(y, s):
    return roc_auc_score(y, s) if len(np.unique(y)) == 2 else np.nan


def signature_stability(df, sig_cols):
    """Per-cohort AUC of each single signature (univariate)."""
    rows = []
    for sig in sig_cols:
        aucs = {}
        for c, g in df.groupby("cohort"):
            if g["response"].nunique() == 2 and min(g["response"].value_counts()) >= 5:
                aucs[c] = _auc(g["response"], g[sig])
        vals = np.array(list(aucs.values()))
        rows.append({"signature": sig, "mean_auc": np.nanmean(vals),
                     "min_auc": np.nanmin(vals), "max_auc": np.nanmax(vals),
                     "n_cohorts": len(vals),
                     "frac_better_than_chance": np.mean(vals > 0.5)})
    return pd.DataFrame(rows).sort_values("mean_auc", ascending=False)


def within_cohort_cv(df, cols, make_model):
    """5-fold CV AUC inside each cohort, averaged over cohorts."""
    aucs = []
    for c, g in df.groupby("cohort"):
        y = g["response"].values
        if y.sum() < 6 or (len(y) - y.sum()) < 6:
            continue
        X = g[cols].fillna(0).values
        oof = np.zeros(len(y))
        skf = StratifiedKFold(5, shuffle=True, random_state=RNG)
        for tr, te in skf.split(X, y):
            m = make_model().fit(X[tr], y[tr])
            oof[te] = m.predict_proba(X[te])[:, 1]
        aucs.append(_auc(y, oof))
    return float(np.nanmean(aucs)), aucs


def loco(df, cols, make_model):
    """Train on all-but-one cohort, test on the held-out cohort."""
    per = {}
    for held in df["cohort"].unique():
        tr = df[df["cohort"] != held]
        te = df[df["cohort"] == held]
        if min(te["response"].value_counts()) < MIN_CLASS:
            continue
        m = make_model().fit(tr[cols].fillna(0).values, tr["response"].values)
        p = m.predict_proba(te[cols].fillna(0).values)[:, 1]
        per[held] = _auc(te["response"].values, p)
    return float(np.nanmean(list(per.values()))), per


def single_signature_loco(df, sig):
    """No training: apply the fixed signature score to each held-out cohort."""
    per = {}
    for held in df["cohort"].unique():
        te = df[df["cohort"] == held]
        if min(te["response"].value_counts()) < MIN_CLASS:
            continue
        per[held] = _auc(te["response"].values, te[sig].values)
    return float(np.nanmean(list(per.values()))), per


def main():
    df, sig_cols, gene_cols = load()
    print(f"cohorts={df['cohort'].nunique()}  n={len(df)}  "
          f"responders={int(df['response'].sum())} ({df['response'].mean()*100:.0f}%)")

    print("\n=== (1) Single-signature stability across cohorts ===")
    stab = signature_stability(df, sig_cols)
    print(stab.round(3).to_string(index=False))
    stab.to_csv(PROC / "signature_stability.csv", index=False)

    logit = lambda: make_pipeline(StandardScaler(),
                                  LogisticRegression(max_iter=2000, C=1.0))
    logit_l2 = lambda: make_pipeline(StandardScaler(),
                                     LogisticRegression(max_iter=2000, C=0.1))
    gbm = lambda: HistGradientBoostingClassifier(
        max_depth=3, learning_rate=0.05, max_iter=300, l2_regularization=1.0,
        random_state=RNG)

    methods = {
        "Best single sig (IFNg_Ayers18)": ("single", "IFNg_Ayers18"),
        "Logistic (11 signatures)": ("model", sig_cols, logit),
        "Logistic-L2 (126 genes)": ("model", gene_cols, logit_l2),
        "GBM (126 genes)": ("model", gene_cols, gbm),
    }

    print("\n=== (2) Within-cohort CV AUC  vs  LOCO AUC (the generalization gap) ===")
    rows = []
    for name, spec in methods.items():
        if spec[0] == "single":
            sig = spec[1]
            wc, _ = within_cohort_cv(df, [sig], lambda: make_pipeline(
                StandardScaler(), LogisticRegression(max_iter=1000)))
            lo, per = single_signature_loco(df, sig)
        else:
            _, cols, mk = spec
            wc, _ = within_cohort_cv(df, cols, mk)
            lo, per = loco(df, cols, mk)
        rows.append({"method": name, "within_cohort_AUC": wc, "LOCO_AUC": lo,
                     "gap": wc - lo})
        print(f"  {name:34}  within={wc:.3f}   LOCO={lo:.3f}   gap={wc-lo:+.3f}")
    res = pd.DataFrame(rows)
    res.to_csv(PROC / "loco_benchmark.csv", index=False)

    # per-held-out-cohort LOCO for the signature-logistic model (for the figure)
    _, per_cohort = loco(df, sig_cols, logit)
    pd.Series(per_cohort, name="LOCO_AUC").to_csv(PROC / "loco_per_cohort.csv")
    print("\nPer-held-out-cohort LOCO AUC (Logistic on signatures):")
    for c, a in sorted(per_cohort.items(), key=lambda x: -x[1]):
        print(f"  {c:34} {a:.3f}")
    return res, stab


if __name__ == "__main__":
    main()
