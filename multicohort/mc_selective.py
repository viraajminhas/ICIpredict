"""
Reliability-aware / selective cross-cohort ICI prediction.

Reframing: transcriptomic ICI prediction fails on *some* cohorts (LOCO AUC
0.47-0.80). Rather than predict everywhere, estimate a per-patient TRANSFER-
RELIABILITY signal and abstain on the unreliable cases.

Reliability = agreement of a jackknife ensemble: for a held-out cohort, we train
one model per (leave-one-more-training-cohort-out) split. If a patient's
prediction flips depending on which training cohorts were used, that patient's
prediction does not transfer -> abstain.

Metric: risk-coverage. Sort held-out patients by reliability, and at each
coverage compute AUC on the retained set. A useful reliability signal makes the
curve rise above the random-abstention baseline.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from mc_data import PROC
from mc_signatures import SIGNATURES

RNG = 20260704
MIN_CLASS = 12


def loco_reliability(df, sig_cols):
    """Per held-out patient: ensemble mean prediction, disagreement, OOD distance."""
    recs = []
    for held in df["cohort"].unique():
        te = df[df["cohort"] == held]
        if min(te["response"].value_counts()) < MIN_CLASS:
            continue
        train_cohorts = [c for c in df["cohort"].unique() if c != held]
        Xte = te[sig_cols].fillna(0).values

        preds = []
        for skip in train_cohorts:                       # jackknife over train cohorts
            tr = df[(df["cohort"] != held) & (df["cohort"] != skip)]
            m = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=2000, C=1.0))
            m.fit(tr[sig_cols].fillna(0).values, tr["response"].values)
            preds.append(m.predict_proba(Xte)[:, 1])
        P = np.vstack(preds)                              # K x n_held
        mean_p, disagree = P.mean(0), P.std(0)

        # OOD: distance to the nearest training-pool centroid per cohort (in sig space)
        tr_all = df[df["cohort"] != held]
        sc = StandardScaler().fit(tr_all[sig_cols].fillna(0).values)
        Zte = sc.transform(Xte)
        Ztr = sc.transform(tr_all[sig_cols].fillna(0).values)
        ood = cdist(Zte, Ztr).min(axis=1)                # nearest-neighbor distance

        for i, (_, row) in enumerate(te.iterrows()):
            recs.append({"cohort": held, "y": int(row["response"]),
                         "p": mean_p[i], "disagree": disagree[i],
                         "conf": abs(mean_p[i] - 0.5), "ood": ood[i]})
    return pd.DataFrame(recs)


def risk_coverage(rel, by, ascending):
    """AUC on the most-reliable fraction at each coverage level."""
    d = rel.sort_values(by, ascending=ascending).reset_index(drop=True)
    covs, aucs = [], []
    for c in np.arange(0.3, 1.01, 0.1):
        k = int(len(d) * c)
        sub = d.iloc[:k]
        if sub["y"].nunique() == 2:
            covs.append(c); aucs.append(roc_auc_score(sub["y"], sub["p"]))
    return np.array(covs), np.array(aucs)


def main():
    df = pd.read_csv(PROC / "multicohort.csv", index_col=0)
    sig_cols = [c for c in SIGNATURES if c in df.columns]
    rel = loco_reliability(df, sig_cols)
    rel.to_csv(PROC / "selective_reliability.csv", index=False)
    full_auc = roc_auc_score(rel["y"], rel["p"])
    print(f"pooled held-out patients: {len(rel)}   full-coverage AUC: {full_auc:.3f}")

    signals = {"disagreement (ours)": ("disagree", True),
               "confidence |p-.5|": ("conf", False),
               "OOD distance": ("ood", True)}
    curves = {}
    rng = np.random.default_rng(RNG)
    print("\nSelective AUC at 50% coverage (abstain on least-reliable half):")
    for name, (col, asc) in signals.items():
        cov, auc = risk_coverage(rel, col, asc)
        curves[name] = (cov, auc)
        at50 = auc[np.argmin(np.abs(cov - 0.5))]
        print(f"  {name:22} AUC@50%={at50:.3f}   (full={full_auc:.3f}, "
              f"gain={at50-full_auc:+.3f})")
    # random-abstention baseline (averaged)
    rand_curve = []
    for c in np.arange(0.3, 1.01, 0.1):
        vals = []
        for _ in range(50):
            sub = rel.sample(frac=min(c, 1.0), random_state=int(rng.integers(1e9)))
            if sub["y"].nunique() == 2:
                vals.append(roc_auc_score(sub["y"], sub["p"]))
        rand_curve.append(np.mean(vals))
    curves["random abstention"] = (np.arange(0.3, 1.01, 0.1), np.array(rand_curve))
    print(f"  {'random abstention':22} AUC@50%={rand_curve[2]:.3f}")

    pd.DataFrame({"coverage": curves["disagreement (ours)"][0],
                  **{k: pd.Series(v[1], index=range(len(v[1]))) for k, v in curves.items()}
                  }).to_csv(PROC / "risk_coverage.csv", index=False)

    # per-cohort transferability: does ensemble disagreement predict cohort AUC?
    print("\nPer-cohort: mean disagreement vs LOCO AUC (is unpredictability foreseeable?)")
    rows = []
    for c, g in rel.groupby("cohort"):
        rows.append({"cohort": c, "auc": roc_auc_score(g["y"], g["p"]),
                     "mean_disagree": g["disagree"].mean(),
                     "mean_ood": g["ood"].mean()})
    pc = pd.DataFrame(rows).sort_values("auc")
    print(pc.round(3).to_string(index=False))
    from scipy.stats import spearmanr
    r, p = spearmanr(pc["mean_disagree"], pc["auc"])
    print(f"\nSpearman(mean_disagree, cohort AUC) = {r:.2f} (p={p:.3f})")
    pc.to_csv(PROC / "per_cohort_transferability.csv", index=False)
    return rel, curves, pc


if __name__ == "__main__":
    main()
