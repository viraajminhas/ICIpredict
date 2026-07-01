"""
Model interpretation + the headline biological analysis.

1. interpret_coxnet : elastic-net Cox hazard ratios (which genes raise/lower risk)
2. rsf_importance   : permutation importance (drop in C-index) for the forest
3. risk_strata_km   : Kaplan-Meier by model-predicted risk tertiles
4. high_tmb_resolution : *the novel bit* -- within TMB-high tumors, do canonical
   resistance-gene mutations (STK11/KEAP1/B2M/JAK1/JAK2/PTEN) identify the
   non-responders that TMB alone cannot? (Directly addresses Rizvi et al.'s
   unexplained high-TMB non-benefiters.)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy.stats import mannwhitneyu
from sklearn.inspection import permutation_importance
from sksurv.ensemble import RandomSurvivalForest
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

from . import config as C
from .build_features import Dataset
from .survival_models import CONTINUOUS, to_surv


def _scale_continuous(X: pd.DataFrame) -> pd.DataFrame:
    Xs = X.copy()
    for c in CONTINUOUS:
        if c in Xs.columns:
            mu, sd = Xs[c].mean(), Xs[c].std() or 1.0
            Xs[c] = (Xs[c] - mu) / sd
    return Xs


def interpret_coxnet(ds: Dataset, alpha: float) -> pd.DataFrame:
    """Hazard ratios from the elastic-net Cox fit on the full cohort."""
    Xs = _scale_continuous(ds.X)
    y = to_surv(ds.event, ds.time)
    model = CoxnetSurvivalAnalysis(l1_ratio=0.5, alphas=[alpha],
                                   fit_baseline_model=False, max_iter=100_000,
                                   normalize=False)
    model.fit(Xs, y)
    coef = np.asarray(model.coef_).ravel()
    out = pd.DataFrame({"feature": Xs.columns, "coef": coef})
    out["hazard_ratio"] = np.exp(out["coef"])
    out["abs_coef"] = out["coef"].abs()
    out = out[out["coef"] != 0].sort_values("abs_coef", ascending=False)
    out["direction"] = np.where(out["coef"] > 0, "worse OS (resistance)",
                                "better OS (benefit)")
    return out.reset_index(drop=True)


def rsf_importance(ds: Dataset, n_repeats: int = 10) -> pd.DataFrame:
    Xs = _scale_continuous(ds.X)
    y = to_surv(ds.event, ds.time)
    rsf = RandomSurvivalForest(n_estimators=300, min_samples_leaf=15,
                               max_features="sqrt", n_jobs=-1,
                               random_state=C.RANDOM_STATE)
    rsf.fit(Xs, y)

    def cindex_scorer(estimator, X, y_true):
        risk = estimator.predict(X)
        return concordance_index_censored(y_true["event"], y_true["time"], risk)[0]

    imp = permutation_importance(rsf, Xs, y, scoring=cindex_scorer,
                                 n_repeats=n_repeats, random_state=C.RANDOM_STATE,
                                 n_jobs=1)
    out = pd.DataFrame({"feature": Xs.columns,
                        "importance": imp.importances_mean,
                        "importance_sd": imp.importances_std})
    return out.sort_values("importance", ascending=False).reset_index(drop=True)


def risk_strata_km(ds: Dataset, oof_risk: np.ndarray):
    """Tertiles of model risk -> KM curves + global log-rank."""
    # rank-transform first so tied/degenerate risk scores can't break qcut
    ranks = pd.Series(oof_risk).rank(method="first")
    q = pd.qcut(ranks, 3, labels=["Low risk", "Intermediate", "High risk"])
    df = pd.DataFrame({"time": ds.time, "event": ds.event, "group": q})
    lr = multivariate_logrank_test(df["time"], df["group"], df["event"])
    km = {}
    for g in ["Low risk", "Intermediate", "High risk"]:
        sub = df[df["group"] == g]
        kmf = KaplanMeierFitter().fit(sub["time"], sub["event"], label=g)
        km[g] = kmf
    return km, df, float(lr.p_value)


def _cox_hr(df: pd.DataFrame, covars: list[str], target: str) -> dict | None:
    """Cox HR of `target`. Unpenalized for valid Wald inference; only if that
    fails to converge (separation in a rare marker) do we fall back to a small
    ridge and flag the estimate as penalized."""
    d = df[["time", "event"] + covars].astype(float)
    for pen in (0.0, 0.1):
        try:
            cph = CoxPHFitter(penalizer=pen).fit(d, duration_col="time",
                                                 event_col="event")
            s = cph.summary.loc[target]
            return {"hr": float(s["exp(coef)"]),
                    "lo": float(s["exp(coef) lower 95%"]),
                    "hi": float(s["exp(coef) upper 95%"]),
                    "p": float(s["p"]), "n": int(len(d)),
                    "penalized": pen > 0}
        except Exception:
            continue
    return None


def resistance_tmb_confounding(ds: Dataset) -> dict:
    """
    *The key analysis.* Canonical resistance-gene mutations correlate with HIGH
    (favorable) TMB, so their harmful effect is MASKED (positive confounding) — a
    naive subgroup split is null. Holding cancer type fixed (within NSCLC, the
    Skoulidis context) and then adding TMB to the model *unmasks* them. We also
    report the pan-cancer effect adjusted for cancer type, which is the dominant
    confounder. This is why an integrated model that sees TMB *and* genes *and*
    cancer type beats TMB alone — and why eyeballing high-TMB cases fails.
    """
    res_genes = [g for g in C.HIGH_TMB_RESISTANCE_SET if g in ds.X.columns]
    df = pd.DataFrame({
        "time": ds.time, "event": ds.event,
        "tmb": ds.meta["tmb"].to_numpy(),
        "log_tmb": ds.meta["log_tmb"].to_numpy(),
        "cancer_type": ds.meta["cancer_type"].to_numpy(),
    })
    for g in res_genes:
        df[g] = ds.X[g].to_numpy()
    df["res"] = (df[res_genes].sum(axis=1) > 0).astype(int)
    ct_dummies = pd.get_dummies(df["cancer_type"], prefix="ct", drop_first=True)
    df = pd.concat([df, ct_dummies.astype(float)], axis=1)
    ct_cols = list(ct_dummies.columns)

    nsclc = df[df["cancer_type"] == C.NSCLC_LABEL].copy()
    assert len(nsclc) > 0, f"no rows for cancer type {C.NSCLC_LABEL!r}"

    # 1) the confound: do resistance-mutated tumors carry higher TMB?
    mut_tmb = df.loc[df["res"] == 1, "tmb"]
    wt_tmb = df.loc[df["res"] == 0, "tmb"]
    mw = mannwhitneyu(mut_tmb, wt_tmb, alternative="greater")

    # 2) HR before vs after adding TMB. NSCLC panels hold cancer type fixed by
    #    construction; pan-cancer panels adjust for cancer type explicitly.
    panels = [
        ("STK11 (NSCLC)", nsclc, "STK11", []),
        ("KEAP1 (NSCLC)", nsclc, "KEAP1", []),
        ("Resistance set (NSCLC)", nsclc, "res", []),
        ("STK11 (pan-cancer, type-adj.)", df, "STK11", ct_cols),
        ("Resistance set (pan-cancer, type-adj.)", df, "res", ct_cols),
    ]
    rows, skipped = [], []
    for name, sub, gene, base in panels:
        n_mut = int(sub[gene].sum())
        n_mut_events = int(sub.loc[sub[gene] == 1, "event"].sum())
        if n_mut < 5 or n_mut_events < C.MIN_EVENTS_PER_ARM:
            skipped.append(f"{name} (mut events={n_mut_events})")
            continue
        before = _cox_hr(sub, base + [gene], gene)            # WITHOUT TMB
        after = _cox_hr(sub, base + ["log_tmb", gene], gene)  # WITH TMB
        if before is None or after is None:
            skipped.append(f"{name} (Cox failed)")
            continue
        rows.append({"marker": name, "n": after["n"], "n_mut": n_mut,
                     "n_mut_events": n_mut_events,
                     "unadj_hr": before["hr"], "unadj_p": before["p"],
                     "adj_hr": after["hr"], "adj_lo": after["lo"],
                     "adj_hi": after["hi"], "adj_p": after["p"],
                     "penalized": bool(before["penalized"] or after["penalized"])})
    if skipped:
        print(f"  [note] resistance markers skipped (underpowered): {', '.join(skipped)}")
    hr_table = pd.DataFrame(rows)

    return {
        "genes": res_genes,
        "nsclc_n": int(len(nsclc)),
        "tmb_median_mut": float(mut_tmb.median()),
        "tmb_median_wt": float(wt_tmb.median()),
        "tmb_mannwhitney_p": float(mw.pvalue),
        "hr_table": hr_table,
        "df": df,
    }
