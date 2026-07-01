"""
Survival models, cross-validated evaluation, and rigorous model comparison.

Endpoint: overall survival (right-censored). All models are scored as *risk*
(higher = worse survival), so Harrell's C, IPCW-C and time-dependent AUC are
directly comparable across models.

Models / baselines
------------------
  FDA_TMB10     : single tissue-agnostic cutoff TMB >= 10 mut/Mb (rule)
  TypeSpecTMB   : Samstein top-20% TMB *within each cancer type* (rule, CV-correct)
  Cox_TMB       : Cox on log(TMB) only                       <- "TMB alone"
  Cox_Clin      : Cox on clinical only (no TMB, no genes)
  Cox_TMB_Clin  : Cox on TMB + clinical (no genes)           <- isolates gene value
  Coxnet_Full   : elastic-net Cox on TMB + clinical + genes  <- interpretable integrated
  RSF_Full      : random survival forest on the full matrix  <- nonlinear integrated
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.linear_model import CoxnetSurvivalAnalysis, CoxPHSurvivalAnalysis
from sksurv.metrics import (concordance_index_censored, concordance_index_ipcw,
                            cumulative_dynamic_auc)
from sksurv.util import Surv

from . import config as C
from .build_features import Dataset

CONTINUOUS = ["log_tmb", "age"]
COXNET_ALPHA_GRID = np.logspace(-3, -1, 10)   # 0.001 .. 0.1


class NestedCoxnetCV(BaseEstimator):
    """Elastic-net Cox that selects its own penalty by *inner* CV on the data it
    is handed. Used inside the outer CV so the penalty is tuned on train folds
    only — removing the non-nested-CV optimism flagged in review."""

    def __init__(self, l1_ratio=0.5, alphas=tuple(COXNET_ALPHA_GRID), cv=4):
        self.l1_ratio = l1_ratio
        self.alphas = alphas
        self.cv = cv

    def fit(self, X, y):
        X = np.asarray(X, float)
        event = y["event"].astype(int)
        skf = StratifiedKFold(self.cv, shuffle=True, random_state=C.RANDOM_STATE)
        scores = np.zeros(len(self.alphas))
        counts = np.zeros(len(self.alphas))
        for tr, te in skf.split(X, event):
            for j, a in enumerate(self.alphas):
                try:
                    m = CoxnetSurvivalAnalysis(l1_ratio=self.l1_ratio, alphas=[a],
                                               fit_baseline_model=False,
                                               max_iter=100_000, normalize=False)
                    m.fit(X[tr], y[tr])
                    risk = m.predict(X[te])
                    scores[j] += concordance_index_censored(
                        y["event"][te], y["time"][te], risk)[0]
                    counts[j] += 1
                except Exception:
                    pass
        mean = np.divide(scores, counts, out=np.full_like(scores, np.nan),
                         where=counts > 0)
        self.best_alpha_ = float(self.alphas[int(np.nanargmax(mean))])
        self.model_ = CoxnetSurvivalAnalysis(
            l1_ratio=self.l1_ratio, alphas=[self.best_alpha_],
            fit_baseline_model=False, max_iter=100_000, normalize=False).fit(X, y)
        return self

    def predict(self, X):
        return self.model_.predict(np.asarray(X, float))


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def to_surv(event: np.ndarray, time: np.ndarray):
    return Surv.from_arrays(event=event.astype(bool), time=time.astype(float))


def make_strata(meta: pd.DataFrame, event: np.ndarray) -> np.ndarray:
    """cancer_type x event strata; rare combos merge into the modal stratum so
    StratifiedKFold never sees a class smaller than n_splits."""
    base = meta["cancer_type"].astype(str).to_numpy() + "|" + event.astype(str)
    s = pd.Series(base)
    counts = s.value_counts()
    rare = counts[counts < C.N_SPLITS].index
    if len(rare):
        s = s.where(~s.isin(rare), other=counts.idxmax())
    return pd.factorize(s)[0]


def _scaler(columns: list[str]) -> ColumnTransformer:
    cont = [c for c in CONTINUOUS if c in columns]
    return ColumnTransformer(
        [("std", StandardScaler(), cont)],
        remainder="passthrough",
    )


# --------------------------------------------------------------------------- #
# model specifications
# --------------------------------------------------------------------------- #
@dataclass
class ModelSpec:
    name: str
    kind: str                       # "estimator" | "rule"
    columns: list[str] | None = None
    builder: object = None          # () -> sklearn-compatible estimator
    rule: object = None             # (train_idx, test_idx, ds) -> risk_test
    label: str = ""


def _rule_fda(train_idx, test_idx, ds: Dataset):
    tmb = ds.meta["tmb"].to_numpy()
    # higher TMB -> lower risk, so risk = -(TMB>=cutoff)
    return -(tmb[test_idx] >= C.FDA_TMB_CUTOFF).astype(float)


def _rule_type_specific(train_idx, test_idx, ds: Dataset):
    """Top-20% TMB within cancer type; thresholds learned on TRAIN only."""
    tmb = ds.meta["tmb"].to_numpy()
    ctype = ds.meta["cancer_type"].to_numpy()
    thresholds = {}
    for t in np.unique(ctype[train_idx]):
        vals = tmb[train_idx][ctype[train_idx] == t]
        thresholds[t] = np.percentile(vals, C.TYPE_SPECIFIC_TOP_PCTL)
    global_thr = np.percentile(tmb[train_idx], C.TYPE_SPECIFIC_TOP_PCTL)
    risk = np.empty(len(test_idx))
    for i, idx in enumerate(test_idx):
        thr = thresholds.get(ctype[idx], global_thr)
        risk[i] = -float(tmb[idx] >= thr)       # high-TMB -> lower risk
    return risk


def build_specs(ds: Dataset) -> list[ModelSpec]:
    all_cols = list(ds.X.columns)
    clin_cols = (ds.feature_groups["demographics"] + ds.feature_groups["drug"]
                 + ds.feature_groups["cancer_type"])
    tmb_clin = ["log_tmb"] + clin_cols

    def cox(alpha=1e-4):
        return lambda: CoxPHSurvivalAnalysis(alpha=alpha)  # tiny ridge for stability

    return [
        ModelSpec("FDA_TMB10", "rule", rule=_rule_fda,
                  label="TMB ≥ 10 (FDA cutoff)"),
        ModelSpec("TypeSpecTMB", "rule", rule=_rule_type_specific,
                  label="Type-specific top-20% TMB (Samstein)"),
        ModelSpec("Cox_TMB", "estimator", columns=["log_tmb"], builder=cox(),
                  label="Cox: TMB alone"),
        ModelSpec("Cox_Clin", "estimator", columns=clin_cols, builder=cox(),
                  label="Cox: clinical only"),
        ModelSpec("Cox_TMB_Clin", "estimator", columns=tmb_clin, builder=cox(),
                  label="Cox: TMB + clinical (no genes)"),
        ModelSpec(
            "Coxnet_Full", "estimator", columns=all_cols,
            builder=lambda: NestedCoxnetCV(l1_ratio=0.5),   # penalty tuned per fold
            label="Elastic-net Cox: TMB + clinical + genes"),
        ModelSpec(
            "RSF_Full", "estimator", columns=all_cols,
            builder=lambda: RandomSurvivalForest(
                n_estimators=300, min_samples_leaf=15, max_features="sqrt",
                n_jobs=-1, random_state=C.RANDOM_STATE),
            label="Random survival forest: full model"),
    ]


def _pipe(spec: ModelSpec) -> Pipeline:
    return Pipeline([("scale", _scaler(spec.columns)), ("model", spec.builder())])


# --------------------------------------------------------------------------- #
# alpha selection for elastic-net Cox (single clean inner CV)
# --------------------------------------------------------------------------- #
def select_coxnet_alpha(ds: Dataset) -> float:
    """Pick the elastic-net penalty for the *interpretation* fit by 5-fold CV.
    (The CV-reported Coxnet_Full uses NestedCoxnetCV, which tunes per fold; this
    single alpha is only for the full-cohort, explicitly in-sample HR figure.)"""
    y = to_surv(ds.event, ds.time)
    X = ds.X
    grid = COXNET_ALPHA_GRID
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=C.RANDOM_STATE)
    strata = make_strata(ds.meta, ds.event)
    scores = np.zeros(len(grid))
    counts = np.zeros(len(grid))
    for tr, te in skf.split(X, strata):
        for j, a in enumerate(grid):
            pipe = Pipeline([("scale", _scaler(list(X.columns))),
                             ("model", CoxnetSurvivalAnalysis(
                                 l1_ratio=0.5, alphas=[a], fit_baseline_model=False,
                                 max_iter=100_000, normalize=False))])
            try:
                pipe.fit(X.iloc[tr], y[tr])
                risk = pipe.predict(X.iloc[te])
                scores[j] += concordance_index_censored(
                    ds.event[te].astype(bool), ds.time[te], risk)[0]
                counts[j] += 1
            except Exception as exc:
                print(f"    [warn] alpha={a:.4f} fold failed: {exc}")
    mean = np.divide(scores, counts, out=np.full_like(scores, np.nan),
                     where=counts > 0)
    best = float(grid[int(np.nanargmax(mean))])
    print(f"  interpretation Coxnet alpha = {best:.4f} "
          f"(CV C-index {np.nanmax(mean):.3f})")
    return best


# --------------------------------------------------------------------------- #
# cross-validated evaluation
# --------------------------------------------------------------------------- #
@dataclass
class CVResult:
    name: str
    label: str
    fold_cindex: list[float] = field(default_factory=list)
    fold_ipcw: list[float] = field(default_factory=list)
    fold_auc: list[float] = field(default_factory=list)   # mean AUC over EVAL_TIMES
    oof_sum: np.ndarray = None
    oof_cnt: np.ndarray = None

    def oof_risk(self) -> np.ndarray:
        return self.oof_sum / np.maximum(self.oof_cnt, 1)

    def summary(self) -> dict:
        return {
            "model": self.name,
            "label": self.label,
            "cindex_mean": float(np.mean(self.fold_cindex)),
            "cindex_sd": float(np.std(self.fold_cindex)),
            "ipcw_mean": float(np.nanmean(self.fold_ipcw)) if self.fold_ipcw else np.nan,
            "auc_mean": float(np.nanmean(self.fold_auc)) if self.fold_auc else np.nan,
        }


def robust_auc(y_tr, y_te, t_te, ev_te, risk, want_times):
    """Time-dependent AUC evaluated PER time point, on the test subset that lies
    within the training follow-up support (cumulative_dynamic_auc's IPCW weights
    are undefined beyond it). Returns (mean_auc_over_valid_times, n_skipped)."""
    tr_max = float(y_tr["time"].max())
    keep = t_te < tr_max
    if keep.sum() < 10 or len(np.unique(risk[keep])) < 2:
        return np.nan, 0
    yk, rk, tk, evk = y_te[keep], risk[keep], t_te[keep], ev_te[keep]
    ev_max = tk[evk].max() if evk.any() else tk.max()
    aucs, n_fail = [], 0
    for tt in want_times:
        if not (tk.min() < tt < min(ev_max, tr_max)):
            continue
        try:
            a, _ = cumulative_dynamic_auc(y_tr, yk, rk, [tt])
            aucs.append(float(a[0]))
        except Exception:
            n_fail += 1
    return (float(np.mean(aucs)) if aucs else np.nan), n_fail


def cross_validate(ds: Dataset, specs: list[ModelSpec]) -> dict[str, CVResult]:
    n = ds.n
    y = to_surv(ds.event, ds.time)
    strata = make_strata(ds.meta, ds.event)
    rskf = RepeatedStratifiedKFold(n_splits=C.N_SPLITS, n_repeats=C.N_REPEATS,
                                   random_state=C.RANDOM_STATE)
    results = {s.name: CVResult(s.name, s.label or s.name,
                                oof_sum=np.zeros(n), oof_cnt=np.zeros(n))
               for s in specs}

    valid_times = np.array(C.EVAL_TIMES, dtype=float)
    n_folds = C.N_SPLITS * C.N_REPEATS
    auc_failures = 0
    for k, (tr, te) in enumerate(rskf.split(ds.X, strata), 1):
        y_tr, y_te = y[tr], y[te]
        ev_te, t_te = ds.event[te].astype(bool), ds.time[te]

        for spec in specs:
            try:
                if spec.kind == "rule":
                    risk = spec.rule(tr, te, ds)
                else:
                    pipe = _pipe(spec)
                    pipe.fit(ds.X[spec.columns].iloc[tr], y_tr)
                    risk = np.asarray(pipe.predict(ds.X[spec.columns].iloc[te]),
                                      dtype=float)
            except Exception as exc:               # numerical edge cases
                print(f"    [warn] {spec.name} fold {k}: {exc}")
                continue

            r = results[spec.name]
            c = concordance_index_censored(ev_te, t_te, risk)[0]
            r.fold_cindex.append(c)
            r.oof_sum[te] += risk
            r.oof_cnt[te] += 1
            # IPCW-C and time-dependent AUC (need train censoring distribution)
            try:
                tau = min(y_tr["time"].max(), t_te.max())
                r.fold_ipcw.append(
                    concordance_index_ipcw(y_tr, y_te, risk, tau=tau)[0])
            except Exception:
                r.fold_ipcw.append(np.nan)
            auc_val, n_fail = robust_auc(y_tr, y_te, t_te, ev_te, risk, valid_times)
            r.fold_auc.append(auc_val)
            auc_failures += n_fail
        if k % C.N_SPLITS == 0:
            print(f"  CV progress: {k}/{n_folds} folds")
    if auc_failures:
        print(f"  [note] {auc_failures} (model x time) AUC points skipped "
              f"(time outside train follow-up support)")
    return results


# --------------------------------------------------------------------------- #
# time-dependent AUC curve (for plotting) for a subset of models
# --------------------------------------------------------------------------- #
def auc_over_time(ds: Dataset, specs: list[ModelSpec], times: list[float],
                  n_repeats: int = 5) -> dict[str, dict]:
    y = to_surv(ds.event, ds.time)
    strata = make_strata(ds.meta, ds.event)
    rskf = RepeatedStratifiedKFold(n_splits=C.N_SPLITS, n_repeats=n_repeats,
                                   random_state=C.RANDOM_STATE)
    times = np.asarray(times, float)
    acc = {s.name: [] for s in specs}     # list of per-fold AUC arrays
    for tr, te in rskf.split(ds.X, strata):
        y_tr, y_te = y[tr], y[te]
        t_te, ev_te = ds.time[te], ds.event[te].astype(bool)
        tr_max = float(y_tr["time"].max())
        keep = t_te < tr_max                      # IPCW support from train
        for spec in specs:
            row = np.full(len(times), np.nan)
            try:
                if spec.kind == "rule":
                    risk = spec.rule(tr, te, ds)
                else:
                    pipe = _pipe(spec)
                    pipe.fit(ds.X[spec.columns].iloc[tr], y_tr)
                    risk = np.asarray(pipe.predict(ds.X[spec.columns].iloc[te]), float)
                rk = risk[keep]
                if keep.sum() >= 10 and len(np.unique(rk)) > 1:
                    yk, tk = y_te[keep], t_te[keep]
                    ev_max = tk[ev_te[keep]].max() if ev_te[keep].any() else tk.max()
                    for i, tt in enumerate(times):
                        if tk.min() < tt < min(ev_max, tr_max):
                            try:
                                a, _ = cumulative_dynamic_auc(y_tr, yk, rk, [tt])
                                row[i] = float(a[0])
                            except Exception:
                                pass
            except Exception:
                pass
            acc[spec.name].append(row)
    out = {}
    for name, rows in acc.items():
        arr = np.vstack(rows)
        out[name] = {"times": times,
                     "auc_mean": np.nanmean(arr, axis=0),
                     "auc_sd": np.nanstd(arr, axis=0)}
    return out


# --------------------------------------------------------------------------- #
# paired bootstrap comparison of two models on pooled OOF risk
# --------------------------------------------------------------------------- #
def bootstrap_cindex_diff(ds: Dataset, risk_a: np.ndarray, risk_b: np.ndarray,
                          n_boot: int = None, seed: int = C.RANDOM_STATE) -> dict:
    """C-index(a) - C-index(b) with a paired patient bootstrap (same resamples)."""
    n_boot = n_boot or C.N_BOOTSTRAP
    rng = np.random.default_rng(seed)
    ev, t = ds.event.astype(bool), ds.time
    n = ds.n

    def cidx(idx, risk):
        return concordance_index_censored(ev[idx], t[idx], risk[idx])[0]

    base = cidx(np.arange(n), risk_a) - cidx(np.arange(n), risk_b)
    diffs = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if ev[idx].sum() < 2:            # need events to define concordance
            diffs[b] = np.nan
            continue
        diffs[b] = cidx(idx, risk_a) - cidx(idx, risk_b)
    diffs = diffs[~np.isnan(diffs)]
    B = len(diffs)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    # two-sided bootstrap p-value with +1 continuity correction (floor 1/(B+1))
    p = 2 * min((np.sum(diffs <= 0) + 1) / (B + 1),
                (np.sum(diffs >= 0) + 1) / (B + 1))
    return {"delta": float(base), "ci_low": float(lo), "ci_high": float(hi),
            "p_value": float(min(p, 1.0)), "n_boot": int(B),
            "c_a": float(cidx(np.arange(n), risk_a)),
            "c_b": float(cidx(np.arange(n), risk_b))}


def paired_fold_test(result_a: "CVResult", result_b: "CVResult") -> dict:
    """Complementary to the bootstrap: a paired test across the identical CV folds
    (both models see the same train/test split each fold), capturing CV/fit
    variability that the pooled-OOF bootstrap ignores. Wilcoxon signed-rank."""
    a = np.asarray(result_a.fold_cindex, float)
    b = np.asarray(result_b.fold_cindex, float)
    m = min(len(a), len(b))
    a, b = a[:m], b[:m]
    d = a - b
    try:
        from scipy.stats import wilcoxon
        p = float(wilcoxon(a, b, zero_method="wilcox").pvalue) if np.any(d != 0) else 1.0
    except Exception:
        p = float("nan")
    return {"delta_mean": float(d.mean()), "delta_sd": float(d.std(ddof=1)),
            "wilcoxon_p": p, "n_folds": int(m)}


def holm_correct(pvals: dict[str, float]) -> dict[str, float]:
    """Holm-Bonferroni step-down adjusted p-values."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    adj, prev = {}, 0.0
    for rank, (k, p) in enumerate(items):
        a = min(1.0, (m - rank) * p)
        prev = max(prev, a)         # enforce monotonicity
        adj[k] = prev
    return adj
