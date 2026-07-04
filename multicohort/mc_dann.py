"""
Domain-adversarial network (DANN; Ganin et al. 2016) for cross-cohort ICI
response prediction. A shared encoder feeds (a) a response classifier and
(b) a cohort discriminator through a gradient-reversal layer, so the learned
representation is predictive of response but invariant to which cohort a sample
came from -> better leave-one-cohort-out (LOCO) transfer.

Compared under identical LOCO splits against:
  * plain MLP (same net, adversary off / lambda=0)
  * Logistic regression on the 11 signatures (classical baseline)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from mc_data import PROC
from mc_signatures import SIGNATURES

torch.set_num_threads(4)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
MIN_CLASS = 12


class GradReverse(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, lamb):
        ctx.lamb = lamb
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad):
        return -ctx.lamb * grad, None


def grad_reverse(x, lamb):
    return GradReverse.apply(x, lamb)


class DANN(nn.Module):
    def __init__(self, d_in, n_domains, hidden=128, z=64, p=0.3):
        super().__init__()
        self.enc = nn.Sequential(
            nn.Linear(d_in, hidden), nn.BatchNorm1d(hidden), nn.ReLU(), nn.Dropout(p),
            nn.Linear(hidden, z), nn.BatchNorm1d(z), nn.ReLU(), nn.Dropout(p))
        self.cls = nn.Linear(z, 1)
        self.dom = nn.Sequential(nn.Linear(z, 64), nn.ReLU(), nn.Linear(64, n_domains))

    def forward(self, x, lamb=0.0):
        z = self.enc(x)
        return self.cls(z).squeeze(1), self.dom(grad_reverse(z, lamb))


def _train(Xtr, ytr, dtr, n_dom, lamb_max=1.0, epochs=200, lr=1e-3, wd=1e-4, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    Xtr = torch.tensor(Xtr, dtype=torch.float32, device=DEVICE)
    ytr = torch.tensor(ytr, dtype=torch.float32, device=DEVICE)
    dtr = torch.tensor(dtr, dtype=torch.long, device=DEVICE)
    pos_w = torch.tensor([(len(ytr) - ytr.sum()) / max(ytr.sum(), 1)], device=DEVICE)
    net = DANN(Xtr.shape[1], n_dom).to(DEVICE)
    opt = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=wd)
    bce = nn.BCEWithLogitsLoss(pos_weight=pos_w)
    ce = nn.CrossEntropyLoss()
    n = len(ytr); bs = min(256, n)
    for ep in range(epochs):
        lamb = lamb_max * (2 / (1 + np.exp(-5 * ep / epochs)) - 1)  # ramp 0->lamb_max
        net.train(); perm = torch.randperm(n, device=DEVICE)
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            if len(idx) < 4:
                continue
            opt.zero_grad()
            yl, dl = net(Xtr[idx], lamb)
            loss = bce(yl, ytr[idx]) + (lamb > 0) * ce(dl, dtr[idx])
            loss.backward(); opt.step()
    net.eval()
    return net


@torch.no_grad()
def _predict(net, X):
    X = torch.tensor(X, dtype=torch.float32, device=DEVICE)
    yl, _ = net(X, 0.0)
    return torch.sigmoid(yl).cpu().numpy()


def loco_dann(df, feat_cols, lamb_max=1.0, seeds=(0, 1, 2), epochs=200):
    per = {}
    for held in df["cohort"].unique():
        tr = df[df["cohort"] != held].copy()
        te = df[df["cohort"] == held]
        if min(te["response"].value_counts()) < MIN_CLASS:
            continue
        sc = StandardScaler().fit(tr[feat_cols].fillna(0).values)
        Xtr = sc.transform(tr[feat_cols].fillna(0).values)
        Xte = sc.transform(te[feat_cols].fillna(0).values)
        doms, dmap = pd.factorize(tr["cohort"])
        preds = []
        for s in seeds:
            net = _train(Xtr, tr["response"].values.astype(float), doms,
                         len(dmap), lamb_max=lamb_max, epochs=epochs, seed=s)
            preds.append(_predict(net, Xte))
        p = np.mean(preds, axis=0)
        per[held] = roc_auc_score(te["response"].values, p)
    return float(np.mean(list(per.values()))), per


def loco_logistic(df, cols):
    per = {}
    for held in df["cohort"].unique():
        tr, te = df[df["cohort"] != held], df[df["cohort"] == held]
        if min(te["response"].value_counts()) < MIN_CLASS:
            continue
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
        m.fit(tr[cols].fillna(0).values, tr["response"].values)
        per[held] = roc_auc_score(te["response"].values,
                                  m.predict_proba(te[cols].fillna(0).values)[:, 1])
    return float(np.mean(list(per.values()))), per


def main(seeds=(0, 1, 2), epochs=150):
    df = pd.read_csv(PROC / "multicohort.csv", index_col=0)
    sig_cols = [c for c in SIGNATURES if c in df.columns]
    gene_cols = [c for c in df.columns if c.startswith("g_")]
    print(f"device={DEVICE}  n={len(df)}  genes={len(gene_cols)}", flush=True)

    results = {}

    def record(name, lo, per):
        results[name] = (lo, per)
        print(f"  {name:32} {lo:.3f}", flush=True)
        pd.DataFrame([{"method": k, "LOCO_AUC": v[0]} for k, v in results.items()]
                     ).to_csv(PROC / "dann_loco.csv", index=False)
        pd.DataFrame({k: v[1] for k, v in results.items()}).to_csv(
            PROC / "dann_per_cohort.csv")

    print("\nLOCO AUC (mean over held-out cohorts with >=12/class):", flush=True)
    record("Logistic (signatures)", *loco_logistic(df, sig_cols))
    record("Logistic (126 genes)", *loco_logistic(df, gene_cols))
    record("MLP (adversary off)", *loco_dann(df, gene_cols, 0.0, seeds, epochs))
    record("DANN (lambda=1.0)", *loco_dann(df, gene_cols, 1.0, seeds, epochs))
    record("DANN (lambda=2.0)", *loco_dann(df, gene_cols, 2.0, seeds, epochs))

    dann = [k for k in results if k.startswith("DANN")]
    best = max(dann, key=lambda k: results[k][0])
    print(f"\nBest DANN: {best} = {results[best][0]:.3f}  vs "
          f"logistic-sigs {results['Logistic (signatures)'][0]:.3f}", flush=True)


if __name__ == "__main__":
    main()
