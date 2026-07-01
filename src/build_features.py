"""
Build the modelling table: survival labels + a clean feature matrix.

Design decisions (all defensible, documented inline):
  * Gene-mutation indicators are restricted to the IMPACT341 *core* gene set
    (genes sequenced on ALL three panels) so that a 0 means "wild-type", never
    "not sequenced". This removes a real panel-version confounder.
  * TMB is panel-normalised (mut/Mb) by MSK, so it is comparable across panels;
    we model log1p(TMB).
  * Outcome is overall survival (OS_MONTHS, OS_STATUS) — the endpoint Samstein
    used. Rows with non-positive follow-up are dropped.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import config as C
from . import cbioportal


@dataclass
class Dataset:
    X: pd.DataFrame                      # feature matrix (numeric / binary), index=sampleId
    time: np.ndarray                     # OS months
    event: np.ndarray                    # 1 = death
    meta: pd.DataFrame                   # sampleId, patientId, CANCER_TYPE, tmb, panel, etc.
    feature_groups: dict[str, list[str]] = field(default_factory=dict)
    gene_cols: list[str] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.X)


_AGE_MIDPOINT = {"<30": 25.0, "31-50": 40.0, "50-60": 55.0, "61-70": 65.0, ">71": 75.0}


def _core_genes() -> set[str]:
    fp = C.DATA_RAW / "impact_core_genes.json"
    if fp.exists():
        with open(fp) as fh:
            return set(json.load(fh))
    # Fallback: derive from the cBioPortal gene-panel API.
    import requests
    panels = {}
    for pid in ("IMPACT341", "IMPACT410", "IMPACT468"):
        r = requests.get(f"{C.API_BASE}/gene-panels/{pid}", timeout=60).json()
        panels[pid] = {g["hugoGeneSymbol"] for g in r.get("genes", [])}
    core = set.intersection(*panels.values())
    with open(fp, "w") as fh:
        json.dump(sorted(core), fh)
    return core


def _gene_indicator_matrix(mutations: pd.DataFrame, samples: list[str],
                           core: set[str]) -> tuple[pd.DataFrame, list[str]]:
    """One binary column per selected gene: 1 if non-silently mutated in sample."""
    m = mutations.copy()
    m = m[m["hugoGeneSymbol"].notna()]
    m = m[~m["mutationType"].isin(C.SILENT_TYPES)]
    m = m[m["hugoGeneSymbol"].isin(core)]              # only reliably-sequenced genes

    # sample x gene presence
    pres = (
        m.assign(val=1)
        .pivot_table(index="sampleId", columns="hugoGeneSymbol", values="val",
                     aggfunc="max", fill_value=0)
    )
    pres = pres.reindex(index=samples, fill_value=0).astype(int)

    freq = pres.mean(axis=0)
    # Curated genes that are in core AND mutated at all.
    curated_in = [g for g in C.CURATED_GENES if g in pres.columns and freq[g] > 0]
    # Data-driven: most recurrently mutated core genes meeting the freq floor.
    recurrent = (
        freq[freq >= C.MIN_GENE_FREQ]
        .sort_values(ascending=False)
        .head(C.TOP_K_RECURRENT)
        .index.tolist()
    )
    keep = sorted(set(curated_in) | set(recurrent))
    return pres[keep], keep


def build(force_download: bool = False) -> Dataset:
    patient, sample, mutations = cbioportal.download_all(force=force_download)
    core = _core_genes()

    df = sample.merge(patient, on="patientId", how="left", suffixes=("", "_pat"))

    # ---- survival labels ----
    # Only accept the coded values "1:DECEASED" / "0:LIVING"; a missing/garbled
    # status must be DROPPED, never silently treated as censored (event=0).
    status = df["OS_STATUS"].astype(str)
    df = df[status.str.match(r"^[01]:")].copy()
    df["event"] = df["OS_STATUS"].astype(str).str.contains("DECEASED").astype(int)
    df["time"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
    df = df[(df["time"].notna()) & (df["time"] > 0)].copy()

    # ---- continuous / categorical clinical features ----
    df["tmb"] = pd.to_numeric(df["TMB_NONSYNONYMOUS"], errors="coerce")
    df = df[df["tmb"].notna()].copy()
    df["log_tmb"] = np.log1p(df["tmb"])

    age_years = pd.to_numeric(df.get("AGE_AT_SEQ_REPORT"), errors="coerce") / 365.25
    age_fallback = df["AGE_GROUP"].map(_AGE_MIDPOINT)
    df["age"] = age_years.where(age_years.between(10, 100), age_fallback)
    df["age"] = df["age"].fillna(df["age"].median())

    df["is_male"] = (df["SEX"].astype(str) == "Male").astype(int)
    df["is_metastasis"] = (df["SAMPLE_TYPE"].astype(str) == "Metastasis").astype(int)

    # cancer type: keep types with >= MIN_TYPE_N, fold the rest into "Other"
    type_counts = df["CANCER_TYPE"].value_counts()
    big_types = type_counts[type_counts >= C.MIN_TYPE_N].index
    df["cancer_type"] = df["CANCER_TYPE"].where(df["CANCER_TYPE"].isin(big_types), "Other")

    # ---- one-hot encodings ----
    drug = pd.get_dummies(df["DRUG_TYPE"].astype(str), prefix="drug").astype(int)
    ctype = pd.get_dummies(df["cancer_type"].astype(str), prefix="ct").astype(int)

    # ---- gene indicators ----
    genes, gene_cols = _gene_indicator_matrix(mutations, df["sampleId"].tolist(), core)
    genes.index = df.index  # align positionally via sampleId reindex above

    # assemble feature matrix
    clinical_cols = ["log_tmb", "age", "is_male", "is_metastasis"]
    X = pd.concat(
        [df[clinical_cols].reset_index(drop=True),
         drug.reset_index(drop=True),
         ctype.reset_index(drop=True),
         genes.reset_index(drop=True)],
        axis=1,
    )
    X.index = df["sampleId"].values

    feature_groups = {
        "tmb": ["log_tmb"],
        "demographics": ["age", "is_male", "is_metastasis"],
        "drug": list(drug.columns),
        "cancer_type": list(ctype.columns),
        "genes": gene_cols,
    }

    meta = df[["sampleId", "patientId", "CANCER_TYPE", "cancer_type", "tmb",
               "log_tmb", "GENE_PANEL", "DRUG_TYPE", "time", "event"]].reset_index(drop=True)

    ds = Dataset(X=X, time=df["time"].to_numpy(float), event=df["event"].to_numpy(int),
                 meta=meta, feature_groups=feature_groups, gene_cols=gene_cols)

    # cache processed
    X.assign(time=ds.time, event=ds.event).to_csv(C.DATA_PROCESSED / "features.csv")
    meta.to_csv(C.DATA_PROCESSED / "meta.csv", index=False)
    with open(C.DATA_PROCESSED / "feature_groups.json", "w") as fh:
        json.dump(feature_groups, fh, indent=2)
    return ds


if __name__ == "__main__":
    ds = build()
    print(f"Final cohort: n={ds.n} | events={ds.event.sum()} "
          f"({ds.event.mean()*100:.1f}%) | features={ds.X.shape[1]}")
    print(f"Gene features ({len(ds.gene_cols)}): {ds.gene_cols}")
    print("Feature groups:", {k: len(v) for k, v in ds.feature_groups.items()})
    print("\nCancer types:")
    print(ds.meta['cancer_type'].value_counts().to_string())
