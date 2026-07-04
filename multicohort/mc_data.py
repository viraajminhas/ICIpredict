"""
Harmonized multi-cohort ICI-RNA dataset from cBioPortal (CRI iAtlas cohorts).

For each cohort: pull expression for a curated immune/TME gene panel + binary
response label, z-score genes WITHIN cohort (removes platform/batch scale), then
compute signature scores. Output one tidy table across all cohorts, LOCO-ready.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from mc_signatures import SIGNATURES, all_genes

API = "https://www.cbioportal.org/api"
ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"; RAW.mkdir(exist_ok=True)
PROC = ROOT / "processed"; PROC.mkdir(exist_ok=True)

# iAtlas-harmonized ICI cohorts with RNA + binary response (RESPONDER/CLINICAL_BENEFIT)
COHORTS = {
    "blca_iatlas_imvigor210_2017": "Bladder",
    "rcc_iatlas_immotion150_2018": "Renal",
    "mel_iatlas_liu_2019": "Melanoma",
    "mel_iatlas_riaz_nivolumab_2017": "Melanoma",
    "mel_iatlas_gide_2019": "Melanoma",
    "mel_iatlas_hugo_ucla_2016": "Melanoma",
    "paad_iatlas_prince_2022": "Pancreatic",
    "brca_iatlas_anders_2022": "TNBC",
    "gbm_iatlas_prins_2019": "Glioma",
    "ccrcc_iatlas_choueiri_2016": "Renal",
}

S = requests.Session(); S.headers.update({"Accept": "application/json"})


def _get(path, **kw):
    for a in range(4):
        try:
            r = S.request(kw.pop("method", "GET"), f"{API}/{path}", timeout=120, **kw)
            if r.status_code == 200:
                return r.json()
        except requests.RequestException:
            pass
        time.sleep(2 ** a)
    raise RuntimeError(f"failed: {path}")


def _entrez(symbols):
    fp = RAW / "entrez_map.json"
    if fp.exists():
        m = json.load(open(fp))
    else:
        g = _get("genes/fetch", method="POST",
                 params={"geneIdType": "HUGO_GENE_SYMBOL"}, json=symbols)
        m = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in g}
        json.dump(m, open(fp, "w"))
    return m


def _expr_profile(sid):
    profs = _get(f"studies/{sid}/molecular-profiles")
    expr = [p for p in profs if p["molecularAlterationType"] == "MRNA_EXPRESSION"]
    # prefer a continuous (non-zscore) profile; else fall back to the first
    cont = [p for p in expr if p.get("datatype") == "CONTINUOUS"]
    return (cont or expr)[0]["molecularProfileId"]


RECIST_POS = {"CR", "PR", "COMPLETE RESPONSE", "PARTIAL RESPONSE"}
RECIST_NEG = {"SD", "PD", "STABLE DISEASE", "PROGRESSIVE DISEASE",
              "CLINICAL PROGRESSIVE DISEASE"}
TRUE_SET = {"TRUE", "YES", "1", "R", "RESPONDER"}
FALSE_SET = {"FALSE", "NO", "0", "NR", "NON-RESPONDER"}


def _labels(sid):
    """Return (sample_attrs, patient_attrs, sample->patient map)."""
    samp, pat, s2p = {}, {}, {}
    for dtype in ("SAMPLE", "PATIENT"):
        rows = _get(f"studies/{sid}/clinical-data",
                    params={"clinicalDataType": dtype, "projection": "DETAILED",
                            "pageSize": 10_000_000})
        for r in rows:
            if dtype == "SAMPLE":
                samp.setdefault(r["sampleId"], {})[r["clinicalAttributeId"]] = r["value"]
                s2p[r["sampleId"]] = r.get("patientId")
            else:
                pat.setdefault(r["patientId"], {})[r["clinicalAttributeId"]] = r["value"]
    return samp, pat, s2p


def _binary_response(sample_id, samp, pat, s2p):
    """Objective response (RECIST) is primary and keyed by sampleId; fall back to
    patient-level RESPONDER / CLINICAL_BENEFIT."""
    srec = samp.get(sample_id, {})
    resp = str(srec.get("RESPONSE", "")).strip().upper()
    if resp in RECIST_POS:
        return 1
    if resp in RECIST_NEG:
        return 0
    prec = pat.get(s2p.get(sample_id), {})
    for key in ("RESPONDER", "CLINICAL_BENEFIT"):
        v = str(prec.get(key, "")).strip().upper()
        if v in TRUE_SET:
            return 1
        if v in FALSE_SET:
            return 0
    return np.nan


def fetch_cohort(sid, force=False):
    fp = RAW / f"{sid}.csv"
    if fp.exists() and not force:
        return pd.read_csv(fp, index_col=0)

    genes = all_genes()
    emap = _entrez(genes)
    eids = [emap[g] for g in genes if g in emap]
    inv = {v: k for k, v in emap.items()}
    pid = _expr_profile(sid)

    data = _get(f"molecular-profiles/{pid}/molecular-data/fetch",
                method="POST", params={"projection": "SUMMARY"},
                json={"entrezGeneIds": eids, "sampleListId": f"{sid}_all"})
    df = pd.DataFrame([{"sampleId": d["sampleId"],
                        "gene": inv.get(d["entrezGeneId"]),
                        "value": d["value"]} for d in data if d.get("value") is not None])
    mat = df.pivot_table(index="sampleId", columns="gene", values="value", aggfunc="mean")

    # labels
    samp, pat, s2p = _labels(sid)
    y = pd.Series({s: _binary_response(s, samp, pat, s2p) for s in mat.index},
                  name="response")
    mat = mat.assign(response=y)
    mat.to_csv(fp)
    return mat


def build(force=False) -> pd.DataFrame:
    frames = []
    for sid, cancer in COHORTS.items():
        raw = fetch_cohort(sid, force=force)
        y = raw["response"]
        expr = raw.drop(columns=["response"])
        keep = y.notna()
        expr, y = expr[keep], y[keep]
        if len(y) < 8 or y.nunique() < 2:
            print(f"  [skip] {sid}: n={len(y)} responders={int(y.sum()) if len(y) else 0}")
            continue
        # within-cohort z-score per gene (batch/scale invariant)
        z = (expr - expr.mean()) / expr.std(ddof=0).replace(0, np.nan)
        # signature scores (signed, mean of available z-scored genes)
        sig = {}
        for name, s in SIGNATURES.items():
            cols = [g for g in s["genes"] if g in z.columns]
            if cols:
                sc = z[cols].mean(axis=1)
                sig[name] = sc * (s["dir"] if s["dir"] != 0 else 1)
        sigdf = pd.DataFrame(sig)
        block = pd.concat([sigdf, z.add_prefix("g_")], axis=1)
        block.insert(0, "cohort", sid)
        block.insert(1, "cancer", cancer)
        block.insert(2, "response", y.astype(int))
        frames.append(block)
        print(f"  {sid:34} n={len(y):>3}  responders={int(y.sum()):>3} "
              f"({y.mean()*100:.0f}%)  genes={z.shape[1]}")
    full = pd.concat(frames, axis=0)
    full.to_csv(PROC / "multicohort.csv")
    return full


if __name__ == "__main__":
    df = build()
    print(f"\nTOTAL: n={len(df)}  cohorts={df['cohort'].nunique()}  "
          f"cancers={df['cancer'].nunique()}  responders={int(df['response'].sum())} "
          f"({df['response'].mean()*100:.0f}%)")
    print(f"signatures={list(SIGNATURES)}  gene features={sum(c.startswith('g_') for c in df.columns)}")
