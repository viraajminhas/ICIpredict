"""
Minimal cBioPortal REST client + on-disk cache.

Pulls the Samstein 2019 ICI cohort (study `tmb_mskcc_2018`):
  * patient-level clinical data   (OS, drug type, sex, age group)
  * sample-level clinical data    (TMB, cancer type, sample type, gene panel)
  * somatic mutation calls        (per-sample, per-gene)

Everything is cached as CSV under data/raw so the rest of the pipeline is fully
reproducible and offline-friendly after the first run.
"""
from __future__ import annotations

import time

import pandas as pd
import requests

from . import config as C

_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json", "User-Agent": "ICIpredict/1.0"})


_RETRYABLE = {429, 500, 502, 503, 504}


def _request(method: str, url: str, **kwargs) -> requests.Response:
    """HTTP with exponential-backoff retry on transient errors only.
    Non-retryable HTTP statuses (e.g. 4xx) raise immediately."""
    last_exc: Exception | None = None
    for attempt in range(5):
        try:
            resp = _SESSION.request(method, url, timeout=120, **kwargs)
        except requests.RequestException as exc:        # network blip -> retry
            last_exc = exc
        else:
            if resp.status_code == 200:
                return resp
            if resp.status_code not in _RETRYABLE:      # fail fast, do not retry
                resp.raise_for_status()
            last_exc = requests.HTTPError(f"{resp.status_code} from {url}")
        wait = 2 ** attempt
        print(f"  [retry {attempt + 1}/5] {url} -> waiting {wait}s")
        time.sleep(wait)
    if last_exc:
        raise last_exc
    raise RuntimeError(f"Failed after retries: {url}")


def _clinical_wide(clinical_data_type: str) -> pd.DataFrame:
    """Fetch PATIENT or SAMPLE clinical data and pivot to a wide table."""
    url = f"{C.API_BASE}/studies/{C.STUDY_ID}/clinical-data"
    params = {
        "clinicalDataType": clinical_data_type,
        "projection": "DETAILED",
        "pageSize": 10_000_000,
    }
    rows = _request("GET", url, params=params).json()
    long = pd.DataFrame(rows)
    id_col = "patientId" if clinical_data_type == "PATIENT" else "sampleId"
    index_cols = list(dict.fromkeys(
        c for c in (id_col, "patientId") if c in long.columns
    ))
    wide = long.pivot_table(
        index=index_cols,
        columns="clinicalAttributeId",
        values="value",
        aggfunc="first",
    ).reset_index()
    wide.columns.name = None
    return wide


def fetch_clinical(force: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    patient_fp = C.DATA_RAW / "clinical_patient.csv"
    sample_fp = C.DATA_RAW / "clinical_sample.csv"
    if not force and patient_fp.exists() and sample_fp.exists():
        return pd.read_csv(patient_fp), pd.read_csv(sample_fp)

    print("Downloading patient-level clinical data ...")
    patient = _clinical_wide("PATIENT")
    print(f"  patients: {len(patient)}")
    print("Downloading sample-level clinical data ...")
    sample = _clinical_wide("SAMPLE")
    print(f"  samples:  {len(sample)}")

    patient.to_csv(patient_fp, index=False)
    sample.to_csv(sample_fp, index=False)
    return patient, sample


def fetch_mutations(force: bool = False) -> pd.DataFrame:
    fp = C.DATA_RAW / "mutations.csv"
    if not force and fp.exists():
        return pd.read_csv(fp)

    print("Downloading somatic mutation calls ...")
    url = f"{C.API_BASE}/molecular-profiles/{C.MUTATION_PROFILE_ID}/mutations/fetch"
    body = {"sampleListId": C.SEQUENCED_SAMPLE_LIST}
    records = _request("POST", url, params={"projection": "DETAILED"},
                       json=body, headers={"Content-Type": "application/json"}).json()
    print(f"  mutation records: {len(records)}")

    def hugo(rec: dict) -> str | None:
        g = rec.get("gene")
        if isinstance(g, dict):
            return g.get("hugoGeneSymbol")
        return rec.get("hugoGeneSymbol")

    muts = pd.DataFrame(
        {
            "sampleId": [r.get("sampleId") for r in records],
            "patientId": [r.get("patientId") for r in records],
            "hugoGeneSymbol": [hugo(r) for r in records],
            "mutationType": [r.get("mutationType") for r in records],
            "proteinChange": [r.get("proteinChange") for r in records],
            "variantType": [r.get("variantType") for r in records],
        }
    )
    muts.to_csv(fp, index=False)
    return muts


def download_all(force: bool = False):
    patient, sample = fetch_clinical(force=force)
    mutations = fetch_mutations(force=force)
    return patient, sample, mutations


if __name__ == "__main__":
    p, s, m = download_all(force=False)
    print("\nPatient cols:", list(p.columns))
    print("Sample cols :", list(s.columns))
    print("Mutation rows:", len(m), "| genes:", m["hugoGeneSymbol"].nunique())
