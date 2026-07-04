"""
Profile candidate public ICI cohorts on cBioPortal: what molecular data
(expression / mutation / CNA) and what outcome labels (response / OS / PFS)
does each have? This determines the multi-cohort study design.
"""
from __future__ import annotations
import json
import requests

API = "https://www.cbioportal.org/api"
S = requests.Session(); S.headers.update({"Accept": "application/json"})

COHORTS = [
    "tmb_mskcc_2018", "blca_iatlas_imvigor210_2017", "rcc_iatlas_immotion150_2018",
    "mixed_allen_2018", "nsclc_pd1_msk_2018", "mel_dfci_2019",
    "mel_iatlas_liu_2019", "mel_iatlas_riaz_nivolumab_2017", "mel_iatlas_gide_2019",
    "mel_iatlas_hugo_ucla_2016", "paad_iatlas_prince_2022", "brca_iatlas_anders_2022",
    "gbm_iatlas_prins_2019", "ccrcc_iatlas_choueiri_2016", "mel_ucla_2016",
    "skcm_mskcc_2014",
]

RESPONSE_KW = ["RESPONSE", "RECIST", "BENEFIT", "BOR", "RESPONDER", "DURABLE", "IRRECIST"]
OS_KW = ["OS_STATUS", "OS_MONTHS"]
PFS_KW = ["PFS_STATUS", "PFS_MONTHS", "PFS"]


def get(path):
    r = S.get(f"{API}/{path}", timeout=60)
    return r.json() if r.status_code == 200 else None


rows = []
for sid in COHORTS:
    profs = get(f"studies/{sid}/molecular-profiles") or []
    types = sorted({p.get("molecularAlterationType") for p in profs})
    attrs = get(f"studies/{sid}/clinical-attributes") or []
    ids = [a.get("clinicalAttributeId", "") for a in attrs]
    has_expr = "MRNA_EXPRESSION" in types
    has_mut = "MUTATION_EXTENDED" in types
    has_cna = "COPY_NUMBER_ALTERATION" in types
    resp = [i for i in ids if any(k in i.upper() for k in RESPONSE_KW)]
    os_ = [i for i in ids if i in OS_KW]
    pfs = [i for i in ids if any(k in i.upper() for k in PFS_KW)]
    study = get(f"studies/{sid}") or {}
    rows.append({
        "study": sid, "n": study.get("allSampleCount"),
        "expr": has_expr, "mut": has_mut, "cna": has_cna,
        "response_fields": resp[:4], "os": bool(os_), "pfs": bool(pfs),
    })

print(f"{'study':32} {'n':>4} {'EXPR':>5} {'MUT':>4} {'CNA':>4} {'OS':>3} {'PFS':>4}  response_fields")
for r in rows:
    print(f"{r['study']:32} {str(r['n']):>4} {str(r['expr']):>5} {str(r['mut']):>4} "
          f"{str(r['cna']):>4} {str(r['os']):>3} {str(r['pfs']):>4}  {', '.join(r['response_fields'])}")

# summary
expr_cohorts = [r for r in rows if r["expr"] and (r["response_fields"] or r["os"])]
print(f"\nCohorts with EXPRESSION + (response or OS): {len(expr_cohorts)}")
print("  ", ", ".join(r["study"] for r in expr_cohorts))
json.dump(rows, open("multicohort/cohort_profile.json", "w"), indent=2)
