"""
Published / canonical transcriptomic signatures for ICI response, plus a curated
immune/TME gene panel. Each signature carries a `direction`: +1 if higher score
is expected to mean BETTER response (so all scores are oriented the same way).

References:
  IFNg_Ayers18  - Ayers et al., J Clin Invest 2017 (T-cell-inflamed GEP)
  IFNg_core     - Ayers 2017 preliminary IFN-gamma set
  CYT           - Rooney et al., Cell 2015 (cytolytic activity: GZMA, PRF1)
  TGFb_EMT      - Mariathasan et al., Nature 2018 (TGF-beta / stroma -> resistance)
  Exhaustion / checkpoints, B-cell/TLS, M2 myeloid, angiogenesis: standard TME sets
"""
from __future__ import annotations

SIGNATURES: dict[str, dict] = {
    "IFNg_Ayers18": {"dir": +1, "genes": [
        "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
        "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
        "PSMB10", "STAT1", "TIGIT"]},
    "IFNg_core": {"dir": +1, "genes": [
        "IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "CXCL11", "HLA-DRA",
        "CCR5", "PRF1", "GZMA"]},
    "CYT": {"dir": +1, "genes": ["GZMA", "PRF1"]},
    "Tcell_CD8": {"dir": +1, "genes": [
        "CD8A", "CD8B", "CD3D", "CD3E", "CD2", "GZMK", "GZMH", "CCL5"]},
    "Checkpoints": {"dir": +1, "genes": [
        "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD274", "PDCD1LG2", "BTLA"]},
    "Bcell_TLS": {"dir": +1, "genes": [
        "MS4A1", "CD79A", "CD79B", "CD19", "TNFRSF13B", "CXCL13"]},
    "NK": {"dir": +1, "genes": ["KLRD1", "KLRK1", "NCR1", "NKG7", "GNLY"]},
    "TGFb_EMT": {"dir": -1, "genes": [
        "TGFB1", "TGFBR2", "ACTA2", "TAGLN", "COL5A1", "COL4A1", "ZEB1", "SNAI2",
        "VIM", "FN1", "FAP", "PDGFRB"]},
    "Myeloid_M2": {"dir": -1, "genes": [
        "CD163", "MRC1", "CSF1R", "IL10", "MSR1", "MARCO"]},
    "Angiogenesis": {"dir": -1, "genes": [
        "VEGFA", "KDR", "PECAM1", "CD34", "ANGPT2", "ESM1"]},
    "Proliferation": {"dir": 0, "genes": ["MKI67", "TOP2A", "BUB1", "CCNB1", "CDK1"]},
}

# Extra immune/TME genes (beyond the signatures) to enrich the learned-model panel.
EXTRA_GENES = [
    "GZMB", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "B2M", "TAP1", "TAP2", "HLA-A",
    "HLA-B", "HLA-C", "HLA-DRB5", "CD4", "FOXP3", "IL2RA", "ENTPD1", "TOX",
    "CXCL10", "CXCL11", "CCL2", "CCL4", "ICOS", "TNFRSF9", "TNFRSF18", "CD28",
    "CD40LG", "SELL", "TCF7", "IL7R", "PRDM1", "BATF", "IRF4", "EOMES", "TBX21",
    "ITGAE", "CD69", "CXCR3", "CCR7", "S100A8", "S100A9", "ARG1", "NOS2", "IDO2",
    "VSIR", "CD276", "VTCN1", "LGALS9", "CEACAM1", "MKI67", "PDGFRA", "THBS1",
    "SPP1", "MMP9", "TIMP1", "IL6", "TNF", "IFNG", "CXCL9",
]


def all_genes() -> list[str]:
    g = set(EXTRA_GENES)
    for s in SIGNATURES.values():
        g.update(s["genes"])
    return sorted(g)
