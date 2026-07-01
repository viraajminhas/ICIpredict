"""
Central configuration: paths, dataset identifiers, gene panels, and modelling
constants for the ICI (immune checkpoint inhibitor) response-prediction project.

Dataset: Samstein et al., "Tumor mutational load predicts survival after
immunotherapy across multiple cancer types", Nature Genetics 2019.
cBioPortal study id: tmb_mskcc_2018 (1,661 ICI-treated patients, pan-cancer,
MSK-IMPACT targeted panel, overall-survival outcomes).

Gene panels below are curated from the project's annotated bibliography:
  - Rizvi 2015 (Science)      -> POLE/POLD1, MMR, KEAP1, smoking/DDR genes
  - Jamieson & Maker 2017     -> acquired resistance: B2M, JAK1/2, antigen presentation
  - Lee/Samstein 2020         -> TMB must be interpreted within cancer type
"""
from __future__ import annotations

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths (project is intentionally outside OneDrive)
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"

for _p in (DATA_RAW, DATA_PROCESSED, RESULTS, FIGURES, TABLES):
    _p.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# cBioPortal API
# ---------------------------------------------------------------------------
API_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "tmb_mskcc_2018"
MUTATION_PROFILE_ID = f"{STUDY_ID}_mutations"
SEQUENCED_SAMPLE_LIST = f"{STUDY_ID}_sequenced"

# ---------------------------------------------------------------------------
# Curated, biology-driven gene panel (grouped, with bibliography rationale)
# Only genes actually present and recurrently mutated in the cohort are kept
# as features (see build_features); absent/ultra-rare genes are dropped.
# ---------------------------------------------------------------------------
GENE_GROUPS: dict[str, list[str]] = {
    # Mismatch-repair deficiency -> hypermutation, highly immunogenic (Rizvi, Jamieson)
    "MMR": ["MSH2", "MSH6", "MLH1", "PMS2", "MSH3", "PMS1", "EPCAM"],
    # Polymerase proofreading -> ultramutated tumors (Rizvi)
    "Proofreading": ["POLE", "POLD1"],
    # Antigen presentation / acquired resistance machinery (Jamieson & Maker)
    "AntigenPresentation": ["B2M", "TAP1", "TAP2", "TAPBP", "CALR", "HLA-A", "HLA-B"],
    # IFN-gamma / JAK-STAT signaling -> resistance when lost (Jamieson & Maker)
    "JAK_IFN": ["JAK1", "JAK2", "STAT1", "IFNGR1", "IFNGR2", "IRF1", "SOCS1"],
    # Driver alterations linked to primary ICI resistance
    "ResistanceDrivers": ["STK11", "KEAP1", "PTEN", "SMARCA4", "NF1", "PBRM1"],
    # Broader DNA-damage-repair genes (immunogenicity)
    "DDR": ["BRCA1", "BRCA2", "ATM", "ATR", "POLQ", "ERCC2", "RAD51", "FANCA", "BAP1"],
    # Common drivers for context / known interactions
    "CommonDrivers": [
        "TP53", "KRAS", "EGFR", "BRAF", "PIK3CA", "APC", "ARID1A", "KMT2D",
        "KMT2C", "RB1", "CDKN2A", "FAT1", "NOTCH1", "LRP1B", "KDM6A", "PIK3R1",
    ],
}
CURATED_GENES: list[str] = sorted({g for genes in GENE_GROUPS.values() for g in genes})

# Genes implicated in *primary resistance despite high TMB* — used in the
# "resolving high-TMB non-responders" analysis (answers Rizvi's open question).
HIGH_TMB_RESISTANCE_SET: list[str] = ["STK11", "KEAP1", "B2M", "JAK1", "JAK2", "PTEN"]

# Genes implicated in enhanced benefit (hypermutation / immunogenic)
BENEFIT_GENE_SET: list[str] = ["POLE", "POLD1", "MSH2", "MSH6", "MLH1", "PMS2"]

# Add the top-K most frequently mutated genes in the cohort as data-driven
# features (union with the curated panel).
TOP_K_RECURRENT = 50
# A gene must be (non-silently) mutated in at least this fraction of samples to
# become a feature column — keeps features informative and avoids noise.
MIN_GENE_FREQ = 0.02

# Mutation types treated as functional (non-silent) for gene indicators.
SILENT_TYPES = {"Silent", "Synonymous", "3'UTR", "5'UTR", "Intron", "IGR",
                "5'Flank", "3'Flank", "RNA"}

# ---------------------------------------------------------------------------
# Clinical / modelling constants
# ---------------------------------------------------------------------------
# FDA tissue-agnostic pembrolizumab TMB cutoff (mut/Mb).
FDA_TMB_CUTOFF = 10.0
# Samstein cancer-type-specific definition of "TMB-high": top quintile within type.
TYPE_SPECIFIC_TOP_PCTL = 80.0
# Minimum samples for a cancer type to be modelled as its own stratum.
MIN_TYPE_N = 30
# Canonical context for the STK11/KEAP1 resistance analysis (Skoulidis et al.).
NSCLC_LABEL = "Non-Small Cell Lung Cancer"
# A marker needs at least this many DEATHS in the mutated arm to be reported.
MIN_EVENTS_PER_ARM = 10

# Cross-validation
N_SPLITS = 5
N_REPEATS = 10
RANDOM_STATE = 20260629
# Time points (months) for time-dependent AUC.
EVAL_TIMES = [6.0, 12.0, 18.0, 24.0]
N_BOOTSTRAP = 2000
