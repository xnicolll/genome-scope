"""Curated tumor-suppressor gene list (Vogelstein 2013 + Sanger CGC, COSMIC fallback).

Confidence: 1.0 = canonical TSG, 0.5 = moderate evidence.
Overridden by backend/data/raw/cosmic/cancer_gene_census.csv if present.
"""

from __future__ import annotations

from .annotations import CosmicCensus

_CANONICAL_TSGS: dict[str, float] = {
    "TP53":    1.0,  # Li-Fraumeni; ~50% of human cancers
    "RB1":     1.0,  # retinoblastoma; first identified TSG
    "APC":     1.0,  # familial adenomatous polyposis, colorectal
    "BRCA1":   1.0,  # hereditary breast/ovarian
    "BRCA2":   1.0,  # hereditary breast/ovarian
    "PTEN":    1.0,  # Cowden syndrome; broad cancer involvement
    "CDKN2A":  1.0,  # p16/INK4a; melanoma + many others
    "VHL":     1.0,  # von Hippel-Lindau; renal
    "WT1":     1.0,  # Wilms tumor
    "NF1":     1.0,  # neurofibromatosis
    "NF2":     1.0,  # schwannomatosis / meningioma
    "STK11":   1.0,  # Peutz-Jeghers; lung adenocarcinoma
    "SMAD4":   1.0,  # juvenile polyposis; pancreatic
    "MLH1":    1.0,  # Lynch syndrome (mismatch repair)
    "MSH2":    1.0,  # Lynch syndrome (mismatch repair)
    "MSH6":    1.0,  # Lynch syndrome
    "PMS2":    1.0,  # Lynch syndrome
    "CDH1":    1.0,  # hereditary diffuse gastric / lobular breast
    "ATM":     1.0,  # ataxia telangiectasia; breast cancer risk
    "TSC1":    1.0,  # tuberous sclerosis
    "TSC2":    1.0,  # tuberous sclerosis
    "PALB2":   1.0,  # Fanconi anemia / breast cancer
    "BAP1":    1.0,  # mesothelioma / uveal melanoma
    "CHEK2":   1.0,  # familial breast cancer
    "RUNX1":   1.0,  # AML; epigenetically silenced
    "ETS2":    0.5,  # context-dependent transcription factor
    "RCAN1":   0.5,  # DSCR1; tumor angiogenesis
    "DYRK1A":  0.5,  # mixed TSG/oncogene evidence
    "TIAM1":   0.5,  # invasion / metastasis
    "DNMT3L":  0.5,  # DNA methyltransferase regulator
    "ERG":     0.5,  # prostate cancer fusion partner
    "TMPRSS2": 0.5,  # prostate cancer fusion partner
    "DCC":     0.5,  # colorectal
    "FBXW7":   0.5,  # broad cancer involvement
    "ARID1A":  0.5,  # SWI/SNF; ovarian, gastric, endometrial
    "KMT2D":   0.5,  # follicular lymphoma
    "KMT2C":   0.5,  # broad
    "CREBBP":  0.5,  # Rubinstein-Taybi; lymphoma
    "EP300":   0.5,  # broad
    "MUTYH":   0.5,  # MUTYH-associated polyposis
}


def builtin_cosmic_census() -> CosmicCensus:
    return CosmicCensus(
        tsg_symbols=set(_CANONICAL_TSGS),
        confidence=dict(_CANONICAL_TSGS),
    )


def n_genes() -> int:
    return len(_CANONICAL_TSGS)
