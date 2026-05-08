"""UCSC knownGene + COSMIC Cancer Gene Census loaders."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


def load_known_genes(
    path: Path | str,
    chrom: str = "chr21",
) -> pd.DataFrame:
    """Per-isoform DataFrame from the UCSC knownGene JSON, filtered to `chrom`."""
    raw = json.loads(Path(path).read_text())
    entries = raw.get("knownGene") or []

    rows: list[dict] = []
    for e in entries:
        if e.get("chrom") != chrom:
            continue
        strand = e.get("strand", "+")
        start = int(e["chromStart"])
        end = int(e["chromEnd"])
        tss = start if strand == "+" else end
        tags = e.get("tag", "") or ""
        rows.append(
            {
                "gene_id": e.get("geneName", ""),
                "isoform_id": e.get("name", ""),
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "tss": tss,
                "gene_type": e.get("geneType", ""),
                "transcript_type": e.get("transcriptType", ""),
                "transcript_class": e.get("transcriptClass", ""),
                "canonical": "Ensembl_canonical" in tags,
            }
        )

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["start", "isoform_id"]).reset_index(drop=True)


@dataclass
class CosmicCensus:
    tsg_symbols: set[str]
    confidence: dict[str, float]

    def __len__(self) -> int:
        return len(self.tsg_symbols)


def load_cosmic_census(
    path: Path | str,
    fallback_to_builtin: bool = True,
) -> CosmicCensus | None:
    """Parse the COSMIC Cancer Gene Census CSV; falls back to the built-in TSG list."""
    p = Path(path)
    if not p.exists():
        if fallback_to_builtin:
            from .builtin_tsg import builtin_cosmic_census
            return builtin_cosmic_census()
        return None

    df = pd.read_csv(p)
    cols = {c.lower().strip().replace(" ", "_"): c for c in df.columns}
    sym_col = cols.get("gene_symbol") or cols.get("gene") or cols.get("name")
    role_col = cols.get("role_in_cancer") or cols.get("role")
    tier_col = cols.get("tier")
    if sym_col is None or role_col is None:
        return None

    role = df[role_col].astype(str).str.lower()
    is_tsg = role.str.contains("tsg") | role.str.contains("tumour suppressor") \
        | role.str.contains("tumor suppressor")
    tsg_df = df[is_tsg].copy()
    symbols = set(tsg_df[sym_col].astype(str).str.strip())

    confidence: dict[str, float] = {}
    if tier_col is not None:
        tier_num = pd.to_numeric(tsg_df[tier_col], errors="coerce").fillna(2)
        for sym, t in zip(tsg_df[sym_col].astype(str).str.strip(), tier_num, strict=True):
            confidence[sym] = 1.0 if t <= 1 else 0.5
    else:
        confidence = {s: 1.0 for s in symbols}

    return CosmicCensus(tsg_symbols=symbols, confidence=confidence)
