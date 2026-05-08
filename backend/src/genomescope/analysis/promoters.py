"""±flank-bp promoter construction + CpG-island overlap (pyranges-backed)."""

from __future__ import annotations

import pandas as pd
import pyranges as pr


DEFAULT_FLANK = 2_000


def build_promoters(
    genes: pd.DataFrame,
    flank: int = DEFAULT_FLANK,
) -> pd.DataFrame:
    """Per-isoform ±flank-bp window around each TSS, start-clamped at 0."""
    if genes.empty:
        return genes.copy()
    g = genes.copy()
    g["start"] = (g["tss"] - flank).clip(lower=0)
    g["end"] = g["tss"] + flank
    return g[
        ["isoform_id", "gene_id", "chrom", "start", "end", "strand", "tss"]
    ].reset_index(drop=True)


def _to_pyranges(df: pd.DataFrame, extra: list[str]) -> pr.PyRanges:
    renamed = df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    keep = ["Chromosome", "Start", "End"] + [c for c in extra if c in renamed.columns]
    return pr.PyRanges(renamed[keep].copy())


def overlap_islands_with_promoters(
    islands: pd.DataFrame,
    promoters: pd.DataFrame,
) -> pd.DataFrame:
    """One row per island × isoform overlap, with overlap_bp."""
    if islands.empty or promoters.empty:
        return pd.DataFrame(
            columns=[
                "island_id", "chrom", "island_start", "island_end",
                "isoform_id", "gene_id", "strand", "tss", "overlap_bp",
            ]
        )

    gr_islands = _to_pyranges(islands, extra=["island_id"])
    gr_promoters = _to_pyranges(
        promoters, extra=["isoform_id", "gene_id", "strand", "tss"]
    )

    joined = gr_islands.join(gr_promoters).df
    if joined.empty:
        return pd.DataFrame(
            columns=[
                "island_id", "chrom", "island_start", "island_end",
                "isoform_id", "gene_id", "strand", "tss", "overlap_bp",
            ]
        )

    ov_start = joined[["Start", "Start_b"]].max(axis=1)
    ov_end = joined[["End", "End_b"]].min(axis=1)
    joined["overlap_bp"] = (ov_end - ov_start).clip(lower=0)

    out = joined.rename(
        columns={
            "Chromosome": "chrom",
            "Start": "island_start",
            "End": "island_end",
        }
    )[
        [
            "island_id",
            "chrom",
            "island_start",
            "island_end",
            "isoform_id",
            "gene_id",
            "strand",
            "tss",
            "overlap_bp",
        ]
    ].sort_values(["island_start", "isoform_id"]).reset_index(drop=True)
    return out
