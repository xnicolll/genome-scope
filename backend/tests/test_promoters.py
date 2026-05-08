"""Tests for promoter construction + CpG island overlap (T-11)."""

from __future__ import annotations

import pandas as pd

from genomescope.analysis.promoters import (
    build_promoters,
    overlap_islands_with_promoters,
)


def _genes() -> pd.DataFrame:
    return pd.DataFrame(
        [
            # + strand: TSS = start
            {"gene_id": "ENSG_A", "isoform_id": "ENST_A1", "chrom": "chr21",
             "start": 10_000, "end": 15_000, "strand": "+", "tss": 10_000,
             "gene_type": "protein_coding", "transcript_type": "protein_coding",
             "transcript_class": "coding", "canonical": True},
            # - strand: TSS = end
            {"gene_id": "ENSG_B", "isoform_id": "ENST_B1", "chrom": "chr21",
             "start": 50_000, "end": 55_000, "strand": "-", "tss": 55_000,
             "gene_type": "lncRNA", "transcript_type": "lncRNA",
             "transcript_class": "nonCoding", "canonical": False},
            # + strand near chrom origin - must be clamped at 0
            {"gene_id": "ENSG_C", "isoform_id": "ENST_C1", "chrom": "chr21",
             "start": 500, "end": 1_500, "strand": "+", "tss": 500,
             "gene_type": "protein_coding", "transcript_type": "protein_coding",
             "transcript_class": "coding", "canonical": True},
        ]
    )


def test_build_promoters_window_size() -> None:
    p = build_promoters(_genes(), flank=2_000)
    assert len(p) == 3
    a = p.set_index("isoform_id")
    assert a.loc["ENST_A1", "start"] == 8_000
    assert a.loc["ENST_A1", "end"] == 12_000
    assert a.loc["ENST_B1", "start"] == 53_000
    assert a.loc["ENST_B1", "end"] == 57_000


def test_build_promoters_clamps_negative_start() -> None:
    p = build_promoters(_genes(), flank=2_000)
    c = p.set_index("isoform_id").loc["ENST_C1"]
    assert c["start"] == 0
    assert c["end"] == 2_500


def test_build_promoters_empty() -> None:
    empty = pd.DataFrame(columns=["gene_id", "isoform_id", "chrom", "start", "end", "strand", "tss"])
    assert build_promoters(empty).empty


def test_overlap_basic_join() -> None:
    promoters = build_promoters(_genes(), flank=2_000)
    islands = pd.DataFrame(
        [
            # overlaps ENST_A1 promoter (8k-12k)
            {"island_id": "isl_1", "chrom": "chr21", "start": 9_500, "end": 10_500},
            # overlaps ENST_B1 promoter (53k-57k)
            {"island_id": "isl_2", "chrom": "chr21", "start": 54_000, "end": 56_500},
            # no overlap with anything
            {"island_id": "isl_3", "chrom": "chr21", "start": 100_000, "end": 101_000},
        ]
    )
    out = overlap_islands_with_promoters(islands, promoters)
    assert len(out) == 2
    out = out.set_index("island_id")
    assert out.loc["isl_1", "isoform_id"] == "ENST_A1"
    assert out.loc["isl_1", "overlap_bp"] == 1_000    # [9500, 10500) ∩ [8000, 12000)
    assert out.loc["isl_2", "overlap_bp"] == 2_500    # [54000, 56500) ∩ [53000, 57000)


def test_overlap_multi_match_island_hits_multiple_isoforms() -> None:
    """One island hitting two nearby promoters must yield two rows."""
    promoters = pd.DataFrame(
        [
            {"isoform_id": "ENST_X", "gene_id": "ENSG_X", "chrom": "chr21",
             "start": 10_000, "end": 14_000, "strand": "+", "tss": 12_000},
            {"isoform_id": "ENST_Y", "gene_id": "ENSG_Y", "chrom": "chr21",
             "start": 13_000, "end": 17_000, "strand": "+", "tss": 15_000},
        ]
    )
    islands = pd.DataFrame(
        [{"island_id": "isl_1", "chrom": "chr21", "start": 13_500, "end": 13_800}]
    )
    out = overlap_islands_with_promoters(islands, promoters)
    assert len(out) == 2
    assert set(out["isoform_id"]) == {"ENST_X", "ENST_Y"}
    assert (out["overlap_bp"] == 300).all()


def test_overlap_empty_inputs_return_empty() -> None:
    empty_islands = pd.DataFrame(columns=["island_id", "chrom", "start", "end"])
    empty_promoters = pd.DataFrame(
        columns=["isoform_id", "gene_id", "chrom", "start", "end", "strand", "tss"]
    )
    assert overlap_islands_with_promoters(empty_islands, empty_promoters).empty
