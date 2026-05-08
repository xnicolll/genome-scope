"""Tests for the isoform-aware methylation report (T-12, T-13, T-14)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from genomescope.analysis.promoters import build_promoters
from genomescope.analysis.report import (
    HYPERMETH_THRESHOLD,
    ReportInputs,
    annotate_islands,
    build_report,
)
from genomescope.data.tcga import MethylationTrack


def _islands() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"island_id": "isl_hyper", "chrom": "chr21", "start": 9_500, "end": 10_500},
            {"island_id": "isl_low",   "chrom": "chr21", "start": 54_000, "end": 55_000},
            {"island_id": "isl_nohit", "chrom": "chr21", "start": 100_000, "end": 100_500},
        ]
    )


def _genes() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"gene_id": "ENSG_APP", "isoform_id": "ENST_APP_1", "chrom": "chr21",
             "start": 10_000, "end": 15_000, "strand": "+", "tss": 10_000,
             "gene_type": "protein_coding", "transcript_type": "protein_coding",
             "transcript_class": "coding", "canonical": True},
            {"gene_id": "ENSG_SOD1", "isoform_id": "ENST_SOD1_1", "chrom": "chr21",
             "start": 50_000, "end": 55_000, "strand": "-", "tss": 55_000,
             "gene_type": "protein_coding", "transcript_type": "protein_coding",
             "transcript_class": "coding", "canonical": True},
        ]
    )


def _brca_track() -> MethylationTrack:
    # isl_hyper (9500-10500) has probes with mean ~0.8
    # isl_low (54000-55000) has probes with mean ~0.1
    # isl_nohit has no probes
    return MethylationTrack(
        chrom="chr21",
        positions=np.array([9_600, 9_800, 10_200, 10_400,
                            54_100, 54_500, 54_900], dtype=np.int64),
        betas=np.array([0.75, 0.82, 0.85, 0.78,
                        0.05, 0.10, 0.15], dtype=np.float64),
    )


def test_annotate_islands_per_cohort_mean_beta() -> None:
    out = annotate_islands(
        _islands(), {"brca": _brca_track()}
    )
    row = out.set_index("island_id")
    assert abs(row.loc["isl_hyper", "mean_beta_brca"] - 0.80) < 0.01
    assert abs(row.loc["isl_low", "mean_beta_brca"] - 0.10) < 0.01
    assert pd.isna(row.loc["isl_nohit", "mean_beta_brca"])


def test_build_report_joins_everything() -> None:
    promoters = build_promoters(_genes(), flank=2_000)
    inputs = ReportInputs(
        islands=_islands(),
        promoters=promoters,
        cohort_tracks={"brca": _brca_track()},
        gene_symbols={"ENSG_APP": "APP", "ENSG_SOD1": "SOD1"},
        tsg_set={"SOD1"},              # pretend SOD1 is a TSG
        cosmic_tier={"SOD1": 0.9, "APP": 0.2},
    )
    report = build_report(inputs)

    # isl_nohit has no promoter overlap → filtered out
    assert set(report["island_id"]) == {"isl_hyper", "isl_low"}

    hyper = report[report["island_id"] == "isl_hyper"].iloc[0]
    low = report[report["island_id"] == "isl_low"].iloc[0]

    assert hyper["gene_symbol"] == "APP"
    assert hyper["isoform_id"] == "ENST_APP_1"
    assert hyper["hypermethylated"]
    assert not hyper["is_tsg"]

    assert low["gene_symbol"] == "SOD1"
    assert low["is_tsg"]
    assert not low["hypermethylated"]


def test_hypermethylated_threshold() -> None:
    promoters = build_promoters(_genes(), flank=2_000)
    inputs = ReportInputs(
        islands=_islands(),
        promoters=promoters,
        cohort_tracks={"brca": _brca_track()},
    )
    report = build_report(inputs)
    for _, row in report.iterrows():
        if not pd.isna(row["max_mean_beta"]):
            assert row["hypermethylated"] == (row["max_mean_beta"] > HYPERMETH_THRESHOLD)


def test_priority_is_beta_times_confidence() -> None:
    promoters = build_promoters(_genes(), flank=2_000)
    inputs = ReportInputs(
        islands=_islands(),
        promoters=promoters,
        cohort_tracks={"brca": _brca_track()},
        gene_symbols={"ENSG_APP": "APP", "ENSG_SOD1": "SOD1"},
        tsg_set={"APP", "SOD1"},
        cosmic_tier={"APP": 0.5, "SOD1": 1.0},
    )
    report = build_report(inputs)
    row = report[report["island_id"] == "isl_hyper"].iloc[0]
    # APP priority = 0.80 (mean) × 0.5 (confidence) ≈ 0.40
    assert abs(row["priority"] - 0.4) < 0.02


def test_multi_cohort_picks_max() -> None:
    """Two cohorts: one hypermethylated, one not. Report should report
    the max and name the winning cohort."""
    islands = pd.DataFrame([
        {"island_id": "isl_1", "chrom": "chr21", "start": 9_500, "end": 10_500},
    ])
    brca = MethylationTrack(
        chrom="chr21",
        positions=np.array([9_600, 10_200]),
        betas=np.array([0.8, 0.9]),
    )
    luad = MethylationTrack(
        chrom="chr21",
        positions=np.array([9_600, 10_200]),
        betas=np.array([0.1, 0.2]),
    )
    promoters = build_promoters(_genes(), flank=2_000)
    report = build_report(
        ReportInputs(
            islands=islands,
            promoters=promoters,
            cohort_tracks={"brca": brca, "luad": luad},
        )
    )
    row = report.iloc[0]
    assert row["max_cohort"] == "brca"
    assert row["max_mean_beta"] > 0.8
    assert row["hypermethylated"]


def test_normal_cohort_drives_delta_vs_normal() -> None:
    """A cohort whose name starts with 'normal' is treated as a healthy
    baseline. delta_vs_normal = max(tumor_betas) - normal_beta per island."""
    islands = pd.DataFrame([
        # cancer-specific: healthy in normal, silenced in tumor
        {"island_id": "isl_cs",  "chrom": "chr21", "start": 9_500,  "end": 10_500},
        # constitutively methylated: methylated in both
        {"island_id": "isl_const", "chrom": "chr21", "start": 9_500, "end": 10_500},
    ])
    promoters = build_promoters(_genes(), flank=2_000)

    # Force both islands to have known per-cohort means by giving each
    # one its own probe stack
    brca = MethylationTrack(
        chrom="chr21",
        positions=np.array([9_700, 9_800], dtype=np.int64),
        betas=np.array([0.85, 0.80]),  # tumor: highly methylated
    )
    normal = MethylationTrack(
        chrom="chr21",
        positions=np.array([9_700, 9_800], dtype=np.int64),
        betas=np.array([0.05, 0.10]),  # normal: unmethylated
    )

    inputs = ReportInputs(
        islands=islands.iloc[:1],     # just isl_cs for this assertion
        promoters=promoters,
        cohort_tracks={"brca": brca, "normal_brca": normal},
    )
    report = build_report(inputs)
    assert not report.empty
    row = report.iloc[0]
    # delta = max_tumor (0.825) - normal_mean (0.075) ≈ +0.75
    assert row["delta_vs_normal"] > 0.6
    assert bool(row["cancer_specific"]) is True
    # max_cohort excludes normal - so it should be 'brca', not 'normal_brca'
    assert row["max_cohort"] == "brca"
    assert "mean_beta_normal" in report.columns


def test_normal_cohort_excluded_from_hypermethylated_flag() -> None:
    """If only the normal cohort has high methylation, an island shouldn't
    be flagged hypermethylated."""
    islands = pd.DataFrame([
        {"island_id": "isl_x", "chrom": "chr21", "start": 9_500, "end": 10_500},
    ])
    promoters = build_promoters(_genes(), flank=2_000)
    brca = MethylationTrack(
        chrom="chr21",
        positions=np.array([9_700]),
        betas=np.array([0.20]),  # tumor: low
    )
    normal = MethylationTrack(
        chrom="chr21",
        positions=np.array([9_700]),
        betas=np.array([0.85]),  # normal: high
    )
    report = build_report(ReportInputs(
        islands=islands, promoters=promoters,
        cohort_tracks={"brca": brca, "normal_brca": normal},
    ))
    row = report.iloc[0]
    assert not bool(row["hypermethylated"])
    assert row["max_cohort"] == "brca"


def test_empty_inputs_return_empty_report() -> None:
    report = build_report(
        ReportInputs(
            islands=pd.DataFrame(columns=["island_id", "chrom", "start", "end"]),
            promoters=pd.DataFrame(
                columns=["isoform_id", "gene_id", "chrom", "start", "end", "strand", "tss"]
            ),
            cohort_tracks={},
        )
    )
    assert report.empty
