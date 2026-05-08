"""Tests for the evaluation helpers."""

from __future__ import annotations

import json
from pathlib import Path

from genomescope.analysis.evaluate import (
    load_ucsc_track_json,
    position_level_metrics,
    segments_to_mask,
)


def test_segments_to_mask_basic() -> None:
    mask = segments_to_mask([(0, 3), (5, 7)], length=10)
    assert mask.tolist() == [True, True, True, False, False, True, True, False, False, False]


def test_perfect_overlap_gives_f1_one() -> None:
    pred = [(0, 5), (10, 15)]
    truth = [(0, 5), (10, 15)]
    m = position_level_metrics(pred, truth, length=20)
    assert m.f1 == 1.0
    assert m.precision == 1.0
    assert m.recall == 1.0


def test_no_overlap_gives_f1_zero() -> None:
    pred = [(0, 5)]
    truth = [(10, 15)]
    m = position_level_metrics(pred, truth, length=20)
    assert m.f1 == 0.0


def test_partial_overlap_f1() -> None:
    pred = [(0, 10)]     # 10 bp
    truth = [(5, 15)]    # 10 bp, 5 bp overlap
    m = position_level_metrics(pred, truth, length=20)
    # tp=5 fp=5 fn=5 → p=0.5 r=0.5 f1=0.5
    assert m.precision == 0.5
    assert m.recall == 0.5
    assert m.f1 == 0.5


def test_load_ucsc_track_json_shape(tmp_path: Path) -> None:
    payload = {
        "downloadTime": "2026",
        "genome": "hg38",
        "cpgIslandExt": [
            {"chrom": "chr21", "chromStart": 100, "chromEnd": 500},
            {"chrom": "chr21", "chromStart": 1000, "chromEnd": 1200},
            {"chrom": "chr20", "chromStart": 50, "chromEnd": 60},  # wrong chrom
        ],
    }
    p = tmp_path / "track.json"
    p.write_text(json.dumps(payload))
    islands = load_ucsc_track_json(p, chrom="chr21")
    assert islands == [(100, 500), (1000, 1200)]
