"""Tests for the TCGA loader + synthetic methylation generator (T-06)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

import pandas as pd
import pytest

from genomescope.data.fasta import Sequence
from genomescope.data.tcga import (
    MethylationTrack,
    _sniff_format,
    load_tcga_cohorts,
    load_tcga_directory,
    synthetic_track_from_fasta,
)


def test_methylation_track_window() -> None:
    track = MethylationTrack(
        chrom="chr21",
        positions=np.array([10, 20, 30, 40, 50], dtype=np.int64),
        betas=np.array([0.1, 0.2, 0.3, 0.4, 0.5]),
    )
    sub = track.window(15, 45)
    assert sub.positions.tolist() == [20, 30, 40]
    assert sub.betas.tolist() == [0.2, 0.3, 0.4]


def test_synthetic_track_produces_sorted_positions() -> None:
    seq = Sequence(chrom="chr21", seq="A" * 5 + "CG" + "T" * 10 + "CG" + "N" * 5)
    track = synthetic_track_from_fasta(seq, islands=[(0, 10)], seed=0)
    assert len(track) == 2
    assert track.positions.tolist() == [5, 17]
    assert np.all(np.diff(track.positions) >= 0)


def test_synthetic_track_respects_island_regions() -> None:
    # sparse CGs outside + a dense CG cluster inside the "island" interval
    left = ("ATCGAT" * 50)            # one CG per 6 bp = ~50 CpGs
    core = ("CGCG" * 40)              # dense cluster = ~80 CpGs
    right = ("ATCGAT" * 50)           # ~50 CpGs
    seq_str = left + core + right
    islands = [(len(left), len(left) + len(core))]
    seq = Sequence(chrom="chr21", seq=seq_str)
    track = synthetic_track_from_fasta(seq, islands=islands, seed=1)
    in_island = (track.positions >= islands[0][0]) & (track.positions < islands[0][1])
    assert in_island.any()
    assert (~in_island).any()
    assert track.betas[in_island].mean() > 0.7
    assert track.betas[~in_island].mean() < 0.3


def test_load_tcga_directory_missing_returns_none(tmp_path: Path) -> None:
    assert load_tcga_directory(tmp_path / "nope", chrom="chr21") is None


def test_load_tcga_directory_averages_replicates(tmp_path: Path) -> None:
    d = tmp_path / "tcga"
    d.mkdir()
    (d / "a.tsv").write_text(
        "chr21\t100\t102\t0.10\nchr21\t200\t202\t0.40\n"
    )
    (d / "b.tsv").write_text(
        "chr21\t100\t102\t0.30\nchr21\t200\t202\t0.60\n"
    )
    track = load_tcga_directory(d, chrom="chr21")
    assert track is not None
    assert track.positions.tolist() == [100, 200]
    assert np.allclose(track.betas, [0.20, 0.50])


def test_load_tcga_directory_filters_other_chromosomes(tmp_path: Path) -> None:
    d = tmp_path / "tcga"
    d.mkdir()
    (d / "a.tsv").write_text(
        "chr21\t100\t102\t0.5\nchr1\t100\t102\t0.9\n"
    )
    track = load_tcga_directory(d, chrom="chr21")
    assert track is not None
    assert track.positions.tolist() == [100]
    assert track.betas.tolist() == [0.5]


def test_sniff_format_detects_sesame(tmp_path: Path) -> None:
    p = tmp_path / "sample.txt"
    p.write_text("cg00000029\t0.71\ncg00000108\t0.05\n")
    assert _sniff_format(p) == "sesame"


def test_sniff_format_detects_bed(tmp_path: Path) -> None:
    p = tmp_path / "sample.bed"
    p.write_text("chr21\t100\t102\t0.5\n")
    assert _sniff_format(p) == "bed"


def test_iter_beta_files_ignores_sidecar_metadata(tmp_path: Path) -> None:
    """gdc-client drops annotations.txt + logs/*.parcel alongside downloads,
    and our own build step writes subset.manifest.txt - skip them all."""
    from genomescope.data.tcga import _iter_beta_files

    d = tmp_path / "tcga"
    (d / "uuid1").mkdir(parents=True)
    (d / "uuid1" / "sample.methylation_array.sesame.level3betas.txt").write_text(
        "cg00000029\t0.5\n"
    )
    (d / "uuid1" / "logs").mkdir()
    (d / "uuid1" / "logs" / "sample.txt.parcel").write_text("binary garbage")
    (d / "uuid2").mkdir()
    (d / "uuid2" / "annotations.txt").write_text("notes\tother\nfoo\tbar\n")
    (d / "subset.manifest.txt").write_text("id\tfilename\tmd5\tsize\tstate\nfoo\n")

    files = _iter_beta_files(d)
    assert len(files) == 1
    assert files[0].name.endswith("level3betas.txt")


def test_load_tcga_cohorts_discovers_subdirectories(tmp_path: Path) -> None:
    root = tmp_path / "tcga"
    (root / "brca_hm450").mkdir(parents=True)
    (root / "luad_hm450").mkdir(parents=True)
    (root / "coad_hm450").mkdir(parents=True)
    (root / "brca_hm450" / "a.bed").write_text("chr21\t100\t102\t0.50\n")
    (root / "luad_hm450" / "b.bed").write_text("chr21\t200\t202\t0.70\n")
    # coad has no files - should be skipped, not raise

    cohorts = load_tcga_cohorts(root, chrom="chr21")
    assert set(cohorts) == {"brca", "luad"}
    assert cohorts["brca"].positions.tolist() == [100]
    assert cohorts["luad"].positions.tolist() == [200]


def test_load_sesame_level3_maps_probes(tmp_path: Path, monkeypatch) -> None:
    """Two SESAMe-format files averaged via an in-memory HM450 probe map."""
    # tiny probe map: 3 chr21 probes, 1 chr1 probe (should be filtered out)
    probe_map = pd.DataFrame(
        {
            "chrom": ["chr21", "chr21", "chr21", "chr1"],
            "start": [1000, 2000, 3000, 500],
            "end":   [1002, 2002, 3002, 502],
        },
        index=pd.Index(["cg_a", "cg_b", "cg_c", "cg_other"], name="probe_id"),
    )
    monkeypatch.setattr(
        "genomescope.data.tcga.load_hm450_probe_map", lambda chrom="chr21": probe_map
    )

    d = tmp_path / "tcga"
    d.mkdir()
    (d / "s1.txt").write_text("cg_a\t0.10\ncg_b\t0.50\ncg_c\t0.90\ncg_other\t0.99\n")
    (d / "s2.txt").write_text("cg_a\t0.30\ncg_b\t0.70\ncg_c\tNA\n")

    track = load_tcga_directory(d, chrom="chr21")
    assert track is not None
    assert track.positions.tolist() == [1000, 2000, 3000]
    # cg_a: (0.10+0.30)/2 = 0.20 ; cg_b: 0.60 ; cg_c: 0.90 (NA dropped)
    assert track.betas.tolist() == pytest.approx([0.20, 0.60, 0.90])
