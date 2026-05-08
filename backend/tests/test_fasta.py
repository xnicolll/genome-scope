"""Tests for genomescope.data.fasta."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from genomescope.data.fasta import (
    NUC_TO_INT,
    Sequence,
    find_cpg_sites,
    load_fasta,
)


def _write_fasta(tmp: Path, name: str, seq: str) -> Path:
    p = tmp / f"{name}.fa"
    p.write_text(f">{name}\n{seq}\n")
    return p


def test_load_fasta_roundtrip(tmp_path: Path) -> None:
    f = _write_fasta(tmp_path, "chr21", "ACGTACGTNN")
    s = load_fasta(f)
    assert s.chrom == "chr21"
    assert s.seq == "ACGTACGTNN"
    assert s.length == 10


def test_load_fasta_uppercases(tmp_path: Path) -> None:
    f = _write_fasta(tmp_path, "chr21", "acgtAcgt")
    assert load_fasta(f).seq == "ACGTACGT"


def test_validate_rejects_invalid_chars(tmp_path: Path) -> None:
    f = _write_fasta(tmp_path, "chr21", "ACGTX")
    with pytest.raises(ValueError, match="unexpected characters"):
        load_fasta(f)


def test_validate_rejects_empty(tmp_path: Path) -> None:
    f = _write_fasta(tmp_path, "chr21", "")
    with pytest.raises(ValueError):
        load_fasta(f)


def test_load_fasta_rejects_multi_record(tmp_path: Path) -> None:
    p = tmp_path / "multi.fa"
    p.write_text(">a\nACGT\n>b\nACGT\n")
    with pytest.raises(ValueError, match="single-record"):
        load_fasta(p)


def test_as_int_array_encoding() -> None:
    s = Sequence(chrom="t", seq="ACGTN")
    arr = s.as_int_array()
    assert arr.dtype == np.uint8
    assert list(arr) == [NUC_TO_INT[c] for c in "ACGTN"]


def test_find_cpg_sites() -> None:
    # CGs at positions 0, 4, 7
    seq = "CGATCGACG"
    assert list(find_cpg_sites(seq)) == [0, 4, 7]


def test_find_cpg_sites_empty() -> None:
    assert list(find_cpg_sites("AAAA")) == []
