"""Tests for the knownGene annotation loader (T-10)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from genomescope.data.annotations import CosmicCensus, load_cosmic_census, load_known_genes


def _write_fake(tmp_path: Path) -> Path:
    payload = {
        "chrom": "chr21",
        "knownGene": [
            {
                "chrom": "chr21",
                "chromStart": 1000,
                "chromEnd": 2000,
                "strand": "+",
                "name": "ENST00000001.1",
                "geneName": "ENSG0000001",
                "geneType": "protein_coding",
                "transcriptType": "protein_coding",
                "transcriptClass": "coding",
                "tag": "Ensembl_canonical,basic",
            },
            {
                "chrom": "chr21",
                "chromStart": 5000,
                "chromEnd": 6000,
                "strand": "-",
                "name": "ENST00000002.1",
                "geneName": "ENSG0000002",
                "geneType": "lncRNA",
                "transcriptType": "lncRNA",
                "transcriptClass": "nonCoding",
                "tag": "basic",
            },
            {
                "chrom": "chr20",   # wrong chrom, must be filtered
                "chromStart": 10,
                "chromEnd": 20,
                "strand": "+",
                "name": "ENST99999.1",
                "geneName": "ENSG9999",
                "tag": "",
            },
        ],
    }
    p = tmp_path / "known.json"
    p.write_text(json.dumps(payload))
    return p


def test_load_known_genes_schema(tmp_path: Path) -> None:
    df = load_known_genes(_write_fake(tmp_path), chrom="chr21")
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert list(df.columns) == [
        "gene_id",
        "isoform_id",
        "chrom",
        "start",
        "end",
        "strand",
        "tss",
        "gene_type",
        "transcript_type",
        "transcript_class",
        "canonical",
    ]


def test_tss_is_strand_aware(tmp_path: Path) -> None:
    df = load_known_genes(_write_fake(tmp_path), chrom="chr21")
    plus = df[df["strand"] == "+"].iloc[0]
    minus = df[df["strand"] == "-"].iloc[0]
    assert plus["tss"] == plus["start"]
    assert minus["tss"] == minus["end"]


def test_canonical_flag(tmp_path: Path) -> None:
    df = load_known_genes(_write_fake(tmp_path), chrom="chr21")
    assert df.loc[df["isoform_id"] == "ENST00000001.1", "canonical"].iloc[0]
    assert not df.loc[df["isoform_id"] == "ENST00000002.1", "canonical"].iloc[0]


def test_other_chromosomes_filtered(tmp_path: Path) -> None:
    df = load_known_genes(_write_fake(tmp_path), chrom="chr21")
    assert (df["chrom"] == "chr21").all()


def test_load_real_chr21_file() -> None:
    """If the real downloaded JSON exists, smoke-test it."""
    real = Path(__file__).resolve().parents[1] / "data" / "raw" / "knownGene_chr21.json"
    if not real.exists():
        return
    df = load_known_genes(real, chrom="chr21")
    assert len(df) > 1000
    assert df["gene_id"].str.startswith("ENSG").any()
    assert df["isoform_id"].str.startswith("ENST").any()
    assert df["canonical"].any()


def test_load_cosmic_missing_returns_builtin_fallback(tmp_path: Path) -> None:
    """When the CSV is absent, default behaviour returns the built-in TSG list."""
    cen = load_cosmic_census(tmp_path / "nope.csv")
    assert cen is not None
    assert len(cen) > 20                           # built-in has ~40 entries
    assert "TP53" in cen.tsg_symbols
    assert "RUNX1" in cen.tsg_symbols              # chr21-relevant TSG
    assert cen.confidence["TP53"] == 1.0


def test_load_cosmic_missing_strict_returns_none(tmp_path: Path) -> None:
    """fallback_to_builtin=False reverts to the original strict behaviour."""
    assert load_cosmic_census(tmp_path / "nope.csv", fallback_to_builtin=False) is None


def test_load_cosmic_parses_tsg_rows(tmp_path: Path) -> None:
    csv = tmp_path / "census.csv"
    csv.write_text(
        "Gene Symbol,Role in Cancer,Tier\n"
        "TP53,TSG,1\n"
        "BRCA1,\"oncogene, TSG\",1\n"
        "KRAS,oncogene,1\n"
        "FOO,TSG,2\n"
    )
    cen = load_cosmic_census(csv)
    assert isinstance(cen, CosmicCensus)
    assert "TP53" in cen.tsg_symbols
    assert "BRCA1" in cen.tsg_symbols
    assert "KRAS" not in cen.tsg_symbols
    assert cen.confidence["TP53"] == 1.0
    assert cen.confidence["FOO"] == 0.5


def test_load_cosmic_handles_legacy_column_name(tmp_path: Path) -> None:
    csv = tmp_path / "legacy.csv"
    csv.write_text("gene,role_in_cancer,tier\nTP53,TSG,1\n")
    cen = load_cosmic_census(csv)
    assert cen is not None
    assert "TP53" in cen.tsg_symbols
