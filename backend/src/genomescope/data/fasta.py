"""FASTA I/O for chromosome-scale sequences."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import SeqIO

VALID_NUCLEOTIDES = frozenset("ACGTN")
NUC_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
INT_TO_NUC = {v: k for k, v in NUC_TO_INT.items()}

_TRANS_TABLE = str.maketrans({"A": "0", "C": "1", "G": "2", "T": "3", "N": "4"})


@dataclass(frozen=True)
class Sequence:
    chrom: str
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)

    def as_int_array(self) -> np.ndarray:
        return np.frombuffer(
            self.seq.translate(_TRANS_TABLE).encode("ascii"), dtype=np.uint8
        ) - ord("0")

    def validate(self) -> None:
        if not self.seq:
            raise ValueError(f"empty sequence for {self.chrom}")
        bad = set(self.seq) - VALID_NUCLEOTIDES
        if bad:
            raise ValueError(f"unexpected characters in {self.chrom}: {sorted(bad)}")


def load_fasta(path: str | Path, uppercase: bool = True) -> Sequence:
    """Load a single-record FASTA file. Raises on multi-record or invalid input."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"no records found in {path}")
    if len(records) > 1:
        raise ValueError(f"expected single-record FASTA, got {len(records)} in {path}")

    rec = records[0]
    seq = str(rec.seq).upper() if uppercase else str(rec.seq)
    sequence = Sequence(chrom=rec.id, seq=seq)
    sequence.validate()
    return sequence


def find_cpg_sites(seq: str) -> np.ndarray:
    """Indices i where seq[i:i+2] == 'CG'."""
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    return np.flatnonzero((arr[:-1] == ord("C")) & (arr[1:] == ord("G")))
