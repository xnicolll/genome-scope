"""Position-level benchmarking against a ground-truth CpG island set."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Metrics:
    precision: float
    recall: float
    f1: float
    tp: int
    fp: int
    fn: int

    def as_dict(self) -> dict[str, float]:
        return {
            "precision": self.precision,
            "recall": self.recall,
            "f1": self.f1,
            "tp": self.tp,
            "fp": self.fp,
            "fn": self.fn,
        }


def segments_to_mask(segments: list[tuple[int, int]], length: int) -> np.ndarray:
    mask = np.zeros(length, dtype=bool)
    for s, e in segments:
        mask[max(0, s) : min(length, e)] = True
    return mask


def position_level_metrics(
    pred: list[tuple[int, int]],
    truth: list[tuple[int, int]],
    length: int,
) -> Metrics:
    p = segments_to_mask(pred, length)
    t = segments_to_mask(truth, length)
    tp = int((p & t).sum())
    fp = int((p & ~t).sum())
    fn = int((~p & t).sum())
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return Metrics(precision=precision, recall=recall, f1=f1, tp=tp, fp=fp, fn=fn)


def load_ucsc_track_json(path: Path, chrom: str = "chr21") -> list[tuple[int, int]]:
    """Parse a UCSC getData/track JSON response into (start, end) intervals."""
    raw = json.loads(Path(path).read_text())
    entries: list[dict] = []
    for v in raw.values():
        if isinstance(v, list) and v and isinstance(v[0], dict):
            entries = v
            break
        if isinstance(v, dict) and chrom in v:
            inner = v[chrom]
            if isinstance(inner, list):
                entries = inner
                break
    return [
        (int(e["chromStart"]), int(e["chromEnd"]))
        for e in entries
        if e.get("chrom", chrom) == chrom
    ]
