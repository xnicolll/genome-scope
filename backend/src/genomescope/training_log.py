"""Append-only CSV log of Baum-Welch training runs."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LOG = ROOT / "data" / "logs" / "training-history.csv"

COLUMNS = [
    "timestamp",
    "model_name",
    "model_type",
    "cohort",
    "window_start",
    "window_end",
    "n_samples",
    "iterations",
    "converged",
    "ll_start",
    "ll_end",
    "ll_delta",
    "state_means",
    "duration_s",
    "note",
]


@dataclass
class TrainingLogEntry:
    model_name: str
    model_type: str
    cohort: str | None
    window_start: int | None
    window_end: int | None
    n_samples: int | None
    iterations: int
    converged: bool
    ll_start: float
    ll_end: float
    state_means: tuple[float, ...]
    duration_s: float
    note: str = ""


def append(entry: TrainingLogEntry, path: Path | None = None) -> Path:
    p = Path(path) if path is not None else DEFAULT_LOG
    p.parent.mkdir(parents=True, exist_ok=True)
    new_file = not p.exists() or p.stat().st_size == 0

    row = {
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "model_name": entry.model_name,
        "model_type": entry.model_type,
        "cohort": entry.cohort or "",
        "window_start": entry.window_start if entry.window_start is not None else "",
        "window_end": entry.window_end if entry.window_end is not None else "",
        "n_samples": entry.n_samples if entry.n_samples is not None else "",
        "iterations": entry.iterations,
        "converged": "yes" if entry.converged else "no",
        "ll_start": f"{entry.ll_start:.4f}",
        "ll_end": f"{entry.ll_end:.4f}",
        "ll_delta": f"{entry.ll_end - entry.ll_start:+.4f}",
        "state_means": ";".join(f"{m:.3f}" for m in entry.state_means),
        "duration_s": f"{entry.duration_s:.2f}",
        "note": entry.note,
    }

    with open(p, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        if new_file:
            w.writeheader()
        w.writerow(row)

    return p


def read_all(path: Path | None = None) -> list[dict]:
    p = Path(path) if path is not None else DEFAULT_LOG
    if not p.exists():
        return []
    with open(p, newline="") as f:
        return list(csv.DictReader(f))
