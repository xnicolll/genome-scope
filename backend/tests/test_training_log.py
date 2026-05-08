"""Tests for the training-history CSV log."""

from __future__ import annotations

from pathlib import Path

from genomescope.training_log import TrainingLogEntry, append, read_all


def _sample(ll_end: float = 45.0) -> TrainingLogEntry:
    return TrainingLogEntry(
        model_name="beta-brca",
        model_type="beta",
        cohort="brca",
        window_start=13_500_000,
        window_end=15_500_000,
        n_samples=3378,
        iterations=5,
        converged=True,
        ll_start=17.1,
        ll_end=ll_end,
        state_means=(0.052, 0.827),
        duration_s=0.12,
        note="",
    )


def test_append_writes_header_on_first_call(tmp_path: Path) -> None:
    path = tmp_path / "history.csv"
    append(_sample(), path)
    lines = path.read_text().splitlines()
    assert lines[0].startswith("timestamp,model_name,model_type")
    assert "beta-brca" in lines[1]


def test_append_is_append_only(tmp_path: Path) -> None:
    path = tmp_path / "history.csv"
    append(_sample(ll_end=40.0), path)
    append(_sample(ll_end=45.0), path)
    append(_sample(ll_end=48.0), path)
    rows = read_all(path)
    assert len(rows) == 3
    assert [float(r["ll_end"]) for r in rows] == [40.0, 45.0, 48.0]


def test_ll_delta_shows_direction(tmp_path: Path) -> None:
    path = tmp_path / "history.csv"
    append(_sample(ll_end=45.0), path)
    rows = read_all(path)
    # sample starts at 17.1 → ll_delta should be positive
    assert rows[0]["ll_delta"].startswith("+")


def test_read_all_on_missing_file(tmp_path: Path) -> None:
    assert read_all(tmp_path / "nope.csv") == []


def test_state_means_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "history.csv"
    append(_sample(), path)
    rows = read_all(path)
    parts = rows[0]["state_means"].split(";")
    assert [float(p) for p in parts] == [0.052, 0.827]
