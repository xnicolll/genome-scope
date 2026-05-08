"""Tests for HMM checkpoint save/load + training-history registry."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from genomescope.hmm.beta_model import beta_cpg_hmm
from genomescope.hmm.checkpoint import (
    Checkpoint,
    TrainingRun,
    load_checkpoint,
    read_registry,
    save_checkpoint,
    upsert_registry,
    utc_now,
)
from genomescope.hmm.model import HMM, standard_cpg_hmm


def test_checkpoint_roundtrip_standard_hmm(tmp_path: Path) -> None:
    original = standard_cpg_hmm()
    path = tmp_path / "std.json"
    save_checkpoint(original, path)

    loaded, history = load_checkpoint(path)
    assert isinstance(loaded, HMM)
    assert loaded.n_states == original.n_states
    assert loaded.n_obs == original.n_obs
    assert np.allclose(loaded.log_init, original.log_init, equal_nan=True)
    # -inf positions survive the JSON round-trip (Python JSON keeps -Infinity)
    np.testing.assert_array_equal(loaded.log_emit, original.log_emit)
    assert loaded.state_labels == original.state_labels
    assert history == []


def test_checkpoint_roundtrip_beta_hmm(tmp_path: Path) -> None:
    original = beta_cpg_hmm()
    path = tmp_path / "beta.json"
    save_checkpoint(original, path)
    loaded, _ = load_checkpoint(path)
    assert loaded.n_states == 2
    assert np.allclose(loaded.alphas, original.alphas)
    assert np.allclose(loaded.betas, original.betas)
    assert loaded.state_labels == original.state_labels


def test_checkpoint_preserves_history(tmp_path: Path) -> None:
    original = beta_cpg_hmm()
    history = [
        TrainingRun(
            timestamp=utc_now(),
            cohort="brca",
            n_samples=3378,
            window_start=13_500_000,
            window_end=15_500_000,
            iterations=5,
            converged=True,
            log_likelihood_start=23541.7,
            log_likelihood_end=23627.1,
            max_iter=15,
            tol=1e-4,
        )
    ]
    path = tmp_path / "beta.json"
    save_checkpoint(original, path, history=history)

    _, loaded_history = load_checkpoint(path)
    assert len(loaded_history) == 1
    assert loaded_history[0].cohort == "brca"
    assert loaded_history[0].log_likelihood_end == 23627.1


def test_checkpoint_rejects_schema_mismatch(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text('{"schema": 999, "model_type": "standard"}')
    try:
        load_checkpoint(path)
    except ValueError as e:
        assert "schema" in str(e)
    else:
        raise AssertionError("expected ValueError for schema mismatch")


def test_registry_insert_and_update(tmp_path: Path) -> None:
    models_dir = tmp_path / "models"
    ckpt_path = models_dir / "brca.json"
    original = beta_cpg_hmm()
    history_v1 = [
        TrainingRun(
            timestamp=utc_now(), cohort="brca", n_samples=3378,
            window_start=0, window_end=100, iterations=5, converged=True,
            log_likelihood_start=100.0, log_likelihood_end=120.0,
            max_iter=15, tol=1e-4,
        )
    ]
    save_checkpoint(original, ckpt_path, history=history_v1)
    upsert_registry(models_dir, "beta-brca", ckpt_path, history_v1)

    entries = read_registry(models_dir)
    assert len(entries) == 1
    assert entries[0].name == "beta-brca"
    assert entries[0].cohorts_seen == ["brca"]
    assert entries[0].total_samples_seen == 3378
    assert entries[0].last_log_likelihood == 120.0

    # Second run against luad - same checkpoint name, extended history
    history_v2 = history_v1 + [
        TrainingRun(
            timestamp=utc_now(), cohort="luad", n_samples=2000,
            window_start=0, window_end=100, iterations=3, converged=True,
            log_likelihood_start=120.0, log_likelihood_end=135.5,
            max_iter=15, tol=1e-4,
        )
    ]
    save_checkpoint(original, ckpt_path, history=history_v2)
    upsert_registry(models_dir, "beta-brca", ckpt_path, history_v2)

    entries = read_registry(models_dir)
    assert len(entries) == 1              # upsert, not duplicate
    assert entries[0].cohorts_seen == ["brca", "luad"]
    assert entries[0].total_samples_seen == 5378
    assert entries[0].last_log_likelihood == 135.5
    assert entries[0].total_runs == 2


def test_registry_handles_multiple_named_models(tmp_path: Path) -> None:
    models_dir = tmp_path / "models"
    for name in ["beta-brca", "beta-luad", "standard-hg38"]:
        p = models_dir / f"{name}.json"
        save_checkpoint(beta_cpg_hmm(), p)
        upsert_registry(models_dir, name, p, [])
    entries = read_registry(models_dir)
    assert {e.name for e in entries} == {"beta-brca", "beta-luad", "standard-hg38"}
