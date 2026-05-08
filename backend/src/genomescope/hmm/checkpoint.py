"""HMM checkpoint serialisation + training-history registry."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .beta_model import BetaHMM
from .model import HMM

CHECKPOINT_SCHEMA = 1


@dataclass
class TrainingRun:
    timestamp: str
    cohort: str | None
    n_samples: int | None
    window_start: int | None
    window_end: int | None
    iterations: int
    converged: bool
    log_likelihood_start: float
    log_likelihood_end: float
    max_iter: int
    tol: float


@dataclass
class Checkpoint:
    schema: int
    model_type: str
    state_labels: list[str]
    log_init: list[float]
    log_trans: list[list[float]]
    log_emit: list[list[float]] | None = None
    n_obs: int | None = None
    alphas: list[float] | None = None
    betas: list[float] | None = None
    history: list[TrainingRun] = field(default_factory=list)

    @classmethod
    def from_hmm(cls, model: HMM | BetaHMM) -> "Checkpoint":
        if isinstance(model, HMM):
            return cls(
                schema=CHECKPOINT_SCHEMA,
                model_type="standard",
                state_labels=list(model.state_labels),
                log_init=model.log_init.tolist(),
                log_trans=model.log_trans.tolist(),
                log_emit=model.log_emit.tolist(),
                n_obs=int(model.n_obs),
            )
        if isinstance(model, BetaHMM):
            return cls(
                schema=CHECKPOINT_SCHEMA,
                model_type="beta",
                state_labels=list(model.state_labels),
                log_init=model.log_init.tolist(),
                log_trans=model.log_trans.tolist(),
                alphas=model.alphas.tolist(),
                betas=model.betas.tolist(),
            )
        raise TypeError(f"unsupported model type: {type(model).__name__}")

    def to_hmm(self) -> HMM | BetaHMM:
        log_init = np.array(self.log_init, dtype=np.float64)
        log_trans = np.array(self.log_trans, dtype=np.float64)
        if self.model_type == "standard":
            if self.log_emit is None or self.n_obs is None:
                raise ValueError("standard checkpoint missing log_emit/n_obs")
            return HMM(
                n_states=log_init.shape[0],
                n_obs=int(self.n_obs),
                log_init=log_init,
                log_trans=log_trans,
                log_emit=np.array(self.log_emit, dtype=np.float64),
                state_labels=list(self.state_labels),
            )
        if self.model_type == "beta":
            if self.alphas is None or self.betas is None:
                raise ValueError("beta checkpoint missing alphas/betas")
            return BetaHMM(
                log_init=log_init,
                log_trans=log_trans,
                alphas=np.array(self.alphas, dtype=np.float64),
                betas=np.array(self.betas, dtype=np.float64),
                state_labels=list(self.state_labels),
            )
        raise ValueError(f"unknown model_type: {self.model_type}")


def save_checkpoint(
    model: HMM | BetaHMM,
    path: Path | str,
    history: list[TrainingRun] | None = None,
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ckpt = Checkpoint.from_hmm(model)
    if history:
        ckpt.history = list(history)
    path.write_text(json.dumps(asdict(ckpt), indent=2))
    return path


def load_checkpoint(path: Path | str) -> tuple[HMM | BetaHMM, list[TrainingRun]]:
    payload = json.loads(Path(path).read_text())
    schema = payload.get("schema", 0)
    if schema != CHECKPOINT_SCHEMA:
        raise ValueError(
            f"checkpoint schema mismatch: got {schema}, expected {CHECKPOINT_SCHEMA}"
        )
    history_raw = payload.pop("history", [])
    ckpt = Checkpoint(**payload)
    return ckpt.to_hmm(), [TrainingRun(**r) for r in history_raw]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ── Registry ─────────────────────────────────────────────────────────────────


@dataclass
class RegistryEntry:
    name: str
    path: str
    model_type: str
    created: str
    updated: str
    total_runs: int
    total_samples_seen: int
    cohorts_seen: list[str]
    last_log_likelihood: float | None


def _registry_path(models_dir: Path) -> Path:
    return models_dir / "registry.json"


def read_registry(models_dir: Path) -> list[RegistryEntry]:
    path = _registry_path(models_dir)
    if not path.exists():
        return []
    return [RegistryEntry(**e) for e in json.loads(path.read_text())]


def write_registry(models_dir: Path, entries: list[RegistryEntry]) -> None:
    path = _registry_path(models_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps([asdict(e) for e in entries], indent=2))


def upsert_registry(
    models_dir: Path,
    name: str,
    checkpoint_path: Path,
    history: list[TrainingRun],
) -> RegistryEntry:
    models_dir = Path(models_dir).resolve()
    checkpoint_path = Path(checkpoint_path).resolve()
    entries = read_registry(models_dir)

    payload = json.loads(checkpoint_path.read_text())
    model_type = payload.get("model_type", "standard")

    try:
        stored_path = str(checkpoint_path.relative_to(models_dir))
    except ValueError:
        stored_path = str(checkpoint_path)

    cohorts_seen: list[str] = []
    samples = 0
    for r in history:
        if r.cohort and r.cohort not in cohorts_seen:
            cohorts_seen.append(r.cohort)
        if r.n_samples:
            samples += r.n_samples

    now = utc_now()
    last_ll = history[-1].log_likelihood_end if history else None

    for i, e in enumerate(entries):
        if e.name == name:
            entries[i] = RegistryEntry(
                name=name,
                path=stored_path,
                model_type=model_type,
                created=e.created,
                updated=now,
                total_runs=len(history),
                total_samples_seen=samples,
                cohorts_seen=cohorts_seen,
                last_log_likelihood=last_ll,
            )
            write_registry(models_dir, entries)
            return entries[i]

    entry = RegistryEntry(
        name=name,
        path=stored_path,
        model_type=model_type,
        created=now,
        updated=now,
        total_runs=len(history),
        total_samples_seen=samples,
        cohorts_seen=cohorts_seen,
        last_log_likelihood=last_ll,
    )
    entries.append(entry)
    write_registry(models_dir, entries)
    return entry
