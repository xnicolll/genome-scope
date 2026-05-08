"""Classical 8-state CpG-island HMM (Durbin et al., 1998).

States = {A, C, G, T} x {depleted, island}; emissions are deterministic on
nucleotide; only transitions are learned via Baum-Welch.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import IntEnum

import numpy as np

N_NUCLEOTIDES = 4
LOG_ZERO = -np.inf


class Region(IntEnum):
    DEPLETED = 0
    ISLAND = 1


def state_index(nuc: int, region: Region) -> int:
    return int(region) * N_NUCLEOTIDES + nuc


def split_state(state: int) -> tuple[int, Region]:
    return state % N_NUCLEOTIDES, Region(state // N_NUCLEOTIDES)


def _log(x: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore"):
        return np.log(x)


@dataclass
class HMM:
    """Log-space HMM parameters with categorical emissions."""

    n_states: int
    n_obs: int
    log_init: np.ndarray
    log_trans: np.ndarray
    log_emit: np.ndarray
    state_labels: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        assert self.log_init.shape == (self.n_states,)
        assert self.log_trans.shape == (self.n_states, self.n_states)
        assert self.log_emit.shape == (self.n_states, self.n_obs)
        if not self.state_labels:
            self.state_labels = [f"S{i}" for i in range(self.n_states)]

    def emit_logp(self, obs: np.ndarray) -> np.ndarray:
        return self.log_emit[:, obs]

    def m_step_emissions(self, obs: np.ndarray, gamma: np.ndarray) -> "HMM":
        total = gamma.sum(axis=0)
        new_emit = np.zeros((self.n_states, self.n_obs), dtype=np.float64)
        for k in range(self.n_obs):
            new_emit[:, k] = gamma[obs == k].sum(axis=0)
        denom = np.where(total > 0, total, 1.0)
        return replace(self, log_emit=_log(new_emit / denom[:, None]))


# Durbin et al. (1998), Table 3.1 - rows are current nucleotide, columns next.
DURBIN_ISLAND = np.array([
    [0.180, 0.274, 0.426, 0.120],  # A+
    [0.171, 0.368, 0.274, 0.188],  # C+ (P(G|C)=0.274 - enriched)
    [0.161, 0.339, 0.375, 0.125],  # G+
    [0.079, 0.355, 0.384, 0.182],  # T+
], dtype=np.float64)

DURBIN_DEPLETED = np.array([
    [0.300, 0.205, 0.285, 0.210],  # A-
    [0.322, 0.298, 0.078, 0.302],  # C- (P(G|C)=0.078 - depleted)
    [0.248, 0.246, 0.298, 0.208],  # G-
    [0.177, 0.239, 0.292, 0.292],  # T-
], dtype=np.float64)


def standard_cpg_hmm(p_switch: float = 1e-3) -> HMM:
    """Textbook 8-state CpG island HMM with Durbin (1998) initialisation."""
    n_states = 2 * N_NUCLEOTIDES
    p_stay = 1.0 - p_switch

    trans = np.zeros((n_states, n_states), dtype=np.float64)
    for i in range(n_states):
        nuc_i, region_i = split_state(i)
        for j in range(n_states):
            nuc_j, region_j = split_state(j)
            target_region = region_i if region_i == region_j else region_j
            row = DURBIN_ISLAND if target_region is Region.ISLAND else DURBIN_DEPLETED
            mass = p_stay if region_i == region_j else p_switch
            trans[i, j] = mass * row[nuc_i, nuc_j]
    trans /= trans.sum(axis=1, keepdims=True)
    assert np.allclose(trans.sum(axis=1), 1.0)

    emit = np.zeros((n_states, N_NUCLEOTIDES), dtype=np.float64)
    for i in range(n_states):
        nuc, _ = split_state(i)
        emit[i, nuc] = 1.0

    init = np.full(n_states, 1.0 / n_states, dtype=np.float64)
    labels = [f"{nuc}_{region.name.lower()}" for region in Region for nuc in "ACGT"]

    return HMM(
        n_states=n_states,
        n_obs=N_NUCLEOTIDES,
        log_init=_log(init),
        log_trans=_log(trans),
        log_emit=_log(emit),
        state_labels=labels,
    )


def region_of_state(state: int) -> Region:
    return split_state(state)[1]


def state_is_island(state: int) -> bool:
    return region_of_state(state) is Region.ISLAND
