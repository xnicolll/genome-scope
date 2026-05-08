"""Beta-emission HMM: each state emits a methylation score in (0, 1).

Hidden states are region labels (depleted, island); per-state Beta(alpha, beta)
parameters are learned via method-of-moments during Baum-Welch.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.stats import beta as beta_dist

from .model import _log

N_BETA_STATES = 2
DEPLETED = 0
ISLAND = 1

_OBS_EPS = 1e-6  # keep observations off the 0/1 boundary where Beta logpdf -> -inf


@dataclass
class BetaHMM:
    """Beta-emission HMM conforming to the HMMLike protocol."""

    log_init: np.ndarray
    log_trans: np.ndarray
    alphas: np.ndarray
    betas: np.ndarray
    state_labels: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.n_states = self.log_init.shape[0]
        assert self.log_trans.shape == (self.n_states, self.n_states)
        assert self.alphas.shape == (self.n_states,)
        assert self.betas.shape == (self.n_states,)
        if not self.state_labels:
            self.state_labels = [f"S{i}" for i in range(self.n_states)]

    def emit_logp(self, obs: np.ndarray) -> np.ndarray:
        obs_c = np.clip(obs, _OBS_EPS, 1.0 - _OBS_EPS)
        return beta_dist.logpdf(
            obs_c[None, :],
            self.alphas[:, None],
            self.betas[:, None],
        )

    def state_means(self) -> np.ndarray:
        return self.alphas / (self.alphas + self.betas)

    def m_step_emissions(self, obs: np.ndarray, gamma: np.ndarray) -> "BetaHMM":
        return update_beta_params(self, obs, gamma)


def beta_cpg_hmm(
    p_switch: float = 1e-2,
    island_alpha: float = 5.0,
    island_beta: float = 0.5,
    depleted_alpha: float = 0.5,
    depleted_beta: float = 5.0,
) -> BetaHMM:
    """2-state Beta HMM with the PRD's hypermethylated-island prior."""
    p_stay = 1.0 - p_switch
    return BetaHMM(
        log_init=_log(np.array([0.5, 0.5])),
        log_trans=_log(np.array([[p_stay, p_switch], [p_switch, p_stay]])),
        alphas=np.array([depleted_alpha, island_alpha], dtype=np.float64),
        betas=np.array([depleted_beta, island_beta], dtype=np.float64),
        state_labels=["depleted", "island"],
    )


def beta_method_of_moments(mean: float, var: float) -> tuple[float, float]:
    """Map (sample mean, variance) to (alpha, beta) via method-of-moments.

    alpha = mean * ((mean*(1-mean)/var) - 1)
    beta  = (1-mean) * ((mean*(1-mean)/var) - 1)
    Returns (1, 1) when the moments are inconsistent with any Beta.
    """
    if not (0.0 < mean < 1.0):
        return 1.0, 1.0
    max_var = mean * (1.0 - mean)
    if var <= 0 or var >= max_var:
        return 1.0, 1.0
    k = max_var / var - 1.0
    return max(mean * k, 1e-3), max((1.0 - mean) * k, 1e-3)


def update_beta_params(
    hmm: BetaHMM,
    obs: np.ndarray,
    gamma: np.ndarray,
) -> BetaHMM:
    new_alphas = hmm.alphas.copy()
    new_betas = hmm.betas.copy()
    total_weight = gamma.sum(axis=0)

    for s in range(hmm.n_states):
        w = total_weight[s]
        if w < 1e-8:
            continue
        mean = float((gamma[:, s] * obs).sum() / w)
        diff = obs - mean
        var = float((gamma[:, s] * diff * diff).sum() / w)
        new_alphas[s], new_betas[s] = beta_method_of_moments(mean, var)

    return replace(hmm, alphas=new_alphas, betas=new_betas)
