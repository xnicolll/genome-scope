"""Log-space Forward-Backward.

Stable across chromosome-scale sequences via scipy.special.logsumexp.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from scipy.special import logsumexp


class HMMLike(Protocol):
    n_states: int
    log_init: np.ndarray
    log_trans: np.ndarray

    def emit_logp(self, obs: np.ndarray) -> np.ndarray: ...


@dataclass
class FBResult:
    log_alpha: np.ndarray
    log_beta: np.ndarray
    log_gamma: np.ndarray
    log_likelihood: float
    log_emit: np.ndarray


def _forward_core(log_init: np.ndarray, log_trans: np.ndarray, emit: np.ndarray) -> np.ndarray:
    T = emit.shape[1]
    log_alpha = np.empty((T, log_init.shape[0]), dtype=np.float64)
    log_alpha[0] = log_init + emit[:, 0]
    for t in range(1, T):
        log_alpha[t] = (
            logsumexp(log_alpha[t - 1][:, None] + log_trans, axis=0) + emit[:, t]
        )
    return log_alpha


def _backward_core(log_trans: np.ndarray, emit: np.ndarray) -> np.ndarray:
    T = emit.shape[1]
    log_beta = np.empty((T, log_trans.shape[0]), dtype=np.float64)
    log_beta[T - 1] = 0.0
    for t in range(T - 2, -1, -1):
        log_beta[t] = logsumexp(
            log_trans + emit[:, t + 1][None, :] + log_beta[t + 1][None, :],
            axis=1,
        )
    return log_beta


def forward(hmm: HMMLike, obs: np.ndarray) -> np.ndarray:
    return _forward_core(hmm.log_init, hmm.log_trans, hmm.emit_logp(obs))


def backward(hmm: HMMLike, obs: np.ndarray) -> np.ndarray:
    return _backward_core(hmm.log_trans, hmm.emit_logp(obs))


def forward_backward(hmm: HMMLike, obs: np.ndarray) -> FBResult:
    emit = hmm.emit_logp(obs)
    log_alpha = _forward_core(hmm.log_init, hmm.log_trans, emit)
    log_beta = _backward_core(hmm.log_trans, emit)
    log_likelihood = float(logsumexp(log_alpha[-1]))
    return FBResult(
        log_alpha=log_alpha,
        log_beta=log_beta,
        log_gamma=log_alpha + log_beta - log_likelihood,
        log_likelihood=log_likelihood,
        log_emit=emit,
    )
