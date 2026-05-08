"""Generic Baum-Welch EM training.

Emission-agnostic: delegates the emission update to `hmm.m_step_emissions(obs, gamma)`
so categorical HMM and BetaHMM share the same loop.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .forward_backward import forward_backward
from .model import _log


@dataclass
class BWResult:
    hmm: object
    log_likelihoods: list[float]
    converged: bool
    iterations: int


def _accumulate_transitions(
    log_alpha: np.ndarray,
    log_beta: np.ndarray,
    log_trans: np.ndarray,
    emit: np.ndarray,
    log_likelihood: float,
    chunk_size: int = 100_000,
) -> np.ndarray:
    """Sum_t xi[t, i, j] without materialising the full (T, N, N) tensor."""
    T = emit.shape[1]
    N = log_trans.shape[0]
    accum = np.zeros((N, N), dtype=np.float64)

    for start in range(0, T - 1, chunk_size):
        end = min(start + chunk_size, T - 1)
        a = log_alpha[start:end, :, None]
        b = emit[:, start + 1 : end + 1].T + log_beta[start + 1 : end + 1, :]
        b = b[:, None, :]
        with np.errstate(invalid="ignore"):
            xi_chunk = np.exp(a + log_trans[None, :, :] + b - log_likelihood)
        np.nan_to_num(xi_chunk, copy=False, nan=0.0)
        accum += xi_chunk.sum(axis=0)
    return accum


def baum_welch(
    hmm,
    obs: np.ndarray,
    max_iter: int = 100,
    tol: float = 1e-4,
    freeze_emissions: bool = True,
    freeze_init: bool = False,
    verbose: bool = False,
) -> BWResult:
    log_likelihoods: list[float] = []
    current = hmm
    prev_ll = -np.inf
    converged = False

    for it in range(1, max_iter + 1):
        fb = forward_backward(current, obs)
        log_likelihoods.append(fb.log_likelihood)
        if verbose:
            print(f"  iter {it:3d}  log-lik = {fb.log_likelihood:.6f}")

        if fb.log_likelihood - prev_ll < tol and it > 1:
            converged = True
            break
        prev_ll = fb.log_likelihood

        gamma = np.exp(fb.log_gamma)
        np.nan_to_num(gamma, copy=False, nan=0.0)

        new_init = gamma[0] if not freeze_init else np.exp(current.log_init)

        xi_sum = _accumulate_transitions(
            fb.log_alpha, fb.log_beta, current.log_trans, fb.log_emit, fb.log_likelihood,
        )
        residence = gamma[:-1].sum(axis=0)
        denom = np.where(residence > 0, residence, 1.0)
        new_trans = xi_sum / denom[:, None]
        zero_rows = residence == 0
        if zero_rows.any():
            new_trans[zero_rows] = np.exp(current.log_trans[zero_rows])
        row_sums = new_trans.sum(axis=1, keepdims=True)
        new_trans = np.where(row_sums > 0, new_trans / row_sums, new_trans)

        current = replace(current, log_init=_log(new_init), log_trans=_log(new_trans))
        if not freeze_emissions:
            current = current.m_step_emissions(obs, gamma)

    return BWResult(
        hmm=current,
        log_likelihoods=log_likelihoods,
        converged=converged,
        iterations=len(log_likelihoods),
    )
