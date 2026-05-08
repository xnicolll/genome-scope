"""Tests for Forward-Backward (T-03).

The canonical sanity checks for any HMM implementation:
  1. sum_s exp(alpha[t, s] * beta[t, s] / P(O)) == 1 for every t
  2. log P(O) from forward == log P(O) from backward
  3. results match a brute-force enumeration on a tiny sequence
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.special import logsumexp

from genomescope.hmm.forward_backward import backward, forward, forward_backward
from genomescope.hmm.model import HMM, standard_cpg_hmm


def _brute_force_loglik(hmm: HMM, obs: np.ndarray) -> float:
    """Sum over every possible state path (only tractable for tiny T)."""
    T = len(obs)
    total = -np.inf
    for path in itertools.product(range(hmm.n_states), repeat=T):
        logp = hmm.log_init[path[0]] + hmm.log_emit[path[0], obs[0]]
        for t in range(1, T):
            logp += hmm.log_trans[path[t - 1], path[t]] + hmm.log_emit[path[t], obs[t]]
        total = np.logaddexp(total, logp)
    return float(total)


def _simple_hmm() -> HMM:
    """Tiny non-degenerate HMM (no -inf emissions) for brute-force tests."""
    log_init = np.log(np.array([0.5, 0.5]))
    log_trans = np.log(np.array([[0.7, 0.3], [0.4, 0.6]]))
    log_emit = np.log(np.array([[0.9, 0.1], [0.2, 0.8]]))
    return HMM(
        n_states=2,
        n_obs=2,
        log_init=log_init,
        log_trans=log_trans,
        log_emit=log_emit,
    )


def test_forward_likelihood_matches_bruteforce() -> None:
    hmm = _simple_hmm()
    obs = np.array([0, 1, 0, 1, 1], dtype=np.int64)
    log_alpha = forward(hmm, obs)
    ll = logsumexp(log_alpha[-1])
    assert np.isclose(ll, _brute_force_loglik(hmm, obs))


def test_backward_matches_forward_likelihood() -> None:
    hmm = _simple_hmm()
    obs = np.array([0, 1, 0, 1, 1], dtype=np.int64)
    log_alpha = forward(hmm, obs)
    log_beta = backward(hmm, obs)
    # log P(O) from backward: logsumexp(log_init + emit[:, 0] + beta[0])
    back_ll = logsumexp(hmm.log_init + hmm.log_emit[:, obs[0]] + log_beta[0])
    fwd_ll = logsumexp(log_alpha[-1])
    assert np.isclose(fwd_ll, back_ll)


def test_posteriors_sum_to_one() -> None:
    hmm = _simple_hmm()
    obs = np.array([0, 1, 0, 1, 1, 0, 0], dtype=np.int64)
    fb = forward_backward(hmm, obs)
    row_sums = np.exp(fb.log_gamma).sum(axis=1)
    assert np.allclose(row_sums, 1.0)


def test_standard_cpg_hmm_on_short_sequence() -> None:
    """Deterministic emissions + -inf log probs must not break FB."""
    hmm = standard_cpg_hmm()
    # ACGTACG → indices [0,1,2,3,0,1,2]
    obs = np.array([0, 1, 2, 3, 0, 1, 2], dtype=np.int64)
    fb = forward_backward(hmm, obs)
    assert np.isfinite(fb.log_likelihood)
    assert np.allclose(np.exp(fb.log_gamma).sum(axis=1), 1.0)


def test_forward_backward_lengths_match() -> None:
    hmm = _simple_hmm()
    obs = np.array([0, 1, 0], dtype=np.int64)
    fb = forward_backward(hmm, obs)
    assert fb.log_alpha.shape == (3, 2)
    assert fb.log_beta.shape == (3, 2)
    assert fb.log_gamma.shape == (3, 2)
