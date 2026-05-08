"""Tests for Baum-Welch training (T-05)."""

from __future__ import annotations

import numpy as np

from genomescope.hmm.baum_welch import baum_welch
from genomescope.hmm.model import HMM, standard_cpg_hmm


def _simple_hmm() -> HMM:
    log_init = np.log(np.array([0.6, 0.4]))
    log_trans = np.log(np.array([[0.7, 0.3], [0.4, 0.6]]))
    log_emit = np.log(np.array([[0.9, 0.1], [0.2, 0.8]]))
    return HMM(n_states=2, n_obs=2, log_init=log_init, log_trans=log_trans, log_emit=log_emit)


def test_log_likelihood_monotone_nondecreasing() -> None:
    """Baum-Welch must not decrease log-likelihood between iterations."""
    hmm = _simple_hmm()
    # biased sequence - should move parameters
    rng = np.random.default_rng(0)
    obs = rng.choice([0, 1], size=200, p=[0.3, 0.7]).astype(np.int64)
    result = baum_welch(hmm, obs, max_iter=30, tol=1e-8, freeze_emissions=False)
    lls = result.log_likelihoods
    # allow tiny numerical wiggle
    for prev, curr in zip(lls, lls[1:], strict=False):
        assert curr >= prev - 1e-8


def test_baum_welch_converges_on_structured_sequence() -> None:
    """On a clearly biased sequence BW should hit the 1e-4 tolerance quickly."""
    hmm = _simple_hmm()
    rng = np.random.default_rng(1)
    # alternating emission regimes → two clear latent states
    obs = np.concatenate([
        rng.choice([0, 1], size=250, p=[0.9, 0.1]),
        rng.choice([0, 1], size=250, p=[0.1, 0.9]),
    ]).astype(np.int64)
    result = baum_welch(hmm, obs, max_iter=100, tol=1e-4, freeze_emissions=False)
    assert result.converged
    assert result.iterations < 100


def test_transitions_remain_valid_probabilities() -> None:
    hmm = _simple_hmm()
    rng = np.random.default_rng(2)
    obs = rng.choice([0, 1], size=150).astype(np.int64)
    result = baum_welch(hmm, obs, max_iter=20, freeze_emissions=False)
    trans = np.exp(result.hmm.log_trans)
    assert np.allclose(trans.sum(axis=1), 1.0)
    init = np.exp(result.hmm.log_init)
    assert np.isclose(init.sum(), 1.0)


def test_standard_cpg_hmm_training_preserves_emissions() -> None:
    """With freeze_emissions=True the emission matrix must be unchanged."""
    hmm = standard_cpg_hmm()
    # small synthetic chr21-like sequence with an island chunk
    seq = "AAAAAATTTTAA" "CGCGCGCGCGCG" "AAAATTTTAAAA"
    obs = np.array([{"A": 0, "C": 1, "G": 2, "T": 3}[c] for c in seq], dtype=np.int64)

    result = baum_welch(hmm, obs, max_iter=10, freeze_emissions=True)
    assert np.allclose(
        np.exp(result.hmm.log_emit), np.exp(hmm.log_emit)
    )


def test_baum_welch_learns_transition_bias() -> None:
    """Train the 8-state CpG HMM on a CpG-rich vs CpG-poor mix -
    it should raise the C_island→G_island transition above its init."""
    hmm = standard_cpg_hmm(p_switch=0.05)
    seq = ("ATAT" * 200) + ("CGCG" * 400) + ("ATAT" * 200)
    obs = np.array([{"A": 0, "C": 1, "G": 2, "T": 3}[c] for c in seq], dtype=np.int64)
    result = baum_welch(hmm, obs, max_iter=30, tol=1e-5, freeze_emissions=True)

    trans = np.exp(result.hmm.log_trans)
    # C_island = state 5, G_island = state 6
    c_island_to_g_island = trans[5, 6]
    init_trans = np.exp(hmm.log_trans)[5, 6]
    assert c_island_to_g_island > init_trans
