"""Tests for the Beta-emission HMM (T-07, T-08)."""

from __future__ import annotations

import numpy as np
from scipy.stats import beta as beta_dist

from genomescope.hmm.baum_welch import baum_welch
from genomescope.hmm.beta_model import (
    BetaHMM,
    beta_cpg_hmm,
    beta_method_of_moments,
    update_beta_params,
)
from genomescope.hmm.forward_backward import forward_backward
from genomescope.hmm.viterbi import viterbi


def test_beta_cpg_hmm_default_shape() -> None:
    hmm = beta_cpg_hmm()
    assert hmm.n_states == 2
    assert hmm.log_init.shape == (2,)
    assert hmm.log_trans.shape == (2, 2)
    assert hmm.alphas.shape == (2,) and hmm.betas.shape == (2,)
    # island mean >> depleted mean
    means = hmm.state_means()
    assert means[1] > 0.8  # island
    assert means[0] < 0.2  # depleted


def test_beta_emit_logp_matches_scipy() -> None:
    hmm = beta_cpg_hmm()
    obs = np.array([0.05, 0.5, 0.95], dtype=np.float64)
    logp = hmm.emit_logp(obs)
    assert logp.shape == (2, 3)
    for s in range(2):
        for t, o in enumerate(obs):
            expected = beta_dist.logpdf(o, hmm.alphas[s], hmm.betas[s])
            assert np.isclose(logp[s, t], expected)


def test_beta_emit_logp_handles_zero_and_one() -> None:
    """Observations exactly at 0 or 1 must not produce NaN."""
    hmm = beta_cpg_hmm()
    obs = np.array([0.0, 1.0], dtype=np.float64)
    logp = hmm.emit_logp(obs)
    assert np.isfinite(logp).all()


def test_beta_hmm_forward_backward_runs() -> None:
    hmm = beta_cpg_hmm()
    rng = np.random.default_rng(0)
    obs = np.clip(rng.beta(1.0, 5.0, size=50), 0.01, 0.99)
    fb = forward_backward(hmm, obs)
    assert fb.log_alpha.shape == (50, 2)
    assert np.isfinite(fb.log_likelihood)
    assert np.allclose(np.exp(fb.log_gamma).sum(axis=1), 1.0)


def test_beta_hmm_viterbi_separates_regimes() -> None:
    """On a sequence with a clear methylation regime switch, Viterbi
    should place the island state where the hypermethylated values are."""
    hmm = beta_cpg_hmm()
    rng = np.random.default_rng(1)
    low = rng.beta(0.5, 5.0, size=80)           # mean ~0.09
    high = rng.beta(5.0, 0.5, size=80)          # mean ~0.91
    obs = np.concatenate([low, high, low])
    result = viterbi(hmm, obs)
    # first and last thirds should be depleted (0), middle should be island (1)
    assert (result.path[:80] == 0).mean() > 0.9
    assert (result.path[80:160] == 1).mean() > 0.9
    assert (result.path[160:] == 0).mean() > 0.9


def test_method_of_moments_roundtrip() -> None:
    """Draw samples from a known Beta, fit params, recover shapes."""
    rng = np.random.default_rng(2)
    true_a, true_b = 4.0, 2.0
    samples = rng.beta(true_a, true_b, size=10_000)
    mean = float(samples.mean())
    var = float(samples.var())
    a, b = beta_method_of_moments(mean, var)
    assert abs(a - true_a) < 0.3
    assert abs(b - true_b) < 0.3


def test_method_of_moments_edge_cases() -> None:
    # mean outside (0, 1) → fall back to (1, 1)
    assert beta_method_of_moments(0.0, 0.1) == (1.0, 1.0)
    assert beta_method_of_moments(1.0, 0.1) == (1.0, 1.0)
    # impossible variance → fall back
    assert beta_method_of_moments(0.5, 10.0) == (1.0, 1.0)


def test_update_beta_params_learns_from_hard_assignments() -> None:
    """Hard posteriors + samples drawn from specific Betas should pull the
    params toward the true values."""
    hmm = beta_cpg_hmm()
    rng = np.random.default_rng(3)
    n = 500
    obs_dep = rng.beta(0.5, 8.0, size=n)
    obs_isl = rng.beta(8.0, 0.5, size=n)
    obs = np.concatenate([obs_dep, obs_isl])
    gamma = np.zeros((2 * n, 2), dtype=np.float64)
    gamma[:n, 0] = 1.0     # first n are depleted
    gamma[n:, 1] = 1.0     # last  n are island

    updated = update_beta_params(hmm, obs, gamma)
    dep_mean = updated.alphas[0] / (updated.alphas[0] + updated.betas[0])
    isl_mean = updated.alphas[1] / (updated.alphas[1] + updated.betas[1])
    assert dep_mean < 0.15
    assert isl_mean > 0.85


def test_baum_welch_trains_beta_hmm_end_to_end() -> None:
    """Full EM: start from a deliberately-bad init, train on mixed data,
    verify the final model recovers the correct regime assignment."""
    # Bad init: mild asymmetry (so EM doesn't land in the symmetric fixed
    # point) but nowhere near the truth yet.
    hmm = beta_cpg_hmm(
        island_alpha=3.0, island_beta=2.0,     # mean 0.60
        depleted_alpha=2.0, depleted_beta=3.0, # mean 0.40
    )
    rng = np.random.default_rng(4)
    obs = np.concatenate(
        [
            rng.beta(0.6, 6.0, size=200),   # depleted regime
            rng.beta(6.0, 0.6, size=200),   # island regime
            rng.beta(0.6, 6.0, size=200),
        ]
    )
    start_ll = float(
        forward_backward(hmm, obs).log_likelihood
    )
    result = baum_welch(
        hmm, obs, max_iter=30, tol=1e-5, freeze_emissions=False, verbose=False
    )
    # Method-of-moments is not a strict MLE update (per PRD T-08), so
    # per-iter monotonicity isn't guaranteed. What matters is that training
    # ends up substantially better than the starting point.
    assert result.log_likelihoods[-1] > start_ll + 50
    # The state means should have separated
    means = result.hmm.state_means()
    assert abs(means[1] - means[0]) > 0.5
    # Viterbi after training should recover the middle segment as one regime.
    # The state labels aren't identified (EM can swap them) so we check
    # consistency within each true regime rather than a specific label.
    path = viterbi(result.hmm, obs).path
    mid_mode = int(np.bincount(path[200:400]).argmax())
    assert (path[200:400] == mid_mode).mean() > 0.85
    other_mode = 1 - mid_mode
    assert (path[:200] == other_mode).mean() > 0.80
    assert (path[400:] == other_mode).mean() > 0.80
