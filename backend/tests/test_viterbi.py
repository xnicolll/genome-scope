"""Tests for Viterbi decoding (T-04)."""

from __future__ import annotations

import itertools

import numpy as np

from genomescope.hmm.model import HMM, standard_cpg_hmm
from genomescope.hmm.viterbi import merge_and_filter, segments, viterbi


def _simple_hmm() -> HMM:
    log_init = np.log(np.array([0.5, 0.5]))
    log_trans = np.log(np.array([[0.7, 0.3], [0.4, 0.6]]))
    log_emit = np.log(np.array([[0.9, 0.1], [0.2, 0.8]]))
    return HMM(n_states=2, n_obs=2, log_init=log_init, log_trans=log_trans, log_emit=log_emit)


def _brute_force_viterbi(hmm: HMM, obs: np.ndarray) -> tuple[list[int], float]:
    T = len(obs)
    best_path: list[int] = []
    best_lp = -np.inf
    for path in itertools.product(range(hmm.n_states), repeat=T):
        lp = hmm.log_init[path[0]] + hmm.log_emit[path[0], obs[0]]
        for t in range(1, T):
            lp += hmm.log_trans[path[t - 1], path[t]] + hmm.log_emit[path[t], obs[t]]
        if lp > best_lp:
            best_lp = lp
            best_path = list(path)
    return best_path, float(best_lp)


def test_viterbi_matches_bruteforce() -> None:
    hmm = _simple_hmm()
    obs = np.array([0, 1, 0, 1, 1, 0], dtype=np.int64)
    result = viterbi(hmm, obs)
    expected_path, expected_lp = _brute_force_viterbi(hmm, obs)
    assert result.path.tolist() == expected_path
    assert np.isclose(result.log_prob, expected_lp)


def test_viterbi_respects_deterministic_emissions() -> None:
    """In the standard CpG HMM, emissions are deterministic -
    the path must emit exactly the observed nucleotides."""
    hmm = standard_cpg_hmm()
    obs = np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.int64)  # ACGTACGT
    result = viterbi(hmm, obs)
    # each state's nucleotide (state % 4) must match the observation
    assert (result.path % 4 == obs).all()


def test_viterbi_log_prob_upper_bounded_by_forward() -> None:
    """Best-path log-prob must be <= total log-likelihood."""
    from genomescope.hmm.forward_backward import forward
    from scipy.special import logsumexp

    hmm = _simple_hmm()
    obs = np.array([0, 1, 0, 1, 1, 0], dtype=np.int64)
    vit = viterbi(hmm, obs)
    fwd_ll = float(logsumexp(forward(hmm, obs)[-1]))
    assert vit.log_prob <= fwd_ll + 1e-9


def test_segments_simple() -> None:
    # island states in standard HMM are 4..7
    path = np.array([0, 4, 5, 6, 1, 2, 7, 7, 3], dtype=np.int32)
    predicate = np.array([False] * 4 + [True] * 4)
    segs = segments(path, predicate)
    # islands at [1,4) and [6,8)
    assert segs == [(1, 4), (6, 8)]


def test_segments_no_matches() -> None:
    path = np.array([0, 1, 2, 3], dtype=np.int32)
    predicate = np.array([False] * 4 + [True] * 4)
    assert segments(path, predicate) == []


def test_segments_full_match() -> None:
    path = np.array([4, 5, 6, 7], dtype=np.int32)
    predicate = np.array([False] * 4 + [True] * 4)
    assert segments(path, predicate) == [(0, 4)]


def test_merge_and_filter_merges_close_segments() -> None:
    segs = [(0, 300), (350, 600), (10_000, 10_500)]
    out = merge_and_filter(segs, min_length=200, merge_gap=100)
    # first two merge (gap=50), third stays separate
    assert out == [(0, 600), (10_000, 10_500)]


def test_merge_and_filter_drops_short_segments() -> None:
    segs = [(0, 150), (5000, 5500)]
    out = merge_and_filter(segs, min_length=200, merge_gap=100)
    assert out == [(5000, 5500)]


def test_merge_and_filter_empty() -> None:
    assert merge_and_filter([], min_length=200, merge_gap=100) == []
