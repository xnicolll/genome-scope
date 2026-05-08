"""Tests for the HMM dataclass and factories (T-02)."""

from __future__ import annotations

import numpy as np

from genomescope.hmm.model import (
    HMM,
    N_NUCLEOTIDES,
    Region,
    region_of_state,
    split_state,
    standard_cpg_hmm,
    state_index,
    state_is_island,
)


def test_state_index_roundtrip() -> None:
    for nuc in range(N_NUCLEOTIDES):
        for region in Region:
            s = state_index(nuc, region)
            assert split_state(s) == (nuc, region)


def test_standard_cpg_hmm_shapes() -> None:
    hmm = standard_cpg_hmm()
    assert isinstance(hmm, HMM)
    assert hmm.n_states == 8
    assert hmm.n_obs == 4
    assert hmm.log_init.shape == (8,)
    assert hmm.log_trans.shape == (8, 8)
    assert hmm.log_emit.shape == (8, 4)
    assert len(hmm.state_labels) == 8


def test_standard_cpg_hmm_rows_are_proper_distributions() -> None:
    hmm = standard_cpg_hmm()
    # rows of exp(log_trans) should sum to 1
    trans = np.exp(hmm.log_trans)
    assert np.allclose(trans.sum(axis=1), 1.0)
    # initial sums to 1
    assert np.isclose(np.exp(hmm.log_init).sum(), 1.0)


def test_standard_cpg_hmm_emissions_are_deterministic() -> None:
    hmm = standard_cpg_hmm()
    emit = np.exp(hmm.log_emit)
    # each state emits exactly one nucleotide with probability 1
    for state in range(hmm.n_states):
        nuc, _ = split_state(state)
        assert emit[state, nuc] == 1.0
        assert emit[state].sum() == 1.0


def test_switch_probability_is_small() -> None:
    hmm = standard_cpg_hmm(p_switch=1e-3)
    trans = np.exp(hmm.log_trans)
    # mass moving from depleted row 0 to any island state
    depleted_to_island = trans[0, 4:].sum()
    assert np.isclose(depleted_to_island, 1e-3)


def test_state_labels_ordered() -> None:
    hmm = standard_cpg_hmm()
    assert hmm.state_labels[:4] == ["A_depleted", "C_depleted", "G_depleted", "T_depleted"]
    assert hmm.state_labels[4:] == ["A_island", "C_island", "G_island", "T_island"]


def test_emit_logp_broadcasts() -> None:
    hmm = standard_cpg_hmm()
    obs = np.array([0, 1, 2, 3, 0], dtype=np.int64)  # A C G T A
    logp = hmm.emit_logp(obs)
    assert logp.shape == (8, 5)
    # state 0 = A_depleted emits A deterministically
    assert logp[0, 0] == 0.0  # log(1)
    assert logp[0, 1] == -np.inf  # log(0)


def test_region_helpers() -> None:
    assert not state_is_island(0)
    assert state_is_island(4)
    assert region_of_state(7) is Region.ISLAND
    assert region_of_state(3) is Region.DEPLETED
