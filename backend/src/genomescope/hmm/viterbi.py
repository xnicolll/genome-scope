"""Log-space Viterbi decoding."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .forward_backward import HMMLike


@dataclass
class ViterbiResult:
    path: np.ndarray
    log_prob: float


def viterbi(hmm: HMMLike, obs: np.ndarray) -> ViterbiResult:
    T = obs.shape[0]
    n = hmm.n_states
    delta = np.empty((T, n), dtype=np.float64)
    psi = np.empty((T, n), dtype=np.int32)

    emit = hmm.emit_logp(obs)
    delta[0] = hmm.log_init + emit[:, 0]
    psi[0] = 0

    for t in range(1, T):
        scores = delta[t - 1][:, None] + hmm.log_trans
        psi[t] = np.argmax(scores, axis=0)
        delta[t] = scores[psi[t], np.arange(n)] + emit[:, t]

    path = np.empty(T, dtype=np.int32)
    path[T - 1] = int(np.argmax(delta[T - 1]))
    for t in range(T - 2, -1, -1):
        path[t] = psi[t + 1, path[t + 1]]

    return ViterbiResult(path=path, log_prob=float(delta[T - 1, path[T - 1]]))


def segments(path: np.ndarray, predicate: np.ndarray) -> list[tuple[int, int]]:
    """Half-open ranges where predicate[state] is True."""
    labels = predicate[path]
    if labels.size == 0:
        return []
    diff = np.diff(labels.astype(np.int8), prepend=0, append=0)
    starts = np.flatnonzero(diff == 1)
    ends = np.flatnonzero(diff == -1)
    return list(zip(starts.tolist(), ends.tolist(), strict=True))


def merge_and_filter(
    segs: list[tuple[int, int]],
    min_length: int = 200,
    merge_gap: int = 100,
) -> list[tuple[int, int]]:
    """Sort, merge any pair within merge_gap, drop anything shorter than min_length.

    Mirrors the Gardiner-Garden & Frommer (1987) curation pipeline.
    """
    if not segs:
        return []
    segs = sorted(segs)
    merged: list[tuple[int, int]] = [segs[0]]
    for s, e in segs[1:]:
        ms, me = merged[-1]
        if s - me <= merge_gap:
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return [(s, e) for s, e in merged if e - s >= min_length]
