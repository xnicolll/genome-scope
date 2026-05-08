"""TCGA methylation β-value loader (SESAMe level-3 + BED) and synthetic fallback."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .fasta import Sequence, find_cpg_sites
from .probes import load_hm450_probe_map


@dataclass
class MethylationTrack:
    chrom: str
    positions: np.ndarray
    betas: np.ndarray

    def __post_init__(self) -> None:
        assert self.positions.shape == self.betas.shape
        assert np.all(np.diff(self.positions) >= 0), "positions must be sorted"

    def __len__(self) -> int:
        return int(self.positions.shape[0])

    def window(self, start: int, end: int) -> "MethylationTrack":
        mask = (self.positions >= start) & (self.positions < end)
        return MethylationTrack(
            chrom=self.chrom,
            positions=self.positions[mask],
            betas=self.betas[mask],
        )


_BETA_SUFFIXES = {".tsv", ".txt", ".bed", ".csv"}
_IGNORE_NAMES = {"annotations.txt", "MANIFEST.txt", "README.md"}
_IGNORE_PATTERNS = ("subset.manifest.", "subset.metadata.", ".manifest.")
_IGNORE_DIRS = {"logs"}


def _sniff_format(path: Path) -> str:
    """'sesame' for probe_id<TAB>beta, 'bed' for chrom<TAB>start<TAB>end<TAB>beta."""
    with open(path) as f:
        for _ in range(5):
            line = f.readline()
            if not line or line.startswith("#"):
                continue
            first = line.split("\t")[0].strip()
            if first.startswith("cg") or first.startswith("rs") or first.startswith("ch."):
                return "sesame"
            if first.startswith("chr"):
                return "bed"
    return "bed"


def _iter_beta_files(directory: Path) -> list[Path]:
    out: list[Path] = []
    for p in directory.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix.lower() not in _BETA_SUFFIXES:
            continue
        if p.name in _IGNORE_NAMES:
            continue
        if any(pat in p.name for pat in _IGNORE_PATTERNS):
            continue
        if any(part in _IGNORE_DIRS for part in p.parts):
            continue
        out.append(p)
    return sorted(out)


def _load_sesame_level3(
    paths: list[Path],
    probe_map: pd.DataFrame,
    chrom: str,
) -> MethylationTrack | None:
    chrom_probes = probe_map[probe_map["chrom"] == chrom] if "chrom" in probe_map.columns else probe_map

    frames: list[pd.Series] = []
    for path in paths:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["probe_id", "beta"],
            comment="#",
            dtype={"probe_id": str, "beta": float},
            na_values=["NA", "nan", "NaN"],
            on_bad_lines="skip",
        )
        df = df.dropna(subset=["beta"])
        df = df[df["probe_id"].isin(chrom_probes.index)]
        if df.empty:
            continue
        frames.append(df.set_index("probe_id")["beta"])

    if not frames:
        return None

    combined = pd.concat(frames, axis=1)
    mean_beta = combined.mean(axis=1, skipna=True).dropna()

    joined = chrom_probes.loc[mean_beta.index].copy()
    joined["beta"] = mean_beta
    joined = joined.sort_values("start")

    return MethylationTrack(
        chrom=chrom,
        positions=joined["start"].to_numpy(dtype=np.int64),
        betas=joined["beta"].to_numpy(dtype=np.float64),
    )


def _load_bed_format(
    paths: list[Path],
    chrom: str,
) -> MethylationTrack | None:
    frames: list[pd.DataFrame] = []
    for path in paths:
        try:
            df = pd.read_csv(
                path,
                sep=None,
                engine="python",
                comment="#",
                header=None,
                usecols=[0, 1, 2, 3],
                names=["chrom", "start", "end", "beta"],
                dtype={"chrom": str, "start": "Int64", "end": "Int64", "beta": float},
                on_bad_lines="skip",
            )
        except (ValueError, pd.errors.ParserError):
            continue
        df = df[df["chrom"] == chrom].dropna()
        if df.empty:
            continue
        frames.append(df)

    if not frames:
        return None

    merged = pd.concat(frames, ignore_index=True)
    averaged = (
        merged.groupby("start", as_index=False)["beta"].mean().sort_values("start")
    )
    return MethylationTrack(
        chrom=chrom,
        positions=averaged["start"].to_numpy(dtype=np.int64),
        betas=averaged["beta"].to_numpy(dtype=np.float64),
    )


def load_tcga_cohorts(
    root: Path,
    chrom: str = "chr21",
) -> dict[str, MethylationTrack]:
    """Map each subdirectory of `root` to a MethylationTrack (platform suffix stripped)."""
    root = Path(root)
    if not root.exists():
        return {}
    tracks: dict[str, MethylationTrack] = {}
    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        track = load_tcga_directory(child, chrom=chrom)
        if track is None:
            continue
        key = child.name
        for suffix in ("_hm450", "_epic", "_hm27"):
            if key.endswith(suffix):
                key = key[: -len(suffix)]
                break
        tracks[key] = track
    return tracks


def load_tcga_directory(
    directory: Path,
    chrom: str = "chr21",
) -> MethylationTrack | None:
    """Load + average TCGA β-files in `directory`; auto-detects SESAMe vs BED format."""
    directory = Path(directory)
    if not directory.exists():
        return None

    files = _iter_beta_files(directory)
    if not files:
        return None

    fmt = _sniff_format(files[0])
    if fmt == "sesame":
        probe_map = load_hm450_probe_map(chrom=chrom)
        return _load_sesame_level3(files, probe_map, chrom)
    return _load_bed_format(files, chrom)


def synthetic_track_from_fasta(
    sequence: Sequence,
    islands: list[tuple[int, int]],
    seed: int = 7,
    island_alpha: float = 5.0,
    island_beta: float = 0.5,
    depleted_alpha: float = 0.5,
    depleted_beta: float = 5.0,
    subset: tuple[int, int] | None = None,
) -> MethylationTrack:
    """Synthetic per-CpG β track: in-island Beta(α,β) hyper, out-of-island depleted."""
    cpg_positions = find_cpg_sites(sequence.seq)
    if subset is not None:
        lo, hi = subset
        mask = (cpg_positions >= lo) & (cpg_positions < hi)
        cpg_positions = cpg_positions[mask]

    rng = np.random.default_rng(seed)

    island_starts = np.array([s for s, _ in islands], dtype=np.int64)
    island_ends = np.array([e for _, e in islands], dtype=np.int64)
    order = np.argsort(island_starts)
    island_starts = island_starts[order]
    island_ends = island_ends[order]
    idx = np.searchsorted(island_starts, cpg_positions, side="right") - 1
    in_island = (idx >= 0) & (cpg_positions < island_ends[np.clip(idx, 0, None)])

    betas = np.empty(cpg_positions.shape[0], dtype=np.float64)
    betas[in_island] = rng.beta(island_alpha, island_beta, size=int(in_island.sum()))
    betas[~in_island] = rng.beta(
        depleted_alpha, depleted_beta, size=int((~in_island).sum())
    )

    return MethylationTrack(
        chrom=sequence.chrom,
        positions=cpg_positions.astype(np.int64),
        betas=betas,
    )
