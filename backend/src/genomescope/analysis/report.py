"""Isoform-aware methylation + cancer report (islands × promoters × cohorts × COSMIC)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd

from ..data.tcga import MethylationTrack
from .promoters import overlap_islands_with_promoters

HYPERMETH_THRESHOLD = 0.6


@dataclass
class ReportInputs:
    islands: pd.DataFrame
    promoters: pd.DataFrame
    cohort_tracks: Mapping[str, MethylationTrack]
    gene_symbols: Mapping[str, str] | None = None
    tsg_set: set[str] | None = None
    cosmic_tier: Mapping[str, float] | None = None


def _mean_beta_for_interval(track: MethylationTrack, start: int, end: int) -> float:
    if len(track) == 0:
        return float("nan")
    mask = (track.positions >= start) & (track.positions < end)
    if not mask.any():
        return float("nan")
    return float(track.betas[mask].mean())


def annotate_islands(
    islands: pd.DataFrame,
    cohort_tracks: Mapping[str, MethylationTrack],
) -> pd.DataFrame:
    out = islands.copy()
    for cohort, track in cohort_tracks.items():
        col = f"mean_beta_{cohort}"
        out[col] = [
            _mean_beta_for_interval(track, int(r["start"]), int(r["end"]))
            for _, r in out.iterrows()
        ]
    return out


def build_report(inputs: ReportInputs) -> pd.DataFrame:
    islands_with_meth = annotate_islands(inputs.islands, inputs.cohort_tracks)

    joined = overlap_islands_with_promoters(
        islands_with_meth.rename(columns={"start": "start", "end": "end"}),
        inputs.promoters,
    )
    if joined.empty:
        return joined

    cohort_cols = [f"mean_beta_{c}" for c in inputs.cohort_tracks]
    joined = joined.merge(
        islands_with_meth[["island_id", *cohort_cols]],
        on="island_id",
        how="left",
    )

    # Cohorts named "normal*" are healthy-tissue baselines; all others are tumor.
    normal_cols = [c for c in cohort_cols if c.startswith("mean_beta_normal")]
    tumor_cols = [c for c in cohort_cols if c not in normal_cols]

    consolidated_cols = tumor_cols if tumor_cols else cohort_cols
    if consolidated_cols:
        joined["max_mean_beta"] = joined[consolidated_cols].max(axis=1, skipna=True)
        # idxmax blows up on all-NA rows; compute per-row so empty-probe islands get "".
        has_any = joined[consolidated_cols].notna().any(axis=1)
        max_cohort = pd.Series("", index=joined.index, dtype=object)
        if has_any.any():
            max_cohort.loc[has_any] = (
                joined.loc[has_any, consolidated_cols]
                .idxmax(axis=1)
                .str.replace("mean_beta_", "", regex=False)
            )
        joined["max_cohort"] = max_cohort
        joined["hypermethylated"] = joined["max_mean_beta"].fillna(0) > HYPERMETH_THRESHOLD
    else:
        joined["max_mean_beta"] = np.nan
        joined["max_cohort"] = ""
        joined["hypermethylated"] = False

    if normal_cols:
        joined["mean_beta_normal"] = joined[normal_cols].mean(axis=1, skipna=True)
        joined["delta_vs_normal"] = (
            joined["max_mean_beta"].fillna(0) - joined["mean_beta_normal"].fillna(0)
        )
        joined["cancer_specific"] = (
            (joined["mean_beta_normal"].fillna(1) < 0.3)
            & (joined["max_mean_beta"].fillna(0) > 0.6)
        )
    else:
        joined["mean_beta_normal"] = np.nan
        joined["delta_vs_normal"] = np.nan
        joined["cancer_specific"] = False

    if inputs.gene_symbols:
        joined["gene_symbol"] = (
            joined["gene_id"].map(inputs.gene_symbols).fillna(joined["gene_id"])
        )
    else:
        joined["gene_symbol"] = joined["gene_id"]

    joined["is_tsg"] = (
        joined["gene_symbol"].isin(inputs.tsg_set) if inputs.tsg_set else False
    )

    if inputs.cosmic_tier:
        joined["cosmic_confidence"] = (
            joined["gene_symbol"].map(inputs.cosmic_tier).fillna(0.0)
        )
    else:
        joined["cosmic_confidence"] = 0.0

    joined["priority"] = (
        joined["max_mean_beta"].fillna(0.0) * joined["cosmic_confidence"]
    )

    if "canonical" in inputs.promoters.columns:
        canonicals = inputs.promoters.set_index("isoform_id")["canonical"]
        joined["canonical"] = joined["isoform_id"].map(canonicals).fillna(False)
    else:
        joined["canonical"] = False

    final_cols = [
        "island_id",
        "chrom",
        "island_start",
        "island_end",
        "gene_id",
        "gene_symbol",
        "isoform_id",
        "canonical",
        "strand",
        "tss",
        "overlap_bp",
        *cohort_cols,
        "max_cohort",
        "max_mean_beta",
        "mean_beta_normal",
        "delta_vs_normal",
        "cancer_specific",
        "hypermethylated",
        "is_tsg",
        "cosmic_confidence",
        "priority",
    ]
    final_cols = [c for c in final_cols if c in joined.columns]
    return (
        joined[final_cols]
        .sort_values(
            ["priority", "max_mean_beta", "overlap_bp"],
            ascending=[False, False, False],
        )
        .reset_index(drop=True)
    )
