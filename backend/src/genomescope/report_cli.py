"""End-to-end isoform methylation report driver."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .analysis.evaluate import load_ucsc_track_json
from .analysis.promoters import build_promoters
from .analysis.report import ReportInputs, build_report
from .data.annotations import load_cosmic_census, load_known_genes
from .data.fasta import load_fasta
from .data.symbols import fetch_symbols
from .data.tcga import load_tcga_cohorts, synthetic_track_from_fasta
from .hmm.model import standard_cpg_hmm
from .hmm.viterbi import merge_and_filter, segments, viterbi
from .pipeline import (
    BOLD,
    CYAN,
    DIM,
    GREEN,
    OUT,
    RAW,
    RESET,
    YELLOW,
    _section,
    _stat,
    encode_sequence,
)

TCGA_ROOT = RAW / "tcga"
COSMIC_PATH = RAW / "cosmic" / "cancer_gene_census.csv"


def _call_islands_standard_hmm(
    seq_slice: str, chrom: str, window_start: int,
) -> pd.DataFrame:
    """Call CpG islands with the standard 8-state nucleotide HMM (Beta HMM would
    only flag the unmethylated subset; the report overlays methylation downstream)."""
    obs = encode_sequence(seq_slice)
    hmm = standard_cpg_hmm()
    vit = viterbi(hmm, obs)
    island_mask = np.array([i >= 4 for i in range(hmm.n_states)])
    raw = segments(vit.path, island_mask)
    segs = merge_and_filter(raw, min_length=200, merge_gap=100)
    return pd.DataFrame(
        [
            {
                "island_id": f"isl_{i:04d}",
                "chrom": chrom,
                "start": s + window_start,
                "end": e + window_start,
            }
            for i, (s, e) in enumerate(segs)
        ]
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=Path, default=RAW / "chr21.fa")
    ap.add_argument("--knownGene", type=Path, default=RAW / "knownGene_chr21.json")
    ap.add_argument("--truth", type=Path, default=RAW / "cpg_islands_chr21.bed")
    ap.add_argument("--cosmic", type=Path, default=COSMIC_PATH)
    ap.add_argument("--offset", type=int, default=13_500_000)
    ap.add_argument("--subset", type=int, default=2_000_000)
    ap.add_argument("--max-iter", type=int, default=15)
    ap.add_argument("--flank", type=int, default=2_000, help="promoter ± flank in bp")
    ap.add_argument("--out", type=Path, default=OUT / "isoform_report.csv")
    ap.add_argument("--out-json", type=Path, default=OUT / "isoform_report.json")
    args = ap.parse_args()

    _section("Loading reference data")
    seq = load_fasta(args.fasta)
    _stat("fasta", f"{seq.chrom}  {seq.length:,} bp")
    start = max(0, args.offset)
    end = min(seq.length, start + args.subset)
    _stat("window", f"[{start:,} - {end:,})")
    seq_slice = seq.seq[start:end]

    genes = load_known_genes(args.knownGene, chrom=seq.chrom)
    _stat("transcripts", f"{len(genes):,}  ({genes['canonical'].sum():,} canonical)")

    promoters = build_promoters(genes, flank=args.flank)

    _section("Loading TCGA cohort methylation")
    cohorts = load_tcga_cohorts(TCGA_ROOT, chrom=seq.chrom)
    if not cohorts:
        _stat("source", f"{YELLOW}synthetic{RESET} (no TCGA data present)")
        truth_all = (
            load_ucsc_track_json(args.truth, chrom=seq.chrom)
            if args.truth.exists() else []
        )
        synthetic = synthetic_track_from_fasta(seq, truth_all, seed=7)
        cohorts = {"synthetic": synthetic}
    else:
        for name, track in cohorts.items():
            _stat(f"cohort {name}", f"{len(track):,} CpG probes")

    _section("Calling islands (standard 8-state nucleotide HMM)")
    islands = _call_islands_standard_hmm(seq_slice, seq.chrom, start)
    _stat("islands called", f"{len(islands):,}")
    primary_name = next(iter(cohorts))

    _section("ENSG → symbol enrichment (MyGene.info, cached)")
    gene_ids = genes["gene_id"].unique().tolist()
    ensg_only = [g for g in gene_ids if g.startswith("ENSG")]
    sym_map = fetch_symbols(ensg_only) if ensg_only else {}
    for g in gene_ids:
        sym_map.setdefault(g, g)
    n_already = sum(1 for g in gene_ids if not g.startswith("ENSG"))
    n_resolved = sum(
        1 for g in gene_ids
        if g.startswith("ENSG") and sym_map.get(g.split(".")[0], "").startswith(("A", "B", "C")) or
           (g.startswith("ENSG") and not sym_map.get(g.split(".")[0], g).startswith("ENSG"))
    )
    _stat("already symbols", f"{n_already:,} / {len(gene_ids):,}")
    _stat("ensg resolved", f"{n_resolved:,} / {len(ensg_only):,}")

    _section("COSMIC Cancer Gene Census")
    cosmic = load_cosmic_census(args.cosmic)
    using_fallback = cosmic is not None and not args.cosmic.exists()
    if cosmic is None:
        _stat("source", f"{YELLOW}missing{RESET} ({args.cosmic})")
        print(f"  {DIM}drop cancer_gene_census.csv at that path to enable full TSG scoring{RESET}")
        tsg_set = None
        conf = None
    elif using_fallback:
        _stat("source", f"{YELLOW}built-in fallback{RESET} (no CSV at {args.cosmic})")
        _stat("tsg genes", f"{len(cosmic):,}  (curated subset - drop full COSMIC CSV for the complete census)")
        tsg_set = cosmic.tsg_symbols
        conf = cosmic.confidence
    else:
        _stat("source", f"{GREEN}COSMIC CSV{RESET}")
        _stat("tsg genes", f"{len(cosmic):,}")
        tsg_set = cosmic.tsg_symbols
        conf = cosmic.confidence

    _section("Building report")
    report = build_report(
        ReportInputs(
            islands=islands,
            promoters=promoters,
            cohort_tracks=cohorts,
            gene_symbols=sym_map,
            tsg_set=tsg_set,
            cosmic_tier=conf,
        )
    )
    _stat("report rows", f"{len(report):,}")
    if not report.empty:
        hyper = int(report["hypermethylated"].sum())
        tsg_hits = int(report["is_tsg"].sum()) if "is_tsg" in report else 0
        _stat("hypermethylated", f"{hyper:,}")
        _stat("tsg hits", f"{tsg_hits:,}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.out, index=False)
    args.out_json.write_text(
        json.dumps(
            {
                "window": {"chrom": seq.chrom, "start": start, "end": end},
                "primary_cohort": primary_name,
                "cohorts": list(cohorts.keys()),
                "n_islands": len(islands),
                "n_report_rows": len(report),
                "rows": report.head(50).to_dict(orient="records"),
            },
            indent=2,
            default=str,
        )
    )
    _section("Saved")
    _stat("csv", str(args.out))
    _stat("json (top 50)", str(args.out_json))

    if not report.empty:
        _section("Top 10 by priority")
        top_cols = [
            c for c in [
                "gene_symbol", "isoform_id", "canonical",
                "max_cohort", "max_mean_beta", "hypermethylated", "is_tsg", "priority",
            ]
            if c in report.columns
        ]
        print(report[top_cols].head(10).to_string(index=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
