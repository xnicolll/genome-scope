"""End-to-end pipeline: FASTA → HMM (standard or Beta) → predicted islands → F1."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from .analysis.evaluate import load_ucsc_track_json, position_level_metrics
from .data.fasta import Sequence, load_fasta
from .data.tcga import (
    MethylationTrack,
    load_tcga_directory,
    synthetic_track_from_fasta,
)
from .hmm.baum_welch import baum_welch
from .hmm.beta_model import BetaHMM, beta_cpg_hmm
from .hmm.checkpoint import (
    TrainingRun,
    load_checkpoint,
    save_checkpoint,
    upsert_registry,
    utc_now,
)
from .hmm.model import HMM, standard_cpg_hmm
from .hmm.viterbi import merge_and_filter, segments, viterbi
from .training_log import TrainingLogEntry, append as append_training_log

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
TCGA_DIR = RAW / "tcga" / "brca_hm450"
OUT = ROOT / "data" / "processed"
MODELS = ROOT / "data" / "models"

BOLD = "\033[1m"
DIM = "\033[2m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
CYAN = "\033[36m"
RESET = "\033[0m"


def _section(title: str) -> None:
    print(f"\n{BOLD}{CYAN}▸ {title}{RESET}")


def _stat(label: str, value: str) -> None:
    print(f"  {DIM}{label:<22}{RESET} {value}")


def encode_sequence(seq: str, seed: int = 42) -> np.ndarray:
    """A/C/G/T → 0/1/2/3 (N replaced with a random nucleotide)."""
    byte_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    lookup = np.full(256, -1, dtype=np.int64)
    lookup[ord("A")] = 0
    lookup[ord("C")] = 1
    lookup[ord("G")] = 2
    lookup[ord("T")] = 3
    encoded = lookup[byte_arr]
    mask = encoded == -1
    if mask.any():
        rng = np.random.default_rng(seed)
        encoded[mask] = rng.integers(0, 4, size=int(mask.sum()))
    return encoded.astype(np.int64)


def _benchmark(
    island_segments: list[tuple[int, int]],
    truth_all: list[tuple[int, int]],
    start: int,
    end: int,
    window_len: int,
) -> dict | None:
    _section("Benchmark vs UCSC ground truth")
    truth: list[tuple[int, int]] = []
    for s, e in truth_all:
        if e <= start or s >= end:
            continue
        truth.append((max(0, s - start), min(window_len, e - start)))
    m = position_level_metrics(island_segments, truth, window_len)
    metrics = m.as_dict()
    _stat("ground truth islands", f"{len(truth):,}")
    _stat("precision", f"{m.precision:.4f}")
    _stat("recall", f"{m.recall:.4f}")
    _stat("f1", f"{m.f1:.4f}")
    badge = f"{GREEN}PASS{RESET}" if m.f1 >= 0.7 else f"{YELLOW}below target{RESET}"
    _stat("target f1 > 0.7", badge)
    return metrics


def _save(out_path: Path, payload: dict) -> None:
    out_path = out_path.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2))
    _section("Saved")
    try:
        rel = out_path.relative_to(ROOT.parent)
    except ValueError:
        rel = out_path
    _stat("output", str(rel))


def run_beta_pipeline(
    fasta_path: Path,
    subset: int | None,
    offset: int,
    max_iter: int,
    train: bool,
    out_path: Path,
    truth_path: Path | None,
    resume_from: Path | None = None,
    checkpoint_out: Path | None = None,
    checkpoint_name: str | None = None,
    cohort: str | None = None,
    tcga_dir: Path | None = None,
) -> dict:
    _section("Loading chr21 FASTA")
    t0 = time.perf_counter()
    seq = load_fasta(fasta_path)
    _stat("chromosome", seq.chrom)
    _stat("length", f"{seq.length:,} bp")
    _stat("elapsed", f"{time.perf_counter() - t0:.2f}s")

    start = max(0, min(offset, seq.length))
    end = seq.length if subset is None else min(seq.length, start + subset)
    _stat("window", f"[{start:,} - {end:,})  ({end - start:,} bp)")

    _section("Methylation track")
    truth_all: list[tuple[int, int]] = []
    if truth_path is not None and truth_path.exists():
        truth_all = load_ucsc_track_json(truth_path, chrom=seq.chrom)

    cohort_dir = tcga_dir if tcga_dir is not None else TCGA_DIR
    track = load_tcga_directory(cohort_dir, chrom=seq.chrom)
    synthetic = track is None
    if synthetic:
        _stat("source", f"{YELLOW}synthetic{RESET} (no TCGA data in {cohort_dir})")
        print(f"  {YELLOW}! F1 on synthetic data is a closed-loop sanity check -{RESET}")
        print(f"  {YELLOW}  the track is generated from the UCSC BED we benchmark against{RESET}")
        print(f"  {YELLOW}  so near-perfect scores are expected. Drop real TCGA files{RESET}")
        print(f"  {YELLOW}  into backend/data/raw/tcga/ for a non-circular comparison.{RESET}")
        track = synthetic_track_from_fasta(seq, truth_all, seed=7)
    else:
        _stat("source", f"{GREEN}TCGA{RESET}  ({len(track):,} CpG sites)")
    track = track.window(start, end)
    _stat("CpGs in window", f"{len(track):,}")
    if len(track) == 0:
        print("  (no CpG sites in window - aborting)")
        return {}

    _section("Beta 2-state HMM")
    prior_history: list[TrainingRun] = []
    if resume_from is not None and resume_from.exists():
        loaded, prior_history = load_checkpoint(resume_from)
        if not isinstance(loaded, BetaHMM):
            raise TypeError(f"{resume_from} is not a Beta HMM checkpoint")
        hmm = loaded
        _stat("init", f"{GREEN}resumed{RESET} from {resume_from.name}")
        _stat("prior runs", f"{len(prior_history)}")
    else:
        hmm = beta_cpg_hmm()
        _stat("init", "fresh (PRD defaults)")
    means = hmm.state_means()
    _stat("init depleted mean", f"{means[0]:.3f}")
    _stat("init island mean",   f"{means[1]:.3f}")

    ll_trace: list[float] = []
    this_run: TrainingRun | None = None
    if train:
        _section(f"Baum-Welch training ({max_iter} max iters)")
        t0 = time.perf_counter()
        result = baum_welch(
            hmm, track.betas, max_iter=max_iter, tol=1e-4,
            freeze_emissions=False, verbose=True,
        )
        hmm = result.hmm
        duration = time.perf_counter() - t0
        _stat("iterations", str(result.iterations))
        _stat("converged", "yes" if result.converged else "no (hit max_iter)")
        _stat("final log-lik", f"{result.log_likelihoods[-1]:.4f}")
        _stat("elapsed", f"{duration:.2f}s")
        ll_trace = result.log_likelihoods
        post_means = hmm.state_means()
        _stat("trained depleted mean", f"{post_means[0]:.3f}")
        _stat("trained island mean",   f"{post_means[1]:.3f}")

        this_run = TrainingRun(
            timestamp=utc_now(),
            cohort=cohort,
            n_samples=int(len(track)),
            window_start=int(start),
            window_end=int(end),
            iterations=int(result.iterations),
            converged=bool(result.converged),
            log_likelihood_start=float(result.log_likelihoods[0]),
            log_likelihood_end=float(result.log_likelihoods[-1]),
            max_iter=int(max_iter),
            tol=1e-4,
        )

        append_training_log(
            TrainingLogEntry(
                model_name=checkpoint_name or "anonymous",
                model_type="beta",
                cohort=cohort,
                window_start=int(start),
                window_end=int(end),
                n_samples=int(len(track)),
                iterations=int(result.iterations),
                converged=bool(result.converged),
                ll_start=float(result.log_likelihoods[0]),
                ll_end=float(result.log_likelihoods[-1]),
                state_means=tuple(float(m) for m in post_means),
                duration_s=float(duration),
                note="synthetic" if synthetic else "",
            )
        )
    else:
        _stat("training", "skipped (--no-train)")

    if checkpoint_out is not None:
        new_history = list(prior_history) + ([this_run] if this_run else [])
        co = Path(checkpoint_out)
        if not co.is_absolute() and co.parent == Path("."):
            co = MODELS / co.name
        saved = save_checkpoint(hmm, co, history=new_history).resolve()
        _section("Checkpoint")
        _stat("saved", str(saved))
        _stat("runs in history", str(len(new_history)))
        if checkpoint_name:
            MODELS.mkdir(parents=True, exist_ok=True)
            entry = upsert_registry(MODELS, checkpoint_name, saved, new_history)
            _stat("registry", f"{entry.name} · {len(entry.cohorts_seen)} cohorts")

    _section("Viterbi decoding")
    t0 = time.perf_counter()
    vit = viterbi(hmm, track.betas)
    # UCSC islands are unmethylated (mean β≈0.1); island state = lower-mean state.
    state_means = hmm.state_means()
    island_state = int(np.argmin(state_means))
    _stat("island state", f"S{island_state}  (mean={state_means[island_state]:.3f})")
    is_island = vit.path == island_state

    raw_segs: list[tuple[int, int]] = []
    in_seg = False
    seg_start = 0
    last_pos = 0
    for pos, flag in zip(track.positions, is_island, strict=True):
        pos = int(pos)
        if flag and not in_seg:
            seg_start = pos
            in_seg = True
        elif not flag and in_seg:
            raw_segs.append((seg_start, last_pos + 2))
            in_seg = False
        last_pos = pos
    if in_seg:
        raw_segs.append((seg_start, last_pos + 2))
    local_segs = [(s - start, e - start) for s, e in raw_segs]
    island_segments = merge_and_filter(local_segs, min_length=200, merge_gap=100)
    total_bp = sum(e - s for s, e in island_segments)

    _stat("raw segments", f"{len(raw_segs):,}")
    _stat("islands (merged)", f"{len(island_segments):,}")
    _stat("total island bp", f"{total_bp:,} ({100 * total_bp / (end - start):.2f}%)")
    _stat("best-path log-p", f"{vit.log_prob:.4f}")
    _stat("elapsed", f"{time.perf_counter() - t0:.2f}s")

    metrics = _benchmark(island_segments, truth_all, start, end, end - start)

    payload = {
        "model": "beta",
        "chrom": seq.chrom,
        "window_start": int(start),
        "window_end": int(end),
        "window_length": int(end - start),
        "cpg_sites": int(len(track)),
        "islands": [
            {"start": int(s) + start, "end": int(e) + start} for s, e in island_segments
        ],
        "total_island_bp": int(total_bp),
        "log_likelihoods": [float(x) for x in ll_trace],
        "viterbi_log_prob": float(vit.log_prob),
        "state_means": hmm.state_means().tolist(),
        "metrics": metrics,
    }
    _save(out_path, payload)
    return payload


def run_pipeline(
    fasta_path: Path,
    subset: int | None,
    offset: int,
    max_iter: int,
    train: bool,
    out_path: Path,
    truth_path: Path | None,
) -> dict:
    _section("Loading chr21 FASTA")
    t0 = time.perf_counter()
    seq = load_fasta(fasta_path)
    _stat("chromosome", seq.chrom)
    _stat("length", f"{seq.length:,} bp")
    _stat("elapsed", f"{time.perf_counter() - t0:.2f}s")

    start = max(0, min(offset, seq.length))
    end = seq.length if subset is None else min(seq.length, start + subset)
    seq_slice = seq.seq[start:end]
    if start or end != seq.length:
        _stat("window", f"[{start:,} - {end:,})  ({end - start:,} bp)")

    _section("Encoding sequence")
    t0 = time.perf_counter()
    obs = encode_sequence(seq_slice)
    n_mask = int((np.frombuffer(seq_slice.encode("ascii"), dtype=np.uint8) == ord("N")).sum())
    _stat("encoded length", f"{obs.shape[0]:,}")
    _stat("N positions", f"{n_mask:,} (filled w/ random)")
    _stat("elapsed", f"{time.perf_counter() - t0:.2f}s")

    _section("Standard 8-state CpG HMM")
    hmm = standard_cpg_hmm()
    _stat("states", str(hmm.n_states))
    _stat("emissions", "deterministic (nucleotide)")

    if train:
        _section(f"Baum-Welch training ({max_iter} max iters)")
        t0 = time.perf_counter()
        result = baum_welch(
            hmm, obs, max_iter=max_iter, tol=1e-4, freeze_emissions=True, verbose=True
        )
        hmm = result.hmm
        _stat("iterations", str(result.iterations))
        _stat("converged", "yes" if result.converged else "no (hit max_iter)")
        _stat("final log-lik", f"{result.log_likelihoods[-1]:.4f}")
        _stat("elapsed", f"{time.perf_counter() - t0:.2f}s")
        ll_trace = result.log_likelihoods
    else:
        _stat("training", "skipped (--no-train)")
        ll_trace = []

    _section("Viterbi decoding")
    t0 = time.perf_counter()
    vit = viterbi(hmm, obs)
    island_mask = np.array([i >= 4 for i in range(hmm.n_states)])
    raw_segments = segments(vit.path, island_mask)
    island_segments = merge_and_filter(raw_segments, min_length=200, merge_gap=100)
    total_bp = sum(e - s for s, e in island_segments)
    _stat("raw segments", f"{len(raw_segments):,}")
    _stat("islands (merged)", f"{len(island_segments):,}")
    _stat("total island bp", f"{total_bp:,} ({100 * total_bp / obs.shape[0]:.2f}%)")
    _stat("best-path log-p", f"{vit.log_prob:.4f}")
    _stat("elapsed", f"{time.perf_counter() - t0:.2f}s")

    metrics = None
    if truth_path is not None and truth_path.exists():
        truth_all = load_ucsc_track_json(truth_path, chrom=seq.chrom)
        metrics = _benchmark(island_segments, truth_all, start, end, obs.shape[0])
    else:
        _stat("benchmark", "skipped (no ground truth file)")

    payload = {
        "model": "standard",
        "chrom": seq.chrom,
        "window_start": int(start),
        "window_end": int(end),
        "window_length": int(obs.shape[0]),
        "subset_used": subset,
        "islands": [
            {"start": int(s) + start, "end": int(e) + start} for s, e in island_segments
        ],
        "total_island_bp": int(total_bp),
        "log_likelihoods": [float(x) for x in ll_trace],
        "viterbi_log_prob": float(vit.log_prob),
        "metrics": metrics,
    }
    _save(out_path, payload)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="GenomeScope HMM pipeline")
    parser.add_argument("--model", choices=["standard", "beta"], default="standard",
                        help="standard = 8-state nucleotide HMM; beta = 2-state Beta methylation HMM")
    parser.add_argument("--chrom", default=None,
                        help="convenience: derive --fasta + --truth from this chromosome name "
                             "(e.g. 'chr19' → data/raw/chr19.fa + data/raw/cpg_islands_chr19.bed)")
    parser.add_argument("--fasta", type=Path, default=None)
    parser.add_argument("--truth", type=Path, default=None)
    parser.add_argument("--subset", type=int, default=None,
                        help="decode only N bases starting at --offset")
    parser.add_argument("--offset", type=int, default=0,
                        help="skip this many bases from the chromosome start")
    parser.add_argument("--max-iter", type=int, default=10)
    parser.add_argument("--no-train", action="store_true",
                        help="skip Baum-Welch, use untrained HMM")
    parser.add_argument("--out", type=Path, default=OUT / "pipeline_run.json")
    parser.add_argument("--resume-from", type=Path, default=None,
                        help="load a prior checkpoint as Baum-Welch starting point")
    parser.add_argument("--checkpoint-out", type=Path, default=None,
                        help="save the trained model + training history to this path")
    parser.add_argument("--checkpoint-name", type=str, default=None,
                        help="register the saved checkpoint under this name")
    parser.add_argument("--cohort", type=str, default=None,
                        help="tag this run with a cohort label (e.g. brca, luad)")
    parser.add_argument("--tcga-dir", type=Path, default=None,
                        help="override the default TCGA directory (e.g. data/raw/tcga/luad_hm450)")
    args = parser.parse_args()

    chrom = args.chrom or "chr21"
    if args.fasta is None:
        args.fasta = RAW / f"{chrom}.fa"
    if args.truth is None:
        args.truth = RAW / f"cpg_islands_{chrom}.bed"

    if not args.fasta.exists():
        print(f"{YELLOW}!{RESET} FASTA not found at {args.fasta}")
        print(f"  run: {BOLD}./run.sh data{RESET}  (or uv run python -m genomescope.scripts.download_data)")
        return 2

    if args.model == "beta":
        run_beta_pipeline(
            fasta_path=args.fasta,
            subset=args.subset,
            offset=args.offset,
            max_iter=args.max_iter,
            train=not args.no_train,
            out_path=args.out,
            truth_path=args.truth if args.truth.exists() else None,
            resume_from=args.resume_from,
            checkpoint_out=args.checkpoint_out,
            checkpoint_name=args.checkpoint_name,
            cohort=args.cohort,
            tcga_dir=args.tcga_dir,
        )
    else:
        if args.resume_from or args.checkpoint_out:
            print(f"{YELLOW}! checkpoint flags currently only supported with --model beta{RESET}")
        run_pipeline(
            fasta_path=args.fasta,
            subset=args.subset,
            offset=args.offset,
            max_iter=args.max_iter,
            train=not args.no_train,
            out_path=args.out,
            truth_path=args.truth if args.truth.exists() else None,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
