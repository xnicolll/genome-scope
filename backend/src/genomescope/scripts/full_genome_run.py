"""Run the standard HMM across multiple chromosomes and aggregate metrics."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"

BOLD = "\033[1m"
DIM = "\033[2m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
CYAN = "\033[36m"
RESET = "\033[0m"


def find_assembled_region(fasta_path: Path) -> tuple[int, int, int]:
    """(start, end, total) bounds of the contiguous non-N region."""
    from ..data.fasta import load_fasta
    seq = load_fasta(fasta_path)
    n_array = np.frombuffer(seq.seq.encode("ascii"), dtype=np.uint8)
    not_n = n_array != ord("N")
    if not not_n.any():
        return 0, 0, len(seq.seq)
    start = int(np.argmax(not_n))
    end = len(seq.seq) - int(np.argmax(not_n[::-1]))
    return start, end, len(seq.seq)


def run_one_chromosome(chrom: str) -> dict:
    print(f"\n{BOLD}{CYAN}━━━ {chrom} ━━━{RESET}")
    fasta = RAW / f"{chrom}.fa"
    truth = RAW / f"cpg_islands_{chrom}.bed"
    if not fasta.exists():
        print(f"  {YELLOW}!{RESET} {fasta.name} missing - skipping")
        return {"chrom": chrom, "skipped": "fasta missing"}
    if not truth.exists():
        print(f"  {YELLOW}!{RESET} {truth.name} missing - skipping")
        return {"chrom": chrom, "skipped": "truth missing"}

    start, end, total = find_assembled_region(fasta)
    assembled = end - start
    print(f"  assembled region: [{start:,} - {end:,})  ({assembled:,} bp)")

    out_path = PROCESSED / f"full_chr_{chrom}.json"
    cmd = [
        sys.executable, "-m", "genomescope.pipeline",
        "--chrom", chrom,
        "--offset", str(start),
        "--subset", str(assembled),
        "--no-train",
        "--out", str(out_path),
    ]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        print(f"  {YELLOW}!{RESET} pipeline failed for {chrom}:")
        print(proc.stderr[-2000:])
        return {"chrom": chrom, "skipped": f"pipeline failed (rc={proc.returncode})"}

    payload = json.loads(out_path.read_text())
    n_islands = len(payload["islands"])
    metrics = payload.get("metrics") or {}
    print(
        f"  islands: {n_islands:,}  "
        f"precision: {metrics.get('precision', 0):.3f}  "
        f"recall: {metrics.get('recall', 0):.3f}  "
        f"F1: {metrics.get('f1', 0):.3f}  "
        f"({elapsed:.1f}s)"
    )
    return {
        "chrom": chrom,
        "assembled_bp": assembled,
        "total_bp": total,
        "n_islands": n_islands,
        "total_island_bp": payload.get("total_island_bp", 0),
        "metrics": metrics,
        "elapsed_s": round(elapsed, 1),
    }


def aggregate(results: list[dict]) -> None:
    valid = [r for r in results if "skipped" not in r]
    skipped = [r for r in results if "skipped" in r]

    summary = {
        "n_chromosomes_attempted": len(results),
        "n_chromosomes_completed": len(valid),
        "skipped": skipped,
        "per_chromosome": valid,
    }
    if valid:
        total_assembled = sum(r["assembled_bp"] for r in valid)
        total_islands = sum(r["n_islands"] for r in valid)
        tp = sum(r["metrics"].get("tp", 0) for r in valid)
        fp = sum(r["metrics"].get("fp", 0) for r in valid)
        fn = sum(r["metrics"].get("fn", 0) for r in valid)
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
        summary["aggregate"] = {
            "total_assembled_bp": total_assembled,
            "total_islands": total_islands,
            "micro_precision": round(precision, 4),
            "micro_recall": round(recall, 4),
            "micro_f1": round(f1, 4),
        }

    out = PROCESSED / "full_genome_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))

    print(f"\n{BOLD}{CYAN}━━━ Full-genome summary ━━━{RESET}")
    if not valid:
        print(f"  {YELLOW}!{RESET} no chromosomes completed - see skipped list in {out}")
        return
    agg = summary["aggregate"]
    print()
    print(f"  Across {len(valid)} chromosomes covering "
          f"{agg['total_assembled_bp'] / 1_000_000:.0f} Mb of human DNA:")
    print()
    print(f"  {BOLD}{agg['total_islands']:,}{RESET} CpG islands flagged")
    print(f"  {BOLD}{agg['micro_precision']:.3f}{RESET} micro-averaged precision  "
          f"(when the model says 'island', it's right this fraction of the time)")
    print(f"  {BOLD}{agg['micro_recall']:.3f}{RESET} micro-averaged recall  "
          f"(the model catches this fraction of the real islands)")
    print(f"  {BOLD}{agg['micro_f1']:.3f}{RESET} F1 score")
    print()
    print(f"  {DIM}per-chromosome breakdown:{RESET}")
    print(f"  {DIM}{'chrom':<8}{'assembled':>15}{'islands':>10}{'F1':>8}{'time':>8}{RESET}")
    for r in valid:
        f1 = r["metrics"].get("f1", 0)
        print(
            f"  {r['chrom']:<8}"
            f"{r['assembled_bp'] / 1_000_000:>12.1f} Mb"
            f"{r['n_islands']:>10,}"
            f"{f1:>8.3f}"
            f"{r['elapsed_s']:>7.1f}s"
        )
    if skipped:
        print()
        print(f"  {YELLOW}skipped:{RESET}")
        for r in skipped:
            print(f"    {r['chrom']}  ({r['skipped']})")
    print()
    print(f"  {DIM}saved → {out.relative_to(ROOT.parent)}{RESET}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--chroms",
        default="chr18,chr19,chr20,chr22,chrY",
        help="comma-separated chromosomes to run",
    )
    args = ap.parse_args()

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    print(f"{BOLD}Running standard HMM across:{RESET} {', '.join(chroms)}")

    results = [run_one_chromosome(c) for c in chroms]
    aggregate(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
