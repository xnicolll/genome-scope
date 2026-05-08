"""Side-by-side benchmark: standard HMM vs Beta-HMM."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

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
    run_beta_pipeline,
    run_pipeline,
)


def _fmt(v: float | None, width: int = 8) -> str:
    if v is None:
        return "-".rjust(width)
    return f"{v:.4f}".rjust(width)


def _delta(a: float, b: float) -> str:
    d = b - a
    colour = GREEN if d > 0 else (YELLOW if d == 0 else "\033[31m")
    sign = "+" if d >= 0 else ""
    return f"{colour}{sign}{d:.4f}{RESET}"


def print_comparison(std: dict, beta: dict) -> None:
    _section("Benchmark comparison")
    std_m = std.get("metrics") or {}
    beta_m = beta.get("metrics") or {}

    header = f"  {DIM}{'model':<12}{'precision':>10}{'recall':>10}{'f1':>10}{'islands':>10}{RESET}"
    print(header)
    print(f"  {'standard':<12}"
          f"{_fmt(std_m.get('precision'), 10)}"
          f"{_fmt(std_m.get('recall'), 10)}"
          f"{_fmt(std_m.get('f1'), 10)}"
          f"{len(std.get('islands', [])):>10}")
    print(f"  {'beta':<12}"
          f"{_fmt(beta_m.get('precision'), 10)}"
          f"{_fmt(beta_m.get('recall'), 10)}"
          f"{_fmt(beta_m.get('f1'), 10)}"
          f"{len(beta.get('islands', [])):>10}")

    if "f1" in std_m and "f1" in beta_m:
        d = beta_m["f1"] - std_m["f1"]
        print()
        print(f"  {BOLD}Δ F1{RESET}  (beta − standard)  {_delta(std_m['f1'], beta_m['f1'])}")
        if d >= 0.02:
            print(f"  {GREEN}✓ PRD target met (Beta-HMM F1 ≥ standard + 0.02){RESET}")
        else:
            print(f"  {YELLOW}! PRD target not met - negative result is also valid per PRD section 5{RESET}")


def main() -> int:
    parser = argparse.ArgumentParser(description="GenomeScope model comparison")
    parser.add_argument("--fasta", type=Path, default=RAW / "chr21.fa")
    parser.add_argument("--truth", type=Path, default=RAW / "cpg_islands_chr21.bed")
    parser.add_argument("--offset", type=int, default=13_500_000)
    parser.add_argument("--subset", type=int, default=2_000_000)
    parser.add_argument("--max-iter", type=int, default=10)
    parser.add_argument("--no-train-standard", action="store_true",
                        help="standard HMM uses Durbin init only (fast)")
    parser.add_argument("--out-dir", type=Path, default=OUT)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    std_out = args.out_dir / "bench_standard.json"
    beta_out = args.out_dir / "bench_beta.json"

    print(f"{BOLD}{CYAN}━━━ STANDARD 8-state CpG HMM ━━━{RESET}")
    std = run_pipeline(
        fasta_path=args.fasta,
        subset=args.subset,
        offset=args.offset,
        max_iter=args.max_iter,
        train=not args.no_train_standard,
        out_path=std_out,
        truth_path=args.truth if args.truth.exists() else None,
    )

    print(f"\n{BOLD}{CYAN}━━━ BETA 2-state methylation HMM ━━━{RESET}")
    beta = run_beta_pipeline(
        fasta_path=args.fasta,
        subset=args.subset,
        offset=args.offset,
        max_iter=args.max_iter,
        train=True,
        out_path=beta_out,
        truth_path=args.truth if args.truth.exists() else None,
    )

    if std and beta:
        print_comparison(std, beta)
        (args.out_dir / "bench_comparison.json").write_text(
            json.dumps({"standard": std, "beta": beta}, indent=2)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
