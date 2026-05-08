"""Standard ∩ Beta HMM ensemble: keep Standard calls that overlap any Beta call."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ..analysis.evaluate import load_ucsc_track_json, position_level_metrics

ROOT = Path(__file__).resolve().parents[3]
PROCESSED = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"

BOLD = "\033[1m"
DIM = "\033[2m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
CYAN = "\033[36m"
RED = "\033[31m"
RESET = "\033[0m"


def _section(title: str) -> None:
    print(f"\n{BOLD}{CYAN}━━━ {title} ━━━{RESET}")


def _stat(label: str, value: str) -> None:
    print(f"  {DIM}{label:<22}{RESET} {value}")


def filter_by_overlap(
    a: list[tuple[int, int]],
    b: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Half-open intervals from `a` that overlap any interval in `b` (linear sweep)."""
    if not a or not b:
        return []
    a_sorted = sorted(a)
    b_sorted = sorted(b)
    out: list[tuple[int, int]] = []
    j = 0
    for ax_start, ax_end in a_sorted:
        while j < len(b_sorted) and b_sorted[j][1] <= ax_start:
            j += 1
        for k in range(j, len(b_sorted)):
            bx_start, bx_end = b_sorted[k]
            if bx_start >= ax_end:
                break
            out.append((ax_start, ax_end))
            break
    return out


def _to_local(islands: list[dict], window_start: int) -> list[tuple[int, int]]:
    return [(int(i["start"]) - window_start, int(i["end"]) - window_start)
            for i in islands]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--standard", type=Path, default=PROCESSED / "full_chr21_standard.json")
    ap.add_argument("--beta", type=Path, default=PROCESSED / "full_chr21_beta.json")
    ap.add_argument("--truth", type=Path, default=RAW / "cpg_islands_chr21.bed")
    ap.add_argument("--out", type=Path, default=PROCESSED / "ensemble.json")
    args = ap.parse_args()

    if not args.standard.exists():
        print(f"{RED}!{RESET} missing {args.standard.name} - run ./run.sh full-chr21-standard")
        return 1
    if not args.beta.exists():
        print(f"{RED}!{RESET} missing {args.beta.name} - run ./run.sh full-chr21-beta")
        return 1

    std = json.loads(args.standard.read_text())
    beta = json.loads(args.beta.read_text())

    if std.get("chrom") != beta.get("chrom"):
        print(f"chrom mismatch: standard={std.get('chrom')}  beta={beta.get('chrom')}")
        return 1
    chrom = std["chrom"]
    start = int(std["window_start"])
    end = int(std["window_end"])
    window_len = end - start

    std_genomic = [(int(i["start"]), int(i["end"])) for i in std["islands"]]
    beta_genomic = [(int(i["start"]), int(i["end"])) for i in beta["islands"]]
    ensemble_genomic = filter_by_overlap(std_genomic, beta_genomic)

    _section(f"Ensemble - {chrom}:{start:,}-{end:,}  ({window_len / 1_000_000:.1f} Mb)")
    _stat("standard islands", f"{len(std_genomic):,}")
    _stat("beta islands", f"{len(beta_genomic):,}")
    _stat("ensemble (intersection)", f"{len(ensemble_genomic):,}")

    truth_all = load_ucsc_track_json(args.truth, chrom=chrom)
    truth_local: list[tuple[int, int]] = []
    for s, e in truth_all:
        if e <= start or s >= end:
            continue
        truth_local.append((max(0, s - start), min(window_len, e - start)))
    _stat("ucsc truth (in window)", f"{len(truth_local):,}")

    std_local = _to_local(std["islands"], start)
    beta_local = _to_local(beta["islands"], start)
    ensemble_local = [(s - start, e - start) for s, e in ensemble_genomic]

    m_std = position_level_metrics(std_local, truth_local, window_len)
    m_beta = position_level_metrics(beta_local, truth_local, window_len)
    m_ens = position_level_metrics(ensemble_local, truth_local, window_len)

    _section("3-way benchmark vs UCSC truth")
    print(
        f"  {DIM}{'model':<14}{'islands':>10}{'precision':>12}{'recall':>10}{'F1':>10}{RESET}"
    )
    for name, m, n in [
        ("standard", m_std, len(std_genomic)),
        ("beta",     m_beta, len(beta_genomic)),
        ("ensemble", m_ens,  len(ensemble_genomic)),
    ]:
        f1_colour = (
            GREEN if m.f1 >= max(m_std.f1, m_beta.f1) + 0.02
            else (YELLOW if m.f1 < min(m_std.f1, m_beta.f1) else RESET)
        )
        print(
            f"  {name:<14}{n:>10,}"
            f"{m.precision:>12.3f}{m.recall:>10.3f}"
            f"  {f1_colour}{m.f1:>7.3f}{RESET}"
        )

    _section("In one sentence")
    print()
    delta_p = m_ens.precision - m_std.precision
    delta_r = m_ens.recall - m_std.recall
    delta_f1 = m_ens.f1 - m_std.f1
    direction = "improves" if delta_f1 > 0.02 else (
        "stays roughly the same as" if abs(delta_f1) <= 0.02 else "underperforms"
    )
    print(
        f"  Filtering the Standard HMM's calls through the Beta HMM "
        f"{direction} F1\n"
        f"  vs the Standard HMM alone "
        f"({m_std.f1:.3f} → {m_ens.f1:.3f}, "
        f"Δ {delta_f1:+.3f}).\n"
    )
    print(f"  Precision: {m_std.precision:.3f} → {m_ens.precision:.3f}  "
          f"(Δ {delta_p:+.3f})")
    print(f"  Recall:    {m_std.recall:.3f} → {m_ens.recall:.3f}  "
          f"(Δ {delta_r:+.3f})")
    print()
    print(f"  {DIM}islands kept: {len(ensemble_genomic):,} of {len(std_genomic):,} "
          f"standard calls = "
          f"{100 * len(ensemble_genomic) / max(1, len(std_genomic)):.1f}%{RESET}")

    payload = {
        "chrom": chrom,
        "window_start": start,
        "window_end": end,
        "window_length": window_len,
        "n_standard": len(std_genomic),
        "n_beta": len(beta_genomic),
        "n_ensemble": len(ensemble_genomic),
        "ensemble_islands": [{"start": s, "end": e} for s, e in ensemble_genomic],
        "standard_metrics": m_std.as_dict(),
        "beta_metrics": m_beta.as_dict(),
        "ensemble_metrics": m_ens.as_dict(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2))
    print(f"\n  {DIM}saved → {args.out.relative_to(ROOT.parent)}{RESET}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
