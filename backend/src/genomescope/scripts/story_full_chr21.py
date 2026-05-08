"""Plain-English summary comparing the full-chr21 run against the 2 Mb sample."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PROCESSED = ROOT / "data" / "processed"

BOLD = "\033[1m"
DIM = "\033[2m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
CYAN = "\033[36m"
RESET = "\033[0m"


def _section(title: str) -> None:
    print(f"\n{BOLD}{CYAN}━━━ {title} ━━━{RESET}")


def _bullet(label: str, body: str) -> None:
    print(f"  {BOLD}{label}{RESET}  {body}")


def _load(name: str) -> dict | None:
    path = PROCESSED / name
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _format_metrics(m: dict | None) -> tuple[str, str, str]:
    if not m:
        return ("-", "-", "-")
    return (
        f"{m['precision']:.3f}",
        f"{m['recall']:.3f}",
        f"{m['f1']:.3f}",
    )


def main() -> int:
    full = _load("full_chr21_standard.json")
    sample = _load("quick.json") or _load("smoke.json")

    if not full:
        print(
            f"{YELLOW}!{RESET} full_chr21_standard.json not found - "
            "run ./run.sh full-chr21-standard first"
        )
        return 1

    n_full = len(full["islands"])
    n_full_bp = full.get("total_island_bp", 0)
    win_full = full["window_length"]
    p_full, r_full, f1_full = _format_metrics(full.get("metrics"))

    _section("Story mode - full chr21 standard HMM")
    print()
    print(f"  Across all {win_full / 1_000_000:.1f} Mb of chr21's real DNA:")
    print()
    _bullet("Islands found:", f"{n_full:,}")
    _bullet("Total island DNA:", f"{n_full_bp:,} bp ({100 * n_full_bp / win_full:.2f}% of the chromosome)")
    _bullet("Precision:", f"{p_full}  (when the model says 'island', it's right this fraction of the time)")
    _bullet("Recall:", f"{r_full}  (the model catches this fraction of the real islands)")
    _bullet("F1 score:", f"{f1_full}  (combined accuracy - higher is better)")

    if sample:
        n_sample = len(sample["islands"])
        win_sample = sample.get("window_length", 0)
        p_s, r_s, f1_s = _format_metrics(sample.get("metrics"))
        _section("How the bigger run compares to the 2 Mb sample")
        print()
        ratio = win_full / max(1, win_sample)
        print(
            f"  We tested on {ratio:.0f}× more DNA than the small {win_sample / 1_000_000:.1f} Mb sample."
        )
        print()
        print(f"  {DIM}{'metric':<14}{'2 Mb sample':>16}{'full chr21':>16}{RESET}")
        print(f"  {'islands':<14}{n_sample:>16,}{n_full:>16,}")
        print(f"  {'precision':<14}{p_s:>16}{p_full:>16}")
        print(f"  {'recall':<14}{r_s:>16}{r_full:>16}")
        print(f"  {'F1':<14}{f1_s:>16}{f1_full:>16}")
        print()
        if abs(float(f1_s) - float(f1_full)) < 0.05:
            print(
                f"  {GREEN}✓{RESET} The numbers are nearly identical - "
                "the small-window result wasn't a fluke."
            )
        else:
            print(
                f"  {YELLOW}!{RESET} The numbers shifted at scale - "
                "this matters for the writeup."
            )

    _section("In one sentence")
    print()
    print(
        f"  The standard HMM, run on the full {win_full / 1_000_000:.0f} Mb of chr21 "
        f"with no training,\n"
        f"  flagged {n_full:,} candidate CpG islands across the chromosome - "
        f"catching {r_full}\n"
        f"  of every real island while being right {p_full} of the time it speaks up."
    )
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
