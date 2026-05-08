"""Download UCSC hg38 FASTA + CpG-island BED + knownGene JSON per chromosome."""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path

import requests
from tqdm import tqdm

RAW = Path(__file__).resolve().parents[3] / "data" / "raw"


def chrom_sources(chrom: str) -> dict[str, dict]:
    return {
        f"{chrom}.fa": {
            "url": (
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/"
                f"{chrom}.fa.gz"
            ),
            "gzip": True,
        },
        f"cpg_islands_{chrom}.bed": {
            "url": (
                "https://api.genome.ucsc.edu/getData/track"
                f"?genome=hg38;track=cpgIslandExt;chrom={chrom}"
            ),
            "gzip": False,
        },
        f"knownGene_{chrom}.json": {
            "url": (
                "https://api.genome.ucsc.edu/getData/track"
                f"?genome=hg38;track=knownGene;chrom={chrom}"
            ),
            "gzip": False,
        },
    }


def _stream_download(url: str, dest: Path) -> None:
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with (
            open(dest, "wb") as f,
            tqdm(total=total, unit="B", unit_scale=True, desc=dest.name) as bar,
        ):
            for chunk in r.iter_content(chunk_size=1 << 16):
                f.write(chunk)
                bar.update(len(chunk))


def _gunzip(src: Path, dest: Path) -> None:
    with gzip.open(src, "rb") as gz, open(dest, "wb") as out:
        shutil.copyfileobj(gz, out)
    src.unlink()


def download_chrom(chrom: str, force: bool = False) -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, info in chrom_sources(chrom).items():
        dest = RAW / name
        if dest.exists() and not force:
            print(f"  [skip] {name} already exists")
            continue
        url = info["url"]
        if info.get("gzip"):
            gz_path = dest.with_suffix(dest.suffix + ".gz")
            _stream_download(url, gz_path)
            _gunzip(gz_path, dest)
        else:
            _stream_download(url, dest)
        print(f"  [ok]   {name}")


def manual_instructions() -> None:
    print(
        "\n"
        "Manual downloads required:\n"
        "  1. TCGA WGBS beta values (open access)\n"
        "     https://portal.gdc.cancer.gov/repository\n"
        "     Filters: Data Category = DNA Methylation,\n"
        "              Data Type = Methylation Beta Value,\n"
        "              Project = TCGA-BRCA, TCGA-LUAD, TCGA-COAD,\n"
        "              Access = open\n"
        "     Place files under: backend/data/raw/tcga/<cohort>/\n"
        "\n"
        "  2. COSMIC Cancer Gene Census (free registration)\n"
        "     https://cancer.sanger.ac.uk/census\n"
        "     Download the CSV, place at: backend/data/raw/cosmic/cancer_gene_census.csv\n"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default="chr21",
                    help="comma-separated list of chromosomes to download")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--quiet-instructions", action="store_true")
    args = ap.parse_args()

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    print(f"Downloading reference data to {RAW}")
    print(f"Chromosomes: {', '.join(chroms)}")
    try:
        for chrom in chroms:
            print(f"\n--- {chrom} ---")
            download_chrom(chrom, force=args.force)
    except requests.HTTPError as e:
        print(f"Download failed: {e}", file=sys.stderr)
        return 1
    if not args.quiet_instructions:
        manual_instructions()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
