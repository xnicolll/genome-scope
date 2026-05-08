"""Illumina HM450 probe → genomic coordinate map (Zhou Lab hg38 manifest)."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

HM450_URL = (
    "https://github.com/zhou-lab/InfiniumAnnotationV1"
    "/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz"
)

ROOT = Path(__file__).resolve().parents[3]
CACHE_DIR = ROOT / "data" / "cache"
HM450_CACHED_GZ = CACHE_DIR / "HM450.hg38.manifest.tsv.gz"
HM450_CACHED_TSV = CACHE_DIR / "HM450.hg38.manifest.tsv"


def _stream_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=180) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with (
            open(dest, "wb") as f,
            tqdm(total=total, unit="B", unit_scale=True, desc=dest.name) as bar,
        ):
            for chunk in r.iter_content(chunk_size=1 << 16):
                f.write(chunk)
                bar.update(len(chunk))


def ensure_hm450_manifest() -> Path:
    if HM450_CACHED_TSV.exists():
        return HM450_CACHED_TSV
    print(f"fetching HM450 probe manifest → {HM450_CACHED_GZ}")
    _stream_download(HM450_URL, HM450_CACHED_GZ)
    with gzip.open(HM450_CACHED_GZ, "rb") as gz, open(HM450_CACHED_TSV, "wb") as out:
        shutil.copyfileobj(gz, out)
    HM450_CACHED_GZ.unlink()
    return HM450_CACHED_TSV


def load_hm450_probe_map(chrom: str | None = "chr21") -> pd.DataFrame:
    """probe_id-indexed DataFrame of (chrom, start, end), optionally chrom-filtered."""
    tsv = ensure_hm450_manifest()
    df = pd.read_csv(
        tsv,
        sep="\t",
        usecols=["CpG_chrm", "CpG_beg", "CpG_end", "Probe_ID"],
        dtype={"CpG_chrm": str, "CpG_beg": "Int64", "CpG_end": "Int64", "Probe_ID": str},
    )
    df = df.rename(
        columns={
            "CpG_chrm": "chrom",
            "CpG_beg": "start",
            "CpG_end": "end",
            "Probe_ID": "probe_id",
        }
    )
    df = df.dropna(subset=["chrom", "start", "probe_id"])
    if chrom is not None:
        df = df[df["chrom"] == chrom].copy()
    df["start"] = df["start"].astype("int64")
    df["end"] = df["end"].astype("int64")
    return df.set_index("probe_id").sort_values("start")
