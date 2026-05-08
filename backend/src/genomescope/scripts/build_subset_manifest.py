"""TCGA subset manifest builder (filter source manifest or query GDC), SHA-256 sampled."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import requests

GDC_FILES = "https://api.gdc.cancer.gov/files"
BATCH = 200
QUERY_PAGE = 500


def fetch_metadata(uuids: list[str]) -> list[dict]:
    """Look up GDC metadata for an explicit list of file UUIDs."""
    out: list[dict] = []
    for i in range(0, len(uuids), BATCH):
        batch = uuids[i : i + BATCH]
        r = requests.post(
            GDC_FILES,
            json={
                "filters": {
                    "op": "in",
                    "content": {"field": "file_id", "value": batch},
                },
                "fields": "file_id,file_name,md5sum,file_size,"
                          "cases.project.project_id,cases.samples.sample_type,"
                          "platform",
                "size": BATCH,
            },
            timeout=60,
        )
        r.raise_for_status()
        out.extend(r.json()["data"]["hits"])
        print(f"  {len(out)}/{len(uuids)} fetched", end="\r", file=sys.stderr)
    print(file=sys.stderr)
    return out


def query_gdc(
    cohort: str,
    platform: str,
    sample_type: str | None = None,
    data_type: str = "Methylation Beta Value",
    access: str = "open",
) -> list[dict]:
    """Page through GDC /files for everything matching the project + platform filters."""
    filters: list[dict] = [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": [cohort]}},
        {"op": "in", "content": {"field": "platform", "value": [platform]}},
        {"op": "in", "content": {"field": "data_type", "value": [data_type]}},
        {"op": "in", "content": {"field": "access", "value": [access]}},
    ]
    if sample_type is not None:
        filters.append({
            "op": "in",
            "content": {"field": "cases.samples.sample_type", "value": [sample_type]},
        })

    body = {
        "filters": {"op": "and", "content": filters},
        "fields": "file_id,file_name,md5sum,file_size,"
                  "cases.project.project_id,cases.samples.sample_type,platform",
        "size": QUERY_PAGE,
        "from": 0,
    }

    out: list[dict] = []
    while True:
        r = requests.post(GDC_FILES, json=body, timeout=60)
        r.raise_for_status()
        page = r.json()["data"]["hits"]
        out.extend(page)
        print(f"  {len(out)} matches so far", end="\r", file=sys.stderr)
        if len(page) < QUERY_PAGE:
            break
        body["from"] += QUERY_PAGE
    print(file=sys.stderr)
    return out


def deterministic_sample(hits: list[dict], n: int) -> list[dict]:
    """Stable N-element subset by sha256(file_id)."""
    def sort_key(h: dict) -> bytes:
        return hashlib.sha256(h["file_id"].encode()).digest()

    return sorted(hits, key=sort_key)[:n]


def write_outputs(
    chosen: list[dict],
    cohort: str,
    platform: str,
    sample_type: str | None,
    out: Path,
    metadata: Path,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = ["id\tfilename\tmd5\tsize\tstate"]
    for h in chosen:
        lines.append(
            f"{h['file_id']}\t{h['file_name']}\t{h['md5sum']}\t{h['file_size']}\treleased"
        )
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote subset manifest → {out}", file=sys.stderr)

    metadata.parent.mkdir(parents=True, exist_ok=True)
    metadata.write_text(
        json.dumps(
            {
                "cohort": cohort,
                "platform": platform,
                "sample_type": sample_type,
                "files": [
                    {
                        "file_id": h["file_id"],
                        "file_name": h["file_name"],
                        "file_size": h["file_size"],
                    }
                    for h in chosen
                ],
            },
            indent=2,
        )
    )
    print(f"wrote metadata       → {metadata}", file=sys.stderr)

    total_mb = sum(h["file_size"] for h in chosen) / 1024 / 1024
    print(f"total download size: {total_mb:.0f} MB", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path, nargs="?",
                    help="path to existing GDC manifest (mode A). "
                         "omit when using --query-gdc (mode B).")
    ap.add_argument("--cohort", default="TCGA-BRCA")
    ap.add_argument("--platform", default="Illumina Human Methylation 450")
    ap.add_argument("--sample-type", default=None,
                    help="e.g. 'Solid Tissue Normal' or 'Primary Tumor'")
    ap.add_argument("--query-gdc", action="store_true",
                    help="fetch matching files directly from GDC (mode B)")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--out", type=Path,
                    default=Path("backend/data/raw/tcga/subset.manifest.txt"))
    ap.add_argument("--metadata", type=Path,
                    default=Path("backend/data/raw/tcga/subset.metadata.json"))
    args = ap.parse_args()

    if args.query_gdc:
        print(f"querying GDC for {args.cohort} + {args.platform}"
              + (f" + sample_type={args.sample_type}" if args.sample_type else ""),
              file=sys.stderr)
        hits = query_gdc(
            cohort=args.cohort,
            platform=args.platform,
            sample_type=args.sample_type,
        )
        if not hits:
            print(f"no matches for {args.cohort} + {args.platform}", file=sys.stderr)
            return 1
        print(f"{len(hits)} files match the filters", file=sys.stderr)
    else:
        if args.source is None:
            ap.error("either provide a source manifest or pass --query-gdc")
        entries = [
            line.split()
            for line in args.source.read_text().splitlines()[1:]
            if line.strip()
        ]
        uuids = [row[0] for row in entries]
        print(f"manifest has {len(uuids)} files total", file=sys.stderr)
        meta = fetch_metadata(uuids)
        def matches(h: dict) -> bool:
            if h["cases"][0]["project"]["project_id"] != args.cohort:
                return False
            if h.get("platform") != args.platform:
                return False
            if args.sample_type is not None:
                sample_types = [
                    s.get("sample_type")
                    for case in h.get("cases", [])
                    for s in case.get("samples", []) or []
                ]
                if args.sample_type not in sample_types:
                    return False
            return True

        hits = [h for h in meta if matches(h)]
        print(f"{len(hits)} files match the filters", file=sys.stderr)

    chosen = deterministic_sample(hits, args.n)
    write_outputs(
        chosen,
        cohort=args.cohort,
        platform=args.platform,
        sample_type=args.sample_type,
        out=args.out,
        metadata=args.metadata,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
