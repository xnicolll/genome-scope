"""ENSG → HGNC symbol resolution via MyGene.info, disk-cached."""

from __future__ import annotations

import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "ensg_to_symbol.json"
MYGENE_URL = "https://mygene.info/v3/gene"
BATCH = 500


def _read_cache() -> dict[str, str]:
    if CACHE.exists():
        return json.loads(CACHE.read_text())
    return {}


def _write_cache(mapping: dict[str, str]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(mapping, indent=2, sort_keys=True))


def fetch_symbols(ensembl_ids: list[str]) -> dict[str, str]:
    """ENSG → HGNC symbol map (version suffixes stripped, unknowns self-map)."""
    cache = _read_cache()

    def norm(i: str) -> str:
        return i.split(".")[0]

    unknown = sorted({norm(i) for i in ensembl_ids if norm(i) not in cache})
    if unknown:
        for i in range(0, len(unknown), BATCH):
            batch = unknown[i : i + BATCH]
            try:
                r = requests.post(
                    f"{MYGENE_URL}/query",
                    data={
                        "q": ",".join(batch),
                        "scopes": "ensembl.gene",
                        "fields": "symbol",
                        "species": "human",
                    },
                    timeout=60,
                )
                r.raise_for_status()
                for hit in r.json():
                    ens = hit.get("query")
                    sym = hit.get("symbol")
                    if ens and sym:
                        cache[ens] = sym
            except requests.RequestException:
                break
        for ens in unknown:
            cache.setdefault(ens, ens)
        _write_cache(cache)

    return {norm(i): cache.get(norm(i), norm(i)) for i in ensembl_ids}
