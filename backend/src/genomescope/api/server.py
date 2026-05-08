"""FastAPI server: reads cached pipeline + report JSON for the frontend dashboard."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from ..hmm.checkpoint import load_checkpoint, read_registry
from ..hmm.model import standard_cpg_hmm
from ..hmm.viterbi import merge_and_filter, segments, viterbi
from ..pipeline import encode_sequence

ROOT = Path(__file__).resolve().parents[3]
PROCESSED = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"
MODELS = ROOT / "data" / "models"

MAX_FASTA_BYTES = 5 * 1024 * 1024
VALID_NUC = set("ACGTNacgtn\n\r\t ")

app = FastAPI(title="GenomeScope API", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_methods=["*"],
    allow_headers=["*"],
)


def _read_json(path: Path) -> dict:
    if not path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"{path.name} not found - run the corresponding ./run.sh target first",
        )
    return json.loads(path.read_text())


@app.get("/api/health")
def health() -> dict:
    return {
        "status": "ok",
        "processed_files": sorted(p.name for p in PROCESSED.glob("*.json")),
    }


@app.get("/api/pipeline/standard")
def standard() -> dict:
    """Most-comprehensive standard-HMM run available (full chr21 → trained → bench)."""
    for name in (
        "full_chr21_standard.json",
        "trained.json",
        "quick.json",
        "smoke.json",
        "bench_standard.json",
    ):
        path = PROCESSED / name
        if path.exists():
            return _read_json(path)
    raise HTTPException(404, "no standard pipeline output found")


@app.get("/api/pipeline/beta")
def beta() -> dict:
    for name in ("beta.json", "bench_beta.json"):
        path = PROCESSED / name
        if path.exists():
            return _read_json(path)
    raise HTTPException(404, "no beta pipeline output found")


@app.get("/api/bench")
def bench() -> dict:
    return _read_json(PROCESSED / "bench_comparison.json")


@app.get("/api/report")
def report() -> dict:
    for name in ("isoform_report_full.json", "isoform_report.json"):
        path = PROCESSED / name
        if path.exists():
            return _read_json(path)
    raise HTTPException(404, "no isoform report found - run ./run.sh report-full")


@app.get("/api/full-genome")
def full_genome() -> dict:
    return _read_json(PROCESSED / "full_genome_summary.json")


@app.get("/api/ensemble")
def ensemble() -> dict:
    return _read_json(PROCESSED / "ensemble.json")


@app.get("/api/models")
def models() -> dict:
    entries = read_registry(MODELS)
    return {
        "n_models": len(entries),
        "models": [
            {
                "name": e.name,
                "model_type": e.model_type,
                "created": e.created,
                "updated": e.updated,
                "total_runs": e.total_runs,
                "total_samples_seen": e.total_samples_seen,
                "cohorts_seen": e.cohorts_seen,
                "last_log_likelihood": e.last_log_likelihood,
            }
            for e in entries
        ],
    }


class PredictRequest(BaseModel):
    fasta: str = Field(..., description="raw FASTA text or just a sequence string")
    model: str = Field("standard", description="'standard' or a registered checkpoint name")
    chrom: str = Field("user", description="label for the input chromosome")


def _parse_fasta_text(text: str) -> tuple[str, str]:
    """(chrom, sequence) from FASTA text or a bare sequence."""
    text = text.strip()
    if not text:
        raise HTTPException(400, "empty input")
    if len(text.encode("utf-8")) > MAX_FASTA_BYTES:
        raise HTTPException(413, f"input exceeds {MAX_FASTA_BYTES // 1024 // 1024} MB limit")

    chrom = "user"
    if text.startswith(">"):
        head, _, body = text.partition("\n")
        chrom = head.lstrip(">").split()[0] or "user"
        text = body
    seq = "".join(c for c in text if c not in "\n\r\t ").upper()
    if not seq:
        raise HTTPException(400, "no sequence content found")
    invalid = set(seq) - set("ACGTN")
    if invalid:
        raise HTTPException(
            400,
            f"sequence contains non-nucleotide characters: {sorted(invalid)[:5]}",
        )
    return chrom, seq


@app.post("/api/predict")
def predict(req: PredictRequest) -> dict:
    """Run the standard nucleotide HMM on user-supplied FASTA."""
    chrom, seq = _parse_fasta_text(req.fasta)

    t0 = time.perf_counter()
    obs = encode_sequence(seq)

    used_checkpoint: str | None = None
    if req.model == "standard":
        hmm = standard_cpg_hmm()
    else:
        entries = read_registry(MODELS)
        match = next((e for e in entries if e.name == req.model), None)
        if match is None:
            raise HTTPException(404, f"unknown model '{req.model}'")
        if match.model_type != "standard":
            raise HTTPException(
                400,
                f"checkpoint '{req.model}' is type '{match.model_type}' - "
                "POST /predict only supports nucleotide HMMs",
            )
        ckpt_path = (MODELS / match.path).resolve()
        loaded, _ = load_checkpoint(ckpt_path)
        hmm = loaded  # type: ignore[assignment]
        used_checkpoint = match.name

    vit = viterbi(hmm, obs)
    island_mask = np.array([i >= 4 for i in range(hmm.n_states)])
    raw = segments(vit.path, island_mask)
    islands = merge_and_filter(raw, min_length=200, merge_gap=100)
    elapsed = time.perf_counter() - t0

    return {
        "chrom": chrom,
        "model": req.model,
        "checkpoint_used": used_checkpoint,
        "window_start": 0,
        "window_end": len(seq),
        "window_length": len(seq),
        "n_islands": len(islands),
        "islands": [{"start": int(s), "end": int(e)} for s, e in islands],
        "total_island_bp": int(sum(e - s for s, e in islands)),
        "viterbi_log_prob": float(vit.log_prob),
        "elapsed_seconds": round(elapsed, 3),
    }


@app.get("/api/truth")
def truth() -> dict:
    """UCSC reference CpG islands for chr21 (track-viewer baseline)."""
    path = RAW / "cpg_islands_chr21.bed"
    if not path.exists():
        raise HTTPException(404, "ground truth BED not present - run ./run.sh data")
    raw = json.loads(path.read_text())
    entries: list[dict] = []
    for v in raw.values():
        if isinstance(v, list) and v and isinstance(v[0], dict):
            entries = v
            break
    islands = [
        {"start": int(e["chromStart"]), "end": int(e["chromEnd"])}
        for e in entries
        if e.get("chrom") == "chr21"
    ]
    return {"chrom": "chr21", "n_islands": len(islands), "islands": islands}
