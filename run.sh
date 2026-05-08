#!/usr/bin/env bash
#
# GenomeScope - one-shot local runner.
#
# Usage:
#   ./run.sh                # env + tests + data + quick pipeline (safe default)
#   ./run.sh env            # sync backend Python env only
#   ./run.sh test           # run the pytest suite
#   ./run.sh data           # download chr21 reference data (UCSC)
#   ./run.sh tcga           # download + ingest real TCGA-BRCA HM450 methylation
#   ./run.sh smoke          # 200kb standard HMM (fastest sanity check)
#   ./run.sh quick          # 2Mb standard HMM, Durbin-init only
#   ./run.sh train          # 2Mb standard HMM with Baum-Welch training
#   ./run.sh beta           # 2Mb Beta methylation HMM, fully trained
#   ./run.sh train-beta     # incremental training - resumes from prior checkpoint
#   ./run.sh daily          # day-of-year windowed training - for launchd
#   ./run.sh schedule-install   # install launchd agent (runs daily at 3am)
#   ./run.sh schedule-remove    # uninstall launchd agent
#   ./run.sh schedule-status    # is the agent loaded? show recent runs
#   ./run.sh bench          # standard vs beta HMM side-by-side F1 comparison
#   ./run.sh report         # phase 3: isoform-aware methylation + TSG report
#   ./run.sh up             # start FastAPI (8000) + Next.js (3000) - background
#   ./run.sh down           # stop both dev servers
#   ./run.sh logs api|web   # tail dev server logs
#   ./run.sh full           # full chr21 pipeline (slow - minutes + ~6GB RAM)
#   ./run.sh full-genome    # run standard HMM across 5 representative chroms
#   ./run.sh clean          # remove generated outputs (keeps raw data)
#
# Later phases (Beta emissions, TCGA overlay, Next.js frontend) will add
# their own targets to this file.

set -euo pipefail

# ─── paths ────────────────────────────────────────────────────────────────────
ROOT="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$ROOT/backend"
UV="${UV:-$HOME/.local/bin/uv}"
if ! [ -x "$UV" ] && command -v uv >/dev/null 2>&1; then
  UV="$(command -v uv)"
fi

# ─── pretty output ────────────────────────────────────────────────────────────
if [ -t 1 ]; then
  BOLD=$'\033[1m'; DIM=$'\033[2m'
  GREEN=$'\033[32m'; YELLOW=$'\033[33m'
  CYAN=$'\033[36m'; RED=$'\033[31m'; RESET=$'\033[0m'
else
  BOLD=""; DIM=""; GREEN=""; YELLOW=""; CYAN=""; RED=""; RESET=""
fi

section() { printf "\n%s%s━━━ %s ━━━%s\n" "$BOLD" "$CYAN" "$1" "$RESET"; }
say()     { printf "%s▸%s %s\n" "$CYAN"   "$RESET" "$1"; }
ok()      { printf "%s✓%s %s\n" "$GREEN"  "$RESET" "$1"; }
warn()    { printf "%s!%s %s\n" "$YELLOW" "$RESET" "$1"; }
err()     { printf "%s✗%s %s\n" "$RED"    "$RESET" "$1" >&2; }

# ─── guards ───────────────────────────────────────────────────────────────────
ensure_uv() {
  if ! [ -x "$UV" ]; then
    err "uv not found (checked $UV)"
    err "install with: curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
  fi
}

ensure_env() {
  section "python environment"
  say "syncing backend with uv"
  (cd "$BACKEND" && "$UV" sync --extra dev >/dev/null)
  ok "env ready"
}

# ─── steps ────────────────────────────────────────────────────────────────────
run_tests() {
  section "tests"
  (cd "$BACKEND" && "$UV" run --extra dev pytest -q)
  ok "all tests passed"
}

run_data() {
  section "reference data"
  say "downloading UCSC chr21 FASTA + CpG islands + knownGene"
  (cd "$BACKEND" && "$UV" run python -m genomescope.scripts.download_data)
  ok "automated downloads complete"
  warn "TCGA + COSMIC are manual - see instructions printed above"
}

#  chr21 p-arm is unassembled (all N) up to ~13.5Mb. Windowed runs start
#  at the first real assembled region with known CpG islands.
WINDOW_OFFSET=13500000

run_smoke() {
  section "pipeline - smoke (200kb, Durbin-init, no training)"
  (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
    --offset "$WINDOW_OFFSET" --subset 200000 --no-train \
    --out data/processed/smoke.json)
}

run_quick() {
  section "pipeline - quick (2Mb, Durbin-init, no training)"
  (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
    --offset "$WINDOW_OFFSET" --subset 2000000 --no-train \
    --out data/processed/quick.json)
}

run_train() {
  section "pipeline - trained (2Mb window, 10 BW iters)"
  warn "Baum-Welch on a small window can overfit - see TRAINING notes"
  (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
    --offset "$WINDOW_OFFSET" --subset 2000000 --max-iter 10 \
    --out data/processed/trained.json)
}

run_beta() {
  section "pipeline - beta methylation hmm (2Mb, trained)"
  (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
    --model beta --offset "$WINDOW_OFFSET" --subset 2000000 --max-iter 15 \
    --out data/processed/beta.json)
}

# Beta HMM on the full assembled chr21 window - needed for the
# standard ⨯ beta ensemble (./run.sh ensemble).
run_full_chr21_beta() {
  section "pipeline - beta methylation hmm on full chr21 (~33 Mb)"
  (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
    --model beta --offset 13500000 --subset 33200000 --max-iter 20 \
    --cohort brca \
    --out data/processed/full_chr21_beta.json)
  ok "full-chr21 beta run complete"
}

# Standard × Beta ensemble: filter the Standard HMM's calls by Beta-HMM
# overlap, then benchmark all 3 vs UCSC. Auto-runs the upstream pipelines
# if their JSONs are missing.
run_ensemble() {
  section "ensemble - standard × beta hmm on full chr21"

  if [ ! -f "$BACKEND/data/processed/full_chr21_standard.json" ]; then
    say "no standard run yet - running ./run.sh full-chr21-standard first"
    run_full_chr21_standard
  fi
  if [ ! -f "$BACKEND/data/processed/full_chr21_beta.json" ]; then
    say "no beta run yet - running ./run.sh full-chr21-beta first"
    run_full_chr21_beta
  fi

  (cd "$BACKEND" && "$UV" run python -m genomescope.scripts.ensemble)
}

# Train a Beta HMM checkpoint that learns over time. The first invocation
# starts from PRD defaults; every subsequent invocation resumes from the
# saved checkpoint, accumulating training history across windows + cohorts.
run_train_beta() {
  section "incremental beta training - checkpoint beta-brca.json"
  local resume_arg=""
  if [ -f "$BACKEND/data/models/beta-brca.json" ]; then
    resume_arg="--resume-from data/models/beta-brca.json"
    say "resuming from existing checkpoint"
  else
    say "fresh checkpoint (PRD defaults)"
  fi
  (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
    --model beta --offset "$WINDOW_OFFSET" --subset 2000000 --max-iter 15 \
    --cohort brca --checkpoint-out beta-brca.json --checkpoint-name beta-brca \
    $resume_arg --out data/processed/beta_trained.json)
}

# Daily incremental training across all available cohorts.
#
# Each invocation runs three back-to-back training sessions (or however
# many cohorts have data on disk):
#   - trains on the FULL assembled chr21 region (~33 Mb)
#   - cycles through every cohort: BRCA tumor → LUAD → BRCA normal
#   - resumes from the same shared checkpoint each step, so by end of day
#     the model has been updated by all three datasets in a row
#
# Cohort directories are checked at runtime - any cohort whose data
# folder is missing is silently skipped.
run_daily() {
  section "daily incremental training (all cohorts)"

  # Available cohorts (cohort label : directory under data/raw/tcga/)
  local pairs=(
    "brca:brca_hm450"
    "luad:luad_hm450"
    "normal_brca:normal_brca_hm450"
  )

  # Filter to cohorts whose data directory actually exists
  local available=()
  for p in "${pairs[@]}"; do
    local dir="${p#*:}"
    if [ -d "$BACKEND/data/raw/tcga/$dir" ] && \
       find "$BACKEND/data/raw/tcga/$dir" -name "*level3betas.txt" -print -quit | grep -q .; then
      available+=("$p")
    fi
  done

  if [ "${#available[@]}" -eq 0 ]; then
    err "no cohorts found under backend/data/raw/tcga/ - run ./run.sh tcga first"
    exit 1
  fi

  say "found ${#available[@]} cohort(s): $(printf '%s ' "${available[@]%%:*}")"
  say "training on full chr21 (~33 Mb) for each"

  local today
  today=$(date +%Y-%m-%d)
  local i=0
  for pair in "${available[@]}"; do
    i=$(( i + 1 ))
    local cohort="${pair%%:*}"
    local dir="${pair#*:}"

    section "[$i/${#available[@]}] cohort: $cohort"

    local resume_arg=""
    if [ -f "$BACKEND/data/models/beta-multi.json" ]; then
      resume_arg="--resume-from data/models/beta-multi.json"
    fi

    (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
      --model beta --offset 13500000 --subset 33200000 --max-iter 15 \
      --cohort "$cohort" --tcga-dir "data/raw/tcga/$dir" \
      --checkpoint-out beta-multi.json --checkpoint-name beta-multi \
      $resume_arg --out "data/processed/daily_${today}_${cohort}.json")
  done

  ok "all ${#available[@]} cohorts trained · appended to backend/data/logs/training-history.csv"
}

# Download + ingest real TCGA methylation data via the GDC manifest.
# Expects gdc-client binary at repo root and gdc_manifest.*.txt alongside.
TCGA_DIR="$BACKEND/data/raw/tcga"
TCGA_SUBSET="$TCGA_DIR/subset.manifest.txt"
TCGA_FILES="$TCGA_DIR/brca_hm450"
GDC_CLIENT="$ROOT/gdc-client"

run_tcga() {
  section "tcga - real methylation data ingestion"

  if ! [ -x "$GDC_CLIENT" ]; then
    err "gdc-client not found at $GDC_CLIENT"
    err "download it from https://gdc.cancer.gov/access-data/gdc-data-transfer-tool"
    exit 1
  fi

  local source_manifest
  source_manifest="$(ls -1t "$ROOT"/gdc_manifest.*.txt 2>/dev/null | head -1 || true)"
  if [ -z "$source_manifest" ]; then
    err "no gdc_manifest.*.txt found at repo root"
    exit 1
  fi
  say "source manifest: $(basename "$source_manifest")"

  mkdir -p "$TCGA_DIR" "$TCGA_FILES"

  if [ -f "$TCGA_SUBSET" ]; then
    ok "subset manifest already exists - reusing $TCGA_SUBSET"
  else
    say "querying GDC API to build 20-file TCGA-BRCA HM450 subset"
    (cd "$BACKEND" && "$UV" run python -m genomescope.scripts.build_subset_manifest \
      "$source_manifest" \
      --cohort "TCGA-BRCA" \
      --platform "Illumina Human Methylation 450" \
      --n 20 \
      --out "$TCGA_SUBSET" \
      --metadata "$TCGA_DIR/subset.metadata.json")
    ok "subset manifest written"
  fi

  say "downloading 20 TCGA-BRCA HM450 files (~250 MB) via gdc-client"
  "$GDC_CLIENT" download -m "$TCGA_SUBSET" -d "$TCGA_FILES" 2>&1 | tail -10 || true

  # gdc-client drops each file under $TCGA_FILES/<uuid>/<filename>.
  # Count how many .txt methylation beta files actually landed.
  local n_files
  n_files=$(find "$TCGA_FILES" -name "*.methylation_array*level3betas.txt" 2>/dev/null | wc -l | tr -d ' ')
  ok "downloaded $n_files level-3 beta files into $TCGA_FILES"

  section "tcga - parsing + caching probe map"
  (cd "$BACKEND" && "$UV" run python -c "
from genomescope.data.tcga import load_tcga_directory
from pathlib import Path
track = load_tcga_directory(Path('data/raw/tcga/brca_hm450'), chrom='chr21')
if track is None:
    raise SystemExit('no usable files found')
print(f'  chr21 CpG probes: {len(track):,}')
print(f'  mean beta:        {track.betas.mean():.3f}')
print(f'  first 5 positions:{track.positions[:5].tolist()}')
")
  ok "tcga ingestion complete - now run: ./run.sh bench"
}

# T-25: pull a fresh TCGA-LUAD subset directly from the GDC API (no
# source manifest required). Drops files into data/raw/tcga/luad_hm450/
# so the existing load_tcga_cohorts auto-discovery picks them up.
TCGA_LUAD_SUBSET="$TCGA_DIR/luad.manifest.txt"
TCGA_LUAD_FILES="$TCGA_DIR/luad_hm450"

run_tcga_luad() {
  section "tcga - TCGA-LUAD second cohort ingestion"

  if ! [ -x "$GDC_CLIENT" ]; then
    err "gdc-client not found at $GDC_CLIENT"
    exit 1
  fi

  mkdir -p "$TCGA_DIR" "$TCGA_LUAD_FILES"

  if [ -f "$TCGA_LUAD_SUBSET" ]; then
    ok "luad manifest already exists - reusing $TCGA_LUAD_SUBSET"
  else
    say "querying GDC API directly for 20 TCGA-LUAD HM450 files"
    (cd "$BACKEND" && "$UV" run python -m genomescope.scripts.build_subset_manifest \
      --query-gdc \
      --cohort TCGA-LUAD \
      --platform "Illumina Human Methylation 450" \
      --sample-type "Primary Tumor" \
      --n 20 \
      --out "$TCGA_LUAD_SUBSET" \
      --metadata "$TCGA_DIR/luad.metadata.json")
    ok "luad manifest written"
  fi

  say "downloading 20 TCGA-LUAD HM450 files via gdc-client"
  "$GDC_CLIENT" download -m "$TCGA_LUAD_SUBSET" -d "$TCGA_LUAD_FILES" 2>&1 | tail -10 || true

  local n_files
  n_files=$(find "$TCGA_LUAD_FILES" -name "*.methylation_array*level3betas.txt" 2>/dev/null | wc -l | tr -d ' ')
  ok "downloaded $n_files level-3 beta files into $TCGA_LUAD_FILES"

  section "tcga - parsing LUAD"
  (cd "$BACKEND" && "$UV" run python -c "
from genomescope.data.tcga import load_tcga_directory
from pathlib import Path
track = load_tcga_directory(Path('data/raw/tcga/luad_hm450'), chrom='chr21')
if track is None:
    raise SystemExit('no usable luad files found')
print(f'  chr21 luad CpG probes: {len(track):,}')
print(f'  mean beta:             {track.betas.mean():.3f}')
")
  ok "LUAD ingestion complete - re-run: ./run.sh report  (now multi-cohort)"
}

# T-24: pull TCGA-BRCA Solid Tissue Normal samples (healthy tissue from
# the same patients as the tumor cohort). Drops files into
# data/raw/tcga/normal_brca_hm450/ - the load_tcga_cohorts helper strips
# the _hm450 suffix so the cohort key becomes "normal_brca", and the
# report module treats any cohort starting with "normal" as a healthy
# baseline (computes delta_vs_normal + cancer_specific flags).
TCGA_NORMAL_SUBSET="$TCGA_DIR/normal_brca.manifest.txt"
TCGA_NORMAL_FILES="$TCGA_DIR/normal_brca_hm450"

run_tcga_normal() {
  section "tcga - TCGA-BRCA solid tissue normal baseline"

  if ! [ -x "$GDC_CLIENT" ]; then
    err "gdc-client not found at $GDC_CLIENT"
    exit 1
  fi

  mkdir -p "$TCGA_DIR" "$TCGA_NORMAL_FILES"

  if [ -f "$TCGA_NORMAL_SUBSET" ]; then
    ok "normal manifest already exists - reusing $TCGA_NORMAL_SUBSET"
  else
    say "querying GDC for 10 TCGA-BRCA HM450 *Solid Tissue Normal* files"
    (cd "$BACKEND" && "$UV" run python -m genomescope.scripts.build_subset_manifest \
      --query-gdc \
      --cohort TCGA-BRCA \
      --platform "Illumina Human Methylation 450" \
      --sample-type "Solid Tissue Normal" \
      --n 10 \
      --out "$TCGA_NORMAL_SUBSET" \
      --metadata "$TCGA_DIR/normal_brca.metadata.json")
    ok "normal manifest written"
  fi

  say "downloading 10 TCGA-BRCA solid-tissue-normal HM450 files"
  "$GDC_CLIENT" download -m "$TCGA_NORMAL_SUBSET" -d "$TCGA_NORMAL_FILES" 2>&1 | tail -10 || true

  local n_files
  n_files=$(find "$TCGA_NORMAL_FILES" -name "*.methylation_array*level3betas.txt" 2>/dev/null | wc -l | tr -d ' ')
  ok "downloaded $n_files normal-tissue beta files into $TCGA_NORMAL_FILES"

  section "tcga - parsing BRCA normal"
  (cd "$BACKEND" && "$UV" run python -c "
from genomescope.data.tcga import load_tcga_directory
from pathlib import Path
track = load_tcga_directory(Path('data/raw/tcga/normal_brca_hm450'), chrom='chr21')
if track is None:
    raise SystemExit('no usable normal files found')
print(f'  chr21 normal CpG probes: {len(track):,}')
print(f'  mean beta:               {track.betas.mean():.3f}')
")
  ok "normal-tissue ingestion complete - re-run: ./run.sh report-full"
}

run_bench() {
  section "benchmark - standard vs beta hmm"
  (cd "$BACKEND" && "$UV" run python -m genomescope.benchmark \
    --offset "$WINDOW_OFFSET" --subset 2000000 --max-iter 10 \
    --no-train-standard)
}

run_report() {
  section "isoform methylation report (phase 3 pipeline)"
  (cd "$BACKEND" && "$UV" run python -m genomescope.report_cli \
    --offset "$WINDOW_OFFSET" --subset 2000000 --max-iter 15)
}

# Full-chr21 isoform-aware methylation report. Combines T-22 + T-23 - runs
# the standard HMM across all 33 Mb of chr21, joins the resulting islands
# with promoters + TCGA-BRCA methylation, and ranks hits using either the
# COSMIC CSV (if present) or the built-in TSG fallback.
run_report_full() {
  section "full chr21 isoform methylation report"
  (cd "$BACKEND" && "$UV" run python -m genomescope.report_cli \
    --offset 13500000 --subset 33200000 --max-iter 15 \
    --out data/processed/isoform_report_full.csv \
    --out-json data/processed/isoform_report_full.json)
}

run_full() {
  section "pipeline - full chr21"
  warn "this can take several minutes and ~6GB RAM"
  (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
    --max-iter 10 --out data/processed/full.json)
}

# T-22: Full chr21 standard HMM run, Durbin-init-only (no Baum-Welch).
# Spans the entire assembled region (offset 13.5M, length 33.2M) so the
# benchmark is across all real chr21 sequence rather than a 2 Mb sample.
run_full_chr21_standard() {
  section "pipeline - full chr21 (standard HMM, Durbin init only)"
  say "running over the full ~33 Mb assembled region of chr21"
  (cd "$BACKEND" && "$UV" run python -m genomescope.pipeline \
    --offset 13500000 --subset 33200000 --no-train \
    --out data/processed/full_chr21_standard.json)
  ok "full-chr21 standard HMM run complete"

  # Plain-English story-mode summary
  (cd "$BACKEND" && "$UV" run python -m genomescope.scripts.story_full_chr21)
}

# T-26: Run the standard HMM across 5 representative chromosomes, beyond
# the chr21 scope. Picks chr18/19/20/22/Y because their assembled sizes
# stay under ~8 GB peak memory for Viterbi on this hardware. Larger
# chromosomes (chr1-3) need a chunked Viterbi which is out of scope.
FULL_GENOME_CHROMS="${FULL_GENOME_CHROMS:-chr18,chr19,chr20,chr22,chrY}"

run_full_genome() {
  section "full-genome scaling (5 representative chromosomes)"

  # Download FASTA + CpG BED + knownGene for every chromosome we plan to run
  say "ensuring reference data is downloaded for: $FULL_GENOME_CHROMS"
  (cd "$BACKEND" && "$UV" run python -m genomescope.scripts.download_data \
    --chroms "$FULL_GENOME_CHROMS" --quiet-instructions)
  ok "reference data ready"

  # Run the HMM on each chromosome and aggregate results
  (cd "$BACKEND" && "$UV" run python -m genomescope.scripts.full_genome_run \
    --chroms "$FULL_GENOME_CHROMS")
}

run_clean() {
  section "clean"
  rm -rf "$BACKEND/data/processed"
  ok "removed backend/data/processed"
}

# ─── launchd scheduling (daily training on macOS) ──────────────────────────────
PLIST_LABEL="com.genomescope.daily"
PLIST_TEMPLATE="$ROOT/scripts/genomescope.daily.plist.template"
PLIST_INSTALLED="$HOME/Library/LaunchAgents/${PLIST_LABEL}.plist"
SCHEDULE_HOUR="${GENOMESCOPE_HOUR:-3}"
SCHEDULE_MINUTE="${GENOMESCOPE_MINUTE:-0}"

run_schedule_install() {
  section "schedule - install launchd agent"
  if [ ! -f "$PLIST_TEMPLATE" ]; then
    err "template missing: $PLIST_TEMPLATE"
    exit 1
  fi
  mkdir -p "$LOGS" "$HOME/Library/LaunchAgents"

  # Render template → installed plist
  sed \
    -e "s|{ROOT}|$ROOT|g" \
    -e "s|{HOME}|$HOME|g" \
    -e "s|{LOGS}|$LOGS|g" \
    -e "s|{HOUR}|$SCHEDULE_HOUR|g" \
    -e "s|{MINUTE}|$SCHEDULE_MINUTE|g" \
    "$PLIST_TEMPLATE" > "$PLIST_INSTALLED"
  ok "wrote $PLIST_INSTALLED"

  # Bootout if already loaded (idempotent reinstall)
  launchctl bootout "gui/$(id -u)/$PLIST_LABEL" 2>/dev/null || true

  launchctl bootstrap "gui/$(id -u)" "$PLIST_INSTALLED"
  ok "loaded agent · will run at ${SCHEDULE_HOUR}:${SCHEDULE_MINUTE} daily"

  say "launchd output logs: $LOGS/launchd.{out,err}.log"
  say "training history:    $BACKEND/data/logs/training-history.csv"
  printf "\n  inspect status: %slaunchctl print gui/\$UID/%s%s\n" "$BOLD" "$PLIST_LABEL" "$RESET"
  printf "  trigger now:    %slaunchctl kickstart -k gui/\$UID/%s%s\n\n" "$BOLD" "$PLIST_LABEL" "$RESET"
}

run_schedule_remove() {
  section "schedule - remove launchd agent"
  launchctl bootout "gui/$(id -u)/$PLIST_LABEL" 2>/dev/null && \
    ok "agent unloaded" || warn "agent was not loaded"
  rm -f "$PLIST_INSTALLED"
  ok "removed $PLIST_INSTALLED"
}

run_schedule_status() {
  section "schedule - status"
  if launchctl print "gui/$(id -u)/$PLIST_LABEL" >/dev/null 2>&1; then
    ok "loaded (runs ${SCHEDULE_HOUR}:${SCHEDULE_MINUTE} daily)"
  else
    warn "not loaded · run ./run.sh schedule-install"
  fi
  if [ -f "$BACKEND/data/logs/training-history.csv" ]; then
    local n_runs
    n_runs=$(($(wc -l < "$BACKEND/data/logs/training-history.csv") - 1))
    say "training-history.csv · $n_runs runs recorded"
    say "last 3:"
    tail -3 "$BACKEND/data/logs/training-history.csv" | column -t -s ,
  else
    say "no training history yet - run ./run.sh train-beta or daily"
  fi
}

# ─── dev orchestrator (FastAPI + Next.js) ──────────────────────────────────────
FRONTEND="$ROOT/frontend"
LOGS="$BACKEND/data/logs"
PIDS="$LOGS/pids"

run_up() {
  section "starting dev servers"
  mkdir -p "$LOGS" "$PIDS"

  if [ -f "$PIDS/api.pid" ] && kill -0 "$(cat "$PIDS/api.pid")" 2>/dev/null; then
    warn "fastapi already running (pid $(cat "$PIDS/api.pid")) - skip"
  else
    say "fastapi   → http://localhost:8000"
    (cd "$BACKEND" && nohup "$UV" run uvicorn genomescope.api.server:app \
      --host 127.0.0.1 --port 8000 \
      > "$LOGS/api.log" 2>&1 & echo $! > "$PIDS/api.pid")
  fi

  if [ -f "$PIDS/web.pid" ] && kill -0 "$(cat "$PIDS/web.pid")" 2>/dev/null; then
    warn "next.js already running (pid $(cat "$PIDS/web.pid")) - skip"
  else
    if [ ! -d "$FRONTEND/node_modules" ]; then
      say "installing frontend dependencies (first run)"
      (cd "$FRONTEND" && npm install --silent)
    fi
    say "next.js   → http://localhost:3000"
    (cd "$FRONTEND" && nohup npm run dev \
      > "$LOGS/web.log" 2>&1 & echo $! > "$PIDS/web.pid")
  fi

  ok "dev servers up - logs at $LOGS/{api,web}.log"
  printf "\n  open: %shttp://localhost:3000%s\n" "$BOLD" "$RESET"
  printf "  stop: %s./run.sh down%s\n\n" "$BOLD" "$RESET"
}

run_down() {
  section "stopping dev servers"
  for name in api web; do
    local pid_file="$PIDS/$name.pid"
    if [ -f "$pid_file" ]; then
      local pid
      pid="$(cat "$pid_file")"
      if kill -0 "$pid" 2>/dev/null; then
        kill "$pid" 2>/dev/null || true
        # next dev spawns child processes - clean them up too
        pkill -P "$pid" 2>/dev/null || true
        ok "stopped $name (pid $pid)"
      else
        warn "$name pid file present but process gone"
      fi
      rm -f "$pid_file"
    fi
  done
}

run_logs() {
  local target="${2:-api}"
  if [ "$target" != "api" ] && [ "$target" != "web" ]; then
    err "logs target must be 'api' or 'web'"
    exit 1
  fi
  tail -f "$LOGS/$target.log"
}

# ─── dispatch ─────────────────────────────────────────────────────────────────
cmd="${1:-all}"
case "$cmd" in
  env)    ensure_uv; ensure_env ;;
  test)   ensure_uv; ensure_env; run_tests ;;
  data)   ensure_uv; ensure_env; run_data ;;
  tcga)   ensure_uv; ensure_env; run_tcga ;;
  tcga-luad) ensure_uv; ensure_env; run_tcga_luad ;;
  tcga-normal) ensure_uv; ensure_env; run_tcga_normal ;;
  smoke)  ensure_uv; ensure_env; run_smoke ;;
  quick)  ensure_uv; ensure_env; run_quick ;;
  train)  ensure_uv; ensure_env; run_train ;;
  beta)   ensure_uv; ensure_env; run_beta ;;
  train-beta) ensure_uv; ensure_env; run_train_beta ;;
  daily)      ensure_uv; ensure_env; run_daily ;;
  schedule-install)  ensure_uv; ensure_env; run_schedule_install ;;
  schedule-remove)   run_schedule_remove ;;
  schedule-status)   run_schedule_status ;;
  bench)  ensure_uv; ensure_env; run_bench ;;
  report) ensure_uv; ensure_env; run_report ;;
  report-full) ensure_uv; ensure_env; run_report_full ;;
  up)     ensure_uv; ensure_env; run_up ;;
  down)   run_down ;;
  logs)   run_logs "$@" ;;
  full)   ensure_uv; ensure_env; run_full ;;
  full-chr21-standard) ensure_uv; ensure_env; run_full_chr21_standard ;;
  full-chr21-beta) ensure_uv; ensure_env; run_full_chr21_beta ;;
  full-genome) ensure_uv; ensure_env; run_full_genome ;;
  ensemble) ensure_uv; ensure_env; run_ensemble ;;
  clean)  run_clean ;;
  all)
    ensure_uv
    ensure_env
    run_tests
    if [ ! -f "$BACKEND/data/raw/chr21.fa" ]; then
      run_data
    else
      section "reference data"
      ok "chr21.fa already present - skipping download"
    fi
    run_quick
    ;;
  -h|--help|help)
    sed -n '2,17p' "$0"
    exit 0
    ;;
  *)
    err "unknown command: $cmd"
    sed -n '2,17p' "$0"
    exit 1
    ;;
esac

section "done"
ok "finished: ./run.sh $cmd"
