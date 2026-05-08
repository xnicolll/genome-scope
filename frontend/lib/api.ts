// Typed wrappers around the FastAPI endpoints.
// Next.js dev server proxies /api/* to localhost:8000 (see next.config.mjs).

export type Island = { start: number; end: number };

export type Metrics = {
  precision: number;
  recall: number;
  f1: number;
  tp: number;
  fp: number;
  fn: number;
};

type WindowedRun = {
  chrom: string;
  window_start: number;
  window_end: number;
  window_length: number;
  islands: Island[];
  total_island_bp: number;
  log_likelihoods: number[];
  viterbi_log_prob: number;
  metrics: Metrics | null;
};

export type StandardRun = WindowedRun & { model: "standard" };

export type BetaRun = WindowedRun & {
  model: "beta";
  cpg_sites: number;
  state_means: [number, number];
};

export type BenchComparison = { standard: StandardRun; beta: BetaRun };

export type IsoformReportRow = {
  island_id: string;
  chrom: string;
  island_start: number;
  island_end: number;
  gene_id: string;
  gene_symbol: string;
  isoform_id: string;
  canonical: boolean;
  strand: "+" | "-";
  tss: number;
  overlap_bp: number;
  max_cohort: string;
  max_mean_beta: number | null;
  mean_beta_normal?: number | null;
  delta_vs_normal?: number | null;
  cancer_specific?: boolean;
  hypermethylated: boolean;
  is_tsg: boolean;
  cosmic_confidence: number;
  priority: number;
  [k: string]: unknown;
};

export type IsoformReport = {
  window: { chrom: string; start: number; end: number };
  primary_cohort: string;
  cohorts: string[];
  n_islands: number;
  n_report_rows: number;
  rows: IsoformReportRow[];
};

export type Truth = { chrom: string; n_islands: number; islands: Island[] };

export type ChromosomeRun = {
  chrom: string;
  assembled_bp: number;
  total_bp: number;
  n_islands: number;
  total_island_bp: number;
  metrics: Metrics;
  elapsed_s: number;
};

export type FullGenome = {
  n_chromosomes_attempted: number;
  n_chromosomes_completed: number;
  skipped: { chrom: string; skipped: string }[];
  per_chromosome: ChromosomeRun[];
  aggregate?: {
    total_assembled_bp: number;
    total_islands: number;
    micro_precision: number;
    micro_recall: number;
    micro_f1: number;
  };
};

// Alias kept for backwards-compat with existing imports.
export type EnsembleMetrics = Metrics;

export type Ensemble = {
  chrom: string;
  window_start: number;
  window_end: number;
  window_length: number;
  n_standard: number;
  n_beta: number;
  n_ensemble: number;
  ensemble_islands: Island[];
  standard_metrics: Metrics;
  beta_metrics: Metrics;
  ensemble_metrics: Metrics;
};

export type Checkpoint = {
  name: string;
  model_type: "standard" | "beta";
  created: string;
  updated: string;
  total_runs: number;
  total_samples_seen: number;
  cohorts_seen: string[];
  last_log_likelihood: number | null;
};

export type ModelList = { n_models: number; models: Checkpoint[] };

export type PredictResult = {
  chrom: string;
  model: string;
  checkpoint_used: string | null;
  window_start: number;
  window_end: number;
  window_length: number;
  n_islands: number;
  islands: Island[];
  total_island_bp: number;
  viterbi_log_prob: number;
  elapsed_seconds: number;
};

async function get<T>(path: string): Promise<T> {
  const r = await fetch(path, { cache: "no-store" });
  if (!r.ok) throw new Error(`${path} > ${r.status} ${r.statusText}`);
  return (await r.json()) as T;
}

async function post<T>(path: string, body: unknown): Promise<T> {
  const r = await fetch(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
    cache: "no-store",
  });
  if (!r.ok) {
    let detail = `${r.status} ${r.statusText}`;
    try {
      const j = (await r.json()) as { detail?: string };
      if (j.detail) detail = j.detail;
    } catch {}
    throw new Error(detail);
  }
  return (await r.json()) as T;
}

export const api = {
  health: () => get<{ status: string; processed_files: string[] }>("/api/health"),
  standard: () => get<StandardRun>("/api/pipeline/standard"),
  beta: () => get<BetaRun>("/api/pipeline/beta"),
  bench: () => get<BenchComparison>("/api/bench"),
  report: () => get<IsoformReport>("/api/report"),
  truth: () => get<Truth>("/api/truth"),
  models: () => get<ModelList>("/api/models"),
  fullGenome: () => get<FullGenome>("/api/full-genome"),
  ensemble: () => get<Ensemble>("/api/ensemble"),
  predict: (fasta: string, model: string = "standard") =>
    post<PredictResult>("/api/predict", { fasta, model }),
};
