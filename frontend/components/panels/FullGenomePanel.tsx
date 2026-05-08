"use client";

import { motion } from "framer-motion";
import { Globe } from "lucide-react";
import { Card, CardBody, CardHeader } from "@/components/ui/Card";
import { api, type FullGenome } from "@/lib/api";
import { useApi } from "@/lib/hooks";
import { cn } from "@/lib/cn";

export function FullGenomePanel() {
  const state = useApi(() => api.fullGenome().catch(() => null), []);

  return (
    <Card>
      <CardHeader
        kicker="Full-genome scaling"
        title="Standard HMM, multi-chromosome"
        right={
          <span className="hidden font-mono text-[10px] uppercase tracking-widest text-ink-400 md:inline">
            ./run.sh full-genome
          </span>
        }
      />
      <CardBody className="overflow-y-auto pt-2">
        {state.status === "loading" && (
          <p className="font-mono text-xs text-ink-400">loading…</p>
        )}
        {state.status === "error" && <EmptyState error={state.error} />}
        {state.status === "ready" &&
          (state.data ? <Body data={state.data} /> : <EmptyState />)}
      </CardBody>
    </Card>
  );
}

function EmptyState({ error }: { error?: string }) {
  return (
    <div className="rounded-chip border border-dashed border-ink-300/40 bg-cream-100/40 px-4 py-6 text-center">
      <Globe className="mx-auto size-5 text-ink-400" strokeWidth={1.5} />
      <p className="mt-2 font-mono text-[11px] text-ink-500">
        no full-genome run found yet
      </p>
      <p className="mt-1 font-mono text-[10px] text-ink-400">
        run <span className="rounded bg-cream-200 px-1">./run.sh full-genome</span>
        {" "}to populate this panel
      </p>
      {error && <p className="mt-2 font-mono text-[10px] text-hyper">{error}</p>}
    </div>
  );
}

function Body({ data }: { data: FullGenome }) {
  const { aggregate: agg, per_chromosome: rows } = data;
  const maxIslands = Math.max(...rows.map((r) => r.n_islands));

  return (
    <div className="space-y-4">
      {agg && (
        <>
          <div className="grid grid-cols-3 gap-2">
            <Metric label="chromosomes" value={data.n_chromosomes_completed.toString()} />
            <Metric
              label="genome covered"
              value={`${(agg.total_assembled_bp / 1_000_000).toFixed(0)} Mb`}
            />
            <Metric label="islands found" value={agg.total_islands.toLocaleString()} />
          </div>

          <div className="rounded-chip border border-ink-300/25 bg-cream-100/60 px-4 py-3">
            <div className="grid grid-cols-3 gap-x-4 gap-y-1 font-mono text-[11px]">
              <Mini label="precision" value={agg.micro_precision.toFixed(3)} />
              <Mini label="recall" value={agg.micro_recall.toFixed(3)} />
              <Mini label="F1" value={agg.micro_f1.toFixed(3)} />
            </div>
            <p className="mt-2 font-mono text-[10px] leading-relaxed text-ink-500">
              recall stays at{" "}
              <span className="text-ink-700">
                {(agg.micro_recall * 100).toFixed(0)}%
              </span>{" "}
              across every chromosome - method generalises beyond chr21
            </p>
          </div>
        </>
      )}

      <ul className="space-y-2">
        {rows.map((r, i) => (
          <motion.li
            key={r.chrom}
            initial={{ opacity: 0, x: -4 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{
              delay: i * 0.04,
              duration: 0.25,
              ease: [0.23, 1, 0.32, 1],
            }}
            className="grid grid-cols-[60px_1fr_auto] items-center gap-3"
          >
            <span className="font-mono text-[11px] font-medium text-ink-700">
              {r.chrom}
            </span>
            <div className="relative h-5 overflow-hidden rounded-full bg-ink-300/15">
              <motion.div
                initial={{ width: 0 }}
                animate={{ width: `${(r.n_islands / maxIslands) * 100}%` }}
                transition={{
                  delay: 0.15 + i * 0.05,
                  duration: 0.5,
                  ease: [0.23, 1, 0.32, 1],
                }}
                className="h-full rounded-full"
                style={{ background: f1Colour(r.metrics.f1) }}
              />
              <span className="absolute inset-0 flex items-center px-2 font-mono text-[10px] text-ink-700">
                {r.n_islands.toLocaleString()} islands
              </span>
            </div>
            <span
              className={cn(
                "w-12 text-right font-mono text-[11px]",
                r.metrics.f1 >= 0.5
                  ? "text-island"
                  : r.metrics.f1 >= 0.4
                    ? "text-ink-700"
                    : "text-hyper",
              )}
            >
              {r.metrics.f1.toFixed(2)}
            </span>
          </motion.li>
        ))}
      </ul>

      <p className="font-mono text-[10px] leading-relaxed text-ink-400">
        bar length = islands flagged · colour = F1 (green high, amber low)
      </p>
    </div>
  );
}

function Metric({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded-chip border border-ink-300/20 bg-cream-50 px-3 py-2.5">
      <p className="font-mono text-[9px] uppercase tracking-widest text-ink-400">
        {label}
      </p>
      <p className="mt-0.5 font-mono text-lg tracking-tight text-ink-900">
        {value}
      </p>
    </div>
  );
}

function Mini({ label, value }: { label: string; value: string }) {
  return (
    <div>
      <span className="block text-[9px] uppercase tracking-widest text-ink-400">
        {label}
      </span>
      <span className="text-ink-900">{value}</span>
    </div>
  );
}

function f1Colour(f1: number): string {
  if (f1 >= 0.55) return "#5BA89A";
  if (f1 >= 0.45) return "#8AB8AB";
  if (f1 >= 0.35) return "#D9C677";
  return "#E89E3D";
}
