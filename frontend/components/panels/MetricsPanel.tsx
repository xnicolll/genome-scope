"use client";

import { motion } from "framer-motion";
import { Sparkles } from "lucide-react";
import { Card, CardBody, CardHeader } from "@/components/ui/Card";
import { StatBlock } from "@/components/Stat";
import { api, type BenchComparison, type EnsembleMetrics } from "@/lib/api";
import { useApi } from "@/lib/hooks";
import { cn } from "@/lib/cn";

type Props = { bench: BenchComparison | null };

export function MetricsPanel({ bench }: Props) {
  const combo = useApi(() => api.ensemble().catch(() => null), []);
  const comboData = combo.status === "ready" ? combo.data : null;

  const std = comboData?.standard_metrics ?? bench?.standard?.metrics;
  const beta = comboData?.beta_metrics ?? bench?.beta?.metrics;
  const ens = comboData?.ensemble_metrics ?? null;

  const rows: {
    label: string;
    sub: string;
    metrics: EnsembleMetrics | null | undefined;
    badge?: "novel" | "best";
  }[] = [
    { label: "Standard 8-state HMM", sub: "nucleotide only", metrics: std },
    {
      label: "Beta methylation HMM",
      sub: "TCGA-BRCA, trained",
      metrics: beta,
      badge: "novel",
    },
  ];
  if (ens) {
    rows.push({
      label: "Combination (Standard × Beta)",
      sub: "intersection",
      metrics: ens,
      badge: "best",
    });
  }

  const maxF1 = Math.max(std?.f1 ?? 0, beta?.f1 ?? 0, ens?.f1 ?? 0);

  return (
    <Card>
      <CardHeader
        kicker="Benchmark"
        title="Model comparison"
        right={
          <span className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
            vs UCSC reference
          </span>
        }
      />
      <CardBody>
        <div className="grid gap-3">
          {rows.map((row, idx) => {
            const isBest =
              row.metrics?.f1 !== undefined &&
              row.metrics?.f1 !== null &&
              Math.abs(row.metrics.f1 - maxF1) < 1e-6 &&
              maxF1 > 0;
            return (
              <motion.div
                key={row.label}
                initial={{ opacity: 0, y: 4 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{
                  delay: idx * 0.04,
                  duration: 0.3,
                  ease: [0.23, 1, 0.32, 1],
                }}
                className={cn(
                  "rounded-chip border p-5 transition-shadow",
                  row.badge === "best" && isBest
                    ? "border-accent bg-accent-50/60 shadow-card"
                    : "border-ink-300/25 bg-cream-100/60",
                )}
              >
                <div className="mb-4 flex items-baseline justify-between gap-3">
                  <div>
                    <p className="flex items-center gap-1.5 text-sm font-medium text-ink-900">
                      {row.badge === "best" && (
                        <Sparkles className="size-3.5 text-accent-500" />
                      )}
                      {row.label}
                    </p>
                    <p className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
                      {row.sub}
                    </p>
                  </div>
                  {row.badge === "novel" && (
                    <span className="rounded-full bg-accent px-2 py-1 font-mono text-[9px] uppercase tracking-widest text-ink-900">
                      novel
                    </span>
                  )}
                  {row.badge === "best" && (
                    <span className="rounded-full bg-ink-900 px-2 py-1 font-mono text-[9px] uppercase tracking-widest text-cream-50">
                      best F1
                    </span>
                  )}
                </div>
                <div className="grid grid-cols-3 gap-2">
                  <MetricCell label="Precision" value={row.metrics?.precision} />
                  <MetricCell label="Recall" value={row.metrics?.recall} />
                  <MetricCell
                    label="F1"
                    value={row.metrics?.f1}
                    emphasised={row.badge === "best"}
                  />
                </div>
              </motion.div>
            );
          })}
        </div>
      </CardBody>
    </Card>
  );
}

function MetricCell({
  label,
  value,
  emphasised,
}: {
  label: string;
  value: number | null | undefined;
  emphasised?: boolean;
}) {
  return (
    <StatBlock
      label={label}
      value={value !== null && value !== undefined ? value.toFixed(3) : "-"}
      tone={emphasised ? "accent" : "default"}
      className={emphasised ? "bg-accent-100/60 ring-1 ring-accent/30" : ""}
    />
  );
}
