"use client";

import { useMemo } from "react";
import { motion } from "framer-motion";
import { ShieldAlert } from "lucide-react";
import { Card, CardBody, CardHeader } from "@/components/ui/Card";
import type { IsoformReport, IsoformReportRow } from "@/lib/api";
import { cn } from "@/lib/cn";

type Props = {
  report: IsoformReport | null;
  onSelect: (id: string | null) => void;
};

export function TSGPanel({ report, onSelect }: Props) {
  const tsgHits = useMemo<IsoformReportRow[]>(() => {
    if (!report) return [];
    const seen = new Map<string, IsoformReportRow>();
    for (const r of report.rows) {
      if (!r.is_tsg) continue;
      const prev = seen.get(r.gene_symbol);
      if (!prev || (r.priority ?? 0) > (prev.priority ?? 0)) {
        seen.set(r.gene_symbol, r);
      }
    }
    return Array.from(seen.values())
      .sort((a, b) => (b.priority ?? 0) - (a.priority ?? 0))
      .slice(0, 5);
  }, [report]);

  const totalTsgRows = report?.rows.filter((r) => r.is_tsg).length ?? 0;

  return (
    <Card>
      <CardHeader
        kicker="Cancer-driver overlap"
        title="TSG hits"
        right={
          <span className="hidden font-mono text-[10px] uppercase tracking-widest text-ink-400 md:inline">
            built-in TSG list
          </span>
        }
      />
      <CardBody className="overflow-y-auto">
        {!report || totalTsgRows === 0 ? (
          <EmptyState />
        ) : (
          <div className="space-y-3">
            <div className="rounded-chip border border-hyper/30 bg-hyper-soft/40 px-4 py-3">
              <p className="font-mono text-[11px] leading-relaxed text-ink-800">
                <span className="text-base font-medium text-ink-900">
                  {totalTsgRows.toLocaleString()}
                </span>{" "}
                <span className="text-ink-600">
                  isoform overlaps land on known cancer-driver genes
                </span>
              </p>
              <p className="mt-1 font-mono text-[10px] text-ink-500">
                top hits ranked by mean β × COSMIC confidence
              </p>
            </div>

            <ol className="space-y-1.5">
              {tsgHits.map((row, i) => (
                <motion.li
                  key={row.gene_symbol}
                  initial={{ opacity: 0, y: 4 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{
                    delay: i * 0.04,
                    duration: 0.25,
                    ease: [0.23, 1, 0.32, 1],
                  }}
                  onClick={() => onSelect(row.island_id)}
                  className={cn(
                    "press flex cursor-pointer items-center justify-between rounded-chip",
                    "border border-ink-300/20 bg-cream-100/60 px-3 py-2.5",
                    "hover:border-ink-300/60 hover:bg-cream-100",
                  )}
                >
                  <div className="flex items-center gap-2">
                    <span className="font-mono text-[10px] tabular-nums text-ink-400">
                      {String(i + 1).padStart(2, "0")}
                    </span>
                    <span className="font-medium text-ink-900">
                      {row.gene_symbol}
                    </span>
                    {row.hypermethylated && (
                      <span className="rounded-full bg-hyper-soft px-1.5 py-0.5 font-mono text-[8px] uppercase tracking-wider text-hyper">
                        hyper
                      </span>
                    )}
                  </div>
                  <div className="flex items-center gap-3 font-mono text-[10px]">
                    <span className="text-ink-500">
                      <span className="text-ink-400">β </span>
                      <span
                        className={
                          row.hypermethylated ? "text-hyper" : "text-ink-700"
                        }
                      >
                        {(row.max_mean_beta ?? 0).toFixed(2)}
                      </span>
                    </span>
                    <span className="text-ink-400">{row.max_cohort || "-"}</span>
                  </div>
                </motion.li>
              ))}
            </ol>

            <p className="font-mono text-[10px] leading-relaxed text-ink-400">
              click any gene {">"} island detail panel opens
            </p>
          </div>
        )}
      </CardBody>
    </Card>
  );
}

function EmptyState() {
  return (
    <div className="rounded-chip border border-dashed border-ink-300/40 bg-cream-100/40 px-4 py-6 text-center">
      <ShieldAlert className="mx-auto size-5 text-ink-400" strokeWidth={1.5} />
      <p className="mt-2 font-mono text-[11px] text-ink-500">no TSG hits yet</p>
      <p className="mt-1 font-mono text-[10px] text-ink-400">
        run <span className="rounded bg-cream-200 px-1">./run.sh report-full</span>
      </p>
    </div>
  );
}
