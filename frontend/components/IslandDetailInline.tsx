"use client";

import { ArrowLeft, Dna } from "lucide-react";
import type { IsoformReportRow } from "@/lib/api";
import { Card, CardBody } from "@/components/ui/Card";

type Props = {
  islandId: string;
  rows: IsoformReportRow[];
  onBack: () => void;
};

function formatBp(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(2)} Mb`;
  if (n >= 1_000) return `${(n / 1_000).toFixed(1)} kb`;
  return `${n} bp`;
}

export function IslandDetailInline({ islandId, rows, onBack }: Props) {
  const matching = rows.filter((r) => r.island_id === islandId);
  const head = matching[0];

  return (
    <Card>
      <div className="flex h-full flex-col">
        <header className="flex shrink-0 items-start gap-3 border-b border-ink-300/20 px-7 py-5">
          <button
            onClick={onBack}
            className="press mt-0.5 inline-flex size-7 shrink-0 items-center justify-center rounded-full border border-ink-300/40 bg-cream-50 text-ink-500 hover:border-ink-400 hover:bg-cream-100 hover:text-ink-900"
            aria-label="Back"
          >
            <ArrowLeft className="size-3.5" />
          </button>
          <div className="min-w-0 flex-1">
            <p className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
              Island detail
            </p>
            <h3 className="mt-0.5 truncate font-mono text-lg font-medium tracking-tight text-ink-900">
              {islandId}
            </h3>
            {head && (
              <p className="mt-1 truncate font-mono text-xs text-ink-500">
                {head.chrom}:{head.island_start.toLocaleString()}-
                {head.island_end.toLocaleString()}
                <span className="ml-2 text-ink-400">
                  ({formatBp(head.island_end - head.island_start)})
                </span>
              </p>
            )}
          </div>
        </header>

        <CardBody className="overflow-y-auto">
          {matching.length === 0 ? (
            <p className="font-mono text-xs text-ink-400">
              No isoform overlap found for this island. It probably sits in
              an intergenic region.
            </p>
          ) : (
            <div className="space-y-4">
              <div className="grid grid-cols-2 gap-2">
                <Stat
                  label="Overlapping isoforms"
                  value={matching.length.toString()}
                />
                {head?.max_mean_beta !== null && head?.max_mean_beta !== undefined && (
                  <Stat
                    label={`Mean β${head.max_cohort ? ` (${head.max_cohort})` : ""}`}
                    value={head.max_mean_beta.toFixed(3)}
                    accent={head.hypermethylated ? "amber" : undefined}
                  />
                )}
              </div>

              <div>
                <p className="mb-2 font-mono text-[10px] uppercase tracking-widest text-ink-400">
                  Isoforms
                </p>
                <ul className="space-y-1.5">
                  {matching.map((row) => (
                    <li
                      key={row.isoform_id}
                      className="press flex items-center justify-between gap-3 rounded-chip border border-ink-300/20 bg-cream-100/70 px-3 py-2.5 hover:border-ink-300/50 hover:bg-cream-100"
                    >
                      <div className="min-w-0">
                        <div className="flex flex-wrap items-baseline gap-1.5">
                          <span className="truncate text-sm font-medium text-ink-900">
                            {row.gene_symbol || row.gene_id}
                          </span>
                          {row.canonical && (
                            <span className="rounded-full bg-accent-100 px-1.5 py-0.5 font-mono text-[8px] uppercase tracking-wider text-ink-700">
                              canonical
                            </span>
                          )}
                          {row.is_tsg && (
                            <span className="rounded-full bg-hyper-soft px-1.5 py-0.5 font-mono text-[8px] uppercase tracking-wider text-hyper">
                              TSG
                            </span>
                          )}
                        </div>
                        <p className="mt-0.5 truncate font-mono text-[10px] text-ink-400">
                          {row.isoform_id} · {row.strand} · TSS{" "}
                          {row.tss.toLocaleString()}
                        </p>
                      </div>
                      <Dna className="size-3.5 shrink-0 text-ink-400" />
                    </li>
                  ))}
                </ul>
              </div>
            </div>
          )}
        </CardBody>
      </div>
    </Card>
  );
}

function Stat({
  label,
  value,
  accent,
}: {
  label: string;
  value: string;
  accent?: "amber";
}) {
  return (
    <div className="rounded-chip border border-ink-300/20 bg-cream-100/70 px-4 py-3">
      <p className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
        {label}
      </p>
      <p
        className={`mt-1 font-mono text-2xl tracking-tight ${
          accent === "amber" ? "text-hyper" : "text-ink-900"
        }`}
      >
        {value}
      </p>
    </div>
  );
}
