"use client";

import { useMemo, useState } from "react";
import { Card, CardBody, CardHeader } from "@/components/ui/Card";
import { cn } from "@/lib/cn";
import type { IsoformReport, IsoformReportRow } from "@/lib/api";

type Props = {
  report: IsoformReport | null;
  onSelect: (islandId: string | null) => void;
  selected: string | null;
};

type FilterMode = "all" | "hyper" | "cancer_specific" | "tsg";

const FILTERS: { id: FilterMode; label: string }[] = [
  { id: "all", label: "all" },
  { id: "hyper", label: "hyper only" },
  { id: "cancer_specific", label: "cancer-specific" },
  { id: "tsg", label: "TSG only" },
];

function applyFilter(rows: IsoformReportRow[], mode: FilterMode): IsoformReportRow[] {
  switch (mode) {
    case "hyper":
      return rows.filter((r) => r.hypermethylated);
    case "cancer_specific":
      return rows.filter((r) => r.cancer_specific === true);
    case "tsg":
      return rows.filter((r) => r.is_tsg);
    default:
      return rows;
  }
}

export function ReportTable({ report, onSelect, selected }: Props) {
  const [mode, setMode] = useState<FilterMode>("all");

  const rows = useMemo(
    () => (report ? applyFilter(report.rows, mode) : []),
    [report, mode],
  );

  const hasNormal = useMemo(
    () =>
      report?.rows.some(
        (r) => r.delta_vs_normal !== null && r.delta_vs_normal !== undefined,
      ) ?? false,
    [report],
  );

  const cancerSpecificCount = report?.rows.filter((r) => r.cancer_specific).length ?? 0;
  const tsgCount = report?.rows.filter((r) => r.is_tsg).length ?? 0;

  const columns = hasNormal
    ? "minmax(0,2fr) minmax(0,1.4fr) 80px 80px 70px 80px"
    : "minmax(0,2fr) minmax(0,1.4fr) 80px 80px";

  return (
    <Card>
      <CardHeader
        kicker="Isoform report"
        title="Gene × cancer overlay"
        right={
          <div className="flex flex-wrap items-center gap-1.5">
            {FILTERS.map((f) => {
              if (f.id === "cancer_specific" && !hasNormal) return null;
              const active = mode === f.id;
              return (
                <button
                  key={f.id}
                  onClick={() => setMode(f.id)}
                  className={cn(
                    "press rounded-full border px-3 py-1.5 font-mono text-[10px] uppercase tracking-widest",
                    active
                      ? "border-ink-700 bg-ink-700 text-cream-50 hover:bg-ink-800"
                      : "border-ink-300/40 text-ink-500 hover:border-ink-400 hover:bg-cream-100 hover:text-ink-700",
                  )}
                >
                  {f.label}
                </button>
              );
            })}
          </div>
        }
      />

      {(cancerSpecificCount > 0 || tsgCount > 0) && (
        <div className="shrink-0 px-7 pb-1 pt-0 font-mono text-[10px] text-ink-500">
          {tsgCount > 0 && (
            <span className="mr-3">
              <span className="text-ink-700">{tsgCount.toLocaleString()}</span>{" "}
              on known cancer drivers
            </span>
          )}
          {cancerSpecificCount > 0 && (
            <span>
              <span className="text-ink-700">
                {cancerSpecificCount.toLocaleString()}
              </span>{" "}
              cancer-specific (healthy {">"} silenced)
            </span>
          )}
        </div>
      )}

      <div
        className="grid shrink-0 items-center gap-x-3 border-b border-ink-300/30 bg-cream-50 px-7 py-3 font-mono text-[10px] uppercase tracking-widest text-ink-400"
        style={{ gridTemplateColumns: columns }}
      >
        <span>Gene</span>
        <span>Isoform</span>
        <span className="text-right">Tumor β</span>
        {hasNormal && (
          <>
            <span className="text-right">Normal β</span>
            <span className="text-right">Δ</span>
          </>
        )}
        <span className="text-right">Cohort</span>
      </div>

      <CardBody className="overflow-y-auto pt-2">
        <ul className="space-y-1">
          {rows.slice(0, 60).map((row) => {
            const isSel = selected === row.island_id;
            const delta = row.delta_vs_normal;
            return (
              <li
                key={`${row.island_id}-${row.isoform_id}`}
                onClick={() => onSelect(isSel ? null : row.island_id)}
                className={cn(
                  "press grid cursor-pointer items-center gap-x-3 rounded-chip px-3 py-2.5 text-sm",
                  isSel ? "bg-accent-50" : "hover:bg-cream-100",
                )}
                style={{ gridTemplateColumns: columns }}
              >
                <div className="flex flex-wrap items-center gap-1.5">
                  <span className="font-medium text-ink-900">
                    {row.gene_symbol}
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
                  {row.cancer_specific && (
                    <span className="rounded-full bg-promoter-soft px-1.5 py-0.5 font-mono text-[8px] uppercase tracking-wider text-promoter">
                      cancer-specific
                    </span>
                  )}
                </div>
                <span className="truncate font-mono text-[11px] text-ink-500">
                  {row.isoform_id}
                </span>
                <span className="text-right font-mono text-[12px]">
                  {row.max_mean_beta !== null && row.max_mean_beta !== undefined ? (
                    <span className={row.hypermethylated ? "text-hyper" : "text-ink-700"}>
                      {row.max_mean_beta.toFixed(3)}
                    </span>
                  ) : (
                    <span className="text-ink-400">-</span>
                  )}
                </span>
                {hasNormal && (
                  <>
                    <span className="text-right font-mono text-[12px]">
                      {row.mean_beta_normal !== null && row.mean_beta_normal !== undefined ? (
                        <span className="text-ink-500">
                          {row.mean_beta_normal.toFixed(3)}
                        </span>
                      ) : (
                        <span className="text-ink-400">-</span>
                      )}
                    </span>
                    <span className="text-right font-mono text-[12px]">
                      {delta !== null && delta !== undefined ? (
                        <span
                          className={
                            delta > 0.2
                              ? "text-hyper"
                              : delta > 0.05
                                ? "text-ink-700"
                                : "text-ink-400"
                          }
                        >
                          {delta >= 0 ? "+" : ""}
                          {delta.toFixed(3)}
                        </span>
                      ) : (
                        <span className="text-ink-400">-</span>
                      )}
                    </span>
                  </>
                )}
                <span className="text-right font-mono text-[10px] uppercase tracking-widest text-ink-500">
                  {row.max_cohort || "-"}
                </span>
              </li>
            );
          })}
        </ul>
        {rows.length === 0 && (
          <p className="py-8 text-center font-mono text-xs text-ink-400">
            {mode === "all"
              ? "No report rows yet."
              : `No rows match the ${FILTERS.find((f) => f.id === mode)?.label} filter.`}
          </p>
        )}
        {report && rows.length > 60 && (
          <p className="pt-2 text-center font-mono text-[10px] text-ink-400">
            showing 60 of {rows.length.toLocaleString()} rows · sorted by priority
          </p>
        )}
      </CardBody>
    </Card>
  );
}
