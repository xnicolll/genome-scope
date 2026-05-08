"use client";

import { useState, type ReactNode } from "react";
import Link from "next/link";
import { motion, AnimatePresence } from "framer-motion";
import { ArrowLeft, Upload } from "lucide-react";
import { api } from "@/lib/api";
import { useApi } from "@/lib/hooks";
import { StatusDot } from "@/components/ui/StatusDot";
import { TrackPanel } from "@/components/panels/TrackPanel";
import { MetricsPanel } from "@/components/panels/MetricsPanel";
import { ReportTable } from "@/components/panels/ReportTable";
import { TSGPanel } from "@/components/panels/TSGPanel";
import { FullGenomePanel } from "@/components/panels/FullGenomePanel";
import { IslandDetailInline } from "@/components/IslandDetailInline";
import { Card, CardBody } from "@/components/ui/Card";

const ROW_1_HEIGHT = "h-[420px]";
const ROW_2_HEIGHT = "h-[520px]";
const ROW_3_HEIGHT = "h-[700px]";

const SWAP_ANIM = {
  initial: { opacity: 0, y: 4 },
  animate: { opacity: 1, y: 0 },
  exit: { opacity: 0, y: -4 },
  transition: { duration: 0.2, ease: [0.23, 1, 0.32, 1] as const },
};

type Source = "track" | "tsg" | "report";
type Selection = { islandId: string; source: Source } | null;

export default function DashboardPage() {
  const [selection, setSelection] = useState<Selection>(null);

  const health = useApi(() => api.health(), []);
  const standard = useApi(() => api.standard(), []);
  const bench = useApi(() => api.bench().catch(() => null), []);
  const report = useApi(() => api.report().catch(() => null), []);
  const truth = useApi(() => api.truth(), []);
  const fullGenome = useApi(() => api.fullGenome().catch(() => null), []);

  const serverOk = health.status === "ready";
  const reportData = report.status === "ready" ? report.data : null;

  const select = (source: Source) => (id: string | null) => {
    setSelection(id ? { islandId: id, source } : null);
  };
  const back = () => setSelection(null);
  const isActive = (source: Source) => selection?.source === source;

  return (
    <main className="min-h-screen px-5 py-8 md:px-10 md:py-12">
      <div className="mx-auto max-w-[1440px]">
        <header className="mb-10 flex flex-wrap items-center justify-between gap-4">
          <div className="flex items-center gap-4">
            <Link
              href="/"
              className="press inline-flex size-9 items-center justify-center rounded-full border border-ink-300/40 bg-cream-50 text-ink-500 hover:border-ink-400 hover:bg-cream-100 hover:text-ink-900"
              aria-label="Back"
            >
              <ArrowLeft className="size-4" />
            </Link>
            <div>
              <p className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
                GenomeScope · dashboard
              </p>
              <h1 className="mt-0.5 text-2xl font-medium tracking-tight text-ink-900">
                HMM-based CpG island detector
              </h1>
              <DataScale
                standardWindowBp={
                  standard.status === "ready" ? standard.data.window_length : null
                }
                fullGenomeMb={
                  fullGenome.status === "ready" && fullGenome.data?.aggregate
                    ? Math.round(fullGenome.data.aggregate.total_assembled_bp / 1_000_000)
                    : null
                }
                cohorts={reportData ? reportData.cohorts : []}
              />
            </div>
          </div>
          <div className="flex items-center gap-3">
            <Link
              href="/upload"
              className="press inline-flex items-center gap-1.5 rounded-full bg-ink px-4 py-2 font-mono text-[10px] uppercase tracking-widest text-cream-50 hover:bg-ink-700"
            >
              <Upload className="size-3" />
              upload fasta
            </Link>
            <div className="flex items-center gap-3 rounded-full border border-ink-300/30 bg-cream-50 px-4 py-2 shadow-card">
              <StatusDot ok={serverOk} />
              <span className="font-mono text-[11px] text-ink-600">
                {serverOk ? "api online" : "api offline"}
              </span>
              {serverOk && health.status === "ready" && (
                <span className="font-mono text-[10px] text-ink-400">
                  · {health.data.processed_files.length} files
                </span>
              )}
            </div>
          </div>
        </header>

        {standard.status === "error" && (
          <Card className="mb-6 border-hyper/40 bg-hyper-soft/40">
            <CardBody>
              <p className="font-mono text-sm text-ink-900">
                couldn&apos;t reach the api · {standard.error}
              </p>
              <p className="mt-1 font-mono text-[11px] text-ink-500">
                run <span className="bg-cream-200 px-1">./run.sh up</span> at the repo root
              </p>
            </CardBody>
          </Card>
        )}

        <div className="space-y-4">
          <div className={`overflow-hidden ${ROW_1_HEIGHT}`}>
            <PanelSwap
              showDetail={isActive("track") && !!selection && !!reportData}
              detail={
                selection && reportData ? (
                  <IslandDetailInline
                    islandId={selection.islandId}
                    rows={reportData.rows}
                    onBack={back}
                  />
                ) : null
              }
              panel={
                standard.status === "ready" && truth.status === "ready" ? (
                  <TrackPanel
                    standard={standard.data}
                    truth={truth.data}
                    report={reportData}
                    selected={null}
                    onSelect={select("track")}
                  />
                ) : (
                  <SkeletonCard height={340} />
                )
              }
            />
          </div>

          <div className={`grid gap-4 overflow-hidden lg:grid-cols-2 ${ROW_2_HEIGHT}`}>
            <PanelSwap
              showDetail={isActive("tsg") && !!selection && !!reportData}
              detail={
                selection && reportData ? (
                  <IslandDetailInline
                    islandId={selection.islandId}
                    rows={reportData.rows}
                    onBack={back}
                  />
                ) : null
              }
              panel={
                <TSGPanel report={reportData} onSelect={select("tsg")} />
              }
            />
            <FullGenomePanel />
          </div>

          <div className={`grid gap-4 overflow-hidden lg:grid-cols-2 ${ROW_3_HEIGHT}`}>
            <MetricsPanel
              bench={bench.status === "ready" ? bench.data : null}
            />
            <PanelSwap
              showDetail={isActive("report") && !!selection && !!reportData}
              detail={
                selection && reportData ? (
                  <IslandDetailInline
                    islandId={selection.islandId}
                    rows={reportData.rows}
                    onBack={back}
                  />
                ) : null
              }
              panel={
                <ReportTable
                  report={reportData}
                  selected={null}
                  onSelect={select("report")}
                />
              }
            />
          </div>
        </div>

        <footer className="relative z-10 mt-12 flex items-center justify-between border-t border-ink-300/30 bg-cream pt-6 font-mono text-[11px] text-ink-400">
          <p>Joe Nicol · CBMFW 4761 Computational Genomics · Spring 2026</p>
          <p>{new Date().toISOString().slice(0, 10)}</p>
        </footer>
      </div>
    </main>
  );
}

function SkeletonCard({ height }: { height: number }) {
  return (
    <div
      className="animate-pulse rounded-card border border-ink-300/30 bg-cream-100"
      style={{ height }}
    />
  );
}

function PanelSwap({
  showDetail,
  detail,
  panel,
}: {
  showDetail: boolean;
  detail: ReactNode;
  panel: ReactNode;
}) {
  return (
    <div className="relative h-full overflow-hidden">
      <AnimatePresence mode="wait" initial={false}>
        <motion.div
          key={showDetail ? "detail" : "panel"}
          {...SWAP_ANIM}
          className="h-full"
        >
          {showDetail ? detail : panel}
        </motion.div>
      </AnimatePresence>
    </div>
  );
}

function DataScale({
  standardWindowBp,
  fullGenomeMb,
  cohorts,
}: {
  standardWindowBp: number | null;
  fullGenomeMb: number | null;
  cohorts: string[];
}) {
  const parts: string[] = [];
  if (standardWindowBp !== null) {
    const mb = standardWindowBp / 1_000_000;
    parts.push(
      mb >= 10 ? `${mb.toFixed(0)} Mb of chr21` : `${mb.toFixed(1)} Mb of chr21`,
    );
  }
  if (cohorts.length > 0) {
    parts.push(
      `${cohorts.length} cohort${cohorts.length > 1 ? "s" : ""} · ${cohorts.join(", ")}`,
    );
  }
  if (fullGenomeMb !== null) {
    parts.push(`+${fullGenomeMb} Mb genome-wide`);
  }
  if (parts.length === 0) return null;
  return (
    <p className="mt-1.5 font-mono text-[11px] text-ink-500">
      {parts.join("  ·  ")}
    </p>
  );
}
