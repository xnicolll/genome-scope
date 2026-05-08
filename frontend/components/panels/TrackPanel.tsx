"use client";

import { Card, CardBody, CardHeader } from "@/components/ui/Card";
import { GenomicTrack, type TrackIsland, type TrackPromoter } from "@/components/GenomicTrack";
import type { IsoformReport, StandardRun, Truth } from "@/lib/api";

type Props = {
  standard: StandardRun;
  truth: Truth;
  report: IsoformReport | null;
  selected: string | null;
  onSelect: (id: string | null) => void;
};

export function TrackPanel({ standard, truth, report, selected, onSelect }: Props) {
  const hyperIds = new Set<string>();
  const reportByWindow = new Map<string, string>();
  if (report) {
    for (const row of report.rows) {
      reportByWindow.set(`${row.island_start}-${row.island_end}`, row.island_id);
      if (row.hypermethylated) {
        hyperIds.add(row.island_id);
      }
    }
  }

  const islands: TrackIsland[] = standard.islands.map((i, idx) => {
    const key = `${i.start}-${i.end}`;
    const id = reportByWindow.get(key) ?? `pred-${idx}`;
    return { id, start: i.start, end: i.end, hypermethylated: hyperIds.has(id) };
  });

  const truthIslands: TrackIsland[] = truth.islands
    .filter((t) => t.end >= standard.window_start && t.start <= standard.window_end)
    .map((t, idx) => ({ id: `truth-${idx}`, start: t.start, end: t.end }));

  const promoters: TrackPromoter[] = report
    ? Array.from(
        new Map(
          report.rows.map((r) => [
            r.isoform_id,
            {
              id: r.isoform_id,
              start: Math.max(0, r.tss - 2000),
              end: r.tss + 2000,
              strand: r.strand,
            } as TrackPromoter,
          ]),
        ).values(),
      )
    : [];

  const hyperCount = islands.filter((i) => i.hypermethylated).length;

  return (
    <Card>
      <CardHeader
        kicker="Phase 4 · chr21 viewer"
        title="Genomic track"
        right={
          <div className="hidden items-center gap-4 md:flex">
            <LegendSwatch color="#5BA89A" label="CpG island" />
            <LegendSwatch color="#E89E3D" label="Hyper-meth" />
            <LegendSwatch color="#7B3F61" label="Promoter" />
          </div>
        }
      />
      <CardBody>
        <GenomicTrack
          windowStart={standard.window_start}
          windowEnd={standard.window_end}
          truth={truthIslands}
          islands={islands}
          promoters={promoters}
          selectedIslandId={selected}
          onSelectIsland={onSelect}
        />
        <div className="mt-3 flex flex-wrap items-center gap-x-6 gap-y-1 font-mono text-[11px] text-ink-500">
          <span>
            {islands.length.toLocaleString()} predicted islands
          </span>
          <span>{truthIslands.length.toLocaleString()} UCSC truth</span>
          <span className="text-hyper">{hyperCount} hyper-methylated</span>
          <span className="ml-auto text-ink-400">
            drag on the track to zoom
          </span>
        </div>
      </CardBody>
    </Card>
  );
}

function LegendSwatch({ color, label }: { color: string; label: string }) {
  return (
    <div className="flex items-center gap-1.5">
      <span
        aria-hidden
        className="block size-2.5 rounded-[2px]"
        style={{ background: color }}
      />
      <span className="font-mono text-[10px] uppercase tracking-widest text-ink-500">
        {label}
      </span>
    </div>
  );
}
