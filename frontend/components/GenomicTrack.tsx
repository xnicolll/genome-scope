"use client";

import * as d3 from "d3";
import { useEffect, useMemo, useRef, useState } from "react";
import { cn } from "@/lib/cn";

export type TrackIsland = {
  id: string;
  start: number;
  end: number;
  hypermethylated?: boolean;
};

export type TrackPromoter = {
  id: string;
  start: number;
  end: number;
  strand?: "+" | "-";
  label?: string;
};

export type TrackProps = {
  windowStart: number;
  windowEnd: number;
  truth: TrackIsland[];
  islands: TrackIsland[];
  promoters: TrackPromoter[];
  selectedIslandId?: string | null;
  onSelectIsland?: (id: string | null) => void;
  height?: number;
};

const LANE_GAP = 14;
const LANE_HEIGHT = 22;
const MARGIN = { top: 24, right: 16, bottom: 40, left: 16 };

const LANES = [
  { key: "truth", label: "UCSC truth", colour: "var(--color-ink-300)" },
  { key: "islands", label: "Predicted islands", colour: "#5BA89A" },
  { key: "hyper", label: "Hyper-methylated", colour: "#E89E3D" },
  { key: "promoters", label: "Gene promoters", colour: "#7B3F61" },
] as const;

type HoverState = {
  island: TrackIsland;
  x: number;
  y: number;
} | null;

export function GenomicTrack({
  windowStart,
  windowEnd,
  truth,
  islands,
  promoters,
  selectedIslandId = null,
  onSelectIsland,
  height = 280,
}: TrackProps) {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const svgRef = useRef<SVGSVGElement | null>(null);
  const [width, setWidth] = useState(800);
  const [zoom, setZoom] = useState<[number, number]>([windowStart, windowEnd]);
  const [hover, setHover] = useState<HoverState>(null);

  useEffect(() => {
    if (!containerRef.current) return;
    const el = containerRef.current;
    const obs = new ResizeObserver(([entry]) => {
      setWidth(Math.max(320, entry.contentRect.width));
    });
    obs.observe(el);
    return () => obs.disconnect();
  }, []);

  useEffect(() => {
    setZoom([windowStart, windowEnd]);
  }, [windowStart, windowEnd]);

  const innerWidth = width - MARGIN.left - MARGIN.right;
  const innerHeight = height - MARGIN.top - MARGIN.bottom;

  const x = useMemo(
    () => d3.scaleLinear().domain(zoom).range([0, innerWidth]).clamp(true),
    [zoom, innerWidth],
  );

  const laneY = (idx: number) =>
    idx * (LANE_HEIGHT + LANE_GAP) + (innerHeight - (LANE_HEIGHT * 4 + LANE_GAP * 3)) / 2;

  useEffect(() => {
    if (!svgRef.current) return;
    const svg = d3.select(svgRef.current);

    const axisGroup = svg.select<SVGGElement>(".axis");
    const axis = d3
      .axisBottom(x)
      .ticks(Math.max(4, Math.min(8, Math.floor(innerWidth / 110))))
      .tickFormat((d) => {
        const mb = Number(d) / 1_000_000;
        return mb >= 1 ? `${mb.toFixed(2)} Mb` : `${(Number(d) / 1000).toFixed(0)} kb`;
      });
    axisGroup
      .attr("transform", `translate(${MARGIN.left}, ${MARGIN.top + innerHeight + 10})`)
      .call(axis as any);
    axisGroup.selectAll("path,line").attr("stroke", "rgba(26,24,20,0.2)");
    axisGroup
      .selectAll("text")
      .attr("fill", "rgba(26,24,20,0.55)")
      .style("font", "11px var(--font-geist-mono, ui-monospace)");

    const brush = d3
      .brushX()
      .extent([
        [0, 0],
        [innerWidth, innerHeight],
      ])
      .on("end", (event) => {
        if (!event.selection) return;
        const [x0, x1] = event.selection as [number, number];
        const range: [number, number] = [x.invert(x0), x.invert(x1)];
        if (range[1] - range[0] > 500) {
          setZoom(range);
        }
        svg.select<SVGGElement>(".brush").call(brush.move as any, null);
      });
    svg
      .select<SVGGElement>(".brush")
      .attr("transform", `translate(${MARGIN.left}, ${MARGIN.top})`)
      .call(brush as any);
    svg
      .select<SVGGElement>(".brush")
      .selectAll(".selection")
      .attr("fill", "rgba(245, 213, 71, 0.25)")
      .attr("stroke", "rgba(232, 158, 61, 0.8)")
      .attr("stroke-width", 1);
  }, [x, innerWidth, innerHeight]);

  const rectFor = (start: number, end: number) => {
    const xStart = Math.max(0, x(start));
    const xEnd = Math.min(innerWidth, x(end));
    return { x: xStart, width: Math.max(1.5, xEnd - xStart) };
  };

  const visibleIslands = islands.filter(
    (i) => i.end >= zoom[0] && i.start <= zoom[1],
  );

  const zoomed = zoom[0] !== windowStart || zoom[1] !== windowEnd;

  return (
    <div ref={containerRef} className="relative w-full">
      <svg
        ref={svgRef}
        width={width}
        height={height}
        className="block select-none"
        style={{ cursor: "crosshair" }}
      >
        {/* brush rendered first so island rects above it receive mouse events */}
        <g className="brush" />

        <g transform={`translate(${MARGIN.left}, ${MARGIN.top})`}>
          {x.ticks(8).map((t) => (
            <line
              key={t}
              x1={x(t)}
              x2={x(t)}
              y1={0}
              y2={innerHeight}
              stroke="rgba(26,24,20,0.04)"
            />
          ))}
        </g>

        <g transform={`translate(${MARGIN.left}, ${MARGIN.top})`}>
          {LANES.map((lane, idx) => (
            <g key={lane.key} transform={`translate(0, ${laneY(idx)})`}>
              <text
                x={0}
                y={-4}
                className="fill-ink-400"
                style={{
                  font: "10px var(--font-geist-mono, ui-monospace)",
                  letterSpacing: "0.08em",
                  textTransform: "uppercase",
                }}
              >
                {lane.label}
              </text>
              <line
                x1={0}
                x2={innerWidth}
                y1={LANE_HEIGHT / 2}
                y2={LANE_HEIGHT / 2}
                stroke="rgba(26,24,20,0.06)"
              />
            </g>
          ))}
        </g>

        <g transform={`translate(${MARGIN.left}, ${MARGIN.top + laneY(0)})`}>
          {truth
            .filter((t) => t.end >= zoom[0] && t.start <= zoom[1])
            .map((t, i) => {
              const r = rectFor(t.start, t.end);
              return (
                <rect
                  key={`truth-${i}`}
                  x={r.x}
                  y={4}
                  width={r.width}
                  height={LANE_HEIGHT - 8}
                  rx={3}
                  fill="rgba(26, 24, 20, 0.14)"
                />
              );
            })}
        </g>

        <g transform={`translate(${MARGIN.left}, ${MARGIN.top + laneY(1)})`}>
          {visibleIslands.map((isl) => {
            const r = rectFor(isl.start, isl.end);
            const selected = isl.id === selectedIslandId;
            const hovered = hover?.island.id === isl.id;
            return (
              <rect
                key={isl.id}
                x={r.x}
                y={2}
                width={r.width}
                height={LANE_HEIGHT - 4}
                rx={4}
                fill={isl.hypermethylated ? "#E89E3D" : "#5BA89A"}
                fillOpacity={selected ? 1 : hovered ? 0.95 : 0.85}
                stroke={
                  selected
                    ? "#1A1814"
                    : hovered
                      ? "rgba(26,24,20,0.4)"
                      : "transparent"
                }
                strokeWidth={selected ? 1.5 : 1}
                style={{
                  cursor: "pointer",
                  transition: "fill-opacity 120ms ease, stroke 120ms ease",
                }}
                onMouseEnter={(e) => {
                  const rect = containerRef.current?.getBoundingClientRect();
                  setHover({
                    island: isl,
                    x: e.clientX - (rect?.left ?? 0),
                    y: e.clientY - (rect?.top ?? 0),
                  });
                }}
                onMouseMove={(e) => {
                  const rect = containerRef.current?.getBoundingClientRect();
                  setHover((h) =>
                    h?.island.id === isl.id
                      ? { ...h, x: e.clientX - (rect?.left ?? 0), y: e.clientY - (rect?.top ?? 0) }
                      : h,
                  );
                }}
                onMouseLeave={() => setHover(null)}
                onClick={(e) => {
                  e.stopPropagation();
                  onSelectIsland?.(selected ? null : isl.id);
                }}
              />
            );
          })}
        </g>

        <g transform={`translate(${MARGIN.left}, ${MARGIN.top + laneY(2)})`}>
          {visibleIslands
            .filter((i) => i.hypermethylated)
            .map((isl) => {
              const r = rectFor(isl.start, isl.end);
              return (
                <rect
                  key={`hy-${isl.id}`}
                  x={r.x}
                  y={2}
                  width={r.width}
                  height={LANE_HEIGHT - 4}
                  rx={4}
                  fill="#E89E3D"
                  fillOpacity={0.9}
                />
              );
            })}
        </g>

        <g transform={`translate(${MARGIN.left}, ${MARGIN.top + laneY(3)})`}>
          {promoters
            .filter((p) => p.end >= zoom[0] && p.start <= zoom[1])
            .slice(0, 300)
            .map((p) => {
              const r = rectFor(p.start, p.end);
              return (
                <rect
                  key={p.id}
                  x={r.x}
                  y={4}
                  width={r.width}
                  height={LANE_HEIGHT - 8}
                  rx={2}
                  fill="#7B3F61"
                  fillOpacity={0.7}
                />
              );
            })}
        </g>

        <g className="axis" />
      </svg>

      {zoomed && (
        <button
          onClick={() => setZoom([windowStart, windowEnd])}
          className={cn(
            "press absolute right-3 top-3 rounded-full border border-ink-300/40",
            "bg-cream-50/95 px-3 py-1 font-mono text-[10px] uppercase tracking-widest",
            "text-ink-600 shadow-card backdrop-blur",
            "hover:border-ink-400 hover:bg-cream-100 hover:text-ink-900",
          )}
        >
          reset zoom
        </button>
      )}

      {hover && (
        <div
          className="pointer-events-none absolute z-30 rounded-chip border border-ink-300/40 bg-cream-50 px-3 py-2 font-mono text-[10px] text-ink-700 shadow-cardHover"
          style={{
            left: Math.min(hover.x + 14, width - 220),
            top: Math.max(hover.y - 60, 8),
          }}
        >
          <div className="flex items-center gap-1.5">
            <span
              className="inline-block size-2 rounded-full"
              style={{ background: hover.island.hypermethylated ? "#E89E3D" : "#5BA89A" }}
            />
            <span className="font-medium text-ink-900">{hover.island.id}</span>
            {hover.island.hypermethylated && (
              <span className="rounded-full bg-hyper-soft px-1.5 py-0.5 text-[8px] uppercase tracking-wider text-hyper">
                hyper
              </span>
            )}
          </div>
          <p className="mt-1 text-ink-500">
            {hover.island.start.toLocaleString()}-
            {hover.island.end.toLocaleString()}
            <span className="ml-2 text-ink-400">
              ({(hover.island.end - hover.island.start).toLocaleString()} bp)
            </span>
          </p>
          <p className="mt-0.5 text-[9px] uppercase tracking-widest text-ink-400">
            click to open detail
          </p>
        </div>
      )}
    </div>
  );
}
