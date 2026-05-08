"use client";

import { useCallback, useEffect, useRef, useState } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { Upload, FileText, Play, AlertCircle, CheckCircle2, History } from "lucide-react";
import { api, type Checkpoint, type PredictResult } from "@/lib/api";
import { useApi } from "@/lib/hooks";
import { Card, CardBody, CardHeader } from "@/components/ui/Card";
import { cn } from "@/lib/cn";

const SAMPLE_FASTA = `>example_island
TTAATATAATTTATTATTTTATTTAAATTTATATTATATTAATAATATATTAAATTATTTTTTTTATATAAATTAATTTT
TTTATATTTTTTATTTTATATTATAATATAATATTTTATTTAATAATAAATTAAAAATAATAATAATTTATTATTAATAA
AATTAATATATAAAATAAATTTTAATAAAATATTAATATAAAATATTTAAATTTATTATTATAATATTAATTATAAAATA
TTTTAATATTAAATAATAATAAAATATAATTTTAATTTATAATATATAATTAAATATAATAATAATAATATAAATAATAA
ATTTAATAATAATAATATAAAATTAATTTTAATATATTTAAAATAATAAAATAATATAATAAATTTAATATTAAATATAT
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
TTAATATATAATATAATATTTAATAATAATATAAATAATAAATTTAATAATAATAATATAAAATTAATTTTAATATATTT
AAAATAATAAAATAATATAATAAATTTAATATTAAATATATTAATTAATAATATAATATTAATAAAAATAATAATAATAA`;

export function UploadPanel() {
  const [text, setText] = useState("");
  const [model, setModel] = useState("standard");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<PredictResult | null>(null);
  const [dragOver, setDragOver] = useState(false);
  const fileRef = useRef<HTMLInputElement | null>(null);

  const modelsState = useApi(() => api.models(), []);
  const standardCheckpoints: Checkpoint[] =
    modelsState.status === "ready"
      ? modelsState.data.models.filter((m) => m.model_type === "standard")
      : [];

  const onSubmit = useCallback(async () => {
    if (!text.trim()) {
      setError("paste a FASTA sequence or load the example");
      return;
    }
    setError(null);
    setBusy(true);
    setResult(null);
    try {
      const r = await api.predict(text, model);
      setResult(r);
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }, [text, model]);

  const onFile = useCallback(async (file: File) => {
    if (file.size > 5 * 1024 * 1024) {
      setError("file exceeds 5MB upload limit");
      return;
    }
    const txt = await file.text();
    setText(txt);
    setError(null);
  }, []);

  const onDrop = useCallback(
    async (e: React.DragEvent<HTMLDivElement>) => {
      e.preventDefault();
      setDragOver(false);
      const f = e.dataTransfer.files?.[0];
      if (f) await onFile(f);
    },
    [onFile],
  );

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.metaKey || e.ctrlKey) && e.key === "Enter") {
        e.preventDefault();
        if (!busy) onSubmit();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [busy, onSubmit]);

  const charCount = text.replace(/[^ACGTNacgtn]/g, "").length;

  return (
    <Card>
      <CardHeader
        kicker="Run analysis"
        title="Upload FASTA"
        right={
          <span className="hidden font-mono text-[10px] uppercase tracking-widest text-ink-400 md:block">
            paste · drop · file pick
          </span>
        }
      />
      <CardBody className="space-y-4">
        <div
          onDragOver={(e) => {
            e.preventDefault();
            setDragOver(true);
          }}
          onDragLeave={() => setDragOver(false)}
          onDrop={onDrop}
          className={cn(
            "relative rounded-chip border border-dashed transition-colors",
            dragOver
              ? "border-accent bg-accent-50/40"
              : "border-ink-300/40 bg-cream-100/40",
          )}
        >
          <textarea
            value={text}
            onChange={(e) => setText(e.target.value)}
            placeholder="> chr_user&#10;ACGTACGT..."
            spellCheck={false}
            className={cn(
              "block min-h-[180px] w-full resize-y rounded-chip bg-transparent",
              "px-4 py-3 font-mono text-[11px] leading-relaxed text-ink-800",
              "outline-none placeholder:text-ink-400/70",
              "focus-visible:ring-2 focus-visible:ring-accent/50",
            )}
          />
          <div className="flex items-center justify-between border-t border-ink-300/20 px-3 py-2">
            <div className="flex items-center gap-3 font-mono text-[10px] uppercase tracking-widest text-ink-400">
              <button
                type="button"
                onClick={() => fileRef.current?.click()}
                className="press inline-flex items-center gap-1.5 hover:text-ink-700"
              >
                <FileText className="size-3" />
                file
              </button>
              <button
                type="button"
                onClick={() => setText(SAMPLE_FASTA)}
                className="press inline-flex items-center gap-1.5 hover:text-ink-700"
              >
                <Upload className="size-3" />
                sample
              </button>
              <button
                type="button"
                onClick={() => {
                  setText("");
                  setResult(null);
                  setError(null);
                }}
                className="press hover:text-ink-700"
              >
                clear
              </button>
            </div>
            <span className="font-mono text-[10px] text-ink-500">
              {charCount.toLocaleString()} bp
            </span>
          </div>
          <input
            ref={fileRef}
            type="file"
            accept=".fa,.fasta,.txt,.fna"
            className="sr-only"
            onChange={(e) => {
              const f = e.target.files?.[0];
              if (f) onFile(f);
              e.target.value = "";
            }}
          />
        </div>

        <div className="flex flex-wrap items-center gap-3">
          <div className="flex flex-1 items-center gap-2">
            <label className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
              Model
            </label>
            <select
              value={model}
              onChange={(e) => setModel(e.target.value)}
              className={cn(
                "press rounded-full border border-ink-300/40 bg-cream-50",
                "px-3 py-1.5 font-mono text-[11px] text-ink-800",
                "outline-none focus-visible:ring-2 focus-visible:ring-accent/50",
              )}
            >
              <option value="standard">standard (PRD defaults)</option>
              {standardCheckpoints.map((c) => (
                <option key={c.name} value={c.name}>
                  {c.name} · {c.total_runs} runs
                </option>
              ))}
            </select>
          </div>
          <button
            onClick={onSubmit}
            disabled={busy || !text.trim()}
            className={cn(
              "press inline-flex items-center gap-2 rounded-full px-5 py-2",
              "text-sm font-medium transition-colors",
              busy || !text.trim()
                ? "cursor-not-allowed bg-ink-300/30 text-ink-400"
                : "bg-ink text-cream-50 hover:bg-ink-700",
            )}
          >
            {busy ? (
              <>
                <Spinner />
                running…
              </>
            ) : (
              <>
                <Play className="size-3.5" />
                run analysis
                <kbd className="ml-1 hidden rounded border border-cream-50/30 px-1 py-px font-mono text-[9px] opacity-70 md:inline">
                  ⌘↵
                </kbd>
              </>
            )}
          </button>
        </div>

        <AnimatePresence>
          {error && (
            <motion.div
              initial={{ opacity: 0, y: -4 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: -4 }}
              transition={{ duration: 0.18, ease: [0.23, 1, 0.32, 1] }}
              className="flex items-start gap-2 rounded-chip border border-hyper/40 bg-hyper-soft/40 px-4 py-3"
            >
              <AlertCircle className="mt-0.5 size-3.5 shrink-0 text-hyper" />
              <p className="font-mono text-[11px] text-ink-800">{error}</p>
            </motion.div>
          )}
        </AnimatePresence>

        <AnimatePresence>
          {result && <ResultBlock result={result} />}
        </AnimatePresence>
      </CardBody>
    </Card>
  );
}

function ResultBlock({ result }: { result: PredictResult }) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 6 }}
      animate={{ opacity: 1, y: 0 }}
      exit={{ opacity: 0, y: 4 }}
      transition={{ duration: 0.25, ease: [0.23, 1, 0.32, 1] }}
      className="space-y-3 rounded-chip border border-island/30 bg-island-soft/40 p-4"
    >
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <CheckCircle2 className="size-4 text-island" />
          <p className="font-mono text-xs font-medium text-ink-900">
            {result.n_islands} island{result.n_islands === 1 ? "" : "s"} found
          </p>
        </div>
        <div className="flex items-center gap-3 font-mono text-[10px] text-ink-500">
          <span>{result.elapsed_seconds.toFixed(2)}s</span>
          <span>{result.window_length.toLocaleString()} bp</span>
        </div>
      </div>

      <MiniTrack length={result.window_length} islands={result.islands} />

      {result.islands.length > 0 && (
        <ul className="space-y-1.5">
          {result.islands.slice(0, 8).map((isl, i) => {
            const len = isl.end - isl.start;
            return (
              <li
                key={i}
                className="flex items-center justify-between rounded-chip bg-cream-50 px-3 py-2 font-mono text-[11px]"
              >
                <span className="text-ink-700">
                  [{isl.start.toLocaleString()} - {isl.end.toLocaleString()})
                </span>
                <span className="text-ink-500">{len.toLocaleString()} bp</span>
              </li>
            );
          })}
          {result.islands.length > 8 && (
            <li className="px-3 font-mono text-[10px] text-ink-400">
              + {result.islands.length - 8} more
            </li>
          )}
        </ul>
      )}

      {result.checkpoint_used && (
        <p className="flex items-center gap-1.5 font-mono text-[10px] text-ink-500">
          <History className="size-3" />
          ran with checkpoint {result.checkpoint_used}
        </p>
      )}
    </motion.div>
  );
}

function MiniTrack({
  length,
  islands,
}: {
  length: number;
  islands: { start: number; end: number }[];
}) {
  return (
    <div className="relative h-3 w-full overflow-hidden rounded-full bg-ink-300/15">
      {islands.map((isl, i) => {
        const x = (isl.start / length) * 100;
        const w = Math.max(0.3, ((isl.end - isl.start) / length) * 100);
        return (
          <div
            key={i}
            className="absolute top-0 h-full rounded-full bg-island"
            style={{ left: `${x}%`, width: `${w}%` }}
          />
        );
      })}
    </div>
  );
}

function Spinner() {
  return (
    <span
      className="size-3.5 animate-spin rounded-full border-2 border-current border-t-transparent"
      aria-hidden
    />
  );
}
