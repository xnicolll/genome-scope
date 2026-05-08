"use client";

import { motion } from "framer-motion";
import { Database, Brain } from "lucide-react";
import { api } from "@/lib/api";
import { useApi } from "@/lib/hooks";
import { Card, CardBody, CardHeader } from "@/components/ui/Card";

function relTime(iso: string): string {
  const d = new Date(iso);
  const diff = (Date.now() - d.getTime()) / 1000;
  if (diff < 60) return `${Math.floor(diff)}s ago`;
  if (diff < 3600) return `${Math.floor(diff / 60)}m ago`;
  if (diff < 86400) return `${Math.floor(diff / 3600)}h ago`;
  return d.toISOString().slice(0, 10);
}

export function ModelsPanel() {
  const state = useApi(() => api.models(), []);

  return (
    <Card>
      <CardHeader
        kicker="Trained checkpoints"
        title="Model registry"
        right={
          <span className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
            backend/data/models/
          </span>
        }
      />
      <CardBody className="pt-2">
        {state.status === "loading" && (
          <p className="font-mono text-xs text-ink-400">loading…</p>
        )}
        {state.status === "error" && (
          <p className="font-mono text-xs text-hyper">
            could not load registry · {state.error}
          </p>
        )}
        {state.status === "ready" && state.data.n_models === 0 && (
          <div className="rounded-chip border border-dashed border-ink-300/40 bg-cream-100/40 px-4 py-6 text-center">
            <Database className="mx-auto size-5 text-ink-400" />
            <p className="mt-2 font-mono text-[11px] text-ink-500">
              no checkpoints trained yet
            </p>
            <p className="mt-1 font-mono text-[10px] text-ink-400">
              run <span className="rounded bg-cream-200 px-1">./run.sh train-beta</span>
            </p>
          </div>
        )}
        {state.status === "ready" && state.data.models.length > 0 && (
          <ul className="space-y-2">
            {state.data.models.map((m, i) => (
              <motion.li
                key={m.name}
                initial={{ opacity: 0, y: 4 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{
                  delay: i * 0.05,
                  duration: 0.25,
                  ease: [0.23, 1, 0.32, 1],
                }}
                className="rounded-chip border border-ink-300/25 bg-cream-100/60 p-4"
              >
                <div className="flex items-start justify-between gap-3">
                  <div className="min-w-0 flex-1">
                    <div className="flex items-center gap-2">
                      <Brain className="size-3.5 text-ink-500" />
                      <h4 className="font-mono text-sm font-medium text-ink-900">
                        {m.name}
                      </h4>
                      <span className="rounded-full bg-ink-300/20 px-1.5 py-0.5 font-mono text-[8px] uppercase tracking-wider text-ink-600">
                        {m.model_type}
                      </span>
                    </div>
                    <div className="mt-2 grid grid-cols-2 gap-x-4 gap-y-1 font-mono text-[10px] text-ink-500 sm:grid-cols-4">
                      <Stat label="runs" value={m.total_runs.toString()} />
                      <Stat label="samples" value={m.total_samples_seen.toLocaleString()} />
                      <Stat label="cohorts" value={m.cohorts_seen.join(", ") || "-"} />
                      <Stat
                        label="updated"
                        value={relTime(m.updated)}
                      />
                    </div>
                    {m.last_log_likelihood !== null && (
                      <p className="mt-2 font-mono text-[10px] text-ink-400">
                        last log-likelihood ·{" "}
                        <span className="text-ink-700">
                          {m.last_log_likelihood.toFixed(3)}
                        </span>
                      </p>
                    )}
                  </div>
                </div>
              </motion.li>
            ))}
          </ul>
        )}
      </CardBody>
    </Card>
  );
}

function Stat({ label, value }: { label: string; value: string }) {
  return (
    <div>
      <span className="block uppercase tracking-widest text-ink-400">{label}</span>
      <span className="text-ink-700">{value}</span>
    </div>
  );
}
