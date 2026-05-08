import { cn } from "@/lib/cn";
import type { ReactNode } from "react";

export function StatBlock({
  label,
  value,
  sub,
  tone,
  className,
}: {
  label: string;
  value: ReactNode;
  sub?: ReactNode;
  tone?: "default" | "accent" | "island" | "hyper" | "promoter";
  className?: string;
}) {
  return (
    <div
      className={cn(
        "rounded-chip border border-ink-300/20 bg-cream-50 px-5 py-4",
        className,
      )}
    >
      <p className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
        {label}
      </p>
      <p
        className={cn(
          "mt-1 font-mono text-3xl tracking-tight",
          tone === "accent" && "text-ink-900",
          tone === "island" && "text-island",
          tone === "hyper" && "text-hyper",
          tone === "promoter" && "text-promoter",
          !tone && "text-ink-900",
        )}
      >
        {value}
      </p>
      {sub && (
        <p className="mt-1 font-mono text-[11px] leading-relaxed text-ink-500">
          {sub}
        </p>
      )}
    </div>
  );
}
