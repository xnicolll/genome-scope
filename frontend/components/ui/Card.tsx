import { cn } from "@/lib/cn";
import type { HTMLAttributes, ReactNode } from "react";

export function Card({ className, ...props }: HTMLAttributes<HTMLDivElement>) {
  return (
    <div
      className={cn(
        "flex h-full flex-col overflow-hidden rounded-card border border-ink-300/30 bg-cream-50 shadow-card",
        "transition-shadow duration-200 ease-out-strong hover:shadow-cardHover",
        className,
      )}
      {...props}
    />
  );
}

export function CardHeader({
  kicker,
  title,
  right,
}: {
  kicker?: string;
  title: ReactNode;
  right?: ReactNode;
}) {
  return (
    <div className="flex shrink-0 items-start justify-between gap-4 px-7 pt-6">
      <div>
        {kicker && (
          <p className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
            {kicker}
          </p>
        )}
        <h3 className="mt-1 text-lg font-medium tracking-tight text-ink-900">
          {title}
        </h3>
      </div>
      {right}
    </div>
  );
}

export function CardBody({
  className,
  ...props
}: HTMLAttributes<HTMLDivElement>) {
  return (
    <div
      className={cn("min-h-0 flex-1 px-7 pb-6 pt-4", className)}
      {...props}
    />
  );
}
