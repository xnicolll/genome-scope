import { cn } from "@/lib/cn";

export function StatusDot({
  ok,
  className,
}: {
  ok: boolean;
  className?: string;
}) {
  return (
    <span
      className={cn(
        "relative inline-block size-1.5 rounded-full",
        ok ? "bg-island" : "bg-ink-400",
        className,
      )}
      aria-hidden
    >
      {ok && (
        <span className="absolute inset-0 animate-ping rounded-full bg-island opacity-60" />
      )}
    </span>
  );
}
