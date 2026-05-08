import Link from "next/link";
import { ArrowRight } from "lucide-react";

export default function LandingPage() {
  return (
    <main className="min-h-screen px-6 py-12 md:px-12 md:py-20">
      <div className="mx-auto max-w-6xl">
        <header className="flex items-center gap-2">
          <span className="size-2.5 rounded-full bg-accent" />
          <span className="font-mono text-xs uppercase tracking-widest text-ink-500">
            GenomeScope
          </span>
        </header>

        <section className="mt-24 max-w-3xl">
          <p className="font-mono text-xs uppercase tracking-widest text-ink-500">
            CBMFW 4761 · Spring 2026 · Joe Nicol
          </p>
          <h1 className="mt-6 text-balance text-5xl font-medium leading-[1.05] tracking-tight text-ink-900 md:text-6xl">
            Finding CpG islands the way{" "}
            <span className="italic text-ink-700">cancer sees them</span>.
          </h1>
          <p className="mt-6 max-w-2xl text-pretty text-lg leading-relaxed text-ink-600">
            An HMM-based CpG island detector that replaces nucleotide emissions
            with Beta-distributed methylation emissions trained directly on TCGA.
          </p>
          <div className="mt-10 flex flex-wrap gap-3">
            <Link
              href="/dashboard"
              className="press inline-flex items-center gap-2 rounded-full bg-ink px-6 py-3 text-sm font-medium text-cream-50 transition-colors duration-200 ease-out-strong hover:bg-ink-700"
            >
              Explore the dashboard
              <ArrowRight className="size-4 transition-transform duration-200 ease-out-strong group-hover:translate-x-0.5" />
            </Link>
            <Link
              href="/upload"
              className="press inline-flex items-center gap-2 rounded-full border border-ink-300/60 bg-cream-50 px-6 py-3 text-sm font-medium text-ink-800 transition-colors duration-200 ease-out-strong hover:border-ink-400 hover:bg-cream-100"
            >
              Upload your own FASTA
            </Link>
          </div>
        </section>
      </div>
    </main>
  );
}
