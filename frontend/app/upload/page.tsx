import Link from "next/link";
import { ArrowLeft } from "lucide-react";
import { UploadPanel } from "@/components/UploadPanel";
import { ModelsPanel } from "@/components/ModelsPanel";

export const metadata = { title: "Upload - GenomeScope" };

export default function UploadPage() {
  return (
    <main className="min-h-screen px-5 py-8 md:px-10 md:py-12">
      <div className="mx-auto max-w-5xl">
        <header className="mb-10 flex flex-wrap items-center justify-between gap-4">
          <div className="flex items-center gap-4">
            <Link
              href="/"
              className="press inline-flex size-9 items-center justify-center rounded-full border border-ink-300/40 bg-cream-50 text-ink-500 hover:border-ink-400 hover:bg-cream-100 hover:text-ink-900"
              aria-label="Back to landing page"
            >
              <ArrowLeft className="size-4" />
            </Link>
            <div>
              <p className="font-mono text-[10px] uppercase tracking-widest text-ink-400">
                GenomeScope · custom analysis
              </p>
              <h1 className="mt-0.5 text-2xl font-medium tracking-tight text-ink-900">
                Bring your own FASTA
              </h1>
            </div>
          </div>
          <Link
            href="/dashboard"
            className="press inline-flex items-center gap-2 rounded-full border border-ink-300/40 bg-cream-50 px-4 py-2 font-mono text-[10px] uppercase tracking-widest text-ink-500 hover:border-ink-400 hover:bg-cream-100 hover:text-ink-900"
          >
            chr21 dashboard
          </Link>
        </header>

        <div className="grid gap-4 lg:grid-cols-5">
          <div className="lg:col-span-3">
            <UploadPanel />
          </div>
          <div className="lg:col-span-2">
            <ModelsPanel />
          </div>
        </div>

        <p className="mt-10 text-center font-mono text-[10px] text-ink-400">
          Sequences are processed in-memory · 5 MB upload limit · standard 8-state HMM only
        </p>
      </div>
    </main>
  );
}
