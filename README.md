# GenomeScope: a Beta-emission HMM for detection of CpG islands and cancer methylation status

A tool to scan human DNA for **CpG islands**. GenomeScope finds
them with a Hidden Markov Model, then overlays real tumour data from
TCGA to flag which ones look cancer-relevant.

Tumours can grow unchecked when cancer-suppressor genes are switched off.
The main off-switch in DNA is methylation, a chemical tag (CH₃) that lands
on **CpG islands** — CG-rich stretches at the start of most genes.
Identifying which islands are silenced in a tumour, and which gene isoforms
they silence, aids in the discovery of candidate driver genes and drug targets.

**Joe Nicol · Columbia CBMFW 4761 · Spring 2026**

## Setup

Requires **macOS or Linux**, **Python 3.11+**, **Node.js 18+**, and **uv** (`brew install uv`).

```bash
git clone https://github.com/xnicolll/hmm-cpg.git
cd hmm-cpg
./run.sh    # installs everything, runs tests, downloads chr21
```

## Quickstart

```bash
./run.sh              # installs everything, runs tests, downloads chr21
./run.sh tcga         # download 20 TCGA-BRCA tumour samples (~250 MB)
./run.sh report-full  # cancer-overlay isoform report
./run.sh up           # launch FastAPI :8000 and Next.js :3000
```

Open **http://localhost:3000**. Stop with `./run.sh down`.

The dashboard is desktop-only (≥ 1280 px). Mobile is future work.

> **TCGA note.** Default cohorts use GDC open-access - no account needed.
> On macOS, `gdc-client` may be Gatekeeper-quarantined on first run; clear it once with
> `xattr -d com.apple.quarantine ./gdc-client`.

## How accurate is it?

Latest run, **full chr21 (33 Mb)**, against the UCSC reference set:

| model                          | precision | recall | F1    | islands |
| ------------------------------ | --------- | ------ | ----- | ------- |
| standard 8-state nucleotide    | 0.226     | 0.940  | 0.364 | 1,337   |
| Beta methylation (TCGA-BRCA)   | 0.110     | 0.408  | 0.173 | -       |
| **combination (intersection)** | **0.801** | 0.489  | **0.607** | -   |

To reproduce: `./run.sh full-chr21-standard` then `./run.sh full-chr21-beta`
then `./run.sh ensemble`.

## Training the model over time

The Beta HMM trains incrementally.

```bash
./run.sh train-beta       # first run: fresh init from Durbin priors
./run.sh train-beta       # second: resumes from the saved checkpoint
./run.sh train-beta       # third, fourth, ...
```

Each run appends one line to
[`backend/data/logs/training-history.csv`](backend/data/logs/training-history.csv)
with timestamp, window, sample count, log-likelihood delta, and learned
state means. The checkpoint itself lives in
`backend/data/models/beta-brca.json` and is exposed at `/api/models`.

### Daily auto-training (macOS)

`./run.sh daily` rotates through three TCGA cohorts (BRCA, LUAD,
BRCA-normal) on a different chr21 window each day. Schedule it via
`launchd`:

```bash
./run.sh schedule-install     # registers a LaunchAgent at 3 am
./run.sh schedule-status      # is it loaded? show recent runs
./run.sh schedule-remove      # unregister
```

If your Mac is asleep at 3 am, `launchd` runs the missed job at the
next wake. Override the time:

```bash
GENOMESCOPE_HOUR=13 GENOMESCOPE_MINUTE=30 ./run.sh schedule-install
```

Trigger immediately:

```bash
launchctl kickstart -k gui/$UID/com.genomescope.daily
```

## All `run.sh` commands

### setup + sanity

| command | what it does |
| --- | --- |
| `./run.sh env` | sync the Python env via `uv` |
| `./run.sh test` | run the pytest suite (97 tests) |
| `./run.sh data` | download chr21 FASTA + UCSC CpG islands + knownGene |
| `./run.sh smoke` | 200 kb pipeline run for fastest sanity check |
| `./run.sh clean` | remove generated outputs, keep raw data |

### TCGA cohorts

| command | what it does |
| --- | --- |
| `./run.sh tcga` | 20 TCGA-BRCA breast tumour samples (HM450, ~250 MB) |
| `./run.sh tcga-luad` | 20 TCGA-LUAD lung tumour samples |
| `./run.sh tcga-normal` | 10 TCGA-BRCA Solid-Tissue-Normal samples |

### model runs

| command | what it does |
| --- | --- |
| `./run.sh quick` | 2 Mb standard HMM, Durbin init only |
| `./run.sh train` | 2 Mb standard HMM with Baum-Welch |
| `./run.sh beta` | 2 Mb Beta methylation HMM, fully trained |
| `./run.sh train-beta` | incremental Beta training, resumes checkpoint |
| `./run.sh daily` | rotating cohort + window schedule (used by launchd) |
| `./run.sh full-chr21-standard` | full 33 Mb chr21 standard HMM |
| `./run.sh full-chr21-beta` | full 33 Mb chr21 Beta HMM |
| `./run.sh full-genome` | scan multiple chromosomes |
| `./run.sh full` | full chr21 pipeline (slow, ~minutes + 6 GB RAM) |

### benchmarks + reports

| command | what it does |
| --- | --- |
| `./run.sh bench` | standard vs Beta HMM, F1 side-by-side |
| `./run.sh ensemble` | three-way comparison incl. the intersection model |
| `./run.sh report` | 2 Mb isoform-aware methylation report |
| `./run.sh report-full` | full chr21 + multi-cohort report (the headline) |

### dev servers + scheduling

| command | what it does |
| --- | --- |
| `./run.sh up` | start FastAPI + Next.js in the background |
| `./run.sh down` | stop both servers |
| `./run.sh logs api\|web` | tail dev-server logs |
| `./run.sh schedule-install` | install daily launchd agent |
| `./run.sh schedule-remove` | uninstall the launchd agent |
| `./run.sh schedule-status` | is the agent loaded + show recent runs |

## Layout

```
hmm-cpg/
├── backend/                  Python HMM engine + FastAPI
│   └── src/genomescope/
│       ├── hmm/              forward-backward, Viterbi, Baum-Welch
│       ├── data/             FASTA, TCGA, HM450 probes, knownGene
│       ├── analysis/         promoter overlap, evaluation, report
│       ├── scripts/          download, ensemble, full-genome, story
│       └── api/              FastAPI server
├── frontend/                 Next.js 14 dashboard
│   ├── app/                  landing, dashboard, upload routes
│   └── components/           D3 track + modular panels
└── run.sh                    one-command orchestrator
```

## Data sources (all open-access)

- **chr21 FASTA** - [UCSC goldenPath hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz)
- **CpG island reference** - UCSC `cpgIslandExt` track
- **Gene annotations** - UCSC `knownGene` (GENCODE V49)
- **TCGA methylation** - GDC Portal, BRCA + LUAD + BRCA-normal HM450 Level-3 β values
- **HM450 probe → genomic coords** - [Zhou Lab InfiniumAnnotation](https://github.com/zhou-lab/InfiniumAnnotationV1)
- **Gene symbols** - [MyGene.info](https://mygene.info), disk-cached on first call
- **Tumour-suppressor list** - built-in 40-gene fallback (Vogelstein 2013 + Sanger CGC)

## Stack

**Backend** - Python 3.11 via `uv` · NumPy · SciPy · BioPython · pandas
· pyranges · scikit-learn · FastAPI · Uvicorn.

**Frontend** - Next.js 14 (App Router, requires Node 18+) · TypeScript · Tailwind
· D3.js v7 · Framer Motion · Geist fonts · lucide-react.

## Troubleshooting

- **Dashboard won't load** — `rm -rf frontend/.next && ./run.sh up`
- **`gdc-client` blocked on macOS** — `xattr -d com.apple.quarantine ./gdc-client`
- **`uv` not found** — `brew install uv` (or see [astral.sh/uv](https://astral.sh/uv))
- **Port 3000 / 8000 already in use** — `./run.sh down` to clean up stale dev servers
