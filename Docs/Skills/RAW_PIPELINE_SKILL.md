# Raw-read processing (`igvfagent raw-pipeline`)

Turns an IGVF accession into an analysed count matrix. This is the skill to
reach for when someone asks to **analyse or process** a dataset, as opposed
to asking what it is — `explain_dataset` answers the second question and
should not be mistaken for an answer to the first.

```bash
# What would it take? Resolves, inventories, pairs reads, infers chemistry.
# Downloads only the (tiny) seqspec.
igvfagent raw-pipeline plan IGVFDS3532MONX

# Do it. Prefers a published matrix; aligns only when there isn't one.
igvfagent raw-pipeline run IGVFDS3532MONX

# Show every command and transfer without performing any.
igvfagent raw-pipeline run IGVFDS3532MONX --dry-run
```

## Long runs: detach

Aligning a real dataset takes tens of minutes to hours — IGVFDS3532MONX is
45.6 GB of reads. The agent runs a tool as a blocking subprocess with no
timeout, inside a web request, so a synchronous run holds the browser open
for the whole job and the result is lost when the socket drops. That is why
a large dataset could only ever be *planned*.

```bash
igvfagent raw-pipeline run IGVFDS3532MONX --detach   # returns a job id at once
igvfagent raw-pipeline status                        # all recent jobs
igvfagent raw-pipeline status <job-id> --tail 30     # one job, more log
```

The job runs in its own session and outlives the conversation, so "did that
analysis finish?" is answerable in a later session. Use `--detach` for
anything beyond a few GB; a small set is fine synchronously.

## Caching — why the second run is instant

Entering the same accession twice is the demo case, and it used to cost the
full hour again: reads went to a per-run timestamped directory, so a repeat
re-downloaded 45.6 GB *and* kept a second copy of it.

Two caches now sit under `Data/RawPipeline/`:

| Path | Holds | Keyed by |
|---|---|---|
| `_fastq/` | downloaded reads | **file** accession (`IGVFFI…`) |
| `_results.json` | completed analyses | accession + technology + workflow + reference |
| `../References/kb/` | kb indices | reference name |

A Portal file's bytes never change, so the read cache needs no invalidation,
and each run hard-links to the single copy — one copy on disk however many
times a dataset is analysed. Downloads print `Cached IGVFFI… (1.87 GB)`
rather than silently skipping, so a fast run is visibly fast instead of
suspicious.

A repeat of an identical analysis returns the existing matrix in **0
seconds** (measured on IGVFDS3532MONX: 730,828,449 reads, 70.9%
pseudoaligned) and still runs the single-cell step, so a demo shows real
UMAPs and markers rather than a cache message. `--no-reuse` forces a
recompute from the reads.

The result index is a hint, not a fact: the matrix it names is checked on
disk before the entry is trusted, so deleting a run directory to reclaim
space makes the entry disappear rather than leaving a confident pointer to
a file that is gone.

**Reclaiming space.** The read cache grows without bound. It is the largest
thing on disk and the cheapest to lose — everything in it re-downloads:

```bash
du -sh Data/RawPipeline/_fastq                  # how much is held
rm -f Data/RawPipeline/_fastq/IGVFFI*.fastq.gz  # drop all cached reads
```

Results (`processed.h5ad`, `markers.csv`, plots) are megabytes, not
gigabytes; keep those.

## Routes

Chosen automatically and printed as `ROUTE:` before anything large moves.

| Route | When | What happens |
|---|---|---|
| `matrix_on_set` | the FileSet carries a count matrix | download it, analyse it |
| `matrix_derived` | an AnalysisSet derived from it publishes one | download that, analyse it |
| `align` | only raw reads exist | download FASTQs → `kb count` → analyse |
| `none` | neither reads nor matrix | say so, do nothing |

A published matrix wins by default because re-aligning reads to reproduce
somebody's published matrix is hours of compute for no new information. Pass
`--force-align` when reprocessing is the actual point.

Where several matrices exist the choice is deliberate, not incidental:
`.h5ad` first (the single-cell pipeline loads it natively), then `.h5`,
`.mtx`, `.tar`, smallest first within a format. On IGVFDS9875NBZW that picks
a 3.7 GB h5ad over a 13 GB tar, and says so.

## The aligner

`kb-python` (kallisto | bustools), installed with the `align` extra:

```bash
pip install 'igvfagent[align]'      # included in [all]
```

It bundles its own kallisto and bustools binaries per platform, so there is
no compiler, no conda, and no separate aligner install — which is what makes
an aligner shippable in a pip-built image at all. Prebuilt human and mouse
indices come from `kb ref -d`, cached under `Data/References/kb/`.

`--workflow kite` assigns feature barcodes (CRISPR sgRNAs). `--workflow nac`
produces nascent/mature matrices.

**What this is not.** kallisto quantifies against a transcriptome. That is
the right tool for the scRNA-seq readout of a CRISPR screen, and wrong for
work that genuinely needs genomic coordinates (structural variants,
allele-specific pileups). `plan` names the route rather than implying the
matrix answers every question.

## Chemistry detection

Inferred from the seqspec, and reported with its evidence. Read **roles are
derived from read length, not from `read_id`**, because the label is not
reliable — real IGVF seqspecs use `R1`/`R2`, `Read1`/`Read2`, and on the
multiome sets the *file accession itself*. A 26 or 28 bp read is the cell
barcode + UMI (10x v2 / v3), a read of ≤12 bp is a sample index, and the long
read is cDNA. A library with a long read and no barcode read is bulk, and is
quantified with `-x BULK --parity single` rather than being reported as
"unknown chemistry".

Override with `--technology` (any kallisto technology string; `kb --list`
enumerates them).

## Read pairing

`kb count` consumes reads strictly positionally — barcode file then
biological file, per pair — so a mispairing does not error, it quantifies the
wrong barcodes against the wrong cDNA. Pairing therefore comes from the
seqspec's per-read file lists where a seqspec exists (**all** of them: a
multiome set carries one per modality), falling back to
`flowcell_id + lane + sequencing_run`. Anything that fails to pair is
reported as `UNPAIRED` rather than appended and hoped for.

## Guards

- `--max-download-gb` (default 100) stops before transferring more.
- Controlled-access files are listed up front; they need an approved DUA.
- Credentials come from `_credentials.py`, so unreleased datasets resolve
  when an access-key pair is configured. Without one, an unreleased
  accession is reported as unresolvable rather than described from a
  free-text search that matched something else.
