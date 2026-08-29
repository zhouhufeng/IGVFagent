# Pathway databases: `pathwaydb`

Pulls the **current** releases of KEGG, Reactome and WikiPathways, normalises
every gene identifier into one namespace, unifies pathways that different
databases describe differently, and loads the result into the local knowledge
graph. Every stored fact records which databases assert it, so agreement
between databases stays queryable rather than collapsing into a duplicate row.

The integration method comes from **IntPath** (Zhou H, Jin J, Zhang H, Yi B,
Wozniak M, Wong L. *BMC Systems Biology* 2012, 6(Suppl 2):S2,
[doi:10.1186/1752-0509-6-S2-S2](https://doi.org/10.1186/1752-0509-6-S2-S2)).
The data is pulled fresh — IntPath's published release is not redistributed
here, and `igvfagent intpath` still serves that 2012 snapshot separately for
comparison.

---

## Why not just query a pathway API

A live API answers "which pathways contain this gene" one database at a time,
and each database names the same biology differently. Asking three databases
gives three overlapping answers to reconcile by hand. Integrating once gives a
single answer that also says *how many independent databases back it* — which
is evidence a single-database answer cannot express.

Pulling in bulk also makes the result reproducible: each download is stamped
with its release identifier and SHA-256, so a figure traces back to the exact
data behind it.

### Current releases vs. the 2012 snapshot

| | IntPath 2012 | pulled 2026-08-29 |
|---|---|---|
| pathways | 582 | **4,016** |
| gene–pathway memberships | 23,873 | **211,120** |
| typed gene–gene relations | 71,116 pairs | **73,925** |
| genes | 5,748 | **14,870** |

IntPath's original pipeline cannot simply be re-run: the WikiPathways SOAP
webservice it consumed has been retired (`listPathways` now returns 404). The
method is reimplemented here against the routes that are supported today.

---

## CLI

```bash
igvfagent pathwaydb sources     # databases, licences, what each contributes
igvfagent pathwaydb pull        # download current releases (~180 MB)
igvfagent pathwaydb build       # normalise, unify, load into the KG
igvfagent pathwaydb status      # what is cached, how old, which release
igvfagent pathwaydb query       # pathways for a gene list
igvfagent pathwaydb evaluate    # score the unification criteria
```

### `pull` — download current releases

```bash
igvfagent pathwaydb pull
igvfagent pathwaydb pull --kgml            # + KEGG maps, needed for typed relations
igvfagent pathwaydb pull --sources kegg    # one database only
igvfagent pathwaydb pull --force           # ignore the 7-day cache freshness window
```

| flag | meaning |
|---|---|
| `--sources` | comma list: `kegg`, `reactome`, `wikipathways` (default all) |
| `--kgml` | also cache the 372 human KGML maps (~370 requests, paced at 3/s) |
| `--force` | re-download even when the cache is fresh |
| `--cache` | cache directory (default `Data/Reference/PathwayDB/`, gitignored) |

Downloads are atomic — an interrupted pull never leaves a truncated file that
would later parse as valid. A partial pull merges into the existing manifest
rather than erasing provenance for sources it did not touch.

### `build` — normalise, unify, ingest

```bash
igvfagent pathwaydb build --relations
igvfagent pathwaydb build --merge lcs             # IntPath's published rule alone
igvfagent pathwaydb build --min-sources 2         # only store corroborated facts
igvfagent pathwaydb build --extra-gmt HumanCyc=~/humancyc.gmt
```

| flag | meaning |
|---|---|
| `--relations` | also build typed gene–gene relations (KEGG KGML + Reactome interactors) |
| `--merge` | unification criterion; see [Unification](#unification-which-pathways-are-the-same-pathway) |
| `--extra-gmt LABEL=PATH` | integrate a local GMT as another source database (repeatable) |
| `--min-sources N` | only ingest facts asserted by at least N databases |
| `--within-db` | also unify names inside one database, as IntPath does (off by default) |
| `--no-hierarchy` | ignore Reactome's own parent/child hierarchy when unifying |
| `--no-kg` | write files without touching the knowledge graph |
| `--no-pull` | fail rather than downloading if the cache is empty |

Outputs land in `Docs/PathwayDB/<timestamp>_<label>/`:

- `memberships.csv` — `gene, pathway, kind, evidence, score, ids`
- `relations.csv` — typed gene–gene edges
- `merges.csv` — **every unification decision with its scores**, so a wrong
  merge is findable and reversible rather than silent
- `summary.json` — counts, per-source totals, release string

### `query` — pathways for a gene list

```bash
igvfagent pathwaydb query --genes CLU,BIN1,PICALM,APOE,TREM2 --top 10 --verbose
```

Reports each pathway, how many query genes hit it, and **which databases**
assert it. Runs against the cache — no network.

### `evaluate` — score the unification criteria

```bash
igvfagent pathwaydb evaluate --agreement --out Docs/PathwayDB/method_evaluation.json
```

See [Choosing a criterion by measurement](#choosing-a-criterion-by-measurement).

### `status` / `sources`

`status` prints the cached files with size, age, SHA-256 prefix and release
string, plus the last build's counts. `sources` prints each database, its
host, its licence, and what it contributes — including what is *not* included
and why.

---

## Agent tools

Both are exposed to the LLM and survive the cross-backend tool cap:

| tool | maps to | use |
|---|---|---|
| `pathway_db_query` | `pathwaydb query` | which pathways a gene set shares, and which databases agree |
| `pathway_db_refresh` | `pathwaydb build` | pull and integrate when pathway data is missing or stale |

The visualiser picks the data up automatically: `pathway-viz network --sources
local` draws from the knowledge graph, which now holds the integrated
membership and relation edges.

---

## Method

### 1. Identifier normalisation

Every identifier resolves through **NCBI `gene_info`** (Entrez → symbol, with
synonyms and Ensembl cross-references), *not* through any pathway database's
own gene list — that keeps the namespace independent of the sources being
merged. 193,813 human records; an unrecognised but symbol-shaped token is kept
rather than dropped, since an unmapped symbol is still a usable fact.

Two source quirks are handled explicitly:

- **Reactome's interactor file ships an entirely empty Entrez column** in the
  current release. Its Ensembl column is the real key, indexed from
  `gene_info`'s `dbXrefs`. Parsing the documented column yields zero rows.
- **WikiPathways' GMT generator flattens `&#39;` to a bare " 39 "**, shipping
  "Alzheimer 39 s disease". Repaired at parse time.

### 2. Unification: which pathways are the same pathway

IntPath compares pathway **names** by longest common subsequence:

```
alignment score = number of aligned characters (the LCS length)
alignment ratio = 2 × score / (len(a) + len(b))

accept if EITHER  (1) score > len(shorter) − 1  AND  ratio ≥ 0.5
              OR  (2) ratio > 0.91
```

then filters with an "error-prone words pair list", groups with a disjoint-set
structure, and names each group by its **shortest** member name. LCS is
computed with a bit-parallel algorithm (verified identical to the dynamic
program over 4,000 random strings) because unification compares millions of
name pairs.

**The published rule needs guards on modern data.** IntPath unified 582 flat,
similarly-named pathways. Reactome — which IntPath did not include — is
explicitly hierarchical and adds three things the 2012 vocabulary lacked:

| problem | example | guard |
|---|---|---|
| single-word generic parents | `disease` is a subsequence of `alzheimer disease` at ratio 0.58, chaining Alzheimer + Chagas + prion into one group | rule (1) requires the shorter name to have ≥ 2 words |
| explicit negative/defective forms | `apoptosis` vs `suppression of apoptosis` | qualifier list — the unary analogue of IntPath's word-pair list |
| family members that differ by a number | `vitamin b12 metabolism` vs `vitamin b6 metabolism` (ratio 0.94) | numeric/roman-numeral token mismatch rejects |
| real parent/child nesting | `cell cycle` vs `cell cycle mitotic` | **Reactome's own hierarchy file**, 9,853 pairs — curated fact, not a heuristic |

The hierarchy guard is the important one: where a source database states the
answer, it is read rather than inferred, and it outranks every similarity
score computed here.

### 3. Choosing a criterion by measurement

There is no gold standard for "these two pathways are the same", so
`pathwaydb evaluate` reports two objective signals instead of an invented
accuracy number:

- **hierarchy violations** — merges Reactome declares to be parent and child.
  Known-wrong. A precision signal.
- **identical-name recall** — of the cross-database pairs whose names are
  already identical (unambiguously the same pathway), how many are found.
  A recall signal, and **biased in favour of the name-based criteria**, which
  get it for free; it is meaningful for comparing the gene-based criteria.

Measured on the 2026-08-29 releases, 4,172 pathways across 3 databases:

| criterion | merges | hierarchy violations | violation rate | identical-name recall | sec |
|---|---|---|---|---|---|
| **consensus** (default) | 171 | **0** | **0.0%** | 100% | 46 |
| `lcs` (IntPath as published) | 207 | 7 | 3.4% | 100% | 38 |
| `tokens` | 156 | 0 | 0.0% | 100% | 7 |
| `jaccard` (gene overlap) | 67 | 2 | 3.0% | 21.5% | 5 |
| `containment` | 4,834 | 126 | 2.6% | 63.3% | 2 |
| `hypergeom` | 46,523 | 318 | 0.7% | 93.7% | 14 |

Findings worth keeping:

- **No single feature separates the two classes.** Gene-set Jaccard was the
  obvious candidate and fails outright: true merges run as low as 0.034 while
  known-wrong nested pairs reach 0.471. Any threshold trades one error for the
  other. Statistical significance (`hypergeom`) is worse — 46,523 merges,
  because two large pathways sharing genes is significant without being the
  same pathway.
- **LCS's false positives are all one kind:** specialisations that read as
  near-identical strings — `cell cycle` / `cell cycle mitotic`,
  `base excision repair` / `PCNA-dependent long patch base excision repair`.
- **Character similarity and word similarity fail differently.** LCS never
  proposes `signaling by hippo` for `hippo signaling pathway` — the words are
  reordered, so the characters do not align. Word overlap never proposes a
  pure spelling variant. Their agreement is only 0.50.

`consensus` (the default) uses that: it takes candidates from **both** the LCS
best hit and the best word-overlap hit — which is what keeps recall — then
merges only when at least 2 of 3 independent signals agree (name alignment
ratio > 0.91, token overlap ≥ 0.5, gene Jaccard ≥ 0.20). Identical names pass
on the name alone; there is nothing to disambiguate. The result is 171 merges
at zero known violations: more merges than `tokens` at the same precision, and
better precision than `lcs` at comparable recall.

Use `--merge lcs` to reproduce IntPath's published rule exactly.

### 4. Typed gene–gene relations

`--relations` builds them from two places, since the membership endpoints
carry no relation type:

- **KEGG KGML** maps, expanded across component families as KEGG renders them
- **Reactome's** curated interactor file

Both map onto IntPath's unified vocabulary, kept distinct rather than
collapsed into `interacts_with` — the difference between a physical
interaction and transcriptional regulation is most of what makes a typed edge
worth storing:

| relation | KG edge | meaning |
|---|---|---|
| `PPrel` | `interacts_with` | protein–protein interaction |
| `ECrel` | `enzyme_relation` | sequential catalysis (shared metabolite) |
| `GPrel` | `in_same_complex` | gene product / complex membership |
| `GErel` | `regulates_expression` | gene-expression regulation |
| `PCrel` | `acts_on_compound` | protein–compound relation |

Current build: 73,925 relations — PPrel 54,002, ECrel 15,996, GErel 3,197,
PCrel 730.

### 5. Provenance and idempotency

- Every download records release id, byte size and SHA-256 in `manifest.json`.
- Every pathway node records the **build that wrote it**. A rebuild prunes
  nodes orphaned by an upstream rename instead of leaving duplicates behind.
  Pathway nodes written by other skills carry no stamp and are never touched.
- `merges.csv` records every unification decision and its scores.

---

## Sources, licences, and what is missing

| database | route | licence | contributes |
|---|---|---|---|
| KEGG | `rest.kegg.jp` | academic use; polite rate limit | membership + typed KGML relations |
| Reactome | `reactome.org/download/current` | CC-BY | membership + curated interactions + hierarchy |
| WikiPathways | `data.wikipathways.org/current/gmt` | CC0 | membership (GMT release) |

Data is fetched at runtime and **never vendored** into the repository. The
cache (`Data/Reference/PathwayDB/`, ~180 MB) and run outputs
(`Docs/PathwayDB/`) are gitignored.

WikiPathways publishes monthly, on the 10th; a refresh shortly after is the
natural cadence.

### BioCyc / HumanCyc — not pulled, and why

IntPath's third source is **not** fetched. Verified 2026-08-29: BioCyc's
download page states data access requires "a BioCyc subscription costing at
least $5,000", and the human pathway listing returns no pathways to an
anonymous client. IGVFagent therefore cannot pull it reproducibly, and
`pathwaydb sources` says so rather than quietly returning two databases where
three are expected.

If **you** hold a licence, integrate your own export — it flows through the
same normalisation, unification and KG ingestion as the fetched sources, and
appears in the per-fact source list like any other database:

```bash
igvfagent pathwaydb build --extra-gmt HumanCyc=~/humancyc.gmt
```

Verified end to end: a HumanCyc "Mismatch repair" set unifies with KEGG
`hsa03430`, Reactome `R-HSA-5358508` and WikiPathways `WP531`, retaining all
four pathway identifiers. Nothing licensed is redistributed — the file stays
on your machine.

`--extra-gmt` is not BioCyc-specific. Any GMT works: MSigDB collections,
BioCarta, PID, or a hand-curated gene set of your own.

---

## Reusing this pattern for other databases

The shape here is general, and is worth reusing wherever IGVFagent needs a
whole reference resource rather than one lookup at a time — **FAVOR** being
the obvious next candidate, since variant annotation has exactly the same
"query one record at a time over HTTP" bottleneck that pathway lookups had.

The reusable parts:

1. **A cache with provenance.** Atomic download, freshness window, release id
   + SHA-256 per file, manifest that merges rather than overwrites.
2. **Normalise into a namespace you do not control.** `GeneIndex` is the
   pathway instance of this; a variant resource needs the same for rsID /
   SPDI / HGVS / CHR:POS:REF:ALT, and the equivalent lesson applies — resolve
   through an authority, not through one source's own identifier list.
3. **Merge with the source's own structure first, heuristics second.**
   Reading Reactome's hierarchy removed every known-wrong merge that string
   similarity produced. Prefer curated fact wherever a database publishes it.
4. **Record every supporting source per fact,** so corroboration is queryable.
5. **Make the merge auditable** (`merges.csv`) and the criterion **measurable**
   (`evaluate`) rather than asserted.
6. **Stamp what you write into the KG,** so a rebuild is idempotent and
   renames do not accumulate orphans.

The two failure modes worth carrying over: a documented column that is empty
in practice (Reactome's Entrez field), and an upstream generator that mangles
its own text (WikiPathways' apostrophes). Both were found by measuring output
counts, not by reading documentation.

---

## Related

- `igvfagent intpath` — the 2012 IntPath release, for comparison
- `igvfagent pathway-viz network` — network figures; reads these KG edges
- `igvfagent enrich pathways` — Enrichr-proxy over-representation testing
