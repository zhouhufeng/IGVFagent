# Reproducing a paper

When a paper publishes its analysis code, reproducing the paper means running
that code. IGVFagent's `paper-code` skill takes the repository named in the
paper's Code Availability statement and runs the authors' analysis
unmodified, in an environment built for it. It checks the output against the
authors' own rendering and reports every figure, section by section.

The design follows [Paper2Agent](https://github.com/jmiao24/Paper2Agent)'s
route for code-backed papers:
- the authors' repository is the scientific source of truth, and nothing in it is reimplemented;
- the source is executed before anything else, and its results are compared with the authors' rendered outputs;
- the environment is isolated and pinned, and only the environment is ever repaired;
- retries are bounded, and blockers are recorded with the exact chunk and error.

**On this page**

- [What it does](#what-it-does)
- [From the browser or a job](#from-the-browser-or-a-job)
- [From the command line](#from-the-command-line)
- [Reading the result](#reading-the-result)
- [Repairs are environment-only](#repairs-are-environment-only)
- [Records, the Reproductions page and the forum](#records-the-reproductions-page-and-the-forum)
- [Reusing a reproduced analysis on new data](#reusing-a-reproduced-analysis-on-new-data)
- [Worked example: Matreyek et al. 2018](#worked-example-matreyek-et-al-2018)

## What it does

| Step | What happens |
|---|---|
| find | The repository comes from the paper's harvest (`bench harvest`). The one named in the Code Availability sentence outranks tools the paper merely used. |
| fetch | It is cloned into `Data/PaperCode/<owner>__<repo>/src` and the commit is recorded. |
| inventory | This finds the main analysis (`.Rmd`, `.qmd`, `.ipynb`, `.R` or `.py`), and the inputs it reads and whether the repository has them. It also finds the packages it loads, the directories it writes to, the authors' own rendered output, the package versions the authors ran (from their comments and rendering) and the analysis sections. |
| env | R uses a micromamba environment with `r-base`, rmarkdown, pandoc, poppler and conda-forge `r-<pkg>`, falling back to CRAN or Bioconductor. Python uses a uv virtual environment. Environments are cached by specification. |
| run | The analysis runs in a work copy, unmodified. An R Markdown file is knit with `keep_md` and figures saved as PNG and PDF. Errors are kept per chunk, so one failing chunk does not hide the rest. A notebook runs with `nbconvert --execute --allow-errors`. |
| compare | Every printed value block is compared with the authors' rendering, as identical, numerically close (1e-6) or the same values in another layout, and so is the figure count. Library warnings are reported separately. |
| replay | This is optional (`--replay`). The analysis runs a second time, and printed values and saved files must come out identical. PDFs are compared without their timestamps. |
| report | `report.md` and `report.html` have one section per top-level heading, holding its figures, errors and printed-value agreement. `summary.json` feeds plan checks. |

**Repositories with many analyses.** Many papers publish one notebook per figure. Pass `--all-entries`, or `all_entries` from a job, to run every notebook and R Markdown file:
- They run in path order in one shared work copy, as the authors ran them, so a later notebook finds what an earlier one wrote.
- Each analysis gets its own sub-report (`entries/NN_<name>/`), with its own agreement, figures and first root-cause error.
- The paper gets one combined summary and report.
- Large files, such as deposited data, are hard-linked into the work copy rather than copied.

**The authors' pipeline, then their notebooks.** Many notebooks read intermediate results that the authors' workflow produced. With `--workflow` (`workflow` from a job), a Snakemake `Snakefile` runs first in the work copy:
- A dry run lists the planned jobs and **every missing input**; those are usually deposited data.
- Then the real run goes with `--keep-going`.
- The report gives the failed rules and the files produced.

After the first pass, **every missing file an analysis names drives more passes**:
- With `--workflow`, the workflow is asked for exactly those files. It dry-runs each and builds every one it can, for example model runs made from the deposited processed objects.
- Any analysis whose missing inputs now exist is re-run, whichever analysis or workflow produced them. This repeats for up to 2 extra passes, which covers notebooks that read another notebook's output.
- The report lists each pass and what is **still missing**. That is usually supplementary tables, or paths on the authors' own machine, which are reported as blockers.

`igvfagent paper-code data <repo> --zenodo <record|DOI> --into <path> [--files GLOB] [--extract]` puts a deposit into the checkout, where later runs find it. It also takes `--url` for a direct link. Checksums are verified and provenance is recorded in `Data/PaperCode/<repo>/data.json`. A job does this itself (`paper_code_fetch_data`) when analyses fail on missing inputs. Missing data that was never deposited is reported as a blocker.

**Compatibility shims for Python notebooks.** Like the R shims, these are loaded into the kernel without editing the notebook, and they are listed in the report whenever they fire. `mpl_missing_style` falls back to the default style when a notebook uses the authors' unpublished matplotlib style, such as `jr`: figures look different, values are unchanged.

**The authors' declared environment comes first.** Python notebooks run on a micromamba environment built from the repository's own `environment.yml`, whenever it pins Python:
- GPU-only packages (`cudatoolkit`, `cudnn` and similar) are dropped on a CPU host.
- If the exact pins don't solve, they are relaxed to major.minor, then to names only, with Python kept at its major.minor.
- `pip:` dependencies install inside the environment. Pinned build tools such as Cython and numpy go first, and source packages are then retried without build isolation. What pip can't build comes from conda-forge/bioconda, and only after that is a pin retried unpinned.
- Every change is recorded in the report.
- Notebooks are driven by IGVF Agent's own nbconvert and nbclient, while their code runs in a kernel launched from the authors' environment (ipykernel, or IRkernel for R). That way the authors' pinned Jupyter packages cannot break the runner, and the kernel name from their machine (such as `jy_anbe_py38`) is not needed. The notebook files are not modified. An analysis that never executed is always counted as a root-cause error, never as clean.

When the authors ran R older than 3.6, the run uses `RNGkind(sample.kind = "Rounding")`. R 3.6 changed `sample()`, and without this the authors' `set.seed()` would give different random draws.

## From the browser or a job

Ask for the reproduction, for example *"Reproduce Matreyek 2018 PTEN VAMP-seq, Nat Genet"*. That becomes a [long-running job](jobs.md), and its protocol makes the code route the reproduction:

1. Resolve and harvest the paper, then `paper_code_find` on the harvest.
2. Run `paper_code_reproduce` (with `replay`) in the background, and wait on `<run_dir>/done.json`.
3. Read `summary.json`. For root-cause errors, re-run with the suggested environment pins, at most 3 repairs.
4. Verify section by section with parallel fresh sub-agents (`delegate_tasks`). Each one says which of the paper's claims its sections' outputs support, contradict or leave untested.
5. Optionally, cross-check against IGVF or MaveDB data. That is a separate stage, never a substitute.

A deterministic harness check backs this up. A reproduction job whose harvested Code Availability names a repository cannot pass verification until a stage has run that code, or has been marked blocked with a reason.

## From the command line

```bash
igvfagent paper-code find --harvest Docs/Benchmark/<run>/harvest.json
igvfagent paper-code inventory FowlerLab/VAMPseq
igvfagent paper-code pipeline FowlerLab/VAMPseq --replay --paper "Matreyek et al. 2018"
igvfagent paper-code pipeline --harvest Docs/Benchmark/<run>/harvest.json --detach
igvfagent paper-code status Docs/PaperCode/<run>
igvfagent paper-code report Docs/PaperCode/<run>        # rewrite the report of a finished run
```

Useful options:
- `--entry` picks the analysis file.
- `--ref` pins a commit.
- `--pin` sets environment pins, repeatable.
- `--strict` stops at the first failing chunk.
- `--timeout-min` limits the run time.

## Reading the result

`Docs/PaperCode/<timestamp>_<owner>_<repo>/`:

| File | Contents |
|---|---|
| `report.md`, `report.html` | per-section report with every figure, root-cause errors with suggested repairs, and the authors' package versions against this run's |
| `summary.json` | `printed.fraction_matched`, `figures.produced` and `figures.reference`, `chunks.root_errors` and `chunks.cascade_errors`, `repair`, `replay.deterministic` |
| `compare.json` | the full printed-value comparison, with the blocks that differ |
| `inventory.json`, `env.json`, `run.json` | what was found, the environment built, the command and new files |
| `work/` | the executed work copy: the knit `.md` and `.html`, figures, and files the analysis saved |
| `figures/` | PNG renderings of the PDFs the analysis saved |

## Repairs are environment-only

The authors' code is never edited. Errors are split into two kinds:
- **root causes:** the chunk, its Rmd line and the message;
- **knock-on errors:** "object not found" in code that needed an object a failing chunk would have created.

For known upstream breakages the report suggests the pin that fixes them. A pin is either `r-<pkg>=<version>` from conda-forge, or `cran:<pkg>@<version>`, built from the CRAN archive with the environment's compilers. A failed pin marks the environment incomplete, so it is never reused.

## Records, the Reproductions page and the forum

Every reproduction or benchmark job that finishes writes one **record per paper**, whatever the outcome. The files are in `Data/Reproductions/<paper_id>/` on the persistent volume:

| File | Contents |
|---|---|
| `record.json` | The best attempt, the full list of attempts, the owner, and the published and forum state. |
| `attempts/<job>.json` | One attempt. It has the paper (DOI, harvested claims, named code repositories), the route (the authors' code or public data), and agreement with the authors' printed values. It also has chunks, figures, replay, shims, benchmark checks, the harness-verified stages, the verifier's verdict, the final report and a selected log. |
| `<job>.html` | One self-contained HTML report, like Paper2Agent's delivered page. It has the headline, the tables, the plan, the final report and the authors'-code report section by section, with **every figure embedded**. |

The outcome comes from harness facts only:
- **reproduced:** the job is done, the verifier passed, and either the authors' code matched at least 90% of their printed values with no root-cause errors, or every benchmark check passed.
- **partial:** anything in between. If the paper publishes its code and an attempt never ran it, that attempt is at best partial.
- **not reproduced:** neither holds.

Attempts on the same paper join one record, matched by DOI or repository.

On the site, go to **📊 Validation → 📑 Reproductions**:
- Search by title, DOI, repository, gene, assay or outcome.
- See every paper with its result, verifier verdict and number of attempts.
- Open any attempt to see its full report, and download the HTML.

Records are private to their owner. **Publish** makes one visible to every signed-in user. A job's panel links to its record, and `/?repro=<paper_id>` opens a record directly.

**Post to discussion.genohub.org** puts a summary in the IGVF Agent › IGVF Agent Help category. The summary includes the headline numbers, the agreement table, benchmark checks, blocked stages, key figures (one per section), the selected run log and a link back to the record.
- The post is authored as the owner's forum account.
- Nothing is sent until the owner has read the preview and confirmed it.
- A later attempt on the same paper is posted as a reply in the same topic.

From the command line:

```bash
igvfagent repro list --query vampseq
igvfagent repro show matreyek2018_pten_vampseq
igvfagent repro record <job_id>                  # refresh one job's record
igvfagent repro record-run Docs/PaperCode/<run> --doi 10.1038/...   # a standalone paper-code run
igvfagent repro backfill                         # every finished job (run by Deploy/redeploy.sh)
igvfagent repro publish <paper_id> [--off]
igvfagent repro post <paper_id>                  # preview; add --yes to post
```

## Reusing a reproduced analysis on new data

When a paper's authors' code has been reproduced, its record becomes a **tool**: `paper_<paper_id>`. It re-runs that code, at the same commit and with the same environment repairs and compatibility shims, using your files in place of the inputs it reads. It reports the result section by section, with every figure. The method stays the authors'; nothing is reimplemented. Paper2Agent's Paper2MCP delivers the same thing, reproduced analyses as tools.

- The agent sees these tools like any other. For example: *"run paper_matreyek2018_pten_vampseq with my DMS table in place of PTEN_variant_data.tsv"*.
- From the command line:
  ```bash
  igvfagent repro tools                                  # every reusable analysis and its inputs
  igvfagent repro apply matreyek2018_pten_vampseq --pten-variant-data-tsv my_scores.tsv [--replay] [--detach]
  ```
- Each input takes a file in the same format as the authors' input; omitted inputs keep theirs.
- The report states the substitution. Its comparison with the authors' rendering then measures how the results changed.

Tested on Matreyek 2018 with every PTEN score multiplied by 0.9: all 67 chunks of the authors' code ran and produced 119 figures, and exactly the printed values that depend on the scores changed (7 blocks).

## Worked example: Matreyek et al. 2018

*Multiplex assessment of protein variant abundance by massively parallel
sequencing* (Nat Genet 2018). The Code Availability statement names
[FowlerLab/VAMPseq](https://github.com/FowlerLab/VAMPseq). That repository has
one R Markdown analysis of 67 chunks in 8 sections, the four supplementary
tables, and the authors' HTML rendering (119 printed-value blocks, 59 images),
all at commit `a960aa9b2034`.

Each run changed only the environment:

| Run | Environment | Result |
|---|---|---|
| 1 | latest conda-forge (R 4.5.3, data.table 1.18.6, reshape2 1.4.5) | Knit ok in 103 s, 115 figures, 104 of 119 printed blocks identical. 2 chunks fail at `melt()` on a data.frame: 3 root errors and 18 knock-on. |
| 2 | adds `r-data.table=1.17.8` | Same failure; data.table 1.17 already raises the error. |
| 3 | adds `cran:data.table@1.16.4 cran:reshape2@1.4.4` | `melt()` now redirects, but gives 3 columns where the authors' next line selects 4. 1 chunk still fails. |
| 4 | adds the `melt_reshape` shim and `RNGkind(sample.kind="Rounding")` (authors' R 3.4.4) | All 67 chunks run; 108 of 110 value blocks match. The 2 permutation-test tables differ. |
| 5 | adds the `jitter_global_rng` shim (authors' ggplot2 2.2.1) | **110 of 110** value blocks match, including both permutation-test tables |
| 6 | latest packages, **no pins**, shims chosen automatically | same as run 5 in one shot: 67 of 67 chunks, 119 figures, 110 of 110 value blocks, deterministic replay |

Final result (run 6, local and on the hosted site):
- **Printed values:** 106 identical, plus 4 with the same values in another layout (correlation matrices wrapped at a different console width). That is all 110 of the authors' value blocks.
- **Library warnings:** the authors' 9 ggplot2 warnings are reported separately. ggplot2 4 words them differently.
- **Replay:** deterministic, with 126 of 126 printed blocks and 65 of 65 saved files identical. The analysis writes its own working directory into its PyMOL scripts, and the comparison normalises that path.

Why the shims were needed. **`melt_reshape`:** In 2018 the analysis attached reshape, then data.table. `melt()` went to data.table, which redirected data.frames to `reshape2::melt`, whose `UseMethod` then found reshape's exported `melt.data.frame` on the search path. R 3.6 stopped searching the search path for S3 methods. The shim is a knitr hook that puts `reshape::melt` first once data.table is attached, which is what `melt()` resolved to in 2018. It is applied by detection: the analysis loads reshape before data.table and calls `melt()`. The Rmd text is unchanged. Paper2Agent's environment manager reached the same diagnosis on this repository.

**`jitter_global_rng`:** ggplot2 before 3.0 jittered points with `base::jitter()` in the global random stream, taking 2n draws per jittered layer. ggplot2 3.0 and later takes one draw for a private seed instead. Every seeded `sample()` after the first jittered plot therefore shifted, including the ClinVar and cancer permutation tests. The shim makes `geom_jitter()` and `position_jitter()` use `seed = NULL`, which runs the same `base::jitter()` in the global stream. It is applied when the authors' ggplot2 is older than 3.0 and the analysis jitters and sets a seed.

---

[← Documentation index](README.md) · [Project README](../../README.md)
