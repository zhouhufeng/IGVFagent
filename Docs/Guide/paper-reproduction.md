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
