# IGVF Agent: Portal data submission skills and CLI

2026-10-16 · zhouhufeng · updated with `create`, the `sheet-*` commands and web submission

## What it is

IGVF Agent's submission skill catches the recurring DACC audit findings in seconds, before you submit, instead of a month later at a meeting. It is a command-line tool (`igvfagent submit`, also runnable alone as `igvf-sub`) plus six agent tools that the IGVF Agent web app and MCP server call.

- It checks objects already on the Portal (prod or staging) and local files before upload, against 24 rules distilled from DACC meetings.
- It reads freely and writes only on request: every write command (`create`, `patch`, `sheet-post`, `sheet-patch`, `sheet-put`, `sheet-upload`) is a dry run unless you add `--apply`, and prod also needs `--yes`.
- It creates objects and uploads files itself (`create`, `sheet-post`, `sheet-upload`), following the official clients' conventions: objects it posted on staging match `iu_register`'s on every submitted field.
- It includes the DACC's Google Sheets metadata-submitter, ported to the command line for TSV, CSV and XLSX sheets (`sheet-*`).
- It needs only Python 3.9 and no packages, so it runs on a laptop, a cluster login node or in a container.

## Before you start

You need your own Portal access keys to see your lab's unreleased objects; released objects are public and need none.

```bash
export IGVF_API_KEY=...      # Access Key ID, from your Portal profile
export IGVF_SECRET_KEY=...   # Access Key Secret
export IGVF_MODE=staging     # staging to rehearse, prod for the real thing
export IGVF_LAB=/labs/<your-lab>/     # optional, as in igvf_utils
export IGVF_AWARD=/awards/<your-award>/
igvfagent submit doctor      # checks the keys and the connection
```

- Staging issues its own keys: a prod key gets HTTP 401 there. Put staging keys in `IGVF_STAGING_API_KEY` / `IGVF_STAGING_SECRET_KEY` (or the generic names with `-m staging`).
- Rehearse on staging (`api.staging.igvf.org`). The old sandbox host now answers HTTP 410, so `--mode sandbox` is mapped to staging automatically.
- Keys come only from the environment. The tool never writes them to files or logs.
- Every object needs a `lab` and an `award`; `doctor` shows which account the keys belong to and the labs it may submit for, and warns if it can read but submit for none.

## CLI reference

Every command reads the Portal and changes nothing, except `patch --apply`. Run them as `igvfagent submit <command>`, or `python3 Scripts/igvf_submission_skill.py <command>` with no install.

| Command | What it does | Example |
| --- | --- | --- |
| `doctor` | Checks keys, mode and connection; shows the account and its labs | `submit doctor` |
| `plan` | Puts a draft list of objects in the order they must be submitted; sends nothing | `submit plan --objects drafts.yaml` |
| `status` | Lists a lab's file sets with release state and open audit counts | `submit status --lab "Jane Doe, Stanford"` |
| `get` | Fetches one object as JSON, optionally only some fields | `submit get IGVFDS6332VCTO --field status files` |
| `audit` | Runs the DACC-findings checklist on file sets or files | `submit audit IGVFDS1234ABCD --recurse --strict` |
| `inputs` | Lists every linked input and flags revoked or archived ones | `submit inputs IGVFFI1234ABCD` |
| `check-revoked` | Walks `derived_from` and flags revoked or archived inputs | `submit check-revoked IGVFFI1234ABCD` |
| `crosscheck` | Says whether corrections in a revised input change your derived file, i.e. whether you must re-upload | `submit crosscheck --old IGVFFI_OLD --new IGVFFI_NEW --mine mine.tsv.gz` |
| `fileset-files` | Member file accessions of a file set, as a paste-ready JSON array | `submit fileset-files IGVFDS1234ABCD --filter-type TabularFile` |
| `validate` | Checks a local file before upload: gzip, md5 and content_md5, ragged rows, BED coordinates | `submit validate peaks.bed.gz --bed3` |
| `template` | The required columns of a Portal profile as a TSV header for `iu_register` | `submit template tabular_file` |
| `enum` | A property's allowed values | `submit enum tabular_file content_type --grep element` |
| `create` | Creates objects from YAML/JSON drafts and/or `iu_register`-style TSV/CSV sheets, in link order; validates files and the live schema; reuses an object that already exists; `--upload` sends file contents | `submit -m staging create --objects drafts.yaml sheet.tsv --profile tabular_file --apply --upload` |
| `sheet-template` | A blank metadata-submitter sheet for a profile, plus a legend of each column's role, type and allowed values | `submit sheet-template tabular_file` |
| `sheet-get` | Fills a sheet's rows from the Portal by accession, uuid or alias | `submit sheet-get tabular_file.tsv` |
| `sheet-post` | POSTs each row and writes the new accession and uuid back into a new sheet | `submit sheet-post curated_set.tsv --apply` |
| `sheet-patch` | PATCHes rows: `--columns a,b`, `--all`, or `--remove x` (GET, drop keys, PUT) | `submit sheet-patch tabular_file.tsv --columns description --apply` |
| `sheet-put` | Replaces whole objects (admin only, `--admin`) | `submit sheet-put lab.tsv --admin --apply` |
| `sheet-upload` | Uploads each row's `#upload_abspath`, or prints the `aws` commands (`--print-cmd`) | `submit sheet-upload tabular_file.tsv --apply` |
| `patch` | Changes one object; a dry run that prints the request and its `iu_register` equivalent unless `--apply` (`--yes` too for prod); `--repoint` swaps one accession inside an array | `submit patch IGVFFI1234ABCD --repoint derived_from=IGVFFI_OLD:IGVFFI_NEW` |

After submission, the companion `igvfagent portal-qc` command checks what is already on the Portal: the Portal's own audit counts (`summary`, `audits`), inputs that were replaced without the derived files being updated (`provenance`), uploads that never validated (`uploads`), file sets missing release-blocking properties (`filesets`) and a per-lab scorecard (`lab`).

## Using it through the agent

In the web app (igvfagent.genohub.org) or the MCP server, ask in plain words; the agent picks one of nine tools. None of them can write to the Portal.

| Tool | Ask it | Runs |
| --- | --- | --- |
| `igvf_submit_audit` | "Will IGVFDS1234ABCD pass the DACC audits?" | `audit` |
| `igvf_submit_inputs` | "Are any inputs of IGVFFI1234ABCD revoked or archived?" | `inputs` |
| `igvf_submit_crosscheck` | "This input was revoked and corrected; do I have to re-upload my file?" | `crosscheck` |
| `igvf_submit_validate` | "Check this BED file before I upload it." (upload it first) | `validate` |
| `igvf_submit_fileset_files` | "Give me the file accessions in IGVFDS1234ABCD as a JSON array." | `fileset-files` |
| `igvf_submit_template` | "What columns does a tabular_file need for iu_register?" | `template` |
| `igvf_submit_sheet_template` | "Give me a sheet to fill in for a curated set." | `sheet-template` |
| `igvf_submit_sheet_get` | "Fill this sheet from what is on the Portal." | `sheet-get` |
| `igvf_submit_create` | "Check my metadata and prepare a staging submission." (dry run only) | `create` |

**Submitting from the web app (staging):**

1. Paste your Portal Access Key ID and Secret into the chat. The app removes them from the message before the model or any transcript sees them, checks which host they work on, and stores them encrypted under your account (the sidebar shows them, with Delete).
2. Upload your data files and metadata sheets to **My data** in the sidebar (your own folder, kept between visits, up to 100 MB per file).
3. Ask the agent to check, audit and prepare the submission. It works with your own keys and can only do a dry run.
4. Review the plan and click **Submit to staging**. Only that click sends it. Submitting to prod is switched off for now.

Once you have given the web app your keys, its submission tools use them, never the deployment's, so it sees exactly what you may see and submit.

## Recommended workflow

Submit each object only after the objects it points at exist; posting out of order is the most common reference error, and `plan` checks it for you.

**Order of objects** (as `plan` enforces it):

1. `software`: the tool itself.
2. `software_version`: must point at a tagged GitHub release, linked by alias (`lab-name:tool`).
3. `workflow`: ties the software versions together.
4. `analysis_step`: a step within the workflow.
5. `analysis_step_version`: what each file links to.
6. The file set: `prediction_set`, `curated_set` or `analysis_set`.
7. `tabular_file`: the data, gzipped, with `derived_from` and `reference_files`.
8. `document`: the file format specification, linked from the file's `file_format_specifications`.
9. For external resources: a tabular file under a curated set (Name, URL, description) instead of uploading the data.

**Where each command fits:**

1. `doctor` on staging, then `plan --objects drafts.yaml` to check the order.
2. `template <profile>` and `enum <profile> <property>` to fill the `iu_register` sheets correctly.
3. `validate <file>` on every local file before upload.
4. Register and upload with the official `iu_register` (igvf_utils).
5. `audit <accession> --recurse` on staging, then fix with `patch` (dry run first) until it is clean.
6. Repeat on prod. Afterwards, `status` and `portal-qc` track release state and the DACC's own audits; `inputs` and `crosscheck` tell you when an upstream file is revoked and whether you must re-upload.

## The checks it encodes

`audit` and `validate` encode 24 rules, each a finding raised at DACC meetings; the test suite checks every one.

| Check | On | Fails or warns when |
| --- | --- | --- |
| `description` | file set | missing, or a phrase rather than 1-2 real sentences |
| `input_file_sets` | prediction set | a file set whose files were inputs is not listed |
| `derived_from` | each output file | an input file is missing; filtering must then be described in the format spec |
| `reference_files` | tabular file | missing (a hard audit error), or pointing at the archived generic GRCh38 file `IGVFFI8743BGYR` |
| `file_format_specifications` | each output file | the spec is attached through `documents` instead |
| `analysis_step_version` | each output file | not linked once the workflow exists |
| gzip | every file | not gzipped (use `gzip -n`) |
| BED type | BED files | no `file_format_type` (`bed3` for a plain coordinate list) |
| Dead inputs | any object | an input is revoked or archived |
| Status | any object | still `in progress` or `submitted` (a warning) |
| Rows | local file | empty, or rows of uneven width |
| BED coordinates | local BED file | too few columns, negative positions, or `end` not greater than `start` (BED is 0-based, half-open) |
| md5 | local file | never fails: prints `md5sum` and `content_md5sum` for the submission sheet |

## Safety and limits

Nothing reaches the Portal unless you run `patch --apply` yourself; everything else only reads.

- `patch` prints the exact request and its `iu_register` equivalent first. It sends only with `--apply`, and prod also needs `--yes`.
- The agent's tools cannot write at all, in the web app or over MCP.
- Keys are read from the environment only and are never logged.

**Not automated yet:**

- Writing the file format specification document for you (planned as `spec-doc`).
- Submitting to prod from the web app (staging only for now).
- Raw data over 100 MB through the web app: use `create --upload` or `sheet-upload` from the command line.

## Where things live and getting help

| What | Where |
| --- | --- |
| Submission CLI | `Scripts/igvf_submission_skill.py` (`igvfagent submit`) |
| Portal-wide QC | `Scripts/igvf_portal_qc_skill.py` (`igvfagent portal-qc`) |
| Submitter's guide | `Docs/Skills/igvf-submission.md` |
| Build spec, from DACC meeting notes | `Docs/BugZilla/Submit-Skill.md` |
| Rule tests (24 rules) | `Benchmarks/test_igvf_submission.py` |
| Agent tool definitions | `Scripts/_tools.py` (`igvf_submit_*`) |

Official DACC tools: [igvf_utils](https://github.com/IGVF-DACC/igvf_utils) (`iu_register` for posting and patching), [igvf-python-client](https://github.com/IGVF-DACC/igvf-python-client), and the read-only [igvf-portal-mcp](https://github.com/IGVF-DACC/igvf-portal-mcp). Profiles at `https://api.data.igvf.org/profiles/<profile>` are the authority on allowed values.

Help: `igvf-portal-help@lists.stanford.edu`. The DACC creates workflow and analysis-step objects on request and performs releases; submission meetings are the third Friday of each month, 10am PT.

