---
name: igvf-submission
description: Submit, audit and fix data on the IGVF data portal (data.igvf.org / sandbox.igvf.org). Use when submitting prediction sets, curated sets, tabular files, BED files or file format specification documents to IGVF; when a DACC audit flags missing derived_from, description, reference_files, file_format_specifications or analysis_step_version; when an input file has been revoked or superseded; or when preparing for the monthly IGVF submission meeting.
---

# IGVF data portal submission

Submission is a slow feedback loop: you submit, the DACC runs audits, and the
next monthly meeting tells you what was missing. The same handful of things are
missing every time. Run the checks yourself first.

## Rehearse on staging, not sandbox

`api.sandbox.igvf.org` is **deprecated** — it answers HTTP 410 with "Please
update your client to use api.staging.igvf.org". Older submission docs and
meeting notes still say sandbox; use **staging** instead. It is a real portal
you can break: do first-time submissions there, then repeat against prod.

Every command takes `-m staging` (`-m sandbox` warns and redirects).

## Credentials

```bash
export IGVF_API_KEY=...        # Access Key ID from the portal profile page
export IGVF_SECRET_KEY=...     # Access Key Secret
export IGVF_MODE=prod          # or staging
export IGVF_LAB=/labs/<lab>/   # optional, used as a default by igvf_utils
export IGVF_AWARD=/awards/<award>/
```

Each host issues its own keys — a prod key returns HTTP 401 on staging, and
vice versa. Get a staging key from the staging profile page.

Your own unreleased objects return HTTP 403 without keys; released objects are
public. Never commit keys — keep them in `~/.igvf/credentials.env` with mode
600 and `source` it.

Verify with `igvfagent submit doctor`. It reports which lab the key may submit for.

## The tool

Absorbed into IGVFagent as `igvfagent submit` (`Scripts/igvf_submission_skill.py`) and `igvfagent portal-qc` (`Scripts/igvf_portal_qc_skill.py`). Both are dependency-free (stdlib Python 3.9+), so they also run standalone.

```bash
igvfagent submit doctor                                   # creds + connectivity
igvfagent submit audit <accession> --recurse              # the recurring-findings checklist
igvfagent submit inputs <accession>                        # flag revoked/archived inputs
igvfagent submit fileset-files <fileset>                   # paste-ready JSON accession array
igvfagent submit validate <path> --bed3                    # local checks before upload
igvfagent submit crosscheck --old A --new B --mine C        # did an input revision hit my file?
igvfagent submit patch <acc> --repoint derived_from=OLD:NEW # safe array element swap
igvfagent submit template <profile> --include ...           # TSV header for iu_register
igvfagent submit enum <profile> <property> --grep <term>    # allowed values
```

`audit` exits non-zero when something will block release, so it works in CI.

## Order of submission

Objects must exist before anything links to them. Bottom-up:

1. **Software** and **software version** — needs a GitHub repo with a *tagged*
   release, e.g. `.../releases/tag/1.0.0`. An untagged repo blocks this.
2. **Workflow** / **analysis step** / **analysis step version** — the DACC
   usually creates these; ask for them and record the returned accession.
3. **Document** with `document_type: file format specification`.
4. **File set** — `prediction_set`, `curated_set`, ... with `description` and
   `input_file_sets` already populated.
5. **Files** — `tabular_file` etc., linking `file_set`, `derived_from`,
   `reference_files`, `file_format_specifications`, `analysis_step_version`.
6. **Upload** the bytes, wait for `upload_status: validated`.
7. **Release** — ask the DACC. Release fails while any input is revoked.

## The checks that actually fail

Every one of these was an open to-do across multiple 2026 meetings.

| What | Where | Note |
|---|---|---|
| `description` | file set | 1-2 real sentences, not a phrase |
| `input_file_sets` | prediction set | every file set you consumed |
| `derived_from` | each output file | every input **file**; filtering must be described in the format spec |
| `reference_files` | tabular file | **omitting it is a hard audit error** |
| `file_format_specifications` | each output file | link here, **never** via `documents` |
| `analysis_step_version` | each output file | link once the workflow exists |
| gzip | every file | `gzip -n`; raised at every single meeting |
| `file_format_type` | BED files | `bed3` for a plain coordinate list |
| 0-based half-open | BED files | `end > start` always |
| `samples` | prediction set | e.g. `IGVFSM0002XTWQ` HepG2, `IGVFSM4381OIBT` K562 |

## Revoked and superseded inputs

An input can be revoked *after* you derive from it, which silently blocks
release of everything downstream. `igvfagent submit inputs <acc>` finds these and reads
`superseded_by` to name the replacement.

When an input is revoked for a data error, there are two very different
outcomes, and guessing wrong wastes a month:

- The revised variants are **not** in your file → only `derived_from` needs
  repointing, then release. No reupload.
- They **are** in your file → fix, re-upload, *and* repoint.

`igvfagent submit crosscheck --old <revoked> --new <corrected> --mine <my file>` decides
this. It compares alleles, not just positions, because a reference-allele
mismatch is corrected in place — the coordinate stays the same while the alleles
change, so position-only comparison misses the dangerous case entirely.

It runs a **control** first: if unrevised input variants do not overlap your
file at all, the coordinate conventions differ and it refuses to give a verdict
rather than reporting a false clean. Trust a clean result only when the control
passed.

## Patching arrays

A REST `PATCH` **replaces** an array. Patching `derived_from` with one
accession deletes the other twenty-four. Either use
`igvfagent submit patch --repoint prop=OLD:NEW`, which rebuilds the array from the live
value, or use `iu_register.py --patch`, which extends arrays by default
(`-w/--overwrite-array-values` to replace instead).

`igvfagent submit patch` is a dry run unless `--apply`, and refuses to write to prod
without `--yes`.

## The official tools

- `iu_register.py -m staging -p <profile> -i <file> -d` — POST/PATCH from TSV,
  JSON or JSONL. `-d` dry run. Add `--patch` and a `record_id` column to modify
  existing objects. Arrays are comma-delimited; `#`-prefixed columns are ignored.
  From [`igvf_utils`](https://github.com/IGVF-DACC/igvf_utils).
- [`igvf-python-client`](https://github.com/IGVF-DACC/igvf-python-client) —
  generated API client.
- [`igvf-portal-mcp`](https://github.com/IGVF-DACC/igvf-portal-mcp) — MCP server,
  **read-only**; useful for exploring, cannot submit.
- Schemas are authoritative and live: `https://api.data.igvf.org/profiles/<profile>`.
  Check an enum there rather than guessing a value.

## Getting help

`igvf-portal-help@lists.stanford.edu`. The DACC creates workflow and analysis
step objects on request, and performs releases. Submission meetings are the
third Friday of the month, 10am PT.
