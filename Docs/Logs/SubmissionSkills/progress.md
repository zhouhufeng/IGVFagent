# IGVF Portal data submission in IGVF Agent: progress

2026-10-17 · zhouhufeng

## Summary

Submitters can now check, audit, create and upload IGVF Portal objects with one tool, from the command line or the web app. Writes rehearse on staging first and send only on an explicit `--apply` or the user's own Submit click.

- **CLI:** `igvfagent submit` has 23 commands covering every stage, from key checks and object ordering to creating, uploading, patching and auditing.
- **The Portal's own audits** now run locally, on Portal objects and on unsubmitted drafts. They reproduce every finding a client can check (327 of 327) and report none the Portal does not.
- **Web app (igvfagent.genohub.org):** users paste their own Portal keys, upload data to a private folder, let the agent dry-run the submission and click **Submit to staging**. Prod submission is switched off.
- **Agent tools:** nine read-only `igvf_submit_*` tools, so the agent can audit and prepare but never write.

Status on 2026-11-20: the submission skills in this repository match igvfagent.genohub.org: the `submit` command with `sheet-*`, `portal-audit`, `send` and the working directory, and the agent's `igvf_submit_*` tools. The hosted site's web-app pieces (pasting keys in the chat, the Submit button) are not in this repository's web app; from your own install you send a plan with `igvfagent submit send`. See [Using it from your own install](#using-it-from-your-own-install-2026-11-20).

## Timeline

| Date | What landed |
|---|---|
| 2026-09-10 | IGVF-Submit absorbed as the `submit` and `portal-qc` skills: `doctor`, `audit` (24 DACC-meeting rules), `inputs`, `crosscheck`, `validate`, `template`, `enum`, `patch` |
| 2026-09-11 | `plan` (dependency order), `check-revoked`, `status` |
| 2026-09-23 | IGVF data model from the lab submission diagrams; lineage context |
| 2026-10-15 | `portal-qc`: provenance by reverse `derived_from` lookup (6 min instead of a full scan), lab slug match, CSV output |
| 2026-10-16 | `submit create`: POSTs drafts in dependency order, validates against the live schema, reuses existing objects, uploads file contents, per-host Portal keys |
| 2026-10-16 | Web submission: keys pasted in chat are removed before the model sees them and stored encrypted per user; "My data" folder; Submit-to-staging button |
| 2026-10-16 | Fix: audit file checks were off after `create` reused the `FILE_TYPES` name; `patch` waits until the Portal shows the write |
| 2026-10-16 | IGVF-DACC igvf-metadata-submitter (Google Sheets) ported as `sheet-template`, `sheet-get`, `sheet-post`, `sheet-patch`, `sheet-put`, `sheet-upload` for TSV/CSV/XLSX |
| 2026-10-16 | Submitter overview, `Docs/Logs/SubmissionSkills/IGVF-Agent-submission-skills.md` |
| 2026-10-17 | igvfd's own audit checkers ported as `submit portal-audit`, with staging and prod fixtures and an agreement report |
| 2026-10-29 | Changing objects already on the Portal from the agent: `sheet-patch` dry runs save a plan; reading staging objects (TST accessions) with your staging key; boto3 for uploads |
| 2026-11-19 | Working directory for local installs (`igvfagent ask/ui --workdir DIR`, or `IGVF_PROJECT_ROOT`): the submission commands find your files there and, inside a sheet, next to the sheet. On a shared deployment a submission may name only the signed-in person's own files |
| 2026-11-20 | `submit send <plan>` sends a saved dry-run plan as checked; every dry run prints it. Each plan records who made it, so the hosted Submit button offers your latest unsent plan, also one made by a deep-analysis job |
| 2026-11-20 | This public repository updated with all of the above |

## What exists now

### Command line: `igvfagent submit`

| Stage | Commands |
|---|---|
| Set up | `doctor` (keys, host, account, labs it may submit for) |
| Prepare | `plan`, `template`, `enum`, `sheet-template`, `sheet-get` |
| Check before upload | `validate` (gzip, md5, ragged rows, BED coordinates), `portal-audit --drafts` |
| Submit | `create` (drafts and `iu_register`-style sheets, `--upload`), `sheet-post`, `sheet-upload` |
| Fix | `patch` (incl. `--repoint`), `sheet-patch`, `sheet-put` (admin) |
| Check after submission | `audit` (24 DACC rules), `portal-audit` (the Portal's own checkers, `--compare` against `@@audit-self`), `inputs`, `check-revoked`, `crosscheck`, `status`, `get`, `fileset-files` |
| Portal-wide QC | `igvfagent portal-qc summary / audits / provenance / uploads / filesets / lab` |

Every write is a dry run unless `--apply` is given, and prod also needs `--yes`. The CLI needs Python 3.9 and no packages.

### The Portal's own audits (`portal-audit`)

The igvfd audit checkers (`src/igvfd/audit`) are ported unchanged into `Scripts/igvfd_audit_rules.py`. `Scripts/igvfd_audits.py` gives them what the server would: `request.embed` becomes a GET on the Portal API, or a lookup of a draft by alias.

Agreement with the Portal's computed audits (`Benchmarks/submission_staging/igvfd_audit_agreement.md`):

| | Objects | Portal findings | Matched | Missed | Extra | Recall | Precision |
|---|---|---|---|---|---|---|---|
| All | 186 | 345 | 327 | 18 | 0 | 94.8% | 100% |
| prod | 127 | 273 | 265 | 8 | 0 | 97.1% | 100% |
| staging | 59 | 72 | 62 | 10 | 0 | 86.1% | 100% |

All 18 misses need data a client cannot have:
- 14 incorrect-file-bucket findings read the server's S3 storage.
- 3 ontology findings need the Portal's ontology database.
- 1 status finding concerns a linked object these keys cannot read; it is reported as `UNVERIFIABLE`, not as a pass.

Tests: `Benchmarks/test_igvfd_audits.py`, with fixtures in `Benchmarks/submission_staging/igvfd_fixtures/`.

### Web app

1. The user pastes their Portal Access Key ID and Secret. `Scripts/portal_keystore.py` removes them from the message before the model or any transcript sees them, checks which host they work on, and stores them encrypted under `IGVF_KEYSTORE_SECRET`.
2. The user uploads data and sheets to **My data**, a private per-user folder of up to 100 MB per file.
3. The agent checks, audits and dry-runs `create` with the user's own keys.
4. The user clicks **Submit to staging**. Only that click writes. Prod needs `IGVF_ALLOW_PROD_SUBMIT=1` on the server and is off.

### Agent tools

`igvf_submit_audit`, `igvf_submit_inputs`, `igvf_submit_crosscheck`, `igvf_submit_validate`, `igvf_submit_fileset_files`, `igvf_submit_template`, `igvf_submit_sheet_template`, `igvf_submit_sheet_get`, `igvf_submit_create` (dry run only). None of them can write to the Portal.

### Verified on staging

Lab `xihong-lin`, award `HG012064`:
- Objects posted by `create`: `TSTDS63011942`, `TSTFI69604689`, and one document. They match what the official `iu_register` posted (`TSTDS45289259`, `TSTFI32661857`) on every submitted field.

## Open items

1. **`portal-audit` is not an agent tool yet.** The CLI has it, but the agent still audits with the 24-rule `igvf_submit_audit`. The next step is an `igvf_submit_portal_audit` tool (read-only), and making the web flow run it on drafts before the Submit button is offered.
2. **File format specification document:** writing it for the user is planned as `spec-doc`.
3. **Prod submission from the web app:** off until staging use by real submitters has been reviewed.
4. **Large raw data through the web app:** uploads are capped at 100 MB per file. Larger files go through `create --upload` or `sheet-upload` from the command line.
5. **Server-only audits:** file bucket and ontology checks cannot run client-side and are reported as skipped. `--compare` shows the Portal's own result for submitted objects.

## Where things live

| What | Where |
|---|---|
| Submission CLI | `Scripts/igvf_submission_skill.py` (`igvfagent submit`) |
| Sheet commands | `Scripts/igvf_metadata_sheet.py` |
| Portal audits | `Scripts/igvfd_audits.py`, `Scripts/igvfd_audit_rules.py` |
| Portal-wide QC | `Scripts/igvf_portal_qc_skill.py` (`igvfagent portal-qc`) |
| Key store | `Scripts/portal_keystore.py` |
| Agent tools | `Scripts/_tools.py` (`igvf_submit_*`) |
| Submitter overview | `Docs/Logs/SubmissionSkills/IGVF-Agent-submission-skills.md` |
| Submitter's guide | `Docs/Skills/igvf-submission.md` |
| Tests | `Benchmarks/test_igvf_submission.py`, `Benchmarks/test_igvfd_audits.py` |

## Using it from your own install (2026-11-20)

```bash
git pull                                   # or clone github.com/zhouhufeng/IGVFagent
pip install -e '.[analysis,ui,llm,submit]' # submit = boto3, for uploading file contents
export IGVF_STAGING_API_KEY=...            # your own staging key pair (Portal profile -> Access keys)
export IGVF_STAGING_SECRET_KEY=...
igvfagent submit -m staging doctor         # keys, account, the labs you may submit for

# 1. Dry run: checks the drafts against the live schema and your files, orders them, saves the exact plan
igvfagent submit create --objects curated_set.json tabular_file.json -m staging
#    -> Submission plan: Data/Submitter/<time>_plan_staging/drafts.json  (2 object(s), staging)
# 2. Send that plan, exactly as checked (prod also needs --yes)
igvfagent submit send Data/Submitter/<time>_plan_staging/drafts.json
# 3. Check the new objects
igvfagent submit -m staging audit <new IGVFDS accession> --recurse   # the DACC-meeting checklist
igvfagent submit -m staging portal-audit <new IGVFDS accession>      # the Portal's own audits
```

With the agent: `igvfagent ask --workdir /path/to/your/data "prepare a staging submission of chr1.tsv.gz as a
tabular file in a curated set"` (or the browser UI started with `igvfagent ui --workdir /path/to/your/data`). The
agent only prepares the dry run; it then gives you the `igvfagent submit send ...` command, which you run yourself
with your keys. To change objects already on the Portal, the agent prepares a `sheet-patch` plan the same way.

Lab and award come from each draft, or from `IGVF_LAB` / `IGVF_AWARD` when a draft has none. Every object needs an
alias with your lab's prefix (`<lab-name>:<id>`); the Portal rejects other prefixes.

