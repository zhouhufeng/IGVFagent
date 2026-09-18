# IGVF Agent — `submit` Skill: Build Specification
 
Build a new first-class skill for IGVFagent (github.com/zhouhufeng/IGVFagent) that
automates IGVF Data Portal submission workflows: validation, dependency-ordered
submission, audit triage, patching, and revoked-input checking. Follow the
existing skill conventions exactly so the skill is dual-drivable (CLI + LLM tool).
 
---
 
## 1. Repo conventions the skill MUST follow
 
- One module: `Scripts/submitter_skill.py`, subcommand-style
  (`python3 Scripts/cli.py submit <subcommand> ...`). Register `"submit"` in the
  `cli.py` dispatch dict.
- Register each subcommand as a typed tool in `Scripts/_tools.py`
  (`Tool(name=..., cli=["submit", "<subcommand>", ...], parameters=<JSON schema>)`).
- ALL endpoints via `Scripts/_endpoints.py` `resolve(name, env_var)` — hex-encoded
  defaults, per-endpoint env override. Add:
  - `portal` (exists): https://data.igvf.org  (override `IGVF_PORTAL_BASE`)
  - `portal_api` (exists): https://api.data.igvf.org
  - `portal_sandbox` (NEW): https://sandbox.igvf.org  (override `IGVF_PORTAL_SANDBOX_BASE`)
- Credentials ONLY from environment; never written to source or logs. Document
  new vars in `.env.example`.
- Outputs under `Data/Submitter/<YYYYMMDD_HHMMSS>_<label>/`; log every run to the
  warehouse `runs` provenance table (skill, subcommand, args_json, started/ended,
  rows_emitted, success) like other skills.
- Prefer wrapping the official clients over raw REST where sensible:
  - igvf_utils (submission): data.igvf.org/help/data-submission/submission-by-igvf-utils/
  - igvf-python-client: github.com/IGVF-DACC/igvf-python-client
  - Portal MCP (reference only): github.com/IGVF-DACC/igvf-portal-mcp
  - API spec: data.igvf.org/help/igvf-api-spec/
  - Schemas/profiles: data.igvf.org/profiles/  (fetch live; do not hardcode schemas)
  - NOTE: verify igvf_utils' expected credential env var names from its docs
    (likely IGVF_API_KEY / IGVF_SECRET_KEY) and adopt the same names.
---
 
## 2. Submission domain model (distilled from DACC meetings — encode as rules)
 
### 2.1 Object graph and REQUIRED submission order
1. `software`  →  2. `software_version` (MUST point to a **tagged GitHub release**;
   link via alias, e.g. `lab-name:tool`)  →  3. `workflow`  →  4. `analysis_step`
   →  5. `analysis_step_version`  →  6. file set (`prediction_set` / `curated_set`
   / `analysis_set`)  →  7. `tabular_file` (the data)  →  8. `document`
   (file format specification)  →  9. curated-set external-source tabular file.
- Files link to workflow via the file's `analysis_step_version` property.
- Every submission needs `lab` and `award` (e.g. /labs/<lab>/, /awards/<award>/).
- `submitter_comment` is free text for caveats (e.g., order of software in an ASV).
### 2.2 Prediction sets
- `file_set_type`: `functional effect`.
- ONE prediction set per cell type; each needs `samples`
  (HepG2: IGVFSM0002XTWQ, K562: IGVFSM4381OIBT — examples, not defaults).
- Needs 1–2 sentence `description`.
- Needs `input_file_sets`: every file set whose files were inputs.
### 2.3 Tabular files
- MUST be gzipped.
- BED outputs: 0-based, `file_format: bed`, `file_format_type: bed3`
  (coordinates only, no scores — model: IGVFFI8092FZKL).
- MUST set `reference_files` for genomic data or audits fire — generic GRCh38
  reference: **IGVFFI8743BGYR**.
- MUST set `derived_from`: ALL input files (accessions/UUIDs). If some inputs
  were filtered/subset, the filtering must be described in the file format
  specification document.
- MUST link a file-format-specification `document`
  (`document_type: "file format specification"`) via the file's
  `file_format_specifications` property — **NOT** via `documents`. The spec
  states 0- vs 1-based, end-inclusive vs -exclusive, and each column's content.
- `content_type` enums seen: `element references` (coordinate lists),
  `external source data` (URL tables under curated sets).
### 2.4 Curated sets / external data
- External resources are NOT uploaded; submit a tabular file under a curated set
  with columns: Name | URL | short description (`content_type: external source data`).
### 2.5 Sequencing data (for completeness)
- Sequencing files require a seqspec document
  (data.igvf.org/help/data-submission/seqspec-submission/).
- Only final VCFs needed for eQTL-type submissions; no intermediates.
### 2.6 Common audits → mechanical fixes (implement as audit triage rules)
| Audit | Fix |
|---|---|
| missing `description` | PATCH description (1–2 sentences) |
| missing `derived_from` | collect input file accessions → PATCH `derived_from: [...]` |
| missing `input_file_sets` | collect input file-set accessions → PATCH on the file set |
| missing genomic reference | PATCH `reference_files: ["IGVFFI8743BGYR"]` |
| unexpected input file sets | reconcile `derived_from` vs `input_file_sets`; document filtering in the format spec |
| missing files (in a set) | upload; ensure gzip |
| revoked upstream input | see revocation workflow (§4, `check-revoked`) |
 
### 2.7 Auth model
- Anonymous portal GET may 403 (e.g., unreleased objects); submission always
  requires credentials. Sandbox (sandbox.igvf.org) is the practice target;
  production is data.igvf.org.
---
 
## 3. Environment variables (add to `.env.example`, never commit values)
 
```
IGVF_PORTAL_BASE=            # prod override
IGVF_PORTAL_SANDBOX_BASE=    # sandbox override
IGVF_API_KEY=                # submission key   (confirm names vs igvf_utils docs)
IGVF_SECRET_KEY=             # submission secret
IGVF_SUBMIT_LAB=             # default lab   e.g. /labs/xihong-lin/
IGVF_SUBMIT_AWARD=           # default award e.g. /awards/HGxxxxxx/
IGVF_SUBMIT_TARGET=sandbox   # sandbox | prod   (DEFAULT: sandbox)
```
 
---
 
## 4. CLI subcommands (each also a registered LLM tool)
 
- `auth-check` — verify credentials against the selected target; print lab/award
  resolution; never print secrets.
- `plan --objects <yaml|sheet>` — given draft metadata, emit the dependency-ordered
  submission plan (§2.1) with per-object required-property checklist; no network
  writes.
- `validate --file <path> --profile <type>` — offline validation: gzip check,
  BED 0-based/bed3 check, column check against the format-spec draft, required
  properties vs live `/profiles/` schema, `reference_files` presence, alias format.
- `spec-doc --file <path> [--based 0|1] [--end exclusive|inclusive]` — scaffold a
  file-format-specification document (markdown + JSON metadata) from the file
  header; user completes column descriptions.
- `submit --objects <yaml|sheet> [--target sandbox|prod] [--dry-run]` — POST in
  dependency order via igvf_utils/python-client. `--dry-run` prints exact payloads.
  Production requires BOTH `--target prod` AND interactive confirmation.
- `audit --id <accession> [--fix-plan]` — fetch audits for an object (and its
  files), map each to the §2.6 fix table, emit a ready-to-apply patch plan.
- `patch --id <accession> --set derived_from=[...],description=...  [--dry-run]` —
  apply targeted PATCHes; always fetch-and-show current values first.
- `check-revoked --id <accession>` — walk `derived_from` upstream; flag any
  revoked/archived/superseded inputs; report the replacement accession from the
  revocation notice; if the revocation ships an invalid-record list (e.g. a
  `.jsonl` of bad alleles), download it and intersect with the local/derived
  file; report which branch applies:
  (a) no overlap → patch `derived_from` to corrected accession only;
  (b) overlap → list affected rows (e.g. swapped ref/alt) for correction,
      then reupload + patch + release.
- `status [--lab <lab>]` — list the lab's submitted objects with release state
  and open audit counts.
All subcommands print artefact paths with the standard announcement prefixes
(`Report:`, `Manifest:`) so the agent runtime captures them.
 
---
 
## 5. Safety rails (non-negotiable)
 
- Default target is SANDBOX; production writes require `--target prod` + confirm.
- Every write subcommand supports `--dry-run`; `plan`/`validate` never write.
- Never log or echo credentials; redact `Authorization` in transcripts.
- All POST/PATCH payloads are saved under the run's output dir before sending
  (auditable, replayable — consistent with the agent's provenance model).
- PATCH never blind-overwrites: read-modify-write with the current value shown.
---
 
## 6. Worked scenario to implement as an integration test (real case)
 
Context: tabular file IGVFFI0856UJHP under prediction set IGVFDS2581EDPS used a
Sherwood input later revoked ("ref allele mismatch to GRCh38, corrected in
IGVFFI1678CDBR"), with 10 invalid alleles listed in `invalid_IGVFFI7160EKDK.jsonl`.
 
`igvfagent submit check-revoked --id IGVFFI0856UJHP` must:
1. detect the revoked input in `derived_from`;
2. surface the replacement accession (IGVFFI1678CDBR);
3. fetch the invalid-allele list and intersect with the file's variants;
4. print branch (a) patch-only or (b) fix-alleles-then-reupload, with the exact
   patch payload for (a).
Use sandbox fixtures / mocked responses for CI; never hit prod in tests.
 
---
 
## 7. Deliverables checklist
 
- [ ] `Scripts/submitter_skill.py` with the subcommands above
- [ ] `cli.py` dispatch entry `"submit"` + `_tools.py` tool registrations
- [ ] `_endpoints.py`: `portal_sandbox` (hex-encoded) + env overrides
- [ ] `.env.example` additions (§3)
- [ ] `Docs/Skills/SUBMITTER_SKILL.md` playbook (usage + the §2 rulebook)
- [ ] Unit tests: validators (gzip/bed3/0-based/spec-linkage/reference_files);
      audit→fix mapping; dependency-order planner
- [ ] Integration test: §6 revocation workflow against mocked portal
- [ ] Dry-run demo transcript in `Docs/Skills/` showing plan → validate →
      submit --dry-run → audit --fix-plan on sandbox
