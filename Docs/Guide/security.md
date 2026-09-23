# Security

How credentials and data are handled. The per-mode data flow is in [`Docs/THREAT_MODEL.md`](../THREAT_MODEL.md).

- The repository ships with an aggressive `.gitignore` that excludes runtime
  outputs, caches, logs, and any `.env` / cookie / token files.
- All credentials are read from environment variables; nothing is hardcoded.
- Logs record request URLs and HTTP status codes only — never credential
  headers.
- Never commit `.env`, browser cookie exports, OAuth tokens, API keys, or
  unreleased / pre-publication data.

---

[← Documentation index](README.md) · [Project README](../../README.md)
