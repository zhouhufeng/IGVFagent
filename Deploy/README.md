# Hosted IGVFagent — `igvfagent.genohub.org`

A shared deployment so researchers can use IGVFagent without installing it
locally or paying for their own API calls. Runs on Jetstream2 allocation
`BIO260320_IU` (Indiana University).

```
internet ──▶ Cloudflare (TLS, DNS) ──▶ cloudflared tunnel (outbound)
                                              │
                                              ▼
                                    nginx gateway (basic auth)
                                              │
                                              ▼
                                    Streamlit app (public mode)
                                              │
                                              ▼
                                    /mnt/igvf-data (2 TB Cinder)
```

Nothing listens on a public port. `cloudflared` dials **out** to Cloudflare,
so the VM's security group only opens SSH and there is no web port to scan.

---

## Provisioned resources

| Resource | Value |
|---|---|
| Instance | `igvfagent-prod`, flavor `m3.medium` (8 vCPU / 30 GB) |
| Image | `Featured-Ubuntu24` (Ubuntu 24.04 LTS) |
| Private IP | assigned on `auto_allocated_network` |
| Floating IP | see `Docs/Secret/DEPLOYMENT.md` — **SSH only** |
| Volume | `igvfagent-data`, 2000 GB, XFS at `/mnt/igvf-data` |
| Keypair | `igvfagent-deploy` (private key in `Docs/Secret/`, gitignored) |
| Security group | `igvfagent` — inbound TCP/22 only |

The LLM runs API-side, so the VM mostly orchestrates HTTP calls and renders
artefacts; 30 GB is comfortable for that plus the occasional scanpy job.

Allocation headroom after this deployment: **185 of 225 cores**, ~722 GB of
877 GB RAM, **1000 GB of 15000 GB volume quota, and 4 of 10 volume slots**.
Volume quota is the binding constraint, not cores — FAVOR already holds
12 TB across 5 volumes (`favor-clickhouse` 4.6 TB, `favor-minio` 3.4 TB,
`favor-expansion` 2 TB, `favor-rocksdb` 1.8 TB, `favor-data-1` 150 GB).

> **Jetstream2 bills core-hours.** An always-on `m3.medium` consumes roughly
> **5.8k SU/month**. Resize up for a heavy analysis and back down after:
>
> ```bash
> openstack server resize --flavor m3.large --wait igvfagent-prod
> openstack server resize confirm igvfagent-prod
> ```
>
> The data volume is mounted by UUID with `nofail`, so it survives a resize.

---

## First-time setup

```bash
ssh -i Docs/Secret/igvfagent-deploy.pem ubuntu@<floating-ip>

# 1. Volume (idempotent; already done on the current host)
sudo bash /srv/igvfagent/Deploy/bootstrap-vm.sh

# 2. Secrets
cd /srv/igvfagent/Deploy
cp .env.prod.example .env.prod && chmod 600 .env.prod
$EDITOR .env.prod          # ANTHROPIC_API_KEY, CF_TUNNEL_TOKEN

# 3. Shared lab password
sudo apt-get install -y apache2-utils
htpasswd -Bc nginx/htpasswd igvf      # prompts for the password
# nginx workers run as uid 101, NOT as ubuntu (uid 1000). A 0600 file owned
# by ubuntu gives every authenticated request a 500 and
#   [crit] open() "/etc/nginx/htpasswd" failed (13: Permission denied)
# in the gateway log — while *unauthenticated* requests still 401 correctly,
# which makes it look like auth works. Hand the file to the nginx uid:
sudo chown 101:101 nginx/htpasswd && sudo chmod 400 nginx/htpasswd

# 4. Up  (--env-file is required: the compose file interpolates from it)
docker compose --env-file .env.prod -f docker-compose.prod.yml up -d
docker compose --env-file .env.prod -f docker-compose.prod.yml ps
```

### Cloudflare Tunnel — already provisioned

Created via the Cloudflare API; recorded here so it can be rebuilt.

| Item | Value |
|---|---|
| Tunnel name | `igvfagent` |
| Tunnel ID | see `Docs/Secret/DEPLOYMENT.md` |
| Config source | `cloudflare` (remotely managed — ingress lives in the API, not on disk) |
| Ingress | `igvfagent.genohub.org` → `http://gateway:80`, fallback `http_status:404` |
| DNS | `CNAME igvfagent.genohub.org → <tunnel-id>.cfargotunnel.com`, proxied |
| Token | `Docs/Secret/cf-tunnel-token.txt` (gitignored) → `CF_TUNNEL_TOKEN` |

Before this, `igvfagent.genohub.org` resolved through the zone's **wildcard**
record to a `searchvity.com` domain-parking page. The explicit CNAME added
here is more specific, so it wins; leave it in place.

To rebuild from scratch, or to re-read the token:

```bash
# token used below is the account-scoped Cloudflare API token
ACCT=<account-id>; TUN=<tunnel-id>
curl -s -H "Authorization: Bearer $CF_API_TOKEN" \
  "https://api.cloudflare.com/client/v4/accounts/$ACCT/cfd_tunnel/$TUN/token"
```

> That API token is account-scoped, so `GET /user/tokens/verify` returns
> **"Invalid API Token"** even though the token is fine. Test it against a
> real endpoint (`/zones?name=genohub.org`) instead.

---

## Operating

Every `docker compose` call needs `--env-file .env.prod`; the compose file
declares `ANTHROPIC_API_KEY` and `CF_TUNNEL_TOKEN` as required (`${VAR:?}`),
so without it even `logs` aborts on a missing-variable error. Shorthand:

```bash
cd /srv/igvfagent/Deploy
alias dc='docker compose --env-file .env.prod -f docker-compose.prod.yml'

dc logs -f app        # agent activity
dc logs -f tunnel     # ingress
dc restart app
dc ps
curl -s localhost/healthz    # unauthenticated probe
```

Update to a new code revision:

```bash
cd /srv/igvfagent && git pull
docker build -t igvfagent:latest .
cd Deploy && docker compose --env-file .env.prod -f docker-compose.prod.yml up -d
```

Rotate the shared password (note the re-`chown` — `htpasswd` rewrites the file
and resets its ownership, which silently breaks the gateway otherwise):

```bash
sudo chown $USER nginx/htpasswd && sudo chmod 600 nginx/htpasswd
htpasswd -B nginx/htpasswd igvf
sudo chown 101:101 nginx/htpasswd && sudo chmod 400 nginx/htpasswd
docker compose --env-file .env.prod -f docker-compose.prod.yml restart gateway
```

### Verifying a deployment

```bash
curl -o /dev/null -w '%{http_code}\n' https://igvfagent.genohub.org/          # 401
curl -o /dev/null -w '%{http_code}\n' -u igvf:wrong https://igvfagent.genohub.org/   # 401
curl -o /dev/null -w '%{http_code}\n' -u igvf:PASS  https://igvfagent.genohub.org/   # 200

# The one that actually matters — Streamlit is a blank page without it.
# Force HTTP/1.1: over HTTP/2 the Upgrade header is meaningless and you
# will get a misleading 200 instead of 101.
curl -sI --http1.1 -u igvf:PASS \
  -H 'Connection: Upgrade' -H 'Upgrade: websocket' \
  -H 'Sec-WebSocket-Version: 13' -H 'Sec-WebSocket-Key: dGhlIHNhbXBsZSBub25jZQ==' \
  https://igvfagent.genohub.org/_stcore/stream | head -1   # 101 Switching Protocols
```

---

## Cost control

Public mode (`IGVF_PUBLIC_MODE=1`) is what keeps a shared key from being an
open wallet:

- the model is **pinned server-side** — no picker in the browser, so nobody
  can select a more expensive model or a local backend that isn't running;
- **max iterations** is capped (`IGVF_PUBLIC_MAX_ITER`). Each iteration is one
  LLM call, so this is the dominant per-visitor cost lever;
- **max tokens/turn** is capped (`IGVF_PUBLIC_MAX_TOKENS`);
- temperature is fixed at 0 for reproducibility.

Visitors may lower these but never raise them. Set a **billing alert on the
Anthropic account** as the real backstop — the caps bound a single run, not
the number of runs.

---

## Security posture

**What is protected**

- No inbound web port; the tunnel is outbound-only.
- Shared-password gate at the nginx layer, ahead of the app.
- The app container drops all capabilities, gets `no-new-privileges`, a PID
  limit, and CPU/memory ceilings.
- Artefact rendering is contained by `Scripts/_pathguard.py`. The UI scrapes
  file paths out of model and tool output and renders them; without the guard
  any absolute path ending in a viewable extension became a download button.
  The guard restricts rendering to the workspace **and** subtracts secret
  material inside it (`Docs/Secret/`, `.env`, `*.pem`, credential-like names).
- Jetstream credentials (`clouds.yaml`) are **never** synced to the VM.

**Known limitation — the workspace is shared, not per-user**

`Scripts/_localstore.py` resolves a single process-global `IGVF_PROJECT_ROOT`
at import time, so every visitor shares one `Data/KG/local_kg.sqlite`, one
DuckDB warehouse, and one `Docs/<skill>/` tree. On this deployment that means:

- users can see each other's run outputs, and
- the local knowledge graph accumulates everyone's queries together.

That is acceptable behind a shared lab password among colleagues; it is **not**
acceptable if the gate is ever opened to the anonymous public. Per-user
isolation needs the module-level root refactored into a per-session value —
tracked as follow-up work, not shipped here.

**Also unaddressed:** the agent executes wrapped CLIs as subprocesses with
model-chosen arguments (`Scripts/_tools.py:3679`) on prompts from whoever is
logged in. The container hardening above bounds the blast radius; it does not
eliminate the class.
