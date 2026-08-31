# Hosted IGVFagent — `igvfagent.genohub.org`

A shared deployment so researchers can use IGVFagent without installing it
locally or paying for their own API calls. Runs on
[Hcloud](https://github.com/zhouhufeng/HCloud).

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
                                    /mnt/igvf-data (data volume)
```

Nothing listens on a public port. `cloudflared` dials **out** to Cloudflare,
so the host firewall only opens SSH and there is no web port to scan.

Host-specific values — addresses, credentials, volume and instance details —
are intentionally **not** in this repo. Operators: see the private runbook.

---

## Requirements

- A Linux host (Ubuntu 24.04 LTS assumed) with Docker Engine + Compose plugin
  and `cloudflared` installed, and a data volume mounted at `/mnt/igvf-data`.
- A Cloudflare zone for the public hostname.
- An Anthropic API key that the deployment will spend against.

The LLM runs API-side, so the host mostly orchestrates HTTP calls and renders
artefacts. 8 vCPU / 30 GB is comfortable for that plus the occasional scanpy
job; the binding constraint is disk for downloaded datasets, not CPU.

---

## First-time setup

```bash
# 1. Data volume (idempotent — formats only an unformatted disk)
sudo bash /srv/igvfagent/Deploy/bootstrap-vm.sh

# 2. Secrets
cd /srv/igvfagent/Deploy
cp .env.prod.example .env.prod && chmod 600 .env.prod
$EDITOR .env.prod          # ANTHROPIC_API_KEY, CF_TUNNEL_TOKEN

# 3. Shared access password
sudo apt-get install -y apache2-utils
htpasswd -Bc nginx/htpasswd igvf      # prompts for the password
# nginx workers run as uid 101, NOT as the login user (uid 1000). A 0600 file
# owned by that user gives every authenticated request a 500 and
#   [crit] open() "/etc/nginx/htpasswd" failed (13: Permission denied)
# in the gateway log — while *unauthenticated* requests still 401 correctly,
# which makes it look like auth works. Hand the file to the nginx uid:
sudo chown 101:101 nginx/htpasswd && sudo chmod 400 nginx/htpasswd

# 4. Up  (--env-file is required: the compose file interpolates from it)
docker compose --env-file .env.prod -f docker-compose.prod.yml up -d
docker compose --env-file .env.prod -f docker-compose.prod.yml ps
```

### Cloudflare Tunnel

Zero Trust → Networks → Tunnels → **Create a tunnel** (type `cloudflared`),
then add a public hostname mapping the site to the `gateway:80` service. Copy
the tunnel token into `CF_TUNNEL_TOKEN`. Creating the public hostname writes
the DNS record automatically.

> An account-scoped Cloudflare API token returns **"Invalid API Token"** from
> `GET /user/tokens/verify` even when it is valid. Test it against a real
> endpoint (`/zones?name=<zone>`) instead.

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
cd /srv/igvfagent && bash Deploy/redeploy.sh
```

`redeploy.sh` pulls, rebuilds with `--no-cache`, recreates the container,
and then **verifies** that the Python actually installed inside the
container matches this checkout. It exits non-zero if it does not.

**Why a script rather than three commands.** The app container does not
mount the source tree — compose mounts only `/mnt/igvf-data/{Data,Docs}`,
and `Scripts/streamlit_app.py` imports `from igvfagent import ...`, i.e.
the package installed into `/opt/venv` by `pip install '.[all]'` at image
**build** time. So `git pull` updates the checkout while the container
keeps serving whatever was baked into the image, and the site looks
unchanged however many times you pull. A bare `docker build` can no-op
too: run it from the wrong directory and it builds the wrong context; let
a cached layer stand or the image id not move and the container is never
recreated. Each of those failures is silent, which is why the script
checks rather than assumes.

To check without changing anything:

```bash
cd /srv/igvfagent && bash Deploy/redeploy.sh --check
```

The manual equivalent, if you would rather run it by hand:

```bash
cd /srv/igvfagent && git pull
cd Deploy
docker compose --env-file .env.prod -f docker-compose.prod.yml build --no-cache app
docker compose --env-file .env.prod -f docker-compose.prod.yml up -d --force-recreate app
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

- the **backend** is fixed server-side, so nobody can point the app at a local
  model that isn't running there and get a confusing connection error;
- the **model menu** is an operator allowlist (`IGVF_PUBLIC_MODELS`), so a
  visitor chooses only from models you are willing to pay for;
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
- Infrastructure credentials are **never** synced to the host.

**Known limitation — the workspace is shared, not per-user**

`Scripts/_localstore.py` resolves a single process-global `IGVF_PROJECT_ROOT`
at import time, so every visitor shares one `Data/KG/local_kg.sqlite`, one
DuckDB warehouse, and one `Docs/<skill>/` tree. On this deployment that means:

- users can see each other's run outputs, and
- the local knowledge graph accumulates everyone's queries together.

That is acceptable behind a shared password among colleagues; it is **not**
acceptable if the gate is ever opened to the anonymous public. Per-user
isolation needs the module-level root refactored into a per-session value —
tracked as follow-up work, not shipped here.

**Also unaddressed:** the agent executes wrapped CLIs as subprocesses with
model-chosen arguments (`Scripts/_tools.py:3679`) on prompts from whoever is
logged in. The container hardening above bounds the blast radius; it does not
eliminate the class.
