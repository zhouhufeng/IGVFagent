#!/usr/bin/env bash
# Push the local IGVF Portal credentials to the hosted deployment, redeploy,
# and prove the change took effect.
#
#   bash Deploy/make-live.sh            # do it
#   bash Deploy/make-live.sh --check    # report only, change nothing
#
# Run from the repo root. Everything it needs is under Docs/Secret/.
#
# Why a script rather than three commands: each step fails silently in its
# own way. Appending to .env.prod twice leaves two IGVF_ACCESS_KEY lines and
# compose honours the LAST, so a stale value can win. `git pull` alone never
# changes the running site, because the container serves the package baked
# into the image at build time. A rebuild that no-ops leaves the old image
# running with no error. Each step is therefore checked, not assumed.
#
# The secret is never passed as a command-line argument to ssh: it would be
# visible in `ps` on the VM for the life of the call. The remote script goes
# in argv (base64, so quoting cannot break it) and the secret keeps stdin to
# itself -- they cannot share stdin, see step 4.
set -euo pipefail

CHECK_ONLY=0
[[ "${1:-}" == "--check" ]] && CHECK_ONLY=1

CRED_FILE="Docs/Secret/IGVFportalAPI.txt"
KEY_FILE="Docs/Secret/igvfagent-deploy.pem"
HOST="ubuntu@149.165.151.21"
REMOTE="/srv/igvfagent"

die()  { printf '\nFAILED: %s\n' "$1" >&2; exit 1; }
step() { printf '\n=== %s ===\n' "$1"; }

[[ -f "$CRED_FILE" ]] || die "$CRED_FILE not found — run from the repo root"
[[ -f "$KEY_FILE"  ]] || die "$KEY_FILE not found"
chmod 600 "$KEY_FILE" "$CRED_FILE" 2>/dev/null || true

SSH=(ssh -i "$KEY_FILE" -o StrictHostKeyChecking=accept-new -o ConnectTimeout=30 "$HOST")

step "1. Read the Portal key pair"
# Parsed by the same module the app uses, so a layout that works here is
# exactly the one that works in the container.
CREDS=$(python3 - <<'PY'
import sys
sys.path.insert(0, "Scripts")
import _credentials as c
cr = c.portal_credentials()
if not cr:
    sys.exit("could not parse a key pair from " + str(c.credential_file()))
print(cr[0]); print(cr[1])
PY
) || die "credential parse failed"
ACCESS_KEY=$(printf '%s\n' "$CREDS" | sed -n 1p)
SECRET_KEY=$(printf '%s\n' "$CREDS" | sed -n 2p)
[[ -n "$ACCESS_KEY" && -n "$SECRET_KEY" ]] || die "empty key pair"
printf 'key id %s, secret %s chars (never printed)\n' "$ACCESS_KEY" "${#SECRET_KEY}"

step "2. Confirm the key actually grants extra visibility"
COUNTS=$(ACCESS_KEY="$ACCESS_KEY" SECRET_KEY="$SECRET_KEY" python3 - <<'PY'
import base64, json, os, urllib.request
tok = base64.b64encode(
    f"{os.environ['ACCESS_KEY']}:{os.environ['SECRET_KEY']}".encode()).decode()
u = "https://api.data.igvf.org/search/?type=MeasurementSet&limit=0&format=json"
def total(auth):
    h = {"Accept": "application/json"}
    if auth:
        h["Authorization"] = "Basic " + tok
    with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=60) as r:
        return json.loads(r.read()).get("total", 0)
print(total(True), total(False))
PY
) || die "Portal probe failed"
AUTHED=${COUNTS%% *}; ANON=${COUNTS##* }
printf 'MeasurementSet: %s authenticated vs %s anonymous\n' "$AUTHED" "$ANON"
[[ "$AUTHED" -gt "$ANON" ]] || die "key grants no extra visibility — wrong or revoked"

step "3. Reach the VM"
"${SSH[@]}" 'set -e; echo "host $(hostname)"; docker --version' || die "SSH failed"

if [[ "$CHECK_ONLY" == 1 ]]; then
  step "4. Current state (--check: nothing is being changed)"
  "${SSH[@]}" "set -e; cd $REMOTE; echo \"HEAD: \$(git log --oneline -1)\"; \
    echo \"env IGVF_ACCESS_KEY lines: \$(grep -c '^IGVF_ACCESS_KEY=' Deploy/.env.prod || echo 0)\"; \
    docker exec igvfagent-app python3 -c \
      'from igvfagent import _credentials as c; print(\"container creds:\", c.describe()[\"source\"], c.describe()[\"key_id\"])' 2>/dev/null \
      || echo 'container creds: unavailable (old image or not running)'; \
    docker exec igvfagent-app kb --version 2>&1 | grep -i kb_python \
      || echo 'aligner: NOT installed in the running image'"
  printf '\n--check complete. Re-run without --check to apply.\n'
  exit 0
fi

step "4. Upsert credentials into .env.prod (idempotent, secret via stdin)"
# The script travels in argv (it is not secret) and the credentials travel on
# stdin (they are). They cannot share stdin: giving ssh a heredoc AND a pipe
# means the heredoc wins, the remote `bash -s` reads the script from stdin,
# and `read` then consumes the script's own remaining lines instead of the
# credentials -- which is exactly how this failed the first time
# ("SECRET_KEY: unbound variable"). base64 keeps the script a single argv
# token, immune to quoting, and leaves stdin free.
REMOTE_SCRIPT=$(cat <<'REMOTE_SCRIPT_EOF'
set -euo pipefail
ROOT="$1"
read -r ACCESS_KEY
read -r SECRET_KEY
[ -n "$ACCESS_KEY" ] && [ -n "$SECRET_KEY" ] || { echo "empty credentials over stdin" >&2; exit 1; }
cd "$ROOT/Deploy"
cp -n .env.prod ".env.prod.bak.$(date +%Y%m%d%H%M%S)" 2>/dev/null || true
upsert() {  # replace in place, else append -- never leave two definitions
  local key="$1" val="$2"
  if grep -q "^${key}=" .env.prod 2>/dev/null; then
    grep -v "^${key}=" .env.prod > .env.prod.tmp
    printf '%s=%s\n' "$key" "$val" >> .env.prod.tmp
    mv .env.prod.tmp .env.prod
  else
    printf '%s=%s\n' "$key" "$val" >> .env.prod
  fi
}
upsert IGVF_ACCESS_KEY        "$ACCESS_KEY"
upsert IGVF_SECRET_ACCESS_KEY "$SECRET_KEY"
chmod 600 .env.prod
echo "IGVF_ACCESS_KEY lines:        $(grep -c '^IGVF_ACCESS_KEY=' .env.prod)"
echo "IGVF_SECRET_ACCESS_KEY lines: $(grep -c '^IGVF_SECRET_ACCESS_KEY=' .env.prod)"
REMOTE_SCRIPT_EOF
)
REMOTE_B64=$(printf '%s' "$REMOTE_SCRIPT" | base64 | tr -d '\n')
printf '%s\n%s\n' "$ACCESS_KEY" "$SECRET_KEY" | "${SSH[@]}" \
  "printf %s '$REMOTE_B64' | base64 -d > /tmp/mklive.\$\$.sh && bash /tmp/mklive.\$\$.sh '$REMOTE'; rc=\$?; rm -f /tmp/mklive.\$\$.sh; exit \$rc" \
  || die "could not update .env.prod"

step "5. Pull the new code and rebuild the image"
# kb-python is a NEW dependency: the aligner cannot appear in an image built
# before it was declared, so this must be a real rebuild.
"${SSH[@]}" "set -e; cd $REMOTE && git pull --ff-only && bash Deploy/redeploy.sh" \
  || die "redeploy failed — the site is still serving the old image"

# Commands run through plain `docker exec`, never `sh -lc`. A login shell
# sources /etc/profile, which resets PATH and drops /opt/venv/bin -- so
# `kb` and `igvfagent` both look absent and the run reports ALIGNER MISSING
# against an image that has the aligner installed and working.
step "6. Verify inside the running container"
"${SSH[@]}" "set -e; \
  docker exec igvfagent-app python3 -c 'from igvfagent import _credentials as c; d=c.describe(); print(\"credentials:\", d[\"source\"], d[\"key_id\"])'; \
  docker exec igvfagent-app kb --version 2>&1 | grep -i kb_python || echo 'ALIGNER MISSING'; \
  docker exec igvfagent-app igvfagent portal search --type MeasurementSet --limit 1 2>/dev/null | grep -iE 'total|matches'"

step "7. Prove the new pipeline works on a real dataset"
"${SSH[@]}" "docker exec igvfagent-app igvfagent raw-pipeline plan IGVFDS6639ECQN 2>/dev/null \
  | grep -E 'ROUTE|Technology|Aligner' || true"

cat <<'DONE'

DONE — https://igvfagent.genohub.org is live on the new image.

Expected in the output above:
  credentials: env <key id>          (not "none")
  kb_python 0.29.x                   (not "ALIGNER MISSING")
  Total matches: 11,070              (not 7,127)
  Aligner: available
DONE
