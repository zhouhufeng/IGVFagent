#!/usr/bin/env bash
#
# Redeploy IGVFagent on the host, and PROVE the new code is running.
#
# The trap this exists to close: the app container does not mount the
# source tree. compose mounts only /mnt/igvf-data/{Data,Docs}, and
# Scripts/streamlit_app.py imports `from igvfagent import ...` — the
# package installed into /opt/venv by `pip install '.[all]'` at image
# BUILD time (Dockerfile). So `git pull` updates the host checkout while
# the container keeps serving whatever was baked into the image, and the
# site looks unchanged no matter how many times you pull.
#
# A plain `docker build` can also quietly do nothing here: run it from
# the wrong directory and it builds the wrong context; let compose reuse
# a cached layer or an unchanged image id and the container is never
# recreated. This script uses compose's own build with --no-cache and
# --force-recreate, then verifies by comparing the Python actually
# installed in the container against the repo checkout.
#
# Usage:  bash Deploy/redeploy.sh          # pull, rebuild, recreate, verify
#         bash Deploy/redeploy.sh --check  # verify only, change nothing
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
COMPOSE=(docker compose --env-file "$HERE/.env.prod" -f "$HERE/docker-compose.prod.yml")
CONTAINER=igvfagent-app
CHECK_ONLY=0
[ "${1:-}" = "--check" ] && CHECK_ONLY=1

# Compare only the TOP-LEVEL modules (Scripts/*.py -> igvfagent/*.py).
# Those are exactly the files `pip install` ships, so the two sides have
# the same file set; hashing the whole tree would false-alarm on any
# subdirectory that is not a package. Hashes are of CONTENT only, sorted,
# so differing paths do not matter.
SHA=""
for c in sha256sum "shasum -a 256"; do
    if command -v ${c%% *} >/dev/null 2>&1; then SHA="$c"; break; fi
done
[ -n "$SHA" ] || { echo "need sha256sum or shasum on the host"; exit 2; }

repo_hash() {
    ( cd "$ROOT/Scripts" && find . -maxdepth 1 -name '*.py' -type f | sort \
        | while read -r f; do $SHA "$f"; done ) \
        | awk '{print $1}' | sort | $SHA | awk '{print $1}'
}

container_hash() {
    docker exec "$CONTAINER" sh -c '
        d=$(ls -d /opt/venv/lib/python*/site-packages/igvfagent 2>/dev/null | head -1)
        [ -n "$d" ] || exit 3
        cd "$d" || exit 3
        find . -maxdepth 1 -name "*.py" -type f | sort \
          | while read -r f; do sha256sum "$f"; done \
          | awk "{print \$1}" | sort | sha256sum
    ' 2>/dev/null | awk '{print $1}'
}

verify() {
    if ! docker ps --format '{{.Names}}' | grep -qx "$CONTAINER"; then
        echo "  container $CONTAINER is not running"
        return 1
    fi
    local r c
    r=$(repo_hash); c=$(container_hash || true)
    echo "  repo      $(cd "$ROOT" && git rev-parse --short HEAD)  code hash ${r:0:12}"
    echo "  container                code hash ${c:0:12}"
    if [ -z "$c" ]; then
        echo "  COULD NOT READ the installed package inside the container"
        return 1
    fi
    if [ "$r" = "$c" ]; then
        echo "  MATCH — the running site is serving this checkout"
        return 0
    fi
    echo "  STALE — the container is running DIFFERENT code from this checkout"
    return 1
}

if [ "$CHECK_ONLY" = 1 ]; then
    echo "Checking deployed code against $ROOT ..."
    verify
    exit $?
fi

echo "==> 1/4  updating the checkout"
git -C "$ROOT" pull --ff-only

echo "==> 2/4  rebuilding the image (no cache — a reused layer is how this silently no-ops)"
"${COMPOSE[@]}" build --no-cache app

echo "==> 3/4  recreating the container"
"${COMPOSE[@]}" up -d --force-recreate app

echo "==> 4/4  verifying the container actually runs the new code"
for i in 1 2 3 4 5 6 7 8 9 10; do
    sleep 3
    if docker ps --format '{{.Names}}' | grep -qx "$CONTAINER"; then break; fi
    echo "  waiting for $CONTAINER to come up ($i/10)"
done

if verify; then
    echo
    echo "Redeploy complete. Hard-refresh the browser (Cmd/Ctrl+Shift+R)."
    exit 0
fi

echo
echo "Redeploy did NOT take effect. The build ran but the container is still"
echo "on old code. Most likely causes, in order:"
echo "  - another container is serving the site:   docker ps | grep igvf"
echo "  - the image tag did not move:              docker images igvfagent"
echo "  - the build used a different context:      check Deploy/docker-compose.prod.yml build.context"
exit 1
