#!/usr/bin/env bash
# Relocate container image storage from the 58 GB root disk to the 2 TB
# data volume. Run ON THE VM as a user with sudo.
#
#   sudo bash Deploy/move-container-storage.sh --check   # report only
#   sudo bash Deploy/move-container-storage.sh           # do it
#
# WHY. Every redeploy builds with --no-cache and leaves its layers behind.
# A dozen deploys in one day took /var/lib/containerd to 31 GB and the root
# filesystem to 73%; the data volume it sits beside is 7% full. Pruning
# recovers the space but the next dozen deploys take it back, so the storage
# belongs on the big disk.
#
# NOT analysis data: results, the FASTQ cache and the run history already
# live on /mnt/igvf-data and are untouched by this. What moves is image
# layers, which cost only rebuild time if lost.
#
# The site is DOWN while this runs -- both daemons stop. Expect a few
# minutes, mostly the copy.
set -euo pipefail

DATA_MNT="/mnt/igvf-data"
NEW_CONTAINERD="$DATA_MNT/containerd"
NEW_DOCKER="$DATA_MNT/docker"
CTR_CFG="/etc/containerd/config.toml"
DOCKER_CFG="/etc/docker/daemon.json"
CHECK_ONLY=0
[[ "${1:-}" == "--check" ]] && CHECK_ONLY=1

die()  { printf '\nFAILED: %s\n' "$1" >&2; exit 1; }
step() { printf '\n=== %s ===\n' "$1"; }

[[ $EUID -eq 0 ]] || die "run with sudo"

step "1. Preconditions"
findmnt -no TARGET "$DATA_MNT" >/dev/null || die "$DATA_MNT is not mounted"
FSTYPE=$(findmnt -no FSTYPE "$DATA_MNT")
echo "  $DATA_MNT: $FSTYPE"
if [[ "$FSTYPE" == "xfs" ]]; then
  # overlay refuses to work on XFS formatted with ftype=0, and the failure
  # is at container start, long after the data has moved. Check first.
  xfs_info "$DATA_MNT" | grep -q "ftype=1" \
    || die "XFS on $DATA_MNT has ftype=0; overlay storage cannot live there"
  echo "  ftype=1 (overlay OK)"
fi
NEED=$(du -xsb /var/lib/containerd /var/lib/docker 2>/dev/null | awk '{s+=$1} END{print s}')
FREE=$(df -B1 --output=avail "$DATA_MNT" | tail -1)
printf '  need %.1f GB, free %.1f GB\n' "$(bc -l <<<"$NEED/1073741824")" \
                                        "$(bc -l <<<"$FREE/1073741824")"
[[ "$FREE" -gt $(( NEED * 2 )) ]] || die "not enough headroom on $DATA_MNT"
df -h / "$DATA_MNT" | awk 'NR>1{printf "  %-16s %5s of %5s (%s)\n", $6,$3,$2,$5}'

if [[ "$CHECK_ONLY" == 1 ]]; then
  printf '\n--check complete. Nothing was changed.\n'
  exit 0
fi

step "2. Stop the stack"
systemctl stop docker.socket 2>/dev/null || true
systemctl stop docker
systemctl stop containerd
sleep 2
pgrep -x containerd >/dev/null && die "containerd still running"
echo "  docker + containerd stopped"

step "3. Copy the data (original kept until verified)"
mkdir -p "$NEW_CONTAINERD" "$NEW_DOCKER"
# -a preserves the hardlinks overlay layers rely on; a plain cp would
# explode the size and break layer sharing.
rsync -aHAX --numeric-ids /var/lib/containerd/ "$NEW_CONTAINERD/"
rsync -aHAX --numeric-ids /var/lib/docker/     "$NEW_DOCKER/"
echo "  copied: $(du -xsh "$NEW_CONTAINERD" | cut -f1) containerd, "\
"$(du -xsh "$NEW_DOCKER" | cut -f1) docker"

step "4. Point the daemons at the new location"
cp -n "$CTR_CFG" "$CTR_CFG.bak.$(date +%Y%m%d%H%M%S)" 2>/dev/null || true
if grep -qE '^\s*root\s*=' "$CTR_CFG"; then
  sed -i "s|^\s*root\s*=.*|root = \"$NEW_CONTAINERD\"|" "$CTR_CFG"
else
  # The key must sit BEFORE the first [table] header or TOML parses it as
  # part of that table and containerd ignores it silently.
  awk -v val="$NEW_CONTAINERD" '
    !done && /^\[/ { print "root = \"" val "\""; print ""; done=1 }
    { print }
    END { if (!done) print "root = \"" val "\"" }
  ' "$CTR_CFG" > "$CTR_CFG.new" && mv "$CTR_CFG.new" "$CTR_CFG"
fi
grep -nE '^root\s*=' "$CTR_CFG" | sed 's/^/  containerd: /'

cp -n "$DOCKER_CFG" "$DOCKER_CFG.bak.$(date +%Y%m%d%H%M%S)" 2>/dev/null || true
python3 - "$DOCKER_CFG" "$NEW_DOCKER" <<'PY'
import json, sys
path, root = sys.argv[1], sys.argv[2]
try:
    with open(path) as fh:
        cfg = json.load(fh)
except (OSError, ValueError):
    cfg = {}
cfg["data-root"] = root
with open(path, "w") as fh:
    json.dump(cfg, fh, indent=4)
print(f"  docker: data-root = {root}")
PY

step "5. Start and verify"
systemctl start containerd
systemctl start docker
sleep 5
systemctl is-active --quiet containerd || die "containerd did not start"
systemctl is-active --quiet docker || die "docker did not start"
ACTUAL=$(docker info --format '{{.DockerRootDir}}')
echo "  docker root: $ACTUAL"
[[ "$ACTUAL" == "$NEW_DOCKER"* ]] || die "docker ignored data-root ($ACTUAL)"
docker image inspect igvfagent:latest >/dev/null 2>&1 \
  || die "igvfagent:latest is missing after the move — do NOT delete the old data"
echo "  igvfagent:latest present"
cd /srv/igvfagent/Deploy
docker compose --env-file .env.prod -f docker-compose.prod.yml up -d
sleep 10
docker exec igvfagent-app curl -fsS -o /dev/null http://localhost:8501/_stcore/health \
  && echo "  app health: OK" || die "app is not healthy"

step "6. Result"
df -h / "$DATA_MNT" | awk 'NR>1{printf "  %-16s %5s of %5s (%s)\n", $6,$3,$2,$5}'
cat <<'DONE'

Verified. The OLD data is still in place, so the root disk has not shrunk
yet. Reclaim it only now that the site is confirmed working:

  sudo rm -rf /var/lib/containerd.old /var/lib/docker.old   # if you renamed
  sudo mv /var/lib/containerd /var/lib/containerd.old
  sudo mv /var/lib/docker     /var/lib/docker.old
  # confirm the site still works, then delete the .old directories

Rollback: remove the "root" key from /etc/containerd/config.toml and
"data-root" from /etc/docker/daemon.json, then restart both daemons — the
original directories are untouched.
DONE
