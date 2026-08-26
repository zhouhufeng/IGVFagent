#!/usr/bin/env bash
# Prepare a freshly-booted Jetstream2 VM to host IGVFagent.
#
# Idempotent: safe to re-run. Formats the attached Cinder volume ONLY if it
# has no filesystem, mounts it at /mnt/igvf-data, and lays out the workspace
# the app container bind-mounts.
#
#   sudo bash Deploy/bootstrap-vm.sh [DEVICE]
#
# DEVICE defaults to the first unmounted, unformatted disk >100G.
set -euo pipefail

MOUNT=/mnt/igvf-data

log() { printf '\n== %s\n' "$*"; }

# ---------------------------------------------------------------- find disk
DEV="${1:-}"
if [[ -z "$DEV" ]]; then
  # Pick the largest block device that has no filesystem and no mountpoint —
  # the root disk always has both, so it can never be selected here.
  DEV=$(lsblk -bdnp -o NAME,SIZE,TYPE,FSTYPE,MOUNTPOINT \
        | awk '$3=="disk" && $4=="" && $5=="" && $2>100*1024^3 {print $1, $2}' \
        | sort -k2 -n | tail -1 | awk '{print $1}')
fi

if [[ -z "$DEV" ]]; then
  if mountpoint -q "$MOUNT"; then
    log "$MOUNT already mounted; nothing to format."
  else
    echo "ERROR: no candidate data disk found and $MOUNT is not mounted." >&2
    lsblk -p -o NAME,SIZE,TYPE,FSTYPE,MOUNTPOINT >&2
    exit 1
  fi
else
  log "Data disk: $DEV ($(lsblk -bdn -o SIZE "$DEV" | numfmt --to=iec))"
  if ! blkid "$DEV" >/dev/null 2>&1; then
    log "No filesystem present — creating XFS (good for many small artefact files)"
    mkfs.xfs -q "$DEV"
  else
    log "Filesystem already present, leaving it alone: $(blkid -o value -s TYPE "$DEV")"
  fi

  UUID=$(blkid -o value -s UUID "$DEV")
  mkdir -p "$MOUNT"
  # Mount by UUID: Cinder device names are not stable across reboots.
  if ! grep -q "$UUID" /etc/fstab; then
    log "Adding to /etc/fstab by UUID=$UUID"
    printf 'UUID=%s  %s  xfs  defaults,noatime,nofail  0  2\n' "$UUID" "$MOUNT" >> /etc/fstab
  fi
  mountpoint -q "$MOUNT" || mount "$MOUNT"
fi

# ------------------------------------------------------------- workspace
log "Preparing workspace under $MOUNT"
mkdir -p "$MOUNT"/Data "$MOUNT"/Docs
# uid/gid 1000 == the non-root `igvf` user inside the container (Dockerfile).
chown -R 1000:1000 "$MOUNT"/Data "$MOUNT"/Docs
chmod 755 "$MOUNT"/Data "$MOUNT"/Docs

log "Disk state"
df -h "$MOUNT"

log "Done. Next: populate Deploy/.env.prod, then"
echo "  cd /srv/igvfagent/Deploy && docker compose -f docker-compose.prod.yml up -d"
