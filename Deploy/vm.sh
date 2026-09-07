#!/usr/bin/env bash
# The one place that knows how to reach the production VM.
#
#   bash Deploy/vm.sh doctor              # why can't I connect? (start here)
#   bash Deploy/vm.sh ssh 'docker ps'     # run a command on the VM
#   bash Deploy/vm.sh install-key PATH    # adopt a key file into ~/.ssh
#   bash Deploy/vm.sh authorize PUBKEY    # let another machine in (needs access)
#
# Also sourceable: `source Deploy/vm.sh` defines vm_key/vm_host/vm_ssh and
# runs no subcommand, so Deploy/make-live.sh resolves the connection exactly
# the way this script does.
#
# Why this file exists. The private half of the OpenStack keypair the VM was
# booted with is a secret, so it is gitignored and does NOT travel with a
# `git clone`. Every checkout on a new machine therefore starts with no way
# in, and the previous hard-coded path (Docs/Secret/igvfagent-deploy.pem)
# turned that into "file not found" — which reads like a bug in the deploy
# script rather than a missing credential. The resolver below searches the
# places the key actually lives, and `doctor` names the one thing that is
# wrong instead of leaving you to bisect ssh flags.
set -euo pipefail

# Public, non-secret facts about the deployment.
VM_HOST_DEFAULT="ubuntu@149.165.151.21"
VM_INSTANCE="igvfagent-prod"           # OpenStack instance name (Jetstream2 / IU)
VM_KEYPAIR="igvfagent-deploy"          # OpenStack keypair it was booted with
VM_KEY_MD5="a5:e9:e4:c0:f3:e0:13:fd:4f:8d:dc:97:3f:9e:a4:32"
VM_REMOTE_ROOT="/srv/igvfagent"

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

# Staged 0600 copy, if the key lives somewhere its mode cannot be fixed
# (an exFAT/NTFS external drive shows every file as 0777 and ssh then
# refuses the key outright — "UNPROTECTED PRIVATE KEY FILE").
# The path is fixed per-process rather than mktemp'd, because vm_key is
# normally called in a command substitution — a subshell, whose assignment to
# a variable the parent's trap reads would be lost, leaking the directory.
_VM_STAGE="${TMPDIR:-/tmp}/igvfagent-vm-$$"
_vm_cleanup() { rm -rf "$_VM_STAGE"; }
trap _vm_cleanup EXIT

vm_host() { printf '%s' "${IGVFAGENT_VM_HOST:-$VM_HOST_DEFAULT}"; }

# Every place the key plausibly lives, best first. $IGVFAGENT_DEPLOY_KEY wins
# so a machine that stores it anywhere at all can say so once and be done.
vm_key_candidates() {
  local c=()
  [[ -n "${IGVFAGENT_DEPLOY_KEY:-}" ]] && c+=("$IGVFAGENT_DEPLOY_KEY")
  c+=("$HOME/.ssh/${VM_KEYPAIR}.pem"
      "$HOME/.ssh/${VM_KEYPAIR}"
      "$HOME/.ssh/${VM_KEYPAIR}.key"
      "$REPO_ROOT/Docs/Secret/${VM_KEYPAIR}.pem"
      "$REPO_ROOT/Docs/Secrets/${VM_KEYPAIR}.pem"
      "$REPO_ROOT/Docs/Secretes/${VM_KEYPAIR}.pem")
  printf '%s\n' "${c[@]}"
}

vm_key_path() {   # first candidate that exists; empty if none
  local p
  while IFS= read -r p; do [[ -f "$p" ]] && { printf '%s' "$p"; return 0; }; done \
    < <(vm_key_candidates)
  return 1
}

vm_key_md5() {    # MD5 fingerprint of a private key, or empty
  ssh-keygen -E md5 -yf "$1" 2>/dev/null \
    | ssh-keygen -E md5 -lf /dev/stdin 2>/dev/null \
    | awk '{print $2}' | sed 's/^MD5://'
}

# Usable path to the key: the original if its mode is already private,
# otherwise a 0600 copy in a private temp dir that dies with this process.
vm_key() {
  local key mode
  key="$(vm_key_path)" || return 1
  mode="$(stat -c '%a' "$key" 2>/dev/null || stat -f '%Lp' "$key" 2>/dev/null || echo 600)"
  if [[ "$mode" != "600" && "$mode" != "400" ]]; then
    chmod 600 "$key" 2>/dev/null || true
    mode="$(stat -c '%a' "$key" 2>/dev/null || stat -f '%Lp' "$key" 2>/dev/null || echo 600)"
  fi
  if [[ "$mode" != "600" && "$mode" != "400" ]]; then
    mkdir -p "$_VM_STAGE"; chmod 700 "$_VM_STAGE"
    cp "$key" "$_VM_STAGE/key"; chmod 600 "$_VM_STAGE/key"
    key="$_VM_STAGE/key"
  fi
  printf '%s' "$key"
}

vm_ssh() {        # vm_ssh [ssh-args...] -- runs against the VM with the right key
  local key; key="$(vm_key)" || {
    printf 'no %s private key on this machine — run: bash Deploy/vm.sh doctor\n' \
      "$VM_KEYPAIR" >&2; return 78; }
  ssh -i "$key" -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new \
      -o ConnectTimeout=20 "$(vm_host)" "$@"
}

# ---------------------------------------------------------------- subcommands

cmd_doctor() {
  local ok=1 key mode fp
  printf '=== IGVFagent VM connectivity ===\n'
  printf 'repo      %s\n' "$REPO_ROOT"
  printf 'host      %s   (OpenStack instance %s, keypair %s)\n' \
    "$(vm_host)" "$VM_INSTANCE" "$VM_KEYPAIR"

  printf '\n--- 1. private key ---\n'
  if key="$(vm_key_path)"; then
    mode="$(stat -c '%a' "$key" 2>/dev/null || stat -f '%Lp' "$key" 2>/dev/null)"
    printf 'found     %s (mode %s)\n' "$key" "$mode"
    fp="$(vm_key_md5 "$key")"
    if [[ -z "$fp" ]]; then
      printf 'FAIL      not a readable private key (passphrase-protected, or corrupt)\n'; ok=0
    elif [[ "$fp" == "$VM_KEY_MD5" ]]; then
      printf 'ok        fingerprint matches OpenStack keypair %s\n' "$VM_KEYPAIR"
    else
      # Not fatal: a machine added later with `authorize` uses its own key,
      # which is the better arrangement anyway (revoke one machine, not all).
      # Step 3 is the real verdict, so fall through to it.
      printf 'note      fingerprint %s\n' "$fp"
      printf '          Not the original boot keypair (%s).\n' "$VM_KEY_MD5"
      printf '          Fine if this machine was added with `vm.sh authorize`;\n'
      printf '          step 3 below settles it either way.\n'
    fi
  else
    printf 'MISSING   no key at any known location. Looked in:\n'
    vm_key_candidates | sed 's/^/            /'
    printf '\n          The key is a secret and is gitignored, so it never arrives\n'
    printf '          with a git clone. Two ways to fix this machine:\n\n'
    printf '          (a) You have the key file somewhere:\n'
    printf '                bash Deploy/vm.sh install-key /path/to/%s.pem\n\n' "$VM_KEYPAIR"
    printf '          (b) You do not, but another machine can already reach the VM.\n'
    printf '              Here:   ssh-keygen -t ed25519 -f ~/.ssh/%s -C "$(whoami)@$(hostname)"\n' "$VM_KEYPAIR"
    printf '              There:  bash Deploy/vm.sh authorize ~/path/to/%s.pub\n' "$VM_KEYPAIR"
    printf '                      (copy this machine'"'"'s .pub over first — it is not secret)\n'
    ok=0
  fi

  printf '\n--- 2. reachability ---\n'
  local hostonly; hostonly="$(vm_host)"; hostonly="${hostonly##*@}"
  if command -v nc >/dev/null 2>&1 && nc -z -w 8 "$hostonly" 22 2>/dev/null; then
    printf 'ok        %s:22 accepts connections\n' "$hostonly"
  elif timeout 10 bash -c "</dev/tcp/$hostonly/22" 2>/dev/null; then
    printf 'ok        %s:22 accepts connections\n' "$hostonly"
  else
    printf 'FAIL      cannot open %s:22 — network, VPN, or the instance is down.\n' "$hostonly"
    printf '          Check state:  openstack server show %s\n' "$VM_INSTANCE"; ok=0
  fi

  printf '\n--- 3. authentication ---\n'
  if [[ "$ok" == 1 ]] && vm_ssh -o BatchMode=yes true 2>/dev/null; then
    printf 'ok        public-key auth accepted\n'
    printf '\n--- 4. the deployment ---\n'
    vm_ssh -o BatchMode=yes "cd $VM_REMOTE_ROOT 2>/dev/null && \
      echo \"checkout  \$(git log --oneline -1)\"; \
      echo \"docker    \$(docker --version 2>&1 | head -1)\"; \
      echo \"running   \$(docker ps --format '{{.Names}}' 2>/dev/null | tr '\n' ' ')\"" \
      || printf 'note      logged in, but %s looks unusual\n' "$VM_REMOTE_ROOT"
    printf '\nAll good — Deploy/make-live.sh will work from this machine.\n'
  elif [[ "$ok" != 1 ]]; then
    printf 'skipped   no usable key — fix step 1 first\n'
  else
    # The key exists and the host answers, so the VM simply does not know
    # this key. Nothing local can fix that; it takes an action on a machine
    # that is already trusted.
    printf 'FAIL      public-key auth refused — the VM does not know this key\n'
    printf '          On a machine that already has access, run:\n'
    printf '            bash Deploy/vm.sh authorize <this-machine>.pub\n'
    printf '          This machine'"'"'s public key (not secret, copy it over):\n'
    local pubkey; pubkey="$(vm_key_path)"
    if [[ -f "${pubkey}.pub" ]]; then
      sed 's/^/            /' "${pubkey}.pub"
    else
      ssh-keygen -yf "$pubkey" 2>/dev/null | sed 's/^/            /' \
        || printf '            (could not derive it from %s)\n' "$pubkey"
    fi
    ok=0
  fi
  [[ "$ok" == 1 ]]
}

cmd_install_key() {
  local src="${1:-}" dest="$HOME/.ssh/${VM_KEYPAIR}.pem" fp
  [[ -f "$src" ]] || { printf 'usage: bash Deploy/vm.sh install-key /path/to/%s.pem\n' "$VM_KEYPAIR" >&2; exit 2; }
  fp="$(vm_key_md5 "$src")"
  [[ -n "$fp" ]] || { printf 'not a usable private key: %s\n' "$src" >&2; exit 1; }
  if [[ "$fp" != "$VM_KEY_MD5" ]]; then
    printf 'note: %s is not the original boot keypair (%s vs %s).\n' \
      "$src" "$fp" "$VM_KEY_MD5" >&2
    printf 'Installing anyway — this is correct if the VM was told about it with\n' >&2
    printf '`vm.sh authorize`. The doctor run below is what actually decides.\n\n' >&2
  fi
  mkdir -p "$HOME/.ssh"; chmod 700 "$HOME/.ssh"
  cp "$src" "$dest"; chmod 600 "$dest"
  printf 'installed %s (0600)\n' "$dest"
  printf 'This is outside the repo on purpose: %s never gets committed by accident.\n' "$dest"
  cmd_doctor
}

cmd_authorize() {
  local pub="${1:-}" yes="${2:-}" body comment
  [[ -f "$pub" ]] || { printf 'usage: bash Deploy/vm.sh authorize /path/to/new-machine.pub [--yes]\n' >&2; exit 2; }
  ssh-keygen -lf "$pub" >/dev/null 2>&1 || { printf 'not a public key: %s\n' "$pub" >&2; exit 1; }
  body="$(awk '{print $1" "$2}' "$pub")"
  # The comment is free text in the .pub file but ends up inside a quoted
  # remote command, so a stray quote would break the argument boundary.
  # Only the characters a comment is ever meaningfully made of survive.
  comment="$(awk '{print $3}' "$pub" | tr -cd 'A-Za-z0-9@._-')"
  printf 'About to grant SSH access to the LIVE deployment (%s):\n' "$(vm_host)"
  printf '  key      %s\n  comment  %s\n' "$(ssh-keygen -lf "$pub" | awk '{print $2}')" "${comment:-<none>}"
  if [[ "$yes" != "--yes" ]]; then
    read -r -p 'Append it to authorized_keys? [y/N] ' a
    [[ "$a" == y || "$a" == Y ]] || { echo 'aborted'; exit 1; }
  fi
  # The key travels as an argument, not on stdin: the heredoc below is itself
  # stdin for the remote bash, so a piped value would be swallowed by it. That
  # is safe here precisely because a PUBLIC key is not secret — unlike the
  # Portal secret in make-live.sh, which must stay off the command line.
  #
  # Idempotent on the key body rather than the whole line, so re-running after
  # the comment changed does not leave two copies.
  vm_ssh "bash -s -- '$body' '${comment:-added-by-vm.sh}'" <<'REMOTE'
set -euo pipefail
BODY="$1"
COMMENT="${2:-added-by-vm.sh}"
mkdir -p ~/.ssh; chmod 700 ~/.ssh
touch ~/.ssh/authorized_keys; chmod 600 ~/.ssh/authorized_keys
if grep -qF "$BODY" ~/.ssh/authorized_keys; then
  echo "already authorized — no change"
else
  # Some keyfiles lack a trailing newline; appending blind would splice this
  # key onto the end of the previous one and silently break both.
  [ -s ~/.ssh/authorized_keys ] && [ "$(tail -c1 ~/.ssh/authorized_keys | wc -l)" -eq 0 ] \
    && printf '\n' >> ~/.ssh/authorized_keys
  printf '%s %s\n' "$BODY" "$COMMENT" >> ~/.ssh/authorized_keys
  echo "authorized"
fi
echo "keys now on file: $(grep -c '^ssh-\|^ecdsa-' ~/.ssh/authorized_keys)"
REMOTE
}

cmd_ssh() { vm_ssh "$@"; }

# Sourced (by make-live.sh) — define the functions and stop.
(return 0 2>/dev/null) && return 0

case "${1:-doctor}" in
  doctor|check|diagnose) shift || true; cmd_doctor "$@" ;;
  install-key)           shift; cmd_install_key "$@" ;;
  authorize)             shift; cmd_authorize "$@" ;;
  ssh)                   shift; cmd_ssh "$@" ;;
  *) printf 'usage: bash Deploy/vm.sh {doctor|ssh <cmd>|install-key <pem>|authorize <pub>}\n' >&2; exit 2 ;;
esac
