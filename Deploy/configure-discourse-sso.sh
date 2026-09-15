#!/usr/bin/env bash
#
# Point the Genohub Discourse community at IGVF Agent as a login provider.
#
#   bash Deploy/configure-discourse-sso.sh --check     # read-only; change nothing
#   bash Deploy/configure-discourse-sso.sh             # apply
#   bash Deploy/configure-discourse-sso.sh --add alice bob
#
# Four changes, all on the forum side:
#
#   1. create the `igvfagent` group, which IS the approval gate
#   2. add the named users to it (defaults to the forum admins)
#   3. set `discourse connect provider secrets` to igvfagent.genohub.org|<secret>
#   4. switch `enable discourse connect provider` on
#
# Neither secret is ever echoed. The API key is read from
# Docs/Secret/discourse-API.txt (gitignored) and the SSO secret straight out of
# Deploy/.env.prod on the deployment host, so nothing has to be pasted and
# nothing lands in a shell history or a terminal scrollback.
#
# Safe to re-run: every step checks the current state first and says "already"
# rather than doing it twice.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

DISCOURSE="${IGVF_DISCOURSE_URL:-https://discussion.genohub.org}"
SITE_HOST="${IGVF_SITE_HOST:-igvfagent.genohub.org}"
GROUP="${IGVF_APPROVAL_GROUP:-igvfagent}"
KEYFILE="${IGVF_DISCOURSE_KEYFILE:-$ROOT/Docs/Secret/discourse-API.txt}"
UA='Mozilla/5.0 (compatible; IGVFagent-setup/1.0)'

CHECK=0
MEMBERS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --check) CHECK=1 ;;
        --add)   shift; while [ $# -gt 0 ] && [ "${1:0:2}" != "--" ]; do MEMBERS+=("$1"); shift; done; continue ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

[ -r "$KEYFILE" ] || { echo "No Discourse API key at $KEYFILE" >&2; exit 2; }
KEY="$(tr -d ' \r\n' < "$KEYFILE")"
[ -n "$KEY" ] || { echo "$KEYFILE is empty" >&2; exit 2; }

api() {  # api METHOD PATH [JSON-BODY]
    local method="$1" path="$2" body="${3:-}"
    if [ -n "$body" ]; then
        curl -sS -m 30 -X "$method" \
             -H "Api-Key: $KEY" -H "Api-Username: system" -H "User-Agent: $UA" \
             -H "Accept: application/json" -H "Content-Type: application/json" \
             -d "$body" "$DISCOURSE$path"
    else
        curl -sS -m 30 -X "$method" \
             -H "Api-Key: $KEY" -H "Api-Username: system" -H "User-Agent: $UA" \
             -H "Accept: application/json" "$DISCOURSE$path"
    fi
}

jq_py() { python3 -c "$1" 2>/dev/null; }

echo "forum:  $DISCOURSE"
echo "site:   $SITE_HOST"
echo "group:  $GROUP"
echo

# --------------------------------------------------------------- state ------
echo "== current state =="
SETTINGS="$(api GET /admin/site_settings.json)"
echo "$SETTINGS" | jq_py "
import json,sys
ss={s['setting']: s for s in json.load(sys.stdin).get('site_settings',[])}
if not ss: print('  ERROR: no settings returned — is the API key an admin key?'); raise SystemExit(1)
p = ss.get('enable_discourse_connect_provider',{}).get('value')
sec = ss.get('discourse_connect_provider_secrets',{}).get('value') or ''
print(f'  enable_discourse_connect_provider : {p}')
print(f'  discourse_connect_provider_secrets: ' + (f'{len(sec)} chars set' if sec else '(empty)'))
" || exit 1

GROUP_ID="$(api GET /groups.json | jq_py "
import json,sys
d=json.load(sys.stdin); gs=d.get('groups', d if isinstance(d,list) else [])
print(next((str(g['id']) for g in gs if g.get('name','').lower()=='${GROUP}'.lower()), ''))
")"
echo "  group '$GROUP'                     : $([ -n "$GROUP_ID" ] && echo "exists (id $GROUP_ID)" || echo 'missing')"

if [ "$CHECK" = 1 ]; then
    echo
    echo "--check: nothing was changed."
    exit 0
fi

# ------------------------------------------------------------ 1. group ------
echo
echo "== 1/4  the approval group =="
if [ -n "$GROUP_ID" ]; then
    echo "  already exists (id $GROUP_ID)"
else
    GROUP_ID="$(api POST /admin/groups.json "$(cat <<JSON
{"group":{"name":"$GROUP","full_name":"IGVF Agent users",
"bio_raw":"Members can sign in to $SITE_HOST. Membership is the approval step — signing up to the community does not grant it on its own.",
"visibility_level":1,"members_visibility_level":1,
"public_admission":false,"public_exit":false,
"allow_membership_requests":false,"mentionable_level":0,"messageable_level":2}}
JSON
)" | jq_py "
import json,sys
d=json.load(sys.stdin); g=d.get('basic_group') or d
print(g.get('id','') if g.get('id') else '')
")"
    [ -n "$GROUP_ID" ] && echo "  created (id $GROUP_ID)" \
        || { echo "  FAILED to create the group"; exit 1; }
fi

# ---------------------------------------------------------- 2. members ------
echo
echo "== 2/4  members =="
if [ "${#MEMBERS[@]:-0}" -eq 0 ]; then
    # Default to the forum's own admins: they can grant themselves access
    # anyway, and starting with an empty group means nobody can get in.
    # A while-read loop rather than mapfile, because macOS ships bash 3.2.
    while IFS= read -r u; do
        [ -n "$u" ] && MEMBERS+=("$u")
    done < <(api GET "/admin/users/list/active.json?limit=100" | jq_py "
import json,sys
for u in json.load(sys.stdin):
    if u.get('admin') and u.get('username') not in ('system','discobot'):
        print(u['username'])
")
    echo "  no --add given; defaulting to forum admins"
fi
if [ "${#MEMBERS[@]:-0}" -eq 0 ]; then
    echo "  no members to add — pass --add <username>" >&2; exit 1
fi
echo "  adding: ${MEMBERS[*]}"
api PUT "/groups/$GROUP_ID/members.json" \
    "{\"usernames\":\"$(IFS=,; echo "${MEMBERS[*]}")\"}" | jq_py "
import json,sys
d=json.load(sys.stdin)
if d.get('success') or d.get('usernames') is not None: print('  ok')
else: print('  response:', json.dumps(d)[:240])
"

# ----------------------------------------------------------- 3. secret ------
echo
echo "== 3/4  the shared secret =="
SSO_SECRET="$(bash "$HERE/vm.sh" ssh "grep '^IGVF_SSO_SECRET=' /srv/igvfagent/Deploy/.env.prod | cut -d= -f2-" 2>/dev/null | tr -d ' \r\n')"
if [ -z "$SSO_SECRET" ]; then
    echo "  could not read IGVF_SSO_SECRET from the host's Deploy/.env.prod" >&2
    exit 1
fi
echo "  read from the host (${#SSO_SECRET} chars, not shown)"
api PUT "/admin/site_settings/discourse_connect_provider_secrets.json" \
    "$(python3 - "$SITE_HOST" "$SSO_SECRET" <<'JSON'
import json, sys
print(json.dumps({"discourse_connect_provider_secrets": f"{sys.argv[1]}|{sys.argv[2]}"}))
JSON
)" >/dev/null && echo "  set for $SITE_HOST"
unset SSO_SECRET

# ------------------------------------------------------------ 4. enable -----
echo
echo "== 4/4  enable the provider =="
api PUT "/admin/site_settings/enable_discourse_connect_provider.json" \
    '{"enable_discourse_connect_provider":"true"}' >/dev/null && echo "  enabled"

# ------------------------------------------------------------- verify -------
echo
echo "== verify =="
sleep 2
CODE="$(curl -sS -m 20 -o /dev/null -w '%{http_code}' -H "User-Agent: $UA" \
        "$DISCOURSE/session/sso_provider")"
echo "  $DISCOURSE/session/sso_provider -> HTTP $CODE"
if [ "$CODE" = "404" ]; then
    echo "  STILL 404 — the provider did not take. Check the setting by hand."
    exit 1
fi
echo "  provider is live (404 would mean still off)"
echo
echo "  Now open https://$SITE_HOST/ and sign in."
echo "  The landing page re-checks every 5 minutes, so its warning banner"
echo "  clears on its own."
