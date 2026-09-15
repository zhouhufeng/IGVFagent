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
#
# The group is both the approval gate AND the place people apply. Discourse
# has a native request-and-approve flow -- allow_membership_requests puts a
# "Request Membership" button on the group page and files each request for a
# group owner to accept or deny -- so there is no application form to build
# and no queue to invent. public_admission stays false or people would simply
# join themselves, which would make the whole gate decorative.
GROUP_BODY="$(python3 - "$GROUP" "$SITE_HOST" <<'JSON'
import json, sys
group, site = sys.argv[1], sys.argv[2]
print(json.dumps({"group": {
    "name": group,
    "full_name": "IGVF Agent Users",
    "bio_raw": (
        f"Members of this group can sign in to **{site}**.\n\n"
        "Membership is the approval step: having an account on this community "
        "does not by itself grant access to the agent, because analyses run on "
        "shared hardware.\n\n"
        "Use **Request Membership** below and say briefly who you are and what "
        "you plan to use it for. An administrator will review it."),
    "visibility_level": 0,           # the group page is publicly linkable
    "members_visibility_level": 2,   # ...but the member list is not public
    "public_admission": False,       # nobody joins themselves
    "public_exit": True,             # leaving needs no ceremony
    "allow_membership_requests": True,
    "membership_request_template": (
        "Who you are, your institution or lab, and what you plan to use "
        "IGVF Agent for."),
    "mentionable_level": 0,
    "messageable_level": 2,
}}))
JSON
)"

echo
echo "== 1/4  the approval group =="
if [ -n "$GROUP_ID" ]; then
    echo "  exists (id $GROUP_ID)"
else
    GROUP_ID="$(api POST /admin/groups.json "$GROUP_BODY" | jq_py "
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
errs=d.get('errors') or []
if d.get('success') or d.get('usernames') is not None:
    print('  ok')
elif errs and all('already members' in e for e in errs):
    print('  already members — nothing to do')
else:
    print('  response:', json.dumps(d)[:240])
"

# Owners BEFORE membership requests. Discourse refuses outright -- "You cannot
# allow membership requests for a group without any owners" -- because a
# request with no owner to notify would sit in a queue nobody can see. The
# route is /groups/<id>/owners.json; the /admin/ one 404s.
echo "  setting owners (a membership request notifies these people)"
api PUT "/groups/$GROUP_ID/owners.json" \
    "{\"usernames\":\"$(IFS=,; echo "${MEMBERS[*]}")\"}" | jq_py "
import json,sys
d=json.load(sys.stdin)
if d.get('success'): print('  owners:', ', '.join(d.get('usernames') or []))
else: print('  owners FAILED:', json.dumps(d)[:200])
"

# Now the settings that need an owner to exist.
echo "  applying group settings"
api PUT "/groups/$GROUP_ID.json" "$GROUP_BODY" | jq_py "
import json,sys
d=json.load(sys.stdin)
if d.get('success'): print('  membership requests enabled')
else: print('  SETTINGS FAILED:', json.dumps(d)[:240]); raise SystemExit(1)
" || exit 1

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
