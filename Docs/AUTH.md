# Accounts, approval and sign-in

The hosted deployment used to be gated by one shared password in an nginx
`htpasswd` file. Everyone who could use it used the same credential, there was
no record of who ran what, and revoking one person meant changing the password
for everyone.

It now uses real accounts. Identity comes from the **Genohub Discourse
community** at `https://discussion.genohub.org`, which already has a signup
flow, email verification and staff moderation. Access to IGVF Agent requires
membership of a Discourse **group** that administrators control.

```
  browser ──► Cloudflare tunnel ──► nginx ──auth_request──► gate (Deploy/auth/gate.py)
                                      │                        │
                                      │                        └──► discussion.genohub.org
                                      └──► Streamlit app  (X-IGVF-User: …)
```

**Signing up and being approved are two separate steps, deliberately.** Anyone
may join the forum. Only members of the approval group may drive an agent that
executes analysis pipelines on a shared machine.

---

## Part 1 — one-time setup on Discourse

Needs a Discourse administrator account. Everything here is in the admin UI.

### 1. Create the approval group

**Admin → Groups → New Group**

| Field | Value |
|---|---|
| Name | `igvfagent` |
| Full name | IGVF Agent users |
| Who can join | Nobody (owners add members) — **not** "anyone" |
| Visibility | Whatever suits; membership visibility does not affect access |

Adding somebody to this group is the approval step. Removing them revokes
access at their next sign-in, and within `IGVF_SESSION_HOURS` for a browser
that is already signed in.

### 2. Enable the SSO provider

**Admin → Settings → Login**

| Setting | Value |
|---|---|
| `enable discourse connect provider` | ✅ on |
| `discourse connect provider secrets` | `igvfagent.genohub.org` \| `<the SSO secret>` |

The secrets setting is a two-column list: the left column is the **hostname**
of the site being logged into, the right column is the shared secret. Discourse
picks the secret by matching the hostname of the `return_sso_url` the gate
sends, so the left column must be exactly `igvfagent.genohub.org`.

Generate the secret with `openssl rand -hex 32`. It goes in two places and
nowhere else: that Discourse setting, and `IGVF_SSO_SECRET` in
`Deploy/.env.prod` on the host.

### 3. Check that the endpoint woke up

```bash
curl -s -o /dev/null -w '%{http_code}\n' https://discussion.genohub.org/session/sso_provider
```

`404` means the provider is still off. Anything else means the setting took.

---

## Part 2 — the host side

In `Deploy/.env.prod`:

```bash
IGVF_SSO_SECRET=<the same value pasted into Discourse>
IGVF_SESSION_SECRET=<a DIFFERENT openssl rand -hex 32>
IGVF_APPROVAL_GROUP=igvfagent
IGVF_SESSION_HOURS=168
# Break-glass, for the cutover. Remove once Discourse sign-in works.
IGVF_BOOTSTRAP_TOKEN=<openssl rand -hex 24>
```

The two secrets **must differ**. The gate refuses to start if they are equal:
one value for both would mean anyone able to mint an SSO payload could mint a
session cookie directly.

`WS_TOKEN` and `Deploy/nginx/htpasswd` are no longer used. Both have been
removed from the live deployment.

Then redeploy. `Deploy/redeploy.sh` refuses while an analysis is running, which
is the behaviour you want — wait for it rather than forcing.

---

### Or do all of it in one command

`Deploy/configure-discourse-sso.sh` performs every step in Parts 1 and 2 above
against the Discourse API, and `--check` reports the current state without
changing anything. Neither secret is echoed or pasted: the API key is read from
the gitignored `Docs/Secret/discourse-API.txt`, the SSO secret straight off the
deployment host. Safe to re-run — each step checks current state first.

```bash
bash Deploy/configure-discourse-sso.sh --check    # read-only
bash Deploy/configure-discourse-sso.sh            # apply
bash Deploy/configure-discourse-sso.sh --add alice bob
```

With no `--add`, it seeds the group with the forum's own administrators: an
empty approval group means nobody can sign in, and admins can grant themselves
access regardless.

## Part 3 — the cutover, in an order that cannot lock you out

Turning a live site over to a new identity provider in one step means that if
anything about the Discourse side is wrong — provider not enabled, secret
mistyped, group not created — **nobody can get in, including the person who
would fix it**. So:

1. **Deploy with `IGVF_BOOTSTRAP_TOKEN` set.** Visit
   `https://igvfagent.genohub.org/_auth/bootstrap?token=<that token>`. You get
   a session without Discourse. If the site loads, the gate, nginx and the app
   are wired correctly.
2. **Now try the real thing**: sign out, then load the site normally. You
   should be bounced to Discourse and back.
3. **When that works, remove `IGVF_BOOTSTRAP_TOKEN`** from `.env.prod` and
   redeploy. The bootstrap URL then returns 404.

**The live deployment has completed this and `IGVF_BOOTSTRAP_TOKEN` is now
unset**, so `/_auth/bootstrap` returns 404. To get back in if Discourse is ever
down or the secret is rotated out from under you, put a fresh token back:

```bash
bash Deploy/vm.sh ssh "cd /srv/igvfagent/Deploy && \
  echo IGVF_BOOTSTRAP_TOKEN=\$(openssl rand -hex 24) >> .env.prod && \
  docker compose --env-file .env.prod -f docker-compose.prod.yml up -d --force-recreate auth"
```

Then read it out of `.env.prod` and visit `/_auth/bootstrap?token=…`. It is
exactly as strong as the shared password it replaced, which is why it is not
left set.

---

## Day to day

**Approve someone**: Admin → Groups → `igvfagent` → add their username. Tell
them to reload. Nothing else is needed on their side.

**Revoke someone**: remove them from the group. They lose access at their next
sign-in; an existing browser session survives until its cookie expires. To cut
every session immediately, rotate `IGVF_SESSION_SECRET` and redeploy — that
invalidates all cookies at once, including your own.

**Someone says they are stuck**: they will have landed on a page naming the
group they need. That page appears only after they have successfully signed in
to Discourse, so if they never reach it the problem is the forum account, not
the approval.

---

## What each person can see

**Results are shared; organisation is private.**

- **Any past answer is readable by anyone signed in**, and is recalled
  automatically when someone asks about the same accession. This is deliberate:
  the expensive case is two people asking the same question a week apart, and
  there is little point paying for that analysis twice to hide a result from a
  colleague who could simply ask for it. Set `IGVF_HISTORY_SHARED=0` for a
  deployment where that trade goes the other way — strict per-user privacy,
  at the cost of re-running work.
- **Projects stay private.** A project is visible only to its owner and the
  people they share it with (`igvfagent project share <username>`, or the
  sidebar). Members can see and add to a shared project; only the owner can
  rename, archive, or change who else is in it.
- Each run still records **who** produced it, so shared does not mean
  anonymous.
- Work recorded before accounts existed belongs to nobody in particular — it
  was produced under the shared password — and stays readable.

`Scripts/test_history_visibility.py` asserts all of the above.

---

## Why it is built this way

**Why not a signup form in IGVF Agent?** It would mean a second credential
store — password hashing, email verification, reset flows, an approval queue —
on a public host, duplicating a forum that already does all of it. A homegrown
credential store is the thing that gets breached.

**Why not Cloudflare Access?** It is a good option and needs no code, but it
has no self-service signup: an administrator adds every email by hand, and user
management lives in the Cloudflare dashboard rather than in the community the
lab already runs.

**Why is auth in nginx rather than in Streamlit?** The app runs analysis CLIs
as subprocesses on prompts from the public internet. It is the thing being
protected; it must not also be the thing deciding who gets in. The app container
publishes no host port and is reachable only through the gateway, which is why
it can trust the `X-IGVF-*` headers it is handed.

**Why did the Safari workaround disappear?** Safari does not send HTTP Basic
Auth credentials on a WebSocket handshake, so the old config had to gate
`/_stcore/stream` on a separate shared token while the page used basic auth.
Safari *does* send cookies. With cookie sessions the page and the socket are
checked the same way, by the same subrequest, and the special case is gone.

---

## Testing

```bash
python3 Deploy/auth/test_gate.py              # the gate, over real HTTP
python3 Scripts/test_history_visibility.py    # who can see whose work
```

The first plays both browser and Discourse against a live gate: forged
signatures, replayed callbacks, spent nonces, tampered cookies, a user outside
the group, and the break-glass path. The second asserts the visibility rules
above, including that sharing a project grants visibility but not control.
