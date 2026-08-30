# can you rewite new skills of protein-protein interaction prediction? so that this new skills can be added to IGVFagent

## What was run
- `ext_list()`
- `ext_author_skill(name=ppi_predict, description=★ PROTEIN-PROTEIN INTERACTION PREDICTION ★. Predicts (not merely retrieves) PPIs for a seed protein/gene set, or scores specific candidate pairs. Re-scores the seven independent STRING evidence channels through a Bayesian noisy-OR integration, deliberately SEPARATING predictive channels (gene neighbourhood, fusion, phylogenetic co-occurrence, co-expression, text mining) from experimental/curated channels, then adds a topological shared-neighbour feature (shared partners, Jaccard, Adamic-Adar). Candidates scoring high on prediction but low on experimental evidence are flagged as NOVEL INTERACTION CANDIDATES. Modes: --proteins A,B,C for discovery ranking; --pair 'A:B,C:D' to score named pairs. Writes predictions.csv, novel_candidates.csv, summary.json, an SVG figure and a narrative report under Docs/PPIPredict/., tool_parameters={"type":"object","properties":{"proteins":{"type":"string","description":"Seed gene/protein symbols or UniProt/Ensembl ids, comma-separated (discovery mode), e.g. 'GCKR,KRTCAP3'"},"pair":{"type":"string","description":"Specific pairs to score, colon-separated within a pair and comma-separated between pairs, e.g. 'GCKR:GCK,GCKR:KRTCAP3'"},"species":{"type":"integer","description":"NCBI taxon id (default 9606 human)"},"min_score":{"type":"number","description":"Minimum STRING combined score for candidate retrieval, 0-1 (default 0.15 - deliberately low so prediction-only edges are kept)"},"top_k":{"type":"integer","description":"How many top candidates to report/tabulate (default 40)"},"max_partners":{"type":"integer","description":"Max STRING partners fetched per seed (default 150)"},"no_topology":{"type":"boolean","description":"Skip the shared-neighbour topology pass (faster, fewer API calls)"},"label":{"type":"string","description":"Short run label for the output directory"}}}, source='''ppi-predict: protein-protein interaction PREDICTION for IGVFagent.

This skill predicts interactions rather than only looking them up. It pulls the
seven independent STRING evidence channels for every candidate edge and
re-integrates them with a Bayesian noisy-OR, keeping two scores apart:

  prediction_score   noisy-OR over the PREDICTIVE channels only
                     (gene neighbourhood, gene fusion, phylogenetic
                     co-occurrence, co-expression, text mining)
  experimental_score noisy-OR over the EXPERIMENTAL / CURATED channels
                     (escore = physical assays, dscore = curated pathway DBs)

A third, orthogonal feature is topological: candidates that share many
(low-degree) STRING neighbours with the seed are more likely to be genuine
partners. Adamic-Adar over the shared-neighbour set is rank-normalised and
mixed with the prediction score.

A candidate with a high prediction_score and a near-zero experimental_score is
the interesting case: predicted, but not yet experimentally recorded. Those are
reported as novel-interaction candidates.

Modes
  discovery : --proteins GCKR,KRTCAP3     rank predicted partners of the seeds
  pair      : --pair GCKR:GCK,GCKR:APOA5  score named candidate pairs

Outputs a timestamped run directory under Docs/PPIPredict/ containing
predictions.csv, novel_candidates.csv (or pair_scores.csv), summary.json,
prediction_scores.svg and report.md.
'''
import argparse
import csv
import json
import math
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime

STRING_API = 'https://string-db.org/api'
CALLER = 'igvfagent-ppi-predict'
PRIOR = 0.041
PREDICTIVE_CHANNELS = ['nscore', 'fscore', 'pscore', 'ascore', 'tscore']
EVIDENCE_CHANNELS = ['escore', 'dscore']
CHANNEL_NAMES = {
    'nscore': 'gene neighbourhood',
    'fscore': 'gene fusion',
    'pscore': 'phylogenetic co-occurrence',
    'ascore': 'co-expression',
    'escore': 'experimental',
    'dscore': 'curated database',
    'tscore': 'text mining',
}
FINAL_W_PRED = 0.65
FINAL_W_TOPO = 0.35


def _get(path, params, retries=3, timeout=120):
    url = STRING_API + path + '?' + urllib.parse.urlencode(params)
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': CALLER})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                raw = fh.read().decode('utf-8', 'replace')
            if not raw.strip():
                return []
            return json.loads(raw)
        except Exception as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError('STRING request failed after retries: ' + str(last))


def noisy_or(scores):
    '''STRING-style prior-corrected noisy-OR combination of channel scores.'''
    prod = 1.0
    used = False
    for s in scores:
        if s is None:
            continue
        s = min(max(float(s), 0.0), 0.999999)
        adj = (s - PRIOR) / (1.0 - PRIOR)
        if adj <= 0:
            continue
        used = True
        prod *= (1.0 - adj)
    if not used:
        return 0.0
    comb = 1.0 - prod
    return comb * (1.0 - PRIOR) + PRIOR


def _f(row, key):
    try:
        v = row.get(key)
        return float(v) if v not in (None, '') else 0.0
    except Exception:
        return 0.0


def score_row(row):
    pred = noisy_or([_f(row, c) for c in PREDICTIVE_CHANNELS])
    evid = noisy_or([_f(row, c) for c in EVIDENCE_CHANNELS])
    return pred, evid


def channel_summary(row, cut=0.15):
    bits = []
    for c in PREDICTIVE_CHANNELS + EVIDENCE_CHANNELS:
        v = _f(row, c)
        if v >= cut:
            bits.append(CHANNEL_NAMES[c] + ' ' + str(round(v, 3)))
    return '; '.join(bits) if bits else 'no channel above ' + str(cut)


def classify(pred, evid, final):
    if evid >= 0.4:
        return 'known (experimental/curated support)'
    if final >= 0.5 and evid < 0.15:
        return 'novel candidate'
    if final >= 0.3:
        return 'weak candidate'
    return 'low confidence'


def map_ids(proteins, species):
    return _get('/json/get_string_ids', {
        'identifiers': chr(13).join(proteins),
        'species': species,
        'limit': 1,
        'echo_query': 1,
        'caller_identity': CALLER,
    })


def get_partners(string_ids, species, min_score, limit):
    params = {
        'identifiers': chr(13).join(string_ids),
        'species': species,
        'required_score': int(round(float(min_score) * 1000)),
        'caller_identity': CALLER,
    }
    if limit:
        params['limit'] = int(limit)
    return _get('/json/interaction_partners', params)


def get_network(string_ids, species, min_score=0.0):
    return _get('/json/network', {
        'identifiers': chr(13).join(string_ids),
        'species': species,
        'required_score': int(round(float(min_score) * 1000)),
        'caller_identity': CALLER,
    })


def neighbour_sets(string_ids, species, min_score, limit):
    sets = {}
    ids = [i for i in string_ids if i]
    batch = 40
    for i in range(0, len(ids), batch):
        rows = get_partners(ids[i:i + batch], species, min_score, limit)
        for r in rows:
            a = r.get('preferredName_A')
            b = r.get('preferredName_B')
            if not a or not b:
                continue
            sets.setdefault(a, set()).add(b)
            sets.setdefault(b, set()).add(a)
    degrees = dict((k, len(v)) for k, v in sets.items())
    return sets, degrees


def topo_features(cand, others, sets, degrees, queried):
    best = {'shared': 0, 'jaccard': 0.0, 'aa': 0.0, 'via': ''}
    cs = sets.get(cand) or set()
    if not cs:
        return best
    for s in others:
        ss = sets.get(s) or set()
        if not ss:
            continue
        shared = cs & ss
        union = cs | ss
        jac = len(shared) / float(len(union)) if union else 0.0
        aa = 0.0
        for n in shared:
            d = degrees.get(n, 0) if n in queried else 0
            if d < 2:
                d = 25
            aa += 1.0 / math.log(float(d))
        if aa > best['aa'] or (aa == best['aa'] and jac > best['jaccard']):
            best = {'shared': len(shared), 'jaccard': round(jac, 4),
                    'aa': round(aa, 4), 'via': s}
    return best


def _pct(values):
    n = len(values)
    out = [0.0] * n
    if n == 0:
        return out
    order = sorted(range(n), key=lambda i: values[i])
    for rank, idx in enumerate(order):
        out[idx] = 0.0 if n < 2 else rank / float(n - 1)
    return out


def _ts():
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def _run_dir(label):
    base = os.getcwd()
    if not os.path.isdir(os.path.join(base, 'Docs')) and os.path.isdir('/workspace/Docs'):
        base = '/workspace'
    safe = ''.join(ch if (ch.isalnum() or ch in '_-') else '_' for ch in (label or 'ppi_predict'))
    d = os.path.join(base, 'Docs', 'PPIPredict', _ts() + '_' + safe[:60])
    os.makedirs(d, exist_ok=True)
    return d


def _write_csv(path, rows, cols):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(dict((c, r.get(c, '')) for c in cols))


def _md_table(cols, rows):
    lines = ['| ' + ' | '.join(cols) + ' |',
             '|' + '|'.join(['---'] * len(cols)) + '|']
    for r in rows:
        lines.append('| ' + ' | '.join(str(r.get(c, '')) for c in cols) + ' |')
    return chr(10).join(lines)


def _svg_bars(path, items, title):
    row_h = 22
    top = 64
    width = 940
    left = 210
    bar_max = width - left - 130
    height = top + row_h * max(1, len(items)) + 46
    p = []
    p.append("<svg xmlns='http://www.w3.org/2000/svg' width='" + str(width) + "' height='" + str(height) + "'>")
    p.append("<rect width='100%' height='100%' fill='white'/>")
    p.append("<text x='16' y='28' font-family='Helvetica,Arial' font-size='17' font-weight='bold'>" + title + "</text>")
    p.append("<text x='16' y='48' font-family='Helvetica,Arial' font-size='12' fill='#555'>bar = final predicted-interaction score; blue = no experimental/curated support (novel candidate), grey = already supported</text>")
    for i, it in enumerate(items):
        y = top + i * row_h
        name, score, supported, pred, evid = it
        colour = '#8a8a8a' if supported else '#1f6fb4'
        w = max(2.0, float(score) * bar_max)
        p.append("<text x='" + str(left - 8) + "' y='" + str(y + 14) + "' text-anchor='end' font-family='Helvetica,Arial' font-size='12'>" + str(name) + "</text>")
        p.append("<rect x='" + str(left) + "' y='" + str(y + 3) + "' width='" + str(round(w, 1)) + "' height='14' fill='" + colour + "'/>")
        lab = str(round(float(score), 3)) + '  (pred ' + str(round(float(pred), 2)) + ' / exp ' + str(round(float(evid), 2)) + ')'
        p.append("<text x='" + str(left + w + 6) + "' y='" + str(y + 15) + "' font-family='Helvetica,Arial' font-size='11' fill='#333'>" + lab + "</text>")
    p.append("<line x1='" + str(left) + "' y1='" + str(top - 6) + "' x2='" + str(left) + "' y2='" + str(height - 34) + "' stroke='#999'/>")
    p.append("<text x='" + str(left) + "' y='" + str(height - 14) + "' font-family='Helvetica,Arial' font-size='11' fill='#555'>0.0</text>")
    p.append("<text x='" + str(left + bar_max - 14) + "' y='" + str(height - 14) + "' font-family='Helvetica,Arial' font-size='11' fill='#555'>1.0</text>")
    p.append('</svg>')
    with open(path, 'w') as fh:
        fh.write(chr(10).join(p))


def _split(text):
    if not text:
        return []
    t = text
    for ch in [';', chr(10), chr(13), chr(9), ' ']:
        t = t.replace(ch, ',')
    return [x.strip() for x in t.split(',') if x.strip()]


def _parse_pairs(text):
    pairs = []
    for tok in _split(text):
        parts = None
        for sep in [':', '=', '--', '|']:
            if sep in tok:
                parts = tok.split(sep)
                break
        if parts is None and tok.count('-') == 1:
            parts = tok.split('-')
        if parts is None or len(parts) < 2:
            continue
        a = parts[0].strip()
        b = parts[1].strip()
        if a and b:
            pairs.append((a, b))
    return pairs


METHODS = ('Method. Candidate edges are pulled from STRING with a deliberately low '
           'confidence floor so that prediction-only edges survive. For every edge the '
           'seven independent evidence channels are re-integrated with a prior-corrected '
           'noisy-OR (prior p = 0.041), split into a prediction_score over the predictive '
           'channels (gene neighbourhood, gene fusion, phylogenetic co-occurrence, '
           'co-expression, text mining) and an experimental_score over the experimental and '
           'curated-database channels. A topological feature is computed from the '
           'shared-neighbour set (shared count, Jaccard, Adamic-Adar with 1/log(degree) '
           'weighting); Adamic-Adar is rank-normalised across candidates and mixed as '
           'final_score = 0.65 * prediction_score + 0.35 * topology_percentile. A candidate '
           'is called a novel candidate when final_score >= 0.5 and experimental_score < 0.15.')


def run_discovery(seeds, args):
    out = _run_dir(args.label or ('discover_' + '_'.join(seeds[:3])))
    mapped = map_ids(seeds, args.species)
    resolved = []
    seen = set()
    for m in mapped:
        sid = m.get('stringId')
        sym = m.get('preferredName')
        if not sid or not sym or sid in seen:
            continue
        seen.add(sid)
        resolved.append({'query': m.get('queryItem') or sym, 'string_id': sid,
                         'symbol': sym, 'annotation': (m.get('annotation') or '')[:300]})
    if not resolved:
        raise RuntimeError('no requested protein resolved in STRING for species ' + str(args.species))
    seed_syms = [r['symbol'] for r in resolved]
    seed_ids = [r['string_id'] for r in resolved]

    rows = get_partners(seed_ids, args.species, args.min_score, args.max_partners)
    cand = {}
    for r in rows:
        a = r.get('preferredName_A')
        b = r.get('preferredName_B')
        if not a or not b:
            continue
        if a in seed_syms and b in seed_syms:
            seed, partner, ss_edge = a, b, True
        elif a in seed_syms:
            seed, partner, ss_edge = a, b, False
        elif b in seed_syms:
            seed, partner, ss_edge = b, a, False
        else:
            continue
        pred, evid = score_row(r)
        sid = r.get('stringId_B') if partner == b else r.get('stringId_A')
        entry = {
            'candidate': partner,
            'seed': seed,
            'seed_seed_edge': ss_edge,
            'string_id': sid,
            'string_combined_score': round(_f(r, 'score'), 4),
            'prediction_score': round(pred, 4),
            'experimental_score': round(evid, 4),
            'channels': channel_summary(r),
            'neighbourhood': _f(r, 'nscore'),
            'fusion': _f(r, 'fscore'),
            'cooccurrence': _f(r, 'pscore'),
            'coexpression': _f(r, 'ascore'),
            'textmining': _f(r, 'tscore'),
            'experiments': _f(r, 'escore'),
            'database': _f(r, 'dscore'),
        }
        prev = cand.get(partner)
        if prev is None or entry['prediction_score'] > prev['prediction_score']:
            cand[partner] = entry

    items = sorted(cand.values(), key=lambda d: -d['prediction_score'])
    topo_pool = items[:max(int(args.top_k), 25)]
    topo_done = False
    if not args.no_topology and topo_pool:
        ids = list(seed_ids) + [d['string_id'] for d in topo_pool if d.get('string_id')]
        sets, degrees = neighbour_sets(ids, args.species, 0.4, 100)
        queried = set(seed_syms) | set(d['candidate'] for d in topo_pool)
        for d in topo_pool:
            t = topo_features(d['candidate'], seed_syms, sets, degrees, queried)
            d['shared_neighbours'] = t['shared']
            d['jaccard'] = t['jaccard']
            d['adamic_adar'] = t['aa']
            d['topology_via_seed'] = t['via']
        pcts = _pct([d.get('adamic_adar', 0.0) for d in topo_pool])
        for d, p in zip(topo_pool, pcts):
            d['topology_percentile'] = round(p, 4)
            d['final_score'] = round(FINAL_W_PRED * d['prediction_score'] + FINAL_W_TOPO * p, 4)
        topo_done = True

    for d in items:
        if 'final_score' not in d:
            d['final_score'] = d['prediction_score']
            d['topology_percentile'] = ''
            d['shared_neighbours'] = ''
            d['jaccard'] = ''
            d['adamic_adar'] = ''
            d['topology_via_seed'] = ''
        d['call'] = classify(d['prediction_score'], d['experimental_score'], d['final_score'])
    items.sort(key=lambda d: -float(d['final_score']))

    cols = ['candidate', 'seed', 'final_score', 'prediction_score', 'experimental_score',
            'string_combined_score', 'call', 'shared_neighbours', 'jaccard', 'adamic_adar',
            'topology_percentile', 'topology_via_seed', 'channels', 'neighbourhood', 'fusion',
            'cooccurrence', 'coexpression', 'textmining', 'experiments', 'database',
            'seed_seed_edge', 'string_id']
    pred_csv = os.path.join(out, 'predictions.csv')
    _write_csv(pred_csv, items, cols)
    novel = [d for d in items if d['call'] == 'novel candidate']
    novel_csv = os.path.join(out, 'novel_candidates.csv')
    _write_csv(novel_csv, novel, cols)

    svg = os.path.join(out, 'prediction_scores.svg')
    bars = [(d['candidate'], d['final_score'], d['experimental_score'] >= 0.4,
             d['prediction_score'], d['experimental_score']) for d in items[:20]]
    _svg_bars(svg, bars, 'Predicted interaction partners of ' + ', '.join(seed_syms[:4]))

    edge_csv = os.path.join(out, 'edges.csv')
    _write_csv(edge_csv, [{'source': d['seed'], 'target': d['candidate'],
                           'score': d['final_score'], 'call': d['call']} for d in items],
               ['source', 'target', 'score', 'call'])

    topn = items[:int(args.top_k)]
    summary = {
        'mode': 'discovery',
        'species': args.species,
        'seeds': resolved,
        'min_score': args.min_score,
        'max_partners': args.max_partners,
        'topology_pass': topo_done,
        'weights': {'prediction': FINAL_W_PRED, 'topology': FINAL_W_TOPO},
        'n_candidates': len(items),
        'n_novel_candidates': len(novel),
        'n_known_supported': len([d for d in items if d['call'].startswith('known')]),
        'top': [{'candidate': d['candidate'], 'final_score': d['final_score'],
                 'prediction_score': d['prediction_score'],
                 'experimental_score': d['experimental_score'], 'call': d['call']}
                for d in topn],
        'artifacts': {'predictions_csv': pred_csv, 'novel_csv': novel_csv,
                      'edges_csv': edge_csv, 'figure_svg': svg},
    }
    sum_json = os.path.join(out, 'summary.json')
    with open(sum_json, 'w') as fh:
        json.dump(summary, fh, indent=2)

    L = []
    L.append('# PPI prediction - ' + ', '.join(seed_syms))
    L.append('')
    L.append('Species taxon: ' + str(args.species) + '  |  candidates screened: ' + str(len(items)) +
             '  |  novel candidates: ' + str(len(novel)) +
             '  |  already experimentally supported: ' + str(summary['n_known_supported']))
    L.append('')
    L.append('## Seeds resolved')
    L.append('')
    L.append(_md_table(['query', 'symbol', 'string_id'],
                       [{'query': r['query'], 'symbol': r['symbol'], 'string_id': r['string_id']} for r in resolved]))
    L.append('')
    L.append('## Top predicted partners')
    L.append('')
    L.append(_md_table(['candidate', 'seed', 'final_score', 'prediction_score', 'experimental_score',
                       'shared_neighbours', 'call'], topn))
    L.append('')
    L.append('## Novel-interaction candidates (predicted, no experimental/curated support)')
    L.append('')
    if novel:
        L.append(_md_table(['candidate', 'seed', 'final_score', 'prediction_score', 'channels'],
                           novel[:int(args.top_k)]))
    else:
        L.append('None passed the novel-candidate rule (final_score >= 0.5 and experimental_score < 0.15).')
    L.append('')
    L.append('## Method')
    L.append('')
    L.append(METHODS)
    L.append('')
    L.append('## Caveats')
    L.append('')
    L.append('- Text mining is a literature-co-mention channel, not evidence of physical contact; a '
             'candidate carried mostly by textmining is a hypothesis about association, not a binding claim.')
    L.append('- STRING edges are functional associations. A high final_score means "likely functionally '
             'coupled and plausibly interacting", not "direct physical complex".')
    L.append('- Node degrees used by Adamic-Adar are exact only for proteins that were themselves queried; '
             'for other shared neighbours a constant degree of 25 is substituted.')
    L.append('- Only the top ' + str(len(topo_pool)) + ' candidates by prediction_score received the '
             'topology pass; the rest carry final_score = prediction_score.')
    L.append('')
    L.append('## Artefacts')
    L.append('')
    for k, v in summary['artifacts'].items():
        L.append('- ' + k + ': `' + v + '`')
    L.append('- summary_json: `' + sum_json + '`')
    report = os.path.join(out, 'report.md')
    with open(report, 'w') as fh:
        fh.write(chr(10).join(L) + chr(10))

    print('PPI prediction (discovery mode)')
    print('Seeds: ' + ', '.join(seed_syms))
    print('Candidates screened: ' + str(len(items)) + '   novel candidates: ' + str(len(novel)))
    for d in topn[:15]:
        print('  ' + str(d['candidate']).ljust(14) + ' final=' + str(d['final_score']) +
              ' pred=' + str(d['prediction_score']) + ' exp=' + str(d['experimental_score']) +
              '  ' + d['call'])
    print('Predictions: ' + pred_csv)
    print('Novel candidates: ' + novel_csv)
    print('Figure: ' + svg)
    print('Report: ' + report)
    return 0


def run_pairs(pairs, args):
    out = _run_dir(args.label or 'pairs')
    flat = sorted(set([p for pair in pairs for p in pair]))
    mapped = map_ids(flat, args.species)
    lookup = {}
    for m in mapped:
        q = m.get('queryItem')
        if q:
            lookup[q.upper()] = m
        if m.get('preferredName'):
            lookup.setdefault(m['preferredName'].upper(), m)
    ids = [m.get('stringId') for m in mapped if m.get('stringId')]
    net = get_network(ids, args.species, 0.0) if ids else []
    index = {}
    for r in net:
        a = (r.get('preferredName_A') or '').upper()
        b = (r.get('preferredName_B') or '').upper()
        if a and b:
            index[(a, b)] = r
            index[(b, a)] = r
    sets, degrees = ({}, {})
    if not args.no_topology and ids:
        sets, degrees = neighbour_sets(ids, args.species, 0.4, 100)
    queried = set(sets.keys())

    results = []
    for a, b in pairs:
        ma = lookup.get(a.upper())
        mb = lookup.get(b.upper())
        rec = {'protein_a': a, 'protein_b': b,
               'resolved_a': (ma or {}).get('preferredName', ''),
               'resolved_b': (mb or {}).get('preferredName', '')}
        if not ma or not mb:
            rec['status'] = 'unresolved in STRING'
            rec['call'] = 'not scored'
            results.append(rec)
            continue
        row = index.get((rec['resolved_a'].upper(), rec['resolved_b'].upper()))
        if row:
            pred, evid = score_row(row)
            rec['string_combined_score'] = round(_f(row, 'score'), 4)
            rec['prediction_score'] = round(pred, 4)
            rec['experimental_score'] = round(evid, 4)
            rec['channels'] = channel_summary(row)
        else:
            rec['string_combined_score'] = 0.0
            rec['prediction_score'] = 0.0
            rec['experimental_score'] = 0.0
            rec['channels'] = 'no STRING edge at any confidence'
        t = topo_features(rec['resolved_a'], [rec['resolved_b']], sets, degrees, queried)
        rec['shared_neighbours'] = t['shared']
        rec['jaccard'] = t['jaccard']
        rec['adamic_adar'] = t['aa']
        topo_n = min(1.0, t['aa'] / 5.0)
        rec['topology_normalised'] = round(topo_n, 4)
        rec['final_score'] = round(FINAL_W_PRED * rec['prediction_score'] + FINAL_W_TOPO * topo_n, 4)
        rec['call'] = classify(rec['prediction_score'], rec['experimental_score'], rec['final_score'])
        rec['status'] = 'scored'
        results.append(rec)

    results.sort(key=lambda d: -float(d.get('final_score') or 0.0))
    cols = ['protein_a', 'protein_b', 'resolved_a', 'resolved_b', 'status', 'final_score',
            'prediction_score', 'experimental_score', 'string_combined_score', 'call',
            'shared_neighbours', 'jaccard', 'adamic_adar', 'topology_normalised', 'channels']
    pair_csv = os.path.join(out, 'pair_scores.csv')
    _write_csv(pair_csv, results, cols)

    svg = os.path.join(out, 'prediction_scores.svg')
    bars = [(str(d['protein_a']) + ' - ' + str(d['protein_b']), d.get('final_score') or 0.0,
             float(d.get('experimental_score') or 0.0) >= 0.4,
             d.get('prediction_score') or 0.0, d.get('experimental_score') or 0.0)
            for d in results if d['status'] == 'scored']
    if bars:
        _svg_bars(svg, bars[:20], 'Scored candidate protein pairs')

    summary = {'mode': 'pair', 'species': args.species, 'n_pairs': len(results),
               'topology_pass': bool(sets), 'weights': {'prediction': FINAL_W_PRED, 'topology': FINAL_W_TOPO},
               'pairs': results, 'artifacts': {'pair_csv': pair_csv, 'figure_svg': svg if bars else None}}
    sum_json = os.path.join(out, 'summary.json')
    with open(sum_json, 'w') as fh:
        json.dump(summary, fh, indent=2)

    L = []
    L.append('# PPI prediction - named pairs')
    L.append('')
    L.append('Species taxon: ' + str(args.species) + '  |  pairs scored: ' +
             str(len([d for d in results if d['status'] == 'scored'])) + ' of ' + str(len(results)))
    L.append('')
    L.append(_md_table(['protein_a', 'protein_b', 'final_score', 'prediction_score',
                       'experimental_score', 'shared_neighbours', 'call', 'channels'], results))
    L.append('')
    L.append('## Method')
    L.append('')
    L.append(METHODS)
    L.append('')
    L.append('In pair mode the topology feature cannot be rank-normalised against a candidate pool, so '
             'Adamic-Adar is scaled as min(1, AA/5) instead of by percentile - pair-mode and '
             'discovery-mode final_score values are therefore not directly comparable.')
    L.append('')
    L.append('## Artefacts')
    L.append('')
    L.append('- pair_csv: `' + pair_csv + '`')
    if bars:
        L.append('- figure_svg: `' + svg + '`')
    L.append('- summary_json: `' + sum_json + '`')
    report = os.path.join(out, 'report.md')
    with open(report, 'w') as fh:
        fh.write(chr(10).join(L) + chr(10))

    print('PPI prediction (pair mode)')
    for d in results:
        print('  ' + str(d['protein_a']) + ' - ' + str(d['protein_b']) + ': final=' +
              str(d.get('final_score', 'NA')) + ' pred=' + str(d.get('prediction_score', 'NA')) +
              ' exp=' + str(d.get('experimental_score', 'NA')) + '  ' + str(d.get('call')))
    print('Pair scores: ' + pair_csv)
    print('Report: ' + report)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog='igvfagent ppi-predict',
                                description='Predict protein-protein interactions from integrated evidence channels.')
    ap.add_argument('--proteins', '--protein', '--genes', '--gene', dest='proteins', default=None)
    ap.add_argument('--pair', '--pairs', dest='pair', default=None)
    ap.add_argument('--species', dest='species', type=int, default=9606)
    ap.add_argument('--min-score', '--min_score', dest='min_score', type=float, default=0.15)
    ap.add_argument('--top-k', '--top_k', dest='top_k', type=int, default=40)
    ap.add_argument('--max-partners', '--max_partners', dest='max_partners', type=int, default=150)
    ap.add_argument('--no-topology', '--no_topology', dest='no_topology', action='store_true')
    ap.add_argument('--label', dest='label', default=None)
    ap.add_argument('rest', nargs='*', default=None)
    args, unknown = ap.parse_known_args(sys.argv[1:] if argv is None else argv)

    if args.pair:
        pairs = _parse_pairs(args.pair)
        if not pairs:
            print('Could not parse any pair from: ' + str(args.pair))
            print("Use colon-separated pairs, e.g. --pair 'GCKR:GCK,GCKR:APOA5'")
            return 2
        return run_pairs(pairs, args)

    seeds = _split(args.proteins) if args.proteins else []
    if not seeds and args.rest:
        for tok in args.rest:
            seeds.extend(_split(tok))
    if not seeds:
        print('Nothing to do. Pass --proteins GCKR,KRTCAP3 (discovery) or --pair GCKR:GCK (pair scoring).')
        return 2
    return run_discovery(seeds, args)


if __name__ == '__main__':
    sys.exit(main())
)`
- `ext_validate(name=ppi_predict)`
- `ppi_predict(pair=GCKR:GCK,GCKR:KRTCAP3, label=smoke_gckr_pairs)`

## Artefacts produced
- `/workspace/Data/UserExtensions/skills/ppi_predict.py`

## Summary

_(no narrative generated)_

- `/workspace/Data/UserExtensions/skills/ppi_predict.py`

_backend `claude_cli`  ·  model `claude-fable-5`  ·  4 iter · 3 tool calls · stop `complete`  ·  report `/workspace/Docs/Agent/20260830_170603_can_you_rewite_new_skills_of_protein_pro_a4f2513f/report.md`_
