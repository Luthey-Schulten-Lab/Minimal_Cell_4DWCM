"""Static first-order knockout impact over the spec.

A lost species breaks every 'all' requirement group it belongs to and, if it was the last member, every 'any' group.
Processes then propagate through depends_on until nothing changes. This is structure, not dynamics: it says which
reactions lose their enzyme and which processes lose a part, not how the cell responds. The model browser runs a port of
this function (gui/app.js: evaluate) and checks it against the per-gene results stored in the spec.
"""

RANK = {'ok': 0, 'reduced': 1, 'blocked': 2}


def gene_species(spec, locus):
    g = spec['genes'][locus]
    out = ['R_' + g['num']]
    if g['type'] == 'protein':
        out.append('P_' + g['num'])
    return out


def _groups(requires, lost, reduced):
    status, why = 'ok', []
    for grp in requires:
        members = grp['species']
        gone = [s for s in members if s in lost]
        weak = [s for s in members if s in reduced and s not in lost]
        if grp['rule'] == 'all' and gone:
            st = grp.get('on_loss', 'blocked')
            why.append('needs all of %s; lost %s' % (', '.join(members), ', '.join(gone)))
        elif grp['rule'] == 'any' and gone and len(gone) == len(members):
            st = grp.get('on_loss', 'blocked')
            why.append('needs one of %s; all lost' % ', '.join(members))
        elif gone or weak:
            st = 'reduced'
            why.append(('%s lost, %s remain' % (', '.join(gone), ', '.join(s for s in members if s not in lost))) if gone
                       else '%s reduced' % ', '.join(weak))
        else:
            continue
        if RANK[st] > RANK[status]:
            status = st
    return status, why


def evaluate(spec, lost, reduced=(), disabled=()):
    lost, reduced = set(lost), set(reduced)
    rx = {rid: {'status': 'blocked', 'why': ['disabled']} for rid in disabled}
    for rid, r in spec['reactions'].items():
        if not r['active'] or rid in rx:
            continue
        st, why = _groups(r['requires'], lost, reduced)
        if st != 'ok':
            rx[rid] = {'status': st, 'why': why}
    pr = {}
    for pid, p in spec['processes'].items():
        st, why = _groups(p['requires'], lost, reduced)
        if st != 'ok':
            pr[pid] = {'status': st, 'why': why}
    changed = True
    while changed:
        changed = False
        for pid, p in spec['processes'].items():
            for dep in p.get('depends_on', []):
                if 'process' in dep:
                    src = pr.get(dep['process'])
                    label = 'process ' + dep['process']
                else:
                    hit = [rid for rid, v in rx.items() if v['status'] == 'blocked' and spec['reactions'][rid]['kind'] == dep['reaction_kind']]
                    src = {'status': 'blocked'} if hit else None
                    label = 'reaction ' + ', '.join(sorted(hit))
                if not src:
                    continue
                st = dep['effect'] if src['status'] == 'blocked' else 'reduced'
                cur = pr.get(pid)
                reason = 'depends on %s (%s)' % (label, src['status'])
                if cur is None or RANK[st] > RANK[cur['status']]:
                    pr[pid] = {'status': st, 'why': (cur['why'] if cur else []) + [reason]}
                    changed = True
                elif reason not in cur['why']:
                    cur['why'].append(reason)
    mets = orphaned_metabolites(spec, rx)
    return {'reactions': rx, 'processes': pr, 'metabolites': mets}


def producers(spec, sid):
    out = []
    for rid, r in spec['reactions'].items():
        if not r['active']:
            continue
        if any(s == sid for s, _ in r['products']) and r.get('kcatF', 1) != 0:
            out.append(rid)
        elif any(s == sid for s, _ in r['substrates']) and r.get('reversible'):
            out.append(rid)
    return out


def orphaned_metabolites(spec, rx):
    """Metabolites that had a producing reaction and now have none left (medium species are fixed, so never orphaned)."""
    blocked = {rid for rid, v in rx.items() if v['status'] == 'blocked'}
    if not blocked:
        return {}
    out = {}
    for sid, sp in spec['species'].items():
        if sp['kind'] not in ('metabolite', 'protein_form'):
            continue
        prods = producers(spec, sid)
        if prods and all(p in blocked for p in prods):
            out[sid] = {'status': 'blocked', 'why': ['every producing reaction is blocked: ' + ', '.join(prods)]}
    return out


def single_knockouts(spec):
    """{locus: {'reactions': {id: status}, 'processes': {id: status}}} for every gene, stored in the spec for the browser's self-check."""
    out = {}
    for locus in spec['genes']:
        res = evaluate(spec, gene_species(spec, locus))
        out[locus] = {k: {i: v['status'] for i, v in res[k].items()} for k in ('reactions', 'processes', 'metabolites')}
    return out
