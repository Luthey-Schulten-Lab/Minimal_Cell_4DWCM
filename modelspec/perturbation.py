"""The perturbation file: what the model browser writes and what the simulator will apply.

    schema: 1
    name: ko_JCVISYN3A_0415
    knockouts:
      - gene: JCVISYN3A_0415
        mode: full                 # full | expression_only | initial_only
    knockdowns:
      - gene: JCVISYN3A_0001
        promoter_scale: 0.25
    initial_protein_counts: {P_0609: 5}
    initial_mrna_means:     {R_0001: 0.0}
    initial_metabolites_mM: {M_atp_c: 1.0}
    medium_mM:              {M_glc__D_e: 0.0}
    reaction_parameters:
      - {reaction: PGI, parameter: kcatF, scale: 0.1}      # or value:
    disabled_reactions: [ATPase]
    rdme_rate_constants:                                  # named jLM rate constants (ribo_bind, deg_bind_rate, secY_on, repOn, ...)
      - {constant: ribo_bind, scale: 0.5}                 # or value:
    gene_rate_scales:                                     # per-gene computed rates; gene "*" = every gene
      - {gene: JCVISYN3A_0001, rate: translation, scale: 0.5}   # translation | rnap_binding | mrna_degradation | secy_insertion | transcription
    diffusion_scales:
      - {constant: diffPtn, scale: 0.5}                   # named diffusion constants, or rna_diff for every RNA / ribosome intermediate
    gip_constants:
      - {constant: riboKcat, value: 10}                   # utility/GIP_rates.py module constants used by the rate formulas

validate() returns (errors, warnings); resolve() turns a valid file into the explicit list of edits, each with the value
it replaces and the simulator function where it takes effect, plus the static impact (modelspec.impact).
"""

import json
import re

from . import SCHEMA_VERSION
from .impact import evaluate, gene_species

KNOCKOUT_MODES = {
    'full': 'no transcription and no protein or RNA at t=0',
    'expression_only': 'no transcription; the t=0 protein and mRNA remain and are only diluted by growth and division '
                       '(the model has no protein degradation)',
    'initial_only': 'protein and mRNA start at 0 but the gene is still transcribed, so the protein is made again',
}
LIST_KEYS = ['knockouts', 'knockdowns', 'reaction_parameters', 'disabled_reactions', 'rdme_rate_constants', 'gene_rate_scales',
             'diffusion_scales', 'gip_constants']
GENE_RATES = {'translation': ('<n>_translat', 'Rxns_RDME.translation (GIP_rates.TranslationRate)'),
              'rnap_binding': ('<n>RNAP_on', 'Rxns_RDME.transcription / transcriptionLong (GIP_rates.RNAP_binding)'),
              'mrna_degradation': ('<n>_RNAdeg', 'Rxns_RDME.degradation_mrna (GIP_rates.mrnaDegradationRate)'),
              'secy_insertion': ('<n>_insertion', 'Rxns_RDME.translocation_secy (GIP_rates.TranslocationRate)'),
              'transcription': ('transcription', 'Rxns_CME.transcription / transcriptionLong (GIP_rates.TranscriptionRate)')}
GIP_EDITABLE = ['rnaPolKcat', 'rnaPolKd', 'riboKcat', 'riboKd', 'ctRNAconc', 'ATPconc', 'UTPconc', 'CTPconc', 'GTPconc']
MAP_KEYS = ['initial_protein_counts', 'initial_mrna_means', 'initial_metabolites_mM', 'medium_mM']
TOP_KEYS = ['schema', 'name', 'description', 'created', 'model'] + LIST_KEYS + MAP_KEYS
NAME = re.compile(r'^[A-Za-z0-9_.\-]{1,80}$')

WHERE = {
    'promoter': 'ImportInitialConditions.initializePromoterStrengths (sim_properties["promoters"])',
    'protein': 'ImportInitialConditions.initializeProteins',
    'mrna': 'ImportInitialConditions.initializeMRNA',
    'trna': 'ImportInitialConditions.initializeTRNA',
    'metabolite': 'ImportInitialConditions.initializeMetabolites',
    'medium': 'ImportInitialConditions.initializeMedium',
    'param_ode': 'Rxns_ODE._random_binding_static / defineNonRandomBindingRxns',
    'param_cme': 'Rxns_CME.tRNAcharging',
    'onoff': 'Rxns_ODE.defineRandomBindingRxns (onoff parameter)',
    'rdme_rc': 'Rxns_RDME / ImportInitialConditions (sim.rateConst)',
    'diffusion': 'Diffusion.generalDiffusionConstants / rnaDiffusion / ribosomeDiffusion (sim.diffusionConst)',
    'gip': 'utility/GIP_rates.py module constant',
}


def empty(name='perturbation'):
    return {'schema': SCHEMA_VERSION, 'name': name, 'description': '', 'knockouts': [], 'knockdowns': [],
            'initial_protein_counts': {}, 'initial_mrna_means': {}, 'initial_metabolites_mM': {}, 'medium_mM': {},
            'reaction_parameters': [], 'disabled_reactions': [], 'rdme_rate_constants': [], 'gene_rate_scales': [],
            'diffusion_scales': [], 'gip_constants': []}


def _num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def param_keys(r):
    return [p['key'] for p in r['params'] if p.get('key')]


def validate(spec, pert):
    err, warn = [], []
    if not isinstance(pert, dict):
        return ['the file is not a mapping'], []
    for k in pert:
        if k not in TOP_KEYS:
            err.append('unknown key %r (allowed: %s)' % (k, ', '.join(TOP_KEYS)))
    if pert.get('schema') != SCHEMA_VERSION:
        err.append('schema must be %d' % SCHEMA_VERSION)
    if not NAME.match(str(pert.get('name', ''))):
        err.append('name must be 1-80 characters of letters, digits, _ . -')
    fp = (pert.get('model') or {}).get('fingerprint')
    if fp and fp != spec['meta']['fingerprint']:
        warn.append('built against inputs %s; the current inputs are %s (input_data/ changed since)' % (fp, spec['meta']['fingerprint']))
    for k in LIST_KEYS:
        if not isinstance(pert.get(k, []) or [], list):
            err.append('%s must be a list' % k)
    for k in MAP_KEYS:
        if not isinstance(pert.get(k, {}) or {}, dict):
            err.append('%s must be a mapping' % k)
    if err:
        return err, warn
    genes, species, rxns = spec['genes'], spec['species'], spec['reactions']

    ko_genes, full_ko = set(), set()
    for i, ko in enumerate(pert.get('knockouts') or []):
        g, mode = ko.get('gene'), ko.get('mode', 'full')
        if g not in genes:
            err.append('knockouts[%d]: unknown gene %r' % (i, g))
            continue
        if mode not in KNOCKOUT_MODES:
            err.append('knockouts[%d]: mode must be one of %s' % (i, ', '.join(KNOCKOUT_MODES)))
        if g in ko_genes:
            err.append('knockouts[%d]: %s listed twice' % (i, g))
        ko_genes.add(g)
        if mode == 'full':
            full_ko.add(g)
        extra = set(ko) - {'gene', 'mode'}
        if extra:
            err.append('knockouts[%d]: unknown keys %s' % (i, ', '.join(sorted(extra))))
    for i, kd in enumerate(pert.get('knockdowns') or []):
        g, s = kd.get('gene'), kd.get('promoter_scale')
        if g not in genes:
            err.append('knockdowns[%d]: unknown gene %r' % (i, g))
        elif g in ko_genes:
            err.append('knockdowns[%d]: %s is also knocked out' % (i, g))
        if not _num(s) or s < 0:
            err.append('knockdowns[%d]: promoter_scale must be a number >= 0' % i)
        elif s > 1:
            warn.append('knockdowns[%d]: promoter_scale %g > 1 is an overexpression' % (i, s))

    for sid, v in (pert.get('initial_protein_counts') or {}).items():
        sp = species.get(sid)
        if not sp or sp['kind'] != 'protein':
            err.append('initial_protein_counts: %r is not a protein id (P_nnnn)' % sid)
        elif not isinstance(v, int) or isinstance(v, bool) or v < 0:
            err.append('initial_protein_counts[%s]: must be an integer >= 0' % sid)
        elif sp['gene'] in full_ko:
            err.append('initial_protein_counts[%s]: gene %s has a full knockout' % (sid, sp['gene']))
    for sid, v in (pert.get('initial_mrna_means') or {}).items():
        sp = species.get(sid)
        if not sp or sp['kind'] != 'mRNA':
            err.append('initial_mrna_means: %r is not an mRNA id (R_nnnn of a protein gene)' % sid)
        elif not _num(v) or v < 0:
            err.append('initial_mrna_means[%s]: must be a number >= 0' % sid)
        elif sp['gene'] in full_ko:
            err.append('initial_mrna_means[%s]: gene %s has a full knockout' % (sid, sp['gene']))
    for sid, v in (pert.get('initial_metabolites_mM') or {}).items():
        sp = species.get(sid)
        if not sp or sp['kind'] != 'metabolite' or sp.get('initial_mM') is None:
            err.append('initial_metabolites_mM: %r is not an intracellular metabolite with an initial concentration' % sid)
        elif not _num(v) or v < 0:
            err.append('initial_metabolites_mM[%s]: must be a number >= 0' % sid)
    for sid, v in (pert.get('medium_mM') or {}).items():
        sp = species.get(sid)
        if not sp or sp['kind'] != 'medium':
            err.append('medium_mM: %r is not a medium species (M_..._e)' % sid)
        elif not _num(v) or v < 0:
            err.append('medium_mM[%s]: must be a number >= 0' % sid)

    seen = set()
    for i, rp in enumerate(pert.get('reaction_parameters') or []):
        rid, key = rp.get('reaction'), rp.get('parameter')
        r = rxns.get(rid)
        if not r:
            err.append('reaction_parameters[%d]: unknown reaction %r' % (i, rid))
            continue
        if not r['active']:
            err.append('reaction_parameters[%d]: %s is not in the model (%s)' % (i, rid, r.get('note', 'inactive')))
        if key not in param_keys(r):
            err.append('reaction_parameters[%d]: %s has no parameter %r (has: %s)' % (i, rid, key, ', '.join(param_keys(r))))
        has_s, has_v = 'scale' in rp, 'value' in rp
        if has_s == has_v:
            err.append('reaction_parameters[%d]: give exactly one of scale, value' % i)
        elif not _num(rp.get('scale', rp.get('value'))) or rp.get('scale', rp.get('value')) < 0:
            err.append('reaction_parameters[%d]: scale/value must be a number >= 0' % i)
        if (rid, key) in seen:
            err.append('reaction_parameters[%d]: %s.%s set twice' % (i, rid, key))
        seen.add((rid, key))
    for rid in pert.get('disabled_reactions') or []:
        r = rxns.get(rid)
        if not r:
            err.append('disabled_reactions: unknown reaction %r' % rid)
        elif not r['active']:
            err.append('disabled_reactions: %s is not in the model' % rid)

    rdme = spec.get('rdme') or {}
    if any(pert.get(k) for k in ('rdme_rate_constants', 'gene_rate_scales', 'diffusion_scales', 'gip_constants')) and not rdme:
        err.append('the spec has no RDME record (build it in the container) so RDME parameters cannot be checked')
        return err, warn

    def scale_or_value(entry, label, i):
        has_s, has_v = 'scale' in entry, 'value' in entry
        if has_s == has_v:
            err.append('%s[%d]: give exactly one of scale, value' % (label, i))
        elif not _num(entry.get('scale', entry.get('value'))) or entry.get('scale', entry.get('value')) < 0:
            err.append('%s[%d]: scale/value must be a number >= 0' % (label, i))
    seen = set()
    for i, e in enumerate(pert.get('rdme_rate_constants') or []):
        c = e.get('constant')
        if c not in rdme.get('rate_constants', {}):
            err.append('rdme_rate_constants[%d]: unknown constant %r (known: %s)' % (i, c, ', '.join(rdme.get('rate_constants', {}))))
        elif not rdme['rate_constants'][c]['used_by']:
            warn.append('rdme_rate_constants[%d]: %s is defined but no reaction uses it' % (i, c))
        scale_or_value(e, 'rdme_rate_constants', i)
        if c in seen:
            err.append('rdme_rate_constants[%d]: %s set twice' % (i, c))
        seen.add(c)
    seen = set()
    for i, e in enumerate(pert.get('gene_rate_scales') or []):
        g, r = e.get('gene'), e.get('rate')
        if g != '*' and g not in genes:
            err.append('gene_rate_scales[%d]: unknown gene %r' % (i, g))
        if r not in GENE_RATES:
            err.append('gene_rate_scales[%d]: rate must be one of %s' % (i, ', '.join(GENE_RATES)))
        elif g in genes and r in ('translation', 'mrna_degradation') and genes[g]['type'] != 'protein':
            err.append('gene_rate_scales[%d]: %s has no %s (not a protein gene)' % (i, g, r))
        elif g in genes and r == 'secy_insertion' and not genes[g].get('membrane_insertion'):
            err.append('gene_rate_scales[%d]: %s is not a SecY client' % (i, g))
        if 'scale' not in e or not _num(e.get('scale')) or e['scale'] < 0:
            err.append('gene_rate_scales[%d]: scale must be a number >= 0' % i)
        if (g, r) in seen:
            err.append('gene_rate_scales[%d]: %s/%s set twice' % (i, g, r))
        seen.add((g, r))
    seen = set()
    dcs = set(rdme.get('diffusion_constants', {})) | {'rna_diff'}
    dcs.discard('zero')
    for i, e in enumerate(pert.get('diffusion_scales') or []):
        c = e.get('constant')
        if c not in dcs:
            err.append('diffusion_scales[%d]: unknown constant %r (known: %s)' % (i, c, ', '.join(sorted(dcs))))
        if 'scale' not in e or not _num(e.get('scale')) or e['scale'] < 0:
            err.append('diffusion_scales[%d]: scale must be a number >= 0' % i)
        if c in seen:
            err.append('diffusion_scales[%d]: %s set twice' % (i, c))
        seen.add(c)
    seen = set()
    for i, e in enumerate(pert.get('gip_constants') or []):
        c = e.get('constant')
        if c not in GIP_EDITABLE:
            err.append('gip_constants[%d]: constant must be one of %s' % (i, ', '.join(GIP_EDITABLE)))
        if 'value' not in e or not _num(e.get('value')) or e['value'] < 0:
            err.append('gip_constants[%d]: value must be a number >= 0' % i)
        if c in seen:
            err.append('gip_constants[%d]: %s set twice' % (i, c))
        seen.add(c)
    return err, warn


def resolve(spec, pert):
    """Explicit edits + static impact. Assumes validate() passed."""
    genes, species, rxns = spec['genes'], spec['species'], spec['reactions']
    edits, lost, reduced = [], set(), set()

    def edit(target, old, new, where, note=''):
        edits.append({'target': target, 'from': old, 'to': new, 'where': WHERE[where], 'note': note})
    for ko in pert.get('knockouts') or []:
        g = genes[ko['gene']]
        mode = ko.get('mode', 'full')
        n = g['num']
        if mode in ('full', 'expression_only'):
            edit('promoter strength %s' % g['locus'], g['promoter'], 0, 'promoter', 'RNAP binding rate becomes 0')
            lost.update(gene_species(spec, g['locus']))
        if mode in ('full', 'initial_only'):
            back = 'still transcribed (promoter %g), so it is made again' % g['promoter'] if mode == 'initial_only' else ''
            if g['type'] == 'protein':
                edit('initial P_%s' % n, g.get('initial_protein'), 0, 'protein', back)
                edit('initial R_%s (Poisson mean)' % n, g.get('mrna_initial_mean'), 0, 'mrna', back)
            elif g['type'] == 'tRNA':
                edit('initial R_%s' % n, g.get('initial_rna'), 0, 'trna', back)
        if mode == 'initial_only':
            reduced.update(gene_species(spec, g['locus']))
    for kd in pert.get('knockdowns') or []:
        g = genes[kd['gene']]
        edit('promoter strength %s' % g['locus'], g['promoter'], g['promoter'] * kd['promoter_scale'], 'promoter',
             'x %g' % kd['promoter_scale'])
        (reduced if kd['promoter_scale'] < 1 else set()).update(gene_species(spec, g['locus']))
    for sid, v in (pert.get('initial_protein_counts') or {}).items():
        old = species[sid].get('initial_count')
        edit('initial %s' % sid, old, v, 'protein', 'placed count; split over compartments as in initializeProteins')
        if old and v < old:
            reduced.add(sid)
    for sid, v in (pert.get('initial_mrna_means') or {}).items():
        edit('initial %s (Poisson mean)' % sid, species[sid].get('initial_mean'), v, 'mrna')
    for sid, v in (pert.get('initial_metabolites_mM') or {}).items():
        edit('initial %s (mM)' % sid, species[sid].get('initial_mM'), v, 'metabolite')
    for sid, v in (pert.get('medium_mM') or {}).items():
        edit('medium %s (mM)' % sid, species[sid].get('medium_mM'), v, 'medium', 'held fixed for the whole run')
    for rp in pert.get('reaction_parameters') or []:
        r = rxns[rp['reaction']]
        p = next(p for p in r['params'] if p.get('key') == rp['parameter'])
        new = p['value'] * rp['scale'] if 'scale' in rp else rp['value']
        edit('%s.%s' % (r['id'], rp['parameter']), p['value'], new, 'param_cme' if r['kind'] == 'cme_trna' else 'param_ode',
             ('x %g' % rp['scale']) if 'scale' in rp else '')
    dis = set()
    for rid in pert.get('disabled_reactions') or []:
        r = rxns[rid]
        dis.add(rid)
        if r['kind'] == 'ode_mm':
            edit('%s.onoff' % rid, 1, 0, 'onoff')
        else:
            for p in r['params']:
                if _num(p['value']) and p['name'] not in ('Radius',):
                    edit('%s.%s' % (rid, p['key']), p['value'], 0, 'param_cme' if r['kind'] == 'cme_trna' else 'param_ode',
                         'disabled: every rate parameter set to 0')
    rdme = spec.get('rdme') or {}
    for e in pert.get('rdme_rate_constants') or []:
        rc = rdme['rate_constants'][e['constant']]
        new = rc['value'] * e['scale'] if 'scale' in e else e['value']
        edit('rate constant %s' % e['constant'], rc['value'], new, 'rdme_rc',
             ('x %g; ' % e['scale'] if 'scale' in e else '') + 'used by ' + ', '.join(rc['used_by']))
    for e in pert.get('gene_rate_scales') or []:
        tmpl, where = GENE_RATES[e['rate']]
        target = ('every gene' if e['gene'] == '*' else e['gene'])
        edits.append({'target': '%s rate of %s' % (e['rate'], target), 'from': 'computed (%s)' % tmpl, 'to': 'x %g' % e['scale'],
                      'where': where, 'note': ''})
        if e['gene'] != '*' and e['scale'] < 1:
            reduced.update(gene_species(spec, e['gene']))
    for e in pert.get('diffusion_scales') or []:
        old = rdme['diffusion_constants'].get(e['constant'], {}).get('value', 'per-RNA value')
        edit('diffusion constant %s' % e['constant'], old, (old * e['scale']) if _num(old) else 'x %g' % e['scale'], 'diffusion', 'x %g' % e['scale'])
    for e in pert.get('gip_constants') or []:
        edit('GIP_rates.%s' % e['constant'], rdme['gip_constants'][e['constant']]['value'], e['value'], 'gip',
             'every rate computed from it changes')
    imp = evaluate(spec, lost, reduced, dis)
    return {'edits': edits, 'impact': imp, 'lost_species': sorted(lost), 'reduced_species': sorted(reduced)}


# ------------------------------------------------------------------------------------------------ YAML
def _scalar(v):
    if v is None:
        return 'null'
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if _num(v):
        return repr(v) if isinstance(v, float) else str(v)
    return json.dumps(str(v))          # a JSON string is a valid YAML double-quoted scalar


def dump(pert):
    """Deterministic YAML that both PyYAML and loads() read back."""
    out = ['# 4DWCM perturbation, modelspec schema %d. Check with: python -m modelspec check <this file>' % SCHEMA_VERSION]
    for k in TOP_KEYS:
        if k not in pert:
            continue
        v = pert[k]
        if isinstance(v, dict):
            if not v:
                out.append('%s: {}' % k)
                continue
            out.append('%s:' % k)
            for kk, vv in v.items():
                out.append('  %s: %s' % (kk, _scalar(vv)))
        elif isinstance(v, list):
            if not v:
                out.append('%s: []' % k)
                continue
            out.append('%s:' % k)
            for item in v:
                if isinstance(item, dict):
                    first = True
                    for kk, vv in item.items():
                        out.append('%s%s: %s' % ('  - ' if first else '    ', kk, _scalar(vv)))
                        first = False
                else:
                    out.append('  - %s' % _scalar(item))
        else:
            out.append('%s: %s' % (k, _scalar(v)))
    return '\n'.join(out) + '\n'


def _parse_scalar(s):
    s = s.strip()
    if s in ('', 'null', '~'):
        return None
    if s == '{}':
        return {}
    if s == '[]':
        return []
    if s in ('true', 'false'):
        return s == 'true'
    if s[0] in '"\'':
        return json.loads(s) if s[0] == '"' else s[1:-1].replace("''", "'")
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return s


def _strip_comment(line):
    out, q = [], None
    for ch in line:
        if q:
            if ch == q:
                q = None
        elif ch in '"\'':
            q = ch
        elif ch == '#':
            break
        out.append(ch)
    return ''.join(out).rstrip()


def loads(text):
    """JSON (the shared model browser saves .json), else PyYAML when present, else a reader for the subset dump() writes."""
    if text.lstrip().startswith('{'):
        return json.loads(text)
    try:
        import yaml
        return yaml.safe_load(text)
    except ImportError:
        pass
    root, key, cur = {}, None, None
    for raw in text.splitlines():
        line = _strip_comment(raw)
        if not line.strip():
            continue
        ind = len(line) - len(line.lstrip())
        body = line.strip()
        if ind == 0:
            k, _, v = body.partition(':')
            key = k.strip()
            root[key] = _parse_scalar(v) if v.strip() else None
            cur = None
        elif body.startswith('- '):
            if root[key] is None:
                root[key] = []
            item = body[2:]
            if re.match(r'^[A-Za-z_][\w]*\s*:', item):
                k, _, v = item.partition(':')
                cur = {k.strip(): _parse_scalar(v)}
                root[key].append(cur)
            else:
                root[key].append(_parse_scalar(item))
                cur = None
        else:
            k, _, v = body.partition(':')
            if cur is not None and ind >= 4:
                cur[k.strip()] = _parse_scalar(v)
            else:
                if root[key] is None:
                    root[key] = {}
                root[key][k.strip()] = _parse_scalar(v)
    return root


def load(path):
    return loads(open(path).read())
