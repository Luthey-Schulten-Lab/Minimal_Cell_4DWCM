"""Compare short 4DWCM runs: key observables per run, and the read-back checks of each perturbation.

    python3 -m modelspec.compare_runs <run_dir> [<run_dir> ...]      (the first run is the reference)

Reads counts_and_fluxes.csv (rows = species / F_<reaction> fluxes, columns = biological seconds) and, when present,
perturbation_applied.json. Standard library only.
"""

import csv
import json
import os
import re
import sys

OBSERVABLES = [
    # (label, how)  how = row name, or ('sum', regex) over row names
    ('RNA polymerase (free RNAP)', 'RNAP'),
    ('RNAP on genes (sum RP_*)', ('sum', r'^RP_\d{4}_C[12]$')),
    ('transcripts made so far (sum RPM_*)', ('sum', r'^RPM_\d{4}$')),
    ('mRNA (sum R_nnnn of protein genes)', ('sum_mrna', None)),
    ('free ribosomes (ribosomeP)', 'ribosomeP'),
    ('mRNA-bound ribosomes (sum RB_*)', ('sum', r'^RB_\d{4}$')),
    ('proteins made so far (sum PM_*)', ('sum', r'^PM_\d{4}$')),
    ('rpoA protein P_0645', 'P_0645'),
    ('pdhC protein P_0227', 'P_0227'),
    ('dnaA protein P_0001', 'P_0001'),
    ('smc protein P_0415', 'P_0415'),
    ('flux PGI', 'F_PGI'),
    ('flux PDH_E3', 'F_PDH_E3'),
    ('flux PDH_acald', 'F_PDH_acald'),
    ('flux GLCpts4 (glucose import)', 'F_GLCpts4'),
    ('cytoplasmic glucose M_glc__D_c', 'M_glc__D_c'),
    ('G6P M_g6p_c', 'M_g6p_c'),
    ('ATP M_atp_c', 'M_atp_c'),
]


def load(run):
    rows, times = {}, []
    with open(os.path.join(run, 'counts_and_fluxes.csv')) as f:
        r = csv.reader(f)
        head = next(r)
        times = [float(t) for t in head[1:]]
        for line in r:
            try:
                rows[line[0]] = [float(x) if x not in ('', 'nan') else float('nan') for x in line[1:]]
            except ValueError:
                continue
    return times, rows


def series(rows, how, mrna_ids):
    if isinstance(how, str):
        return rows.get(how)
    kind, pat = how
    keys = [k for k in rows if (re.match(pat, k) if kind == 'sum' else k in mrna_ids)]
    if not keys:
        return None
    n = len(rows[keys[0]])
    return [sum(rows[k][i] for k in keys) for i in range(n)]


def fmt(v):
    if v is None:
        return '—'
    if v != v:
        return 'nan'
    if abs(v) >= 1e4 or (abs(v) < 1e-3 and v != 0):
        return '%.3g' % v
    return ('%.4g' % v)


def main(runs):
    spec_p = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'modelspec', 'out', 'model_spec.json')
    mrna_ids = set()
    if os.path.exists(spec_p):
        spec = json.load(open(spec_p))
        mrna_ids = {k for k, v in spec['species'].items() if v['kind'] == 'mRNA'}
    data = {}
    for run in runs:
        name = os.path.basename(run.rstrip('/'))
        try:
            data[name] = load(run)
        except FileNotFoundError:
            print('%s: no counts_and_fluxes.csv (run not finished or crashed)' % name)
    names = list(data)
    if not names:
        return
    print('%-40s' % 'observable (t=0 -> last second, mean)' + ''.join('%26s' % n for n in names))
    for label, how in OBSERVABLES:
        cells = []
        for n in names:
            times, rows = data[n]
            s = series(rows, how, mrna_ids)
            if s is None:
                cells.append('—')
                continue
            vals = [v for v in s if v == v]
            mean = sum(vals) / len(vals) if vals else float('nan')
            cells.append('%s -> %s (%s)' % (fmt(s[0]), fmt(s[-1]), fmt(mean)))
        print('%-40s' % label + ''.join('%26s' % c for c in cells))
    print('\nlast biological second: ' + ', '.join('%s %g' % (n, data[n][0][-1]) for n in names))
    for run in runs:
        p = os.path.join(run, 'perturbation_applied.json')
        if os.path.exists(p):
            a = json.load(open(p))
            print('\n%s: perturbation %r, %d/%d read-back checks passed, %d edits applied' % (
                os.path.basename(run.rstrip('/')), a['name'], a['checks_passed'], a['checks_passed'] + a['checks_failed'], len(a['log'])))
            for c in a['checks']:
                if not c['ok']:
                    print('   FAILED %s: expected %s, got %s' % (c['what'], c['expected'], c['got']))


if __name__ == '__main__':
    main(sys.argv[1:])
