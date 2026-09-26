#!/usr/bin/env python3
"""Regenerate counts_fluxes_temp/counts_fluxes_<t>.csv (the per-time backups FileSaving keeps) from counts_and_fluxes.csv:
file <t> = 'name,value' lines of column <t> (header row 'Time,<t>'). usage: cftemp_regen.py <run_dir> <out_dir> [--verify]

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import os, sys
from concurrent.futures import ProcessPoolExecutor
def table(p):
    rows = open(p, 'rb').read().split(b'\n')
    if rows and rows[-1] == b'': rows = rows[:-1]
    return [r.split(b',') for r in rows]
def main(run, out, verify):
    tab = table(os.path.join(run, 'counts_and_fluxes.csv')); os.makedirs(out, exist_ok=True); bad = []; n = 0
    names = [r[0] for r in tab]
    for c, t in enumerate(tab[0][1:], 1):
        b = b''.join(nm + b',' + r[c] + b'\n' for nm, r in zip(names, tab))
        fn = 'counts_fluxes_%s.csv' % t.decode(); open(os.path.join(out, fn), 'wb').write(b); n += 1
        if verify:
            o = os.path.join(run, 'counts_fluxes_temp', fn)
            if not os.path.exists(o) or open(o, 'rb').read() != b: bad.append(fn)
    extra = sorted(set(os.listdir(os.path.join(run, 'counts_fluxes_temp'))) - set(os.listdir(out))) if verify else []
    print({'run': run, 'files': n, 'mismatch': len(bad), 'examples': bad[:5], 'originals_not_regenerated': extra[:5]})
if __name__ == '__main__': main(sys.argv[1], sys.argv[2], '--verify' in sys.argv)
