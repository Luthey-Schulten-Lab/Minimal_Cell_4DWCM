#!/usr/bin/env python3
"""extract_run.sh helper. pre <out>: remember directory mtimes. post <out> [--keep-sidecar]: check that everything the lean run did not
write is back, give regenerated files a mode and an mtime (0644; data.lammps_<s> = mtime of dna_monomers_<s>.bin, chromosome.lammpstrj =
the last monomer file's, counts_fluxes_temp/* = counts_and_fluxes.csv's), restore directory mtimes, drop wcm_sidecar/.

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import os, sys, json, re, shutil
out = os.path.abspath(sys.argv[2]); sc = os.path.join(out, 'wcm_sidecar'); D = os.path.join(out, 'DNA')
if sys.argv[1] == 'pre':
    json.dump({r: os.stat(os.path.join(out, r)).st_mtime_ns for r in ('.', 'DNA', 'counts_fluxes_temp') if os.path.isdir(os.path.join(out, r))},
              open(os.path.join(sc, '.dir_mtimes.json'), 'w')); sys.exit(0)
steps = sorted(int(m.group(1)) for m in (re.match(r'dna_monomers_(\d+)\.bin$', f) for f in os.listdir(D)) if m)
problems = []
def stamp(p, ref):
    os.chmod(p, 0o644); t = os.stat(ref).st_mtime_ns; os.utime(p, ns=(t, t))
for s in steps:
    p = os.path.join(D, 'data.lammps_%d' % s)
    if os.path.exists(p): stamp(p, os.path.join(D, 'dna_monomers_%d.bin' % s))
    elif '--skip-restart-files' not in sys.argv: problems.append('missing DNA/data.lammps_%d' % s)
tr = os.path.join(D, 'chromosome.lammpstrj')
if os.path.exists(tr):
    stamp(tr, os.path.join(D, 'dna_monomers_%d.bin' % steps[-1]))
    with open(tr, 'rb') as f: nfr = f.read().count(b'ITEM: TIMESTEP\n')
    if nfr != len(steps) + 1: problems.append('trajectory has %d frames, expected %d' % (nfr, len(steps) + 1))
else: problems.append('missing DNA/chromosome.lammpstrj')
cf = os.path.join(out, 'counts_fluxes_temp'); csv = os.path.join(out, 'counts_and_fluxes.csv'); ncf = 0
if os.path.isdir(cf):
    for f in os.listdir(cf): stamp(os.path.join(cf, f), csv); ncf += 1
dm = json.load(open(os.path.join(sc, '.dir_mtimes.json')))
if '--keep-sidecar' not in sys.argv: shutil.rmtree(sc)
for r, t in sorted(dm.items(), key=lambda x: -x[0].count('/')): os.utime(os.path.join(out, r), ns=(t, t))
res = {'out': out, 'steps': len(steps), 'restart_files': sum(os.path.exists(os.path.join(D, 'data.lammps_%d' % s)) for s in steps),
       'trajectory': os.path.exists(tr), 'counts_temp_files': ncf, 'problems': problems}
print(json.dumps(res)); sys.exit(1 if problems else 0)
