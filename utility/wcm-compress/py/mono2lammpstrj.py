#!/usr/bin/env python3
"""mono2lammpstrj — rebuild a 4DWCM DNA trajectory (LAMMPS custom dump, as btree's write_dump writes chromosome.lammpstrj)
from files every run keeps:
  DNA/dna_monomers_<step>.bin        coordinates of frame k = hook steps[k-1], printed with %g
  DNA/chromo_topo_<step>.dat         leaves → ter/ori/fork types (btree::prepare_types), from the NEXT step's file
  DNA/loops/loops_<step>.txt         SMC anchor (7) / hinge (8) beads of hook steps[k-1] (override mono and fork beads)
  DNA/chromosome_operations_<step>.inp   the hook's boundary directive → membrane shell + box (bdry.py)
  DNA/chromosome_operations_<step>.inp + the run's BD lengths (r_bdry) → membrane shell (btree's own boundary_surface.cpp,
                                     bdry/bdry_gen) and box header (btree calc_bbox) — bdry.py
Frames that no kept file describes come from a literal sidecar (--literal DIR with chromosome_literal.lammpstrj + .idx):
  frame 0 (the minimised start inside the first hook), and frames where a loop leg indexes past the DNA into the shell
  (rng50 frame 164: two shell particles retyped 7/8 and moved by BD).
usage: mono2lammpstrj.py <run_dir> -o out.lammpstrj [--literal DIR] [--stride N] [--first K] [--last K] [--no-boundary] [--jobs N]
Default output (all frames, boundary on, sidecar given) is byte-identical to the dump.

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import os, re, sys, argparse
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

HDR = 'ITEM: ATOMS id type x y z c_id_track c_type_track\n'
RUN = None; STEPS = None; NOBD = False

def steps_of(run):
    D = os.path.join(run, 'DNA')
    return sorted(int(m.group(1)) for m in (re.match(r'dna_monomers_(\d+)\.bin$', f) for f in os.listdir(D)) if m)

def topo(run, step):
    lv = []
    for l in open(os.path.join(run, 'DNA', 'chromo_topo_%d.dat' % step)).read().split('\n')[1:]:
        if l.strip(): lv.append(tuple(int(x) for x in l.split(',')[1:6]))   # start_link, start, mid, end, end_link (1-based)
    return lv

LOOP_BLOCK = int(os.environ.get('M2L_LOOP_BLOCK', '-1'))   # which 'Number of loops' block when a file holds several

def loops(run, step):
    blocks = []
    for l in open(os.path.join(run, 'DNA', 'loops', 'loops_%d.txt' % step)).read().split('\n'):
        if l.startswith('Number'): blocks.append([])
        elif l and not l.startswith('Replication'): blocks[-1].append(tuple(int(x) for x in l.split('\t')[:2]))
    return blocks[LOOP_BLOCK] if blocks else []

def types(n, lv, an):
    t = [3] * (n + 1)
    for sl, s, m, e, el in lv:                         # btree::prepare_types (base 3: mono, 4 ori, 5 ter, 6 fork)
        for i in range(s, e + 1): t[i] = 3
        t[m] = 4
        if sl == e and el == s: t[s] = 5
        else: t[sl] = 6; t[el] = 6
    for a, h in an:                                    # loop anchors / hinges win over mono and fork, not over ori / ter;
        if a <= n and t[a] in (3, 6): t[a] = 7         # a loop leg beyond this hook's chromosome length is not in the dump
        if h <= n and t[h] in (3, 6): t[h] = 8
    return t

def irregular(run, steps, k):
    """True when frame k is not regenerable: a loop leg of hook steps[k-1] lies beyond the DNA (it retypes and moves a shell bead)."""
    if k == 0: return True
    st = steps[k - 1]; n = os.path.getsize(os.path.join(run, 'DNA', 'dna_monomers_%d.bin' % st)) // 24
    return any(a > n or h > n for a, h in loops(run, st))

_SH = {}
def _shell(d, rb):
    import bdry
    if d not in _SH: _SH[d] = bdry.shell(d, rb)[0]
    return _SH[d]

def bdry_directive(run, step):
    for l in open(os.path.join(run, 'DNA', 'chromosome_operations_%d.inp' % step)):
        if '_bdry:' in l and not l.startswith('#'): return l.strip()
    raise ValueError('no boundary directive for step %d' % step)

def frame(k):
    import bdry
    run, steps = RUN, STEPS; st = steps[k - 1]
    pre = os.path.join(run, 'DNA', 'dna_monomers_%d_prerotation.bin' % st)   # division: the dump precedes rotateChromosome
    src = pre if os.path.exists(pre) else os.path.join(run, 'DNA', 'dna_monomers_%d.bin' % st)
    x = np.fromfile(src, '<f8').reshape(-1, 3).tolist(); n = len(x)
    lv = topo(run, steps[k] if k < len(steps) else st); t = types(n, lv, loops(run, st))
    rb = bdry.r_bdry(run); d = bdry_directive(run, st); sh = [] if NOBD else _shell(d, rb).tolist()
    box = '\n'.join(bdry.box_lines(run, st, rb=rb)) + '\n'
    out = ['ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS ff ff ff\n' % (10000 if k == 1 else 0, n + len(sh)), box, HDR]
    out.append(''.join('%d %d %g %g %g %d %d\n' % (i + 1, t[i + 1], p[0], p[1], p[2], i + 1, t[i + 1]) for i, p in enumerate(x)))
    out.append(''.join('%d 1 %g %g %g %d 1\n' % (n + j + 1, p[0], p[1], p[2], n + j + 1) for j, p in enumerate(sh)))
    return ''.join(out).encode()

def literal_frames(ldir):
    """{frame index: bytes} from <ldir>/chromosome_literal.lammpstrj + .idx (one frame index per line, file order)."""
    if not ldir: return {}
    idx = [int(x) for x in open(os.path.join(ldir, 'chromosome_literal.idx')).read().split()]
    b = open(os.path.join(ldir, 'chromosome_literal.lammpstrj'), 'rb').read(); o = []; q = 0
    while True:
        i = b.find(b'ITEM: TIMESTEP\n', q)
        if i < 0: break
        o.append(i); q = i + 1
    return {k: b[a:e] for k, a, e in zip(idx, o, o[1:] + [len(b)])}

def init(run, steps, nobd):
    global RUN, STEPS, NOBD
    RUN, STEPS, NOBD = run, steps, nobd

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--literal'); ap.add_argument('--stride', type=int, default=1); ap.add_argument('--first', type=int, default=0)
    ap.add_argument('--last', type=int); ap.add_argument('--no-boundary', action='store_true'); ap.add_argument('--jobs', type=int, default=28)
    a = ap.parse_args(); steps = steps_of(a.run); last = a.last if a.last is not None else len(steps)
    lit = literal_frames(a.literal)
    ks = [k for k in range(a.first, last + 1, a.stride) if k > 0 or 0 in lit]    # frame 0 exists only as a literal
    gen = [k for k in ks if k not in lit]; warn = [k for k in gen if irregular(a.run, steps, k)]
    if warn: print('WARNING: frames %s are not regenerable from kept files (loop leg past the DNA) and have no literal copy; '
                   'written best-effort' % warn[:10], file=sys.stderr)
    with open(a.out, 'wb') as f, ProcessPoolExecutor(a.jobs, initializer=init, initargs=(a.run, steps, a.no_boundary)) as ex:
        it = ex.map(frame, gen, chunksize=2)
        for k in ks: f.write(lit[k] if k in lit else next(it))
    print('wrote %s: %d frames (%d literal)%s' % (a.out, len(ks), sum(k in lit for k in ks), ' (no boundary)' if a.no_boundary else ''))

if __name__ == '__main__': main()
