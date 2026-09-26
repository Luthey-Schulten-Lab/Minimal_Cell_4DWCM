#!/usr/bin/env python3
"""Regenerate DNA/data.lammps_<step> by replaying the hook's own btree directive file (chromosome_operations_<step>.inp) up to and
including its last sys_write_sim_read_LAMMPS_data before 'translocate' — btree_chromo's own code writes the file. Needs the same
btree_chromo build the run used (a lean run carries it in wcm_sidecar/btree) and a GPU (btree_chromo links libcuda).
usage: data_replay.py <run_dir> <run_name> <scratch_dir> <all|part:i/n|step,...> [--jobs N] [--compare] [--emit DIR] [--skip s,...]
env:   BTREE_DIR (default /Software/btree_chromo), INPUT_DIR (default /inp), BTREE_BUILD_DIR (default build_fused)

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import os, re, sys, subprocess, argparse, hashlib, shutil, json, time
from concurrent.futures import ThreadPoolExecutor
BTREE_DIR = os.environ.get('BTREE_DIR', '/Software/btree_chromo').rstrip('/')        # btree_chromo tree (binary + DNA model)
INPUT_DIR = os.environ.get('INPUT_DIR', '/inp').rstrip('/')                          # folder holding the BD lengths file
BT = BTREE_DIR + '/' + os.environ.get('BTREE_BUILD_DIR', 'build_fused') + '/apps/btree_chromo'
IN_CMDS = ('input_state', 'load_mono_coords', 'load_loops', 'simulator_load_loop_params')
OUT_CMDS = ('prepare_simulator', 'dump_topology', 'output_state', 'write_loops', 'sys_write_sim_read_LAMMPS_data', 'write_mono_coords')

def ops_for(run, name, step, wd):
    """The hook's directive file cut after its last data write before 'translocate', with the run's own paths (whatever machine or
    container wrote them) mapped: <data dir>/ → run dir, BD lengths → INPUT_DIR, DNA model → BTREE_DIR, every output → wd."""
    src = open(os.path.join(run, 'DNA', 'chromosome_operations_%d.inp' % step)).read().split('\n'); last = None
    for i, l in enumerate(src):
        if l.startswith('translocate'): break
        if l.startswith('sys_write_sim_read_LAMMPS_data'): last = i
    if last is None: return None
    prefix = None                                                  # '<...>/<run name>/' as the run saw it
    for l in src:
        cmd, _, arg = l.partition(':')
        if cmd in IN_CMDS and '/DNA/' in arg: prefix = arg[:arg.rindex('/DNA/') + 1]; break
    out = []
    for l in src[:last + 1]:
        if not l.strip(): continue
        cmd, _, arg = l.partition(':')
        if cmd in IN_CMDS and prefix: arg = arg.replace(prefix, run.rstrip('/') + '/')
        elif cmd == 'load_BD_lengths': arg = INPUT_DIR + '/' + os.path.basename(arg)
        elif cmd == 'simulator_set_DNA_model': arg = BTREE_DIR + '/' + os.path.basename(arg.rstrip('/'))
        elif cmd in OUT_CMDS:
            base = arg.split(',', 1); fn = os.path.basename(base[0])
            if cmd == 'sys_write_sim_read_LAMMPS_data': fn = 'data.lammps_%d' % step     # lean runs name a rolling file here
            arg = os.path.join(wd, fn) + (',' + base[1] if len(base) > 1 else '')
        elif cmd == 'simulator_set_output_details': arg = wd.rstrip('/') + '/,' + arg.split(',', 1)[1]
        elif cmd == 'simulator_set_nProc': arg = '1'
        out.append(cmd + (':' + arg if arg else ''))
    return '\n'.join(out) + '\n'

def replay(args):
    run, name, outd, step, compare, emit = args
    wd = os.path.join(outd, 'w_%d' % step); os.makedirs(wd, exist_ok=True)
    ops = ops_for(run, name, step, wd)
    if ops is None: return {'step': step, 'ok': False, 'why': 'no data write in ops'}
    open(os.path.join(wd, 'ops.inp'), 'w').write(ops); t = time.time()
    r = subprocess.run([BT, os.path.join(wd, 'ops.inp')], cwd=wd, stdout=open(os.path.join(wd, 'btree.out'), 'w'), stderr=subprocess.STDOUT)
    res = {'step': step, 'rc': r.returncode, 's': round(time.time() - t, 1)}
    got = os.path.join(wd, 'data.lammps_%d' % step); res['ok'] = r.returncode == 0 and os.path.exists(got)
    if emit and res['ok']:                                       # extraction: move the file into the run's DNA/ and drop the scratch
        shutil.move(got, os.path.join(emit, 'data.lammps_%d' % step)); shutil.rmtree(wd); return res
    if compare and os.path.exists(got):
        a = open(got, 'rb').read(); b = open(os.path.join(run, 'DNA', 'data.lammps_%d' % step), 'rb').read()
        res['identical'] = a == b
        if a != b:
            A = a.split(b'\n'); B = b.split(b'\n'); i = next((j for j in range(min(len(A), len(B))) if A[j] != B[j]), None)
            res['first_diff'] = (i, B[i][:80].decode() if i is not None else None, A[i][:80].decode() if i is not None else None, len(A), len(B))
        else: shutil.rmtree(wd)
    return res
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('name'); ap.add_argument('out'); ap.add_argument('steps')
    ap.add_argument('--jobs', type=int, default=4); ap.add_argument('--compare', action='store_true'); ap.add_argument('--emit')
    ap.add_argument('--skip', default=''); a = ap.parse_args()
    allst = sorted(int(m.group(1)) for m in (re.match(r'dna_monomers_(\d+)\.bin$', f) for f in os.listdir(os.path.join(a.run, 'DNA'))) if m)
    if a.steps == 'all': st = allst
    elif a.steps.startswith('part:'):                            # part:i/n → every n-th step starting at i (split across GPUs)
        i, n = (int(x) for x in a.steps[5:].split('/')); st = allst[i::n]
    else: st = [allst[int(x[1:])] if x.startswith('@') else int(x) for x in a.steps.split(',')]
    skip = {int(x) for x in a.skip.split(',') if x}; st = [s for s in st if s not in skip]
    with ThreadPoolExecutor(a.jobs) as ex:
        R = list(ex.map(replay, [(a.run, a.name, a.out, s, a.compare, a.emit) for s in st]))
    for r in R: print(json.dumps(r))
    print(json.dumps({'steps': len(R), 'identical': sum(1 for r in R if r.get('identical')), 'failed': [r['step'] for r in R if not r['ok']][:20]}))
