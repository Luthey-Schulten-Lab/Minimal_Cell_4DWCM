"""
Chromosome configuration via btree_chromo (Benjamin Gilbert), including replication.

Authors
-------
Ron Acda — resident btree_chromo --serve process, array-based chromosome lattice update,
    lean output (wcm_sidecar)
    (using an iterative LLM-guided workflow: https://github.com/quarkron/iterative-hillclimber/tree/main)
Alfia Parvez — SMC looping hooks, BD walltime scale, in-place LAMMPS update,
    CPU isolation (``DNA_CPU_CORES``), topo/loops rescue paths
Zane Thornburg — original btree_chromo coupling and DNA particle remap
"""

import numpy as np

import os
import re

import subprocess

import shutil

from subprocess import Popen, PIPE

import time as timepy

from Bio import SeqIO
from Bio.Seq import Seq

from utility.LatticeFunctions import *

import utility.GIP_rates as GIP



# Biological time between DNA hook updates (seconds); must match Hook.next_DNA_time increment.
DNA_HOOK_INTERVAL_S = 4.0

# template_replicate.inp advances loops / BD once per 2 s batch inside repeat:30
TEMPLATE_HOOK_BATCH_S = 2.0

# Syn3A chromosome length in 10 bp RDME/LAMMPS beads (loop_params N=)
SYN3A_CHROMO_BEADS = 54338


_wcm_coords_cache = [None]   # (DNA_lattice_coords dict, its (n, 3) int64 array) of the last placeNewChromosome


def _dna_work_dir(sim_properties):
    return sim_properties['working_directory'] + 'DNA/'


def _loops_dir(sim_properties):
    loops_dir = _dna_work_dir(sim_properties) + 'loops/'
    os.makedirs(loops_dir, exist_ok=True)
    return loops_dir


def _loops_file(sim_properties, timestep):
    return _loops_dir(sim_properties) + 'loops_{:d}.txt'.format(timestep)


def _loop_prng_seed(sim_properties, timestep):
    return int(timestep / 10000 + sim_properties['dna_rng_seed'])


def _translocate_speed(sim_properties):
    """Per template_replicate.inp / run_btree_chromo_replicate.py (v_bps/10 per 2 s batch)."""
    v_bps = sim_properties.get('dna_loop_translocate_bps', 200)
    return max(1, int(round(v_bps / 10.0)))


def _translocate_steps_for_hook(sim_properties):
    """translocate:N,T steps for one DNA hook, scaled from the template 2 s batch size."""
    speed = _translocate_speed(sim_properties)
    hook_s = sim_properties.get('dna_hook_interval_s', DNA_HOOK_INTERVAL_S)
    return max(1, int(round(speed * hook_s / TEMPLATE_HOOK_BATCH_S)))



# the run does not write what the kept files regenerate byte for byte (utility/wcm-compress):
#   DNA/chromosome.lammpstrj   -> btree dumps each hook's frame into wcm_sidecar/hookdump.lammpstrj; the frame is kept in
#                                 wcm_sidecar/literal/chromosome_literal.lammpstrj (+ .idx) only for the first hook (frames 0, 1)
#                                 and for irregular frames (a shell atom retyped by a loop leg past the DNA), then the file is dropped
#   DNA/data.lammps_<step>     -> one rolling wcm_sidecar/data.lammps_current (the btree -> LAMMPS hand-off); the first hook's
#                                 file is written straight into wcm_sidecar/literal/DNA/ (it holds the equilibrated start)
#   counts_fluxes_temp/*.csv   -> pruned after the final concatenation when each file equals its column (FileSaving)
# WCM_LEAN_OUTPUT=0 = the original layout. WCM_LEAN_SHADOW=1 = lean AND the original files (gate runs).
def _lean():
    return os.environ.get('WCM_LEAN_OUTPUT', '1') != '0'


def _lean_shadow():
    return _lean() and os.environ.get('WCM_LEAN_SHADOW', '0') == '1'


def _sidecar_dir(workDir):
    d = os.path.dirname(workDir.rstrip('/')) + '/wcm_sidecar/'
    os.makedirs(d + 'literal/DNA', exist_ok=True)
    return d


def _data_path(workDir, timestep, sim_properties):
    if not _lean():
        return workDir + 'data.lammps_{:d}'.format(timestep)
    if _is_first_dna_step(sim_properties):
        return _sidecar_dir(workDir) + 'literal/DNA/data.lammps_{:d}'.format(timestep)
    return _sidecar_dir(workDir) + 'data.lammps_current'


def _output_details(workDir):
    if not _lean():
        return 'simulator_set_output_details:' + workDir + ',chromosome\n'
    return 'simulator_set_output_details:' + _sidecar_dir(workDir) + ',hookdump\n'


def _split_frames(b):
    o = [m.start() for m in re.finditer(rb'ITEM: TIMESTEP\n', b)]
    return [b[a:e] for a, e in zip(o, o[1:] + [len(b)])]


def _frame_irregular(fr, n_dna):
    """A shell atom (id > n_dna) with a type other than 1: a loop leg past the DNA retyped it; the frame is not regenerable."""
    body = fr.split(b'ITEM: ATOMS', 1)[1].split(b'\n', 1)[1]
    for line in body.split(b'\n'):
        if line:
            i, t = line.split(b' ', 2)[:2]
            if int(i) > n_dna and t != b'1':
                return True
    return False


def _lean_setup_sidecar(sim_properties, workDir, step):
    """First hook: LEAN.json and the exact btree build + BD lengths the regenerators replay with (~8 MB)."""
    import json, shutil
    sc = _sidecar_dir(workDir)
    if os.path.exists(sc + 'LEAN.json'):
        return
    headDir = sim_properties['dna_software_directory']; bdir = headDir + 'btree_chromo/'
    build = os.environ.get('BTREE_BUILD_DIR', 'build')
    try:
        os.makedirs(sc + 'btree/' + build + '/apps', exist_ok=True)
        shutil.copy2(bdir + build + '/apps/btree_chromo', sc + 'btree/' + build + '/apps/btree_chromo')
        model = os.environ.get('DNA_MODEL_DIR_NAME', 'LAMMPS_DNA_model_kk')
        if not os.path.exists(sc + 'btree/' + model):
            shutil.copytree(bdir + model, sc + 'btree/' + model)
        ops = open(workDir + 'chromosome_operations_{:d}.inp'.format(step)).read()
        m = re.search(r'^load_BD_lengths:(\S+)$', ops, re.M)
        if m:
            os.makedirs(sc + 'input_data', exist_ok=True); shutil.copy2(m.group(1), sc + 'input_data/' + os.path.basename(m.group(1)))
    except Exception as e:
        print('could not copy the btree build into the sidecar:', e)
    json.dump({'format': 1, 'writer': 'lean output', 'run_dir_name': os.path.basename(os.path.dirname(workDir.rstrip('/'))),
               'not_written': ['DNA/chromosome.lammpstrj', 'DNA/data.lammps_<step>', 'counts_fluxes_temp/counts_fluxes_<t>.csv'],
               'regenerate_with': 'utility/wcm-compress/extract_run.sh',
               'btree_build_dir': build, 'shadow': _lean_shadow()}, open(sc + 'LEAN.json', 'w'), indent=1)


def _lean_collect(sim_properties, workDir, step):
    """After hook <step>'s outputs are complete: keep its dumped frame(s) only if not regenerable, then drop the scratch dump."""
    if not _lean():
        return
    sc = _sidecar_dir(workDir); hd = sc + 'hookdump.lammpstrj'
    if not os.path.exists(hd):
        return
    b = open(hd, 'rb').read(); frames = _split_frames(b)
    steps = sorted(int(m.group(1)) for m in (re.match(r'dna_monomers_(\d+)\.bin$', f) for f in os.listdir(workDir)) if m)
    k = steps.index(step) + 1 if step in steps else len(steps)
    first = (steps and step == steps[0])
    idx = list(range(0, len(frames))) if first else list(range(k - len(frames) + 1, k + 1))
    n_dna = os.path.getsize(workDir + 'dna_monomers_{:d}.bin'.format(step)) // 24
    keep = [(i, fr) for i, fr in zip(idx, frames) if first or _frame_irregular(fr, n_dna)]
    if keep:
        with open(sc + 'literal/chromosome_literal.lammpstrj', 'ab') as f:
            for _, fr in keep: f.write(fr)
        with open(sc + 'literal/chromosome_literal.idx', 'a') as f:
            for i, _ in keep: f.write('%d\n' % i)
    if _lean_shadow():
        with open(workDir + 'chromosome.lammpstrj', 'ab') as f: f.write(b)
        import shutil
        src = sc + ('literal/DNA/data.lammps_{:d}'.format(step) if first else 'data.lammps_current')
        if os.path.exists(src): shutil.copy2(src, workDir + 'data.lammps_{:d}'.format(step))
    if first:
        _lean_setup_sidecar(sim_properties, workDir, step)
    os.remove(hd)


def leanFinalize(sim_properties):
    """End of the run: wait for the last DNA hook and collect its frame (the next hook would have done it)."""
    if _lean() and sim_properties.get('last_DNA_step') is not None:
        workDir = sim_properties['working_directory'] + 'DNA/'
        if os.path.exists(_sidecar_dir(workDir) + 'hookdump.lammpstrj'):
            checkLastChromosome(sim_properties)
        cur = _sidecar_dir(workDir) + 'data.lammps_current'       # the rolling hand-off file is scratch once the run ends
        if os.path.exists(cur):
            os.remove(cur)


def _is_first_dna_step(sim_properties):
    return sim_properties['last_DNA_step'] is None


def _parse_loops_snapshot(loops_path):
    """Return (loop_line_strings, fork_metadata_line) from the last snapshot in a loops file."""
    if not os.path.isfile(loops_path):
        return [], None

    with open(loops_path, 'r') as f:
        lines = f.readlines()

    header_idx = -1
    for i, line in enumerate(lines):
        if line.startswith('Number of loops:'):
            header_idx = i

    if header_idx < 0:
        return [], None

    fork_line = None
    if header_idx + 1 < len(lines) and lines[header_idx + 1].startswith('Replication forks:'):
        fork_line = lines[header_idx + 1].rstrip('\n')

    records = []
    for line in lines[header_idx + 2:]:
        stripped = line.strip()
        if not stripped:
            break
        parts = stripped.split()
        if len(parts) >= 2:
            records.append(line if line.endswith('\n') else line + '\n')

    return records, fork_line


def _write_loops_snapshot(loops_path, records, fork_line=None):
    with open(loops_path, 'w') as f:
        f.write('Number of loops: {:d}\n'.format(len(records)))
        if fork_line is not None:
            f.write(fork_line + '\n')
        else:
            f.write('Replication forks: 0, 0\n')
        for line in records:
            f.write(line)


def _prepare_loops_file_for_hook(sim_properties, from_timestep, m_target):
    """
    Trim the previous hook's loops snapshot to m_target rows so it matches numSmc
    written to loop_params.txt for this hook.
    """
    loops_path = _loops_file(sim_properties, from_timestep)
    records, fork_line = _parse_loops_snapshot(loops_path)

    if len(records) > m_target:
        print('Trimming loops file from {:d} to {:d} SMCs'.format(len(records), m_target))
        records = records[:m_target]
    elif len(records) < m_target:
        print('WARNING: loops file has {:d} entries but RDME requests {:d} SMCs; '
              'btree_chromo will birth the remainder after load'.format(len(records), m_target))

    _write_loops_snapshot(loops_path, records, fork_line)
    return loops_path


# Wall-clock seconds to wait for a background btree_chromo step before rescue.
# The first hook equilibrates with long BD + translocate and often exceeds 5 min.
DNA_BTREE_WAIT_TIMEOUT_S = 1800


def _chromo_topo_path(sim_properties, step):
    return _dna_work_dir(sim_properties) + 'chromo_topo_{:d}.dat'.format(step)


def _write_prereplication_chromo_topo(sim_properties, step, n_beads=SYN3A_CHROMO_BEADS):
    """Write a non-replicated topology file matching btree dump_topology format."""
    topo_path = _chromo_topo_path(sim_properties, step)
    ori_idx = n_beads // 2
    with open(topo_path, 'w') as f:
        f.write('size={}\n'.format(n_beads))
        f.write('m({}),{},1,{},{},1\n'.format(n_beads, n_beads, ori_idx, n_beads))
    return topo_path


def _ensure_chromo_topo(sim_properties, step=None):
    """Return chromo_topo path, synthesizing pre-replication topology if btree omitted it."""
    if step is None:
        step = sim_properties['last_DNA_step']
    topo_path = _chromo_topo_path(sim_properties, step)
    if os.path.isfile(topo_path) and os.path.getsize(topo_path) > 0:
        return topo_path
    n_beads = sim_properties.get('counts', {}).get('chromosome', SYN3A_CHROMO_BEADS)
    if not n_beads or n_beads < SYN3A_CHROMO_BEADS:
        n_beads = SYN3A_CHROMO_BEADS
    print('WARNING: missing {}; writing pre-replication topology ({} beads)'.format(
        topo_path, n_beads))
    return _write_prereplication_chromo_topo(sim_properties, step, n_beads)


# the btree_chromo process started for each DNA step, so checkLastChromosome can wait on its exit instead of polling
# for its output files at 0.5-2 s intervals (on average a quarter to a full second of extra RDME idling per DNA hook whenever
# btree is the longer arm; a size check can also see a file btree is still writing). {step: Popen}
_btree_procs = {}


# one resident `btree_chromo --serve` process runs the DNA steps' directive files one after another (each through a
# fresh btree driver, as a separate process would), so the ~1.6 s per DNA step of process start, CUDA context, MPI and Kokkos/LAMMPS
# initialisation is paid once per run. A step handle waits on the server's "BTREE_SERVE_DONE <rc> <s> <file>" line the way
# checkLastChromosome waited on the step's process exit. Any failure (a binary without --serve, a server that exits or reports a
# nonzero rc) falls back to the per-step process for that step and every later one. WCM_BTREE_SERVE=0 disables the server.
import atexit, select


class _BtreeServer:
    def __init__(self, cmd, env):
        self.cmd, self.env, self.proc, self.dead = list(cmd), dict(env), None, False

    def start(self, timeout=120.0):
        try:
            self.proc = Popen(self.cmd + ['--serve'], stdin=PIPE, stdout=PIPE, stderr=None, env=self.env, bufsize=0)   # unbuffered: select sees every line
        except Exception as exc:
            print('btree server: start failed ({}), per-step processes'.format(exc)); self.dead = True; return False
        line = self._readline(timeout)
        if line is None or not line.startswith('BTREE_SERVE_READY'):
            print('btree server: no READY line ({!r}), per-step processes'.format(line)); self.stop(); self.dead = True; return False
        print('btree server: started, pid {}'.format(self.proc.pid))
        return True

    def _readline(self, timeout):
        """Next protocol line, None on timeout or EOF (EOF marks the server dead)."""
        end = timepy.time() + timeout
        while True:
            left = end - timepy.time()
            if left <= 0:
                return None
            r, _, _ = select.select([self.proc.stdout], [], [], left)
            if not r:
                return None
            raw = self.proc.stdout.readline()
            if not raw:
                self.dead = True
                return None
            line = raw.decode(errors='replace').strip()
            if line.startswith('BTREE_SERVE_'):
                return line

    def submit(self, path):
        try:
            self.proc.stdin.write((path + '\n').encode()); self.proc.stdin.flush(); return True
        except Exception as exc:
            print('btree server: submit failed ({}), per-step processes'.format(exc)); self.dead = True; return False

    def stop(self):
        if self.proc is None:
            return
        try:
            if self.proc.poll() is None:
                self.proc.stdin.write(b'EXIT\n'); self.proc.stdin.flush(); self.proc.wait(timeout=30)
        except Exception:
            try: self.proc.kill()
            except Exception: pass
        self.proc = None


class _ServedStep:
    """Stands in for the step's Popen in _btree_procs: wait(timeout) returns when the server reports the file done."""
    def __init__(self, server, path, fallback_cmd, env):
        self.server, self.path, self.fallback_cmd, self.env, self.done, self.rc = server, path, fallback_cmd, env, False, None

    def wait(self, timeout=None):
        if self.done:
            return self.rc
        end = None if timeout is None else timepy.time() + timeout
        while not self.server.dead:
            left = 1e9 if end is None else end - timepy.time()
            if left <= 0:
                raise subprocess.TimeoutExpired(self.path, timeout)
            line = self.server._readline(min(left, 3600.0))
            if line is None:
                continue
            parts = line.split(' ', 3)
            if parts[0] == 'BTREE_SERVE_DONE' and len(parts) == 4 and parts[3] == self.path:
                self.rc = int(parts[1])
                print('btree server: {} done rc={} in {} s'.format(self.path, parts[1], parts[2]))
                self.done = True
                if self.rc != 0:   # as a failed per-step process: no rerun; later steps get their own processes
                    print('btree server: nonzero rc, per-step processes from now on'); self.server.stop(); self.server.dead = True
                return self.rc
        # the server died with this step unfinished: run the step in its own process, as without the server
        print('btree server: unavailable, running {} in its own process'.format(self.path))
        p = Popen(self.fallback_cmd, stdin=None, stdout=subprocess.DEVNULL, stderr=None, env=self.env)
        self.rc = p.wait(timeout=None if end is None else max(1.0, end - timepy.time()))
        self.done = True
        return self.rc


_btree_server = None


def _launch_btree(launch_cmd, directives, env):
    """The step's handle: a served request when the resident server runs, else the step's own process."""
    global _btree_server
    if os.environ.get('WCM_BTREE_SERVE', '1') != '0':
        base = launch_cmd[:-1]
        if _btree_server is None:
            _btree_server = _BtreeServer(base, env)
            if _btree_server.start():
                atexit.register(_btree_server.stop)
        if not _btree_server.dead and _btree_server.cmd == list(base) and _btree_server.env == dict(env) and _btree_server.submit(directives):
            return _ServedStep(_btree_server, directives, launch_cmd, env)
    return Popen(launch_cmd, stdin=None, stdout=subprocess.DEVNULL, stderr=None, env=env)


def _dna_outputs_ready(workDir, step):
    mono = workDir + 'dna_monomers_{:d}.bin'.format(step)
    topo = workDir + 'chromo_topo_{:d}.dat'.format(step)
    return (
        os.path.isfile(mono) and os.path.getsize(mono) > 0
        and os.path.isfile(topo) and os.path.getsize(topo) > 0
    )


def _load_loops_directive(sim_properties, loops_path):
    """template_replicate.inp line 23 — omitted on the first DNA hook only."""
    if _is_first_dna_step(sim_properties):
        return ''
    return 'load_loops:' + loops_path


def _equilibrate_loops_directive(sim_properties):
    """template_replicate.inp line 27 — first hook only."""
    if not _is_first_dna_step(sim_properties):
        return ''
    steps = sim_properties.get('dna_loop_equilibrate_steps', 360000)
    return 'translocate:{:d},F'.format(steps)


def _append_string(sim_properties):
    """template_replicate.inp {append_string} on line 24."""
    if _is_first_dna_step(sim_properties):
        return 'noappend,first'
    return 'append,nofirst'


def _initial_soft_harmonic_directive(sim_properties):
    """template_replicate.inp line 24 — initial frame before the hook repeat body."""
    run_steps = sim_properties.get('dna_initial_soft_harmonic_steps', 10000)
    thermo_freq = 1000
    steps_before_output = sim_properties.get('dna_initial_soft_harmonic_output', 20000)
    return 'simulator_run_soft_harmonic:{:d},{:d},{:d},{}'.format(
        run_steps, thermo_freq, steps_before_output, _append_string(sim_properties))


def _run_dynamics_directive(sim_properties):
    """template_replicate.inp line 44 — BD during the hook, scaled from 2 s template defaults."""
    hook_s = sim_properties.get('dna_hook_interval_s', DNA_HOOK_INTERVAL_S)
    scale = hook_s / TEMPLATE_HOOK_BATCH_S
    wall_scale = float(sim_properties.get('dna_bd_walltime_scale', 1.0))
    run_steps = max(
        1000,
        int(sim_properties.get('dna_bd_run_steps', 20000) * scale * wall_scale),
    )
    # dump_freq (arg 3) == run_steps → one lammpstrj frame at the end of each hook BD batch
    dump_freq = run_steps
    return 'simulator_run_soft_harmonic:{:d},1000,{:d},append,nofirst'.format(
        run_steps, dump_freq)


def _num_smc(sim_properties):
    """
    Active SMC loopers for btree_chromo (numSmc in loop_params.txt).

    Uses the same basis as the legacy simulator_run_loops loop count:
    P_0415 is the SMC protein count on the RDME lattice; each complex uses two
    copies, so complexes = P_0415 / 2. Only the bound fraction participates in
    loop extrusion (default 1.0 = all complexes active).
    """
    bound_fraction = float(sim_properties.get('dna_smc_bound_fraction', 1.0))
    n_complexes = sim_properties['counts']['P_0415'] / 2.0
    return max(1, int(n_complexes * bound_fraction))


def _write_loop_params_file(sim_properties, num_smc):
    """Write loop_params.txt for the current hook, overriding numSmc for btree_chromo."""
    work_dir = _dna_work_dir(sim_properties)
    template_path = sim_properties['head_directory'] + 'input_data/loop_params.txt'
    out_path = work_dir + 'loop_params.txt'

    overrides = {'numSmc': str(num_smc)}
    with open(template_path, 'r') as template_file:
        lines = template_file.readlines()

    with open(out_path, 'w') as out_file:
        for line in lines:
            stripped = line.strip()
            if stripped and not stripped.startswith('#') and '=' in stripped:
                key = stripped.split('=', 1)[0].strip()
                if key in overrides:
                    out_file.write('{}={}\n'.format(key, overrides[key]))
                    continue
            out_file.write(line)

    print('DNA hook loop_params: numSmc={:d} (P_0415={:d}, bound_fraction={})'.format(
        num_smc,
        int(sim_properties['counts']['P_0415']),
        sim_properties.get('dna_smc_bound_fraction', 1.0),
    ))

    return out_path


def _write_replicate_hook_protocol(f, sim_properties, timestep, workDir, rep_started):
    """
    SMC looping + BD for one 4 s DNA hook, following template_replicate.inp and
    run_btree_chromo_replicate.py (load_loops, equilibrate_loops, repeat body, run_dynamics).
    """
    m_target = _num_smc(sim_properties)

    loops_path = None
    if not _is_first_dna_step(sim_properties):
        loops_path = _prepare_loops_file_for_hook(
            sim_properties, sim_properties['last_DNA_step'], m_target)

    loop_params_path = _write_loop_params_file(sim_properties, m_target)
    f.write('simulator_load_loop_params:' + loop_params_path + '\n')

    load_loops = _load_loops_directive(sim_properties, loops_path)
    if load_loops:
        f.write(load_loops + '\n')

    if _is_first_dna_step(sim_properties):
        f.write(_initial_soft_harmonic_directive(sim_properties) + '\n')

    equilibrate_loops = _equilibrate_loops_directive(sim_properties)
    if equilibrate_loops:
        f.write(equilibrate_loops + '\n')

    # One hook-interval block (template_replicate.inp repeat body, lines 31-44)
    f.write('sync_simulator_and_system\n')

    if rep_started:
        if _is_first_dna_step(sim_properties):
            repCw = 40
            repCcw = 40
        else:
            repCw, repCcw = getReplicatedSegments(sim_properties)
        f.write('set_initial_state\n')
        f.write('transform:m_cw' + str(repCw) + '_ccw' + str(repCcw) + '\n')
        f.write('set_final_state\n')
        f.write('output_state:' + workDir + 'rep_state_{:d}.txt\n'.format(timestep))
        f.write('map_replication\n')

    f.write('write_loops:' + _loops_file(sim_properties, timestep) + '\n')
    f.write('sys_write_sim_read_LAMMPS_data:' + _data_path(workDir, timestep, sim_properties) + '\n')

    use_fork_repulsion = bool(sim_properties.get('dna_fork_partition_repulsion', False))
    if use_fork_repulsion and not checkDaughtersFullyPartitioned(sim_properties):
        print('DNA not fully partitioned yet; fork partition repulsion ON during looping')
        f.write('switch_fork_partition_repulsion:T\n')

    f.write('translocate:{:d},T\n'.format(_translocate_steps_for_hook(sim_properties)))
    # Prefer in-place LAMMPS coord update vs a second full write/read (dna_second_roundtrip_inplace).
    if bool(sim_properties.get('dna_second_roundtrip_inplace', True)):
        f.write('sys_update_lammps_inplace:' + _data_path(workDir, timestep, sim_properties) + '\n')
    else:
        f.write('sys_write_sim_read_LAMMPS_data:' + _data_path(workDir, timestep, sim_properties) + '\n')
    f.write('simulator_form_loops:F\n')
    f.write('simulator_minimize_topoDNA_harmonic:1000\n')
    f.write('simulator_set_delta_t:2.5E+7\n')
    f.write(_run_dynamics_directive(sim_properties) + '\n')


#########################################################################################
def updateChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict, updateRegions):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    if sim_properties['last_DNA_step'] is None:
        
        DNA_lattice_coords = sim_properties['DNAcoords']
        
        genome = sim_properties['genome']
        
    elif sim_properties['last_DNA_step'] is not None:
    
        checkLastChromosome(sim_properties)
        
        if ('rotated_DNA' not in sim_properties.keys()) and sim_properties['division_started']:
            
            rotateChromosome(time, sim_properties)

#         deleteOldChromosome(lattice)
    
        region_dict, DNA_lattice_coords = placeNewChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict)

        genome = remapDNA(sim_properties)

        sim_properties['genome'] = genome
        
        moveDnaParticles(sim_properties, lattice, DNA_lattice_coords)
        
    if np.rint(time) < sim_properties['total_time']:

        writeRiboObstacleFile(time, ribo_site_dict, sim_properties)

    #     writeMembraneBoundaryFile(time, region_dict, sim_properties)

        writeChromosomeInputFile(time, sim_properties, updateRegions)

        runNewChromosome(time, sim_properties)
        
        print('Started running new chromosome in background')
        
    return region_dict, DNA_lattice_coords, genome
#########################################################################################


#########################################################################################
def updateChromosomeDivision(time, lattice, sim_properties, region_dict, ribo_site_dict):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers

    Returns:
        region_dict, DNA_lattice_coords, genome
    Called by:
        Hook.run() when (division_started and updateRegions) is True.
    Description:
        Division-phase chromosome update for the protein_science (persistent-SMC)
        branch. writeChromosomeInputFile() already folds the dividing-membrane
        boundary (overlapping_spheres_bdry) into the normal replicate-hook protocol
        when division_started is True ("only the membrane boundary shape changes"),
        so a division DNA step uses the same persistent-SMC btree_chromo path as
        growth, with updateRegions forced True.

        This mirrors the control flow of the upstream updateChromosomeDivision
        (rotate -> place -> remap -> move -> optional late-division partition ->
        write + run) but routes btree_chromo through writeChromosomeInputFile /
        runNewChromosome instead of the legacy writeDivisionChromosomeInputFile /
        runDivChromosome (which do not exist on this branch and emit the old
        non-persistent-SMC directives).
    """
    print('DNA Division Step')

    if sim_properties['last_DNA_step'] is None:

        DNA_lattice_coords = sim_properties['DNAcoords']

        genome = sim_properties['genome']

    elif sim_properties['last_DNA_step'] is not None:

        checkLastChromosome(sim_properties)

        if 'rotated_DNA' not in sim_properties.keys():

            rotateChromosome(time, sim_properties)

        region_dict, DNA_lattice_coords = placeNewChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict)

        genome = remapDNA(sim_properties)

        sim_properties['genome'] = genome

        moveDnaParticles(sim_properties, lattice, DNA_lattice_coords)

        sim_properties['DNAcoords'] = DNA_lattice_coords

    if np.rint(time) < sim_properties['total_time']:

        # Late-division daughter separation (matches the upstream gating). Runs
        # btree_chromo synchronously once to push the two daughters into opposite
        # lobes. Disable with sim_properties['dna_partition_at_division'] = False.
        if sim_properties.get('dna_partition_at_division', True) \
                and (sim_properties['divH'] > 190) \
                and not checkDaughtersFullyPartitioned(sim_properties):

            if 'partitionedDNA' not in sim_properties.keys():

                partitionChromosomes(sim_properties)

                sim_properties['partitionedDNA'] = True

        writeRiboObstacleFile(time, ribo_site_dict, sim_properties)

        writeChromosomeInputFile(time, sim_properties, True)

        runNewChromosome(time, sim_properties)

        print('Started running new dividing chromosome in background')

    return region_dict, DNA_lattice_coords, genome
#########################################################################################


#########################################################################################
def checkLastChromosome(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    print("Checking for configuration from previous DNA step")
    
    workDir = sim_properties['working_directory']+'DNA/'
    
    step = sim_properties['last_DNA_step']
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(step)
    topo_file = workDir + 'chromo_topo_{:d}.dat'.format(step)
    
    print(DNAfile)
    
    last_DNA_complete = _dna_outputs_ready(workDir, step)
    
    if not last_DNA_complete:
        
        print("Waiting on BRGDNA to complete configuration (monomers + chromo_topo)")
        
    DNA_wait = 0.0
    sleep_interval = 0.5  # Start with 0.5 second checks for faster detection
    wait_timeout = float(sim_properties.get('dna_btree_wait_timeout_s', DNA_BTREE_WAIT_TIMEOUT_S))

    # when this process launched btree for this step, wait for it to exit (returns as soon as it does), then check the
    # outputs as before; with no handle (e.g. after a restart) or no outputs after the exit, the polling loop below runs as before.
    proc = _btree_procs.pop(step, None)
    if proc is not None and not last_DNA_complete:
        t_wait0 = timepy.time()
        try:
            proc.wait(timeout=wait_timeout)
        except subprocess.TimeoutExpired:
            pass
        DNA_wait = timepy.time() - t_wait0
        last_DNA_complete = _dna_outputs_ready(workDir, step)
    
    while not last_DNA_complete:
        
        if DNA_wait >= wait_timeout:
            
            rescueDNA(sim_properties)
            
            print('Created Rescue Files')
            
            return None
        
        timepy.sleep(sleep_interval)
        DNA_wait = DNA_wait + sleep_interval
        
        last_DNA_complete = _dna_outputs_ready(workDir, step)
        
        # Exponential backoff: increase interval gradually, but cap at 2 seconds
        # This allows fast detection when file appears soon, but doesn't waste CPU on long waits
        if DNA_wait < 5.0:
            sleep_interval = 0.5  # First 5 seconds: check every 0.5s
        elif DNA_wait < 15.0:
            sleep_interval = 1.0  # 5-15 seconds: check every 1s
        else:
            sleep_interval = 2.0  # After 15 seconds: check every 2s
        
    print("Waited seconds: "+str(DNA_wait))

    if _lean():
        if proc is not None:
            try:
                proc.wait(timeout=wait_timeout)            # the whole directive list finished (dump, monomers, loops, state)
            except subprocess.TimeoutExpired:
                pass
        _lean_collect(sim_properties, workDir, step)
        
    print("Chromosome configuration ready to load")
        
    return None
#########################################################################################


#########################################################################################
def placeNewChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    timestep = int(time/sim_properties['timestep'])
    
    workDir = sim_properties['working_directory']+'DNA/'
    
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
    fileType = DNAfile.split('.')[1] #'bin'
    
    print(fileType)
    
    N_edges = sim_properties['lattice_edges']
    
    N_edges_x = N_edges[0] # Number of subvolumes making up and edge of the simulation space N x N x N
    N_edges_y = N_edges[1]
    N_edges_z = N_edges[2]

    N_2_x = N_edges_x/2
    N_2_y = N_edges_y/2
    N_2_z = N_edges_z/2
    
    DNAsites = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    DNA_lattice_coords = {}
    
    if fileType=='bin':
        
        with open(DNAfile,'rb') as f:
            
            DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
        
        # Validate DNA file content
        if len(DNAbin) == 0:
            print(f"ERROR: DNA file {DNAfile} is empty (0 bytes)! btree_chromo may have failed.")
            print("Attempting to use rescue DNA configuration...")
            rescueDNA(sim_properties)
            rescueDNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
            if os.path.isfile(rescueDNAfile) and os.path.getsize(rescueDNAfile) > 0:
                print(f"Using rescue DNA file: {rescueDNAfile}")
                with open(rescueDNAfile,'rb') as f:
                    DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
            else:
                raise ValueError(f"DNA file is empty and rescue file also unavailable. Cannot proceed.")
        
        if len(DNAbin) % 3 != 0:
            print(f"ERROR: DNA file {DNAfile} has invalid size: {len(DNAbin)} floats (not divisible by 3)")
            print(f"This indicates file corruption. Expected size divisible by 3 (each particle needs x,y,z coordinates).")
            print("Attempting to use rescue DNA configuration...")
            rescueDNA(sim_properties)
            rescueDNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
            if os.path.isfile(rescueDNAfile) and os.path.getsize(rescueDNAfile) > 0:
                print(f"Using rescue DNA file: {rescueDNAfile}")
                with open(rescueDNAfile,'rb') as f:
                    DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
                # Validate rescue file too
                if len(DNAbin) % 3 != 0:
                    raise ValueError(f"Rescue DNA file also has invalid size: {len(DNAbin)} floats (not divisible by 3)")
            else:
                raise ValueError(f"DNA file has invalid size and rescue file also unavailable. Cannot proceed.")
        
        # Final safety check right before reshape (in case validation somehow didn't catch it)
        if len(DNAbin) % 3 != 0:
            print(f"FATAL: DNA file {DNAfile} has invalid size: {len(DNAbin)} floats (not divisible by 3) - validation was bypassed!")
            print("Attempting to use rescue DNA configuration...")
            rescueDNA(sim_properties)
            rescueDNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
            if os.path.isfile(rescueDNAfile) and os.path.getsize(rescueDNAfile) > 0:
                print(f"Using rescue DNA file: {rescueDNAfile}")
                with open(rescueDNAfile,'rb') as f:
                    DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
                if len(DNAbin) % 3 != 0:
                    raise ValueError(f"Rescue DNA file also has invalid size: {len(DNAbin)} floats (not divisible by 3)")
            else:
                raise ValueError(f"DNA file has invalid size ({len(DNAbin)} floats, not divisible by 3) and rescue file also unavailable. Cannot proceed.")
        
        # Attempt reshape with error handling
        try:
            DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
        except ValueError as e:
            if "cannot reshape" in str(e):
                print(f"ERROR: Failed to reshape DNA file {DNAfile} - {str(e)}")
                print(f"File size: {len(DNAbin)} floats (should be divisible by 3 for x,y,z coordinates)")
                print("Attempting to use rescue DNA configuration...")
                rescueDNA(sim_properties)
                rescueDNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
                if os.path.isfile(rescueDNAfile) and os.path.getsize(rescueDNAfile) > 0:
                    print(f"Using rescue DNA file: {rescueDNAfile}")
                    with open(rescueDNAfile,'rb') as f:
                        DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
                    if len(DNAbin) % 3 != 0:
                        raise ValueError(f"Rescue DNA file also has invalid size: {len(DNAbin)} floats (not divisible by 3)")
                    DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
                else:
                    raise ValueError(f"DNA file reshape failed and rescue file also unavailable. Cannot proceed.")
            else:
                raise  # Re-raise if it's a different ValueError
        
        if len(DNAcoords) == 0:
            print(f"ERROR: DNA file {DNAfile} loaded but contains 0 particles!")
            print("Attempting to use rescue DNA configuration...")
            rescueDNA(sim_properties)
            rescueDNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
            if os.path.isfile(rescueDNAfile) and os.path.getsize(rescueDNAfile) > 0:
                print(f"Using rescue DNA file: {rescueDNAfile}")
                with open(rescueDNAfile,'rb') as f:
                    DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
                if len(DNAbin) % 3 != 0:
                    raise ValueError(f"Rescue DNA file has invalid size: {len(DNAbin)} floats (not divisible by 3)")
                DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
            else:
                raise ValueError(f"DNA file has 0 particles and rescue file also unavailable. Cannot proceed.")
        
        print(DNAcoords.shape)
        # the per-monomer loop as array operations (the same float64 operations element by element; int() truncation =
        # astype(int64)); DNA_lattice_coords gets the same keys and lists of Python ints. The array is kept for moveDnaParticles.
        scale = 10*sim_properties['lattice_spacing']
        lat = np.empty((len(DNAcoords), 3), dtype=np.int64)
        lat[:, 0] = (np.ceil((DNAcoords[:, 0]*1e-9)/scale)+N_2_x).astype(np.int64)
        lat[:, 1] = (np.ceil((DNAcoords[:, 1]*1e-9)/scale)+N_2_y).astype(np.int64)
        lat[:, 2] = (np.ceil((DNAcoords[:, 2]*1e-9)/scale)+N_2_z).astype(np.int64)
        DNAsites[lat[:, 0], lat[:, 1], lat[:, 2]] = True
        DNA_lattice_coords = dict(zip(range(1, len(DNAcoords)+1), lat.tolist()))
        _wcm_coords_cache[0] = (DNA_lattice_coords, lat)

    elif fileType=='xyz':
        
        with open(DNAfile,'rb') as f:
            
            for line_number,line in enumerate(f):
                
                if line_number == 0:
                    
                    DNAparticleCount=int(line)
                    
                elif line_number == 1:
                    
                    continue
                    
                else:
                    
                    particleIdx = int(line_number-2)
                    
                    atomic_symbol, x, y, z = line.split()
                    
                    x_lattice = np.ceil((float(x)*1e-9)/(10*sim_properties['lattice_spacing']))+N_2_x
                    y_lattice = np.ceil((float(y)*1e-9)/(10*sim_properties['lattice_spacing']))+N_2_y
                    z_lattice = np.ceil((float(z)*1e-9)/(10*sim_properties['lattice_spacing']))+N_2_z

                    DNAsites[int(x_lattice),int(y_lattice),int(z_lattice)] = True
                    
                    DNA_lattice_coords[particleIdx+1] = [int(x_lattice),int(y_lattice),int(z_lattice)]
        
#     sim_properties['DNAcoords'] = DNA_lattice_coords

    # the site-type updates of the former per-site loops, written through the site-lattice view (indexed [a3, a2, a1]
    # for API calls (a1, a2, a3); the calls here pass (site[2], site[1], site[0]) of a mask index, so the view is indexed like the
    # masks), with the same precedence per former-DNA site (cytoplasm, outer cytoplasm, membrane, otherwise each ribosome type's
    # cross then centre in dict order, later ones winning) and then the new DNA sites. One real setSiteType call (value unchanged)
    # clears the lattice's GPU-synced flag whenever the former loops made at least one call.
    sv = lattice.getSiteLatticeView()
    A = region_dict['DNA']['shape'] == True
    called = np.zeros(A.shape, dtype=bool)
    if A.any():
        c = A & (region_dict['cytoplasm']['full_shape'] == True)
        o = A & ~c & (region_dict['outer_cytoplasm']['full_shape'] == True)
        m = A & ~c & ~o & (region_dict['membrane']['shape'] == True)
        r = A & ~c & ~o & ~m
        for ribo_type, type_dict in ribo_site_dict.items():
            crossID = type_dict['cross_idx']; centerID = type_dict['center_idx']
            hit = r & (region_dict[crossID]['shape'] == True); sv[hit] = region_dict[crossID]['index']; called |= hit
            hit = r & (region_dict[centerID]['shape'] == True); sv[hit] = region_dict[centerID]['index']; called |= hit
        sv[m] = region_dict['membrane']['index']
        sv[o] = region_dict['outer_cytoplasm']['index']
        sv[c] = region_dict['cytoplasm']['index']
        called |= c | o | m

    region_dict['DNA']['shape'] = DNAsites

    region_dict['cytoplasm']['shape'] = region_dict['cytoplasm']['full_shape'] & ~region_dict["DNA"]["shape"]
    region_dict['outer_cytoplasm']['shape'] = region_dict['outer_cytoplasm']['full_shape'] & ~region_dict["DNA"]["shape"]

    newDNA = region_dict['DNA']['shape'] == True
    sv[newDNA] = region_dict['DNA']['index']
    called |= newDNA
    first = np.argwhere(called)
    if len(first) > 0:
        site = first[0]
        lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), int(sv[site[0], site[1], site[2]]))

    print("Updated DNA Sites")
    
    return region_dict, DNA_lattice_coords
#########################################################################################
    
    
#########################################################################################
def _safe_add_dna_particle(lattice, a, b, c, particleID, dna_sites=None, max_occ=15):
    """Place particleID at lattice site (a, b, c) with occupancy-overflow handling.

    moveDnaParticles was historically the only particle-placement path that called
    lattice.addParticle() with no occupancy check, so a saturated voxel (16 particles)
    raised InvalidParticleException and killed the whole run. This mirrors the
    relocation fallback already used in RibosomesRDME/Growth/Division: if the target
    site is full, place at the nearest DNA lattice site that still has capacity.
    Returns True if the particle was placed.
    """
    a, b, c, particleID = int(a), int(b), int(c), int(particleID)
    try:
        if lattice.getOccupancy(a, b, c) < max_occ:
            lattice.addParticle(a, b, c, particleID)
            return True
    except Exception:
        pass
    if dna_sites is not None and len(dna_sites) > 0:
        d2 = ((dna_sites[:, 0] - a) ** 2 + (dna_sites[:, 1] - b) ** 2 + (dna_sites[:, 2] - c) ** 2)
        for idx in np.argsort(d2)[:1024]:
            ta, tb, tc = int(dna_sites[idx, 0]), int(dna_sites[idx, 1]), int(dna_sites[idx, 2])
            try:
                if lattice.getOccupancy(ta, tb, tc) < max_occ:
                    lattice.addParticle(ta, tb, tc, particleID)
                    return True
            except Exception:
                continue
    print("Warning: moveDnaParticles could not place particle", particleID,
          "near", (a, b, c), "- all nearby DNA sites full")
    return False


def moveDnaParticles(sim_properties, lattice, DNA_lattice_coords):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
                        
    genome = sim_properties['genome']
    OldDnaParticleCoords = sim_properties['DNAcoords']
    
    plattice = lattice.getParticleLatticeView()

    # Candidate DNA sites (in lattice.addParticle arg order = [z, y, x]) for
    # relocating a particle when its target voxel is already at max occupancy.
    _cached = _wcm_coords_cache[0]
    if _cached is not None and _cached[0] is DNA_lattice_coords:   # the same values, from placeNewChromosome's array
        _dna_site_arr = np.ascontiguousarray(_cached[1][:, ::-1]) if len(_cached[1]) else None
    else:
        _dna_site_arr = (np.array([[v[2], v[1], v[0]] for v in DNA_lattice_coords.values()],
                                  dtype=np.int64)
                         if DNA_lattice_coords else None)

    for locusTag, locusDict in genome.items():
        
        locusNum = locusTag.split('_')[1]
        
#         trsc2ID = sim_properties['name_to_index']['RP_' + locusNum + '_t']
        
        starts = locusDict['startIndex']
        ends = locusDict['endIndex']
        
        old_starts = locusDict['prevStartIndex']
        old_ends = locusDict['prevEndIndex']
        
#         print(OldDnaParticleCoords)

        for i in range(len(starts)):
        
            chromo = i + 1
            
            geneID = sim_properties['name_to_index']['G_' + locusNum + '_C' + str(int(chromo))]
        
            start = starts[i]
            
            if (i+1)>len(old_starts):
                
                newStartXYZ = DNA_lattice_coords[start]

                _safe_add_dna_particle(lattice, newStartXYZ[2], newStartXYZ[1], newStartXYZ[0], geneID, _dna_site_arr)

                print('Placed new gene for ', locusTag)
                    
            else:
                
                old_start = old_starts[i]

                oldStartXYZ = OldDnaParticleCoords[old_start]

                newStartXYZ = DNA_lattice_coords[start]
                
                rnasequence = locusDict["RNAsequence"]
                
                if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
                    
                    if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID):

                        deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID)

                        _safe_add_dna_particle(lattice, newStartXYZ[2], newStartXYZ[1], newStartXYZ[0], geneID, _dna_site_arr)

                        continue

                    for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):

                        trscID = sim_properties['name_to_index']['RP_' + locusNum + '_' + str(i) + '_C' + str(int(chromo))]

                        if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID):

                            deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID)

                            _safe_add_dna_particle(lattice, newStartXYZ[2], newStartXYZ[1], newStartXYZ[0], trscID, _dna_site_arr)

                            break

                        trsc2ID = sim_properties['name_to_index']['RP_' + locusNum + '_t_' + str(i) + '_C' + str(int(chromo))]

                        if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID):

                            deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID)

                            _safe_add_dna_particle(lattice, newStartXYZ[2], newStartXYZ[1], newStartXYZ[0], trsc2ID, _dna_site_arr)

                            break
                    
                else:
                    
                    trscID = sim_properties['name_to_index']['RP_' + locusNum + '_C' + str(int(chromo))]

                    if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID):

                        deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID)

                        _safe_add_dna_particle(lattice, newStartXYZ[2], newStartXYZ[1], newStartXYZ[0], geneID, _dna_site_arr)
                        
                        continue

                    if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID):

                        deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID)

                        _safe_add_dna_particle(lattice, newStartXYZ[2], newStartXYZ[1], newStartXYZ[0], trscID, _dna_site_arr)
                        
                        continue
                    
#                 if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID):

#                     deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID)

#                     lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), trsc2ID)

#                 else:

#                     print('Something is wrong with gene: ', locusTag)
        
#     sim_properties['DNAcoords'] = DNA_lattice_coords

    for feature, fdict in sim_properties['chromosome_features'].items():
        
        oldPartXYZ = OldDnaParticleCoords[int(fdict['index'])]
        
        partIdx = sim_properties['name_to_index'][feature]
        
        if checkParticle(plattice, int(oldPartXYZ[2]), int(oldPartXYZ[1]), int(oldPartXYZ[0]), partIdx):
        
            newPartXYZ = DNA_lattice_coords[int(fdict['index'])]

            deleteParticle(plattice, int(oldPartXYZ[2]), int(oldPartXYZ[1]), int(oldPartXYZ[0]), partIdx)

            _safe_add_dna_particle(lattice, newPartXYZ[2], newPartXYZ[1], newPartXYZ[0], partIdx, _dna_site_arr)
        
    print('Placed chemically relevant DNA particles')
        
    return None
#########################################################################################
    
    
#########################################################################################   
def runNewChromosome(time, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    timestep = int(time/sim_properties['timestep'])
    
#     os.system('LAMMPS_NEW_LOAD twistable_BD_OMP')
    
#     headDir = sim_properties['head_directory'] + 'external_programs/'
    
    headDir = sim_properties['dna_software_directory']
    workDir = sim_properties['working_directory']+'DNA/'
    
    DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)
    
#     DNA_executable = headDir + 'btree_chromo/build/apps/program ' + DirectivesFname
#     print(DNA_executable)
    
#     os.system(DNA_executable)

#     if sim_properties['last_DNA_step'] is not None:

#         LAMMPStraj= workDir + 'chromosome_{:d}.lammpstrj'.format(sim_properties['last_DNA_step'])

#         os.system('rm ' + LAMMPStraj)
    
    _btree_build = os.environ.get('BTREE_BUILD_DIR', 'build')
    DNAargs = [headDir + 'btree_chromo/' + _btree_build + '/apps/btree_chromo', DirectivesFname]
    
    # Get DNA GPU ID from environment (set by sbatch script as container GPU 1)
    # Inside Docker, GPUs are renumbered: first GPU = 0 (RDME), second GPU = 1 (DNA)
    dna_gpu_id = os.environ.get("DNA_GPU_ID", "1")
    # Create environment for subprocess with only the DNA GPU visible
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(dna_gpu_id)

    # Optional CPU isolation for btree_chromo via DNA_CPU_CORES + taskset.
    launch_cmd = DNAargs
    dna_cores = os.environ.get("DNA_CPU_CORES", "").strip()
    if dna_cores and shutil.which("taskset"):
        launch_cmd = ["taskset", "-c", dna_cores] + DNAargs
        print("btree_chromo pinned to CPUs {} (taskset)".format(dna_cores))

    _btree_procs[timestep] = _launch_btree(launch_cmd, DirectivesFname, env)
    
    lastStep = sim_properties['last_DNA_step']
    
    sim_properties['last_last_DNA_step'] = lastStep
    
    sim_properties['last_DNA_step'] = timestep
    
    print("Set next DNA configuration to run in background")
    
    return
#########################################################################################

    
#########################################################################################
def writeChromosomeInputFile(time, sim_properties, updateRegions):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
#     headDir = sim_properties['head_directory'] + 'external_programs/'
    headDir = sim_properties['dna_software_directory']
    workDir = sim_properties['working_directory']+'DNA/'
    
    rep_started = sim_properties['rep_started']
    
    cyto_radius_angstroms = int((sim_properties['cyto_radius'])*sim_properties['lattice_spacing']*10/1e-9)
    timestep = int(time/sim_properties['timestep'])
    
    rng_number = int(timestep/10000+sim_properties['dna_rng_seed'])
    
    processor_number = 8 #25
    
    DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)
    
    DnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(timestep)
    
    try:
        PrevDnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    except:
        PrevDnaBinFname = workDir + 'x_chain_Syn3A_chromosome_init_rep00001.bin'
    
    RiboFname = workDir + 'ribo_obstacles_{:d}.bin'.format(timestep)
    
    MemBoundaryFname = workDir + 'mem_boundary_{:d}.bin'.format(timestep)
    
    DnaXyzFname = workDir + 'dna_monomers_{:d}.xyz'.format(timestep)
    
    DnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(timestep)
    
    try:
        PrevDnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(sim_properties['last_DNA_step'])
    except:
        PrevDnaQuatFname = None
    
    with open(DirectivesFname, 'w') as f:

        loop_seed = _loop_prng_seed(sim_properties, timestep)
        f.write('btree_prng_seed:10\n')
        f.write('replicator_prng_seed:10\n')
        f.write('loop_prng_seed:{:d}\n'.format(loop_seed))

#         f.write('new_chromo:54338\n')
        
        if not rep_started:
            f.write('new_chromo:54338\n')
#             sim_properties['rep_started'] = True
        else:
            f.write('input_state:' + workDir + 'rep_state_{:d}.txt\n'.format(sim_properties['last_DNA_step']))

#         f.write('load_BD_lengths:' + headDir + 'btree_chromo/test_case/in_BD_lengths_LAMMPS_test.txt\n')
        f.write('load_BD_lengths:' + sim_properties['head_directory'] + 'input_data/in_BD_lengths_LAMMPS_test.txt\n')
        
        f.write('load_mono_coords:' + PrevDnaBinFname + ',row\n')
        
#         if PrevDnaQuatFname is not None:
#             f.write('load_mono_quats:' + PrevDnaQuatFname + ',row\n')
        
#         if not sim_properties['division_started']:
            
        if not sim_properties['division_started']:
            f.write('spherical_bdry:{:d},0,0,0\n'.format(int(cyto_radius_angstroms - sim_properties['lattice_spacing']*10/1e-9)))
        else:
            # Same replicate hook as pre-division; only the membrane boundary shape changes.
            cyto_radius_angstroms = int((sim_properties['divR']) * 10 - sim_properties['lattice_spacing']*10/1e-9)
            cyto_height_angstroms = int((sim_properties['divH']) * 10)
            f.write('overlapping_spheres_bdry:{:d},{:d},0,0,0,0,0,1\n'.format(cyto_height_angstroms, cyto_radius_angstroms))

#         if (sim_properties['gamma_V'] == 1.0) and sim_properties['division_started']:
#             f.write('spherical_bdry:{:d},0,0,0\n'.format(cyto_radius_angstroms))
            
#         if (sim_properties['gamma_V'] < 1.0) and sim_properties['division_started']:
#             memFname = sim_properties['DNA_membrane_file']
#             f.write('load_bdry_coords:' + memFname + ',row\n')
        
#         f.write('load_bdry_coords:' + MemBoundaryFname + ',row\n')

#         f.write('write_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
            
        f.write('prepare_simulator:' + workDir + 'log_{:d}.log\n'.format(timestep))
        f.write('simulator_set_prng_seed:{:d}\n'.format(rng_number))
        f.write('simulator_set_nProc:{:d}\n'.format(processor_number))
        f.write('simulator_set_DNA_model:' + headDir + 'btree_chromo/' + os.environ.get('DNA_MODEL_DIR_NAME', 'LAMMPS_DNA_model_kk') + '\n')
#         f.write('simulator_set_output_details:' + workDir + ',chromosome_{:d}\n'.format(timestep))
        f.write(_output_details(workDir))

        f.write('simulator_set_delta_t:1.0E+5\n')
        
        # TURN OFF TWISTING #
        f.write('switch_twisting_angles:F\n')
        f.write('switch_ellipsoids:F\n')
        
#         f.write('switch_Ori_bdry_attraction:T\n')
#         f.write('switch_Ori_pair_repulsion:T\n')

        f.write('dump_topology:'+ workDir +'chromo_topo_{:d}.dat,1\n'.format(timestep))

        f.write('sys_write_sim_read_LAMMPS_data:' + _data_path(workDir, timestep, sim_properties) + '\n')

        _write_replicate_hook_protocol(f, sim_properties, timestep, workDir, rep_started)

        f.write('sync_simulator_and_system\n')

        #if not updateRegions:
        #    f.write('load_ribo_coords:' + RiboFname + ',row\n')
        
        #f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
        
        #f.write('simulator_minimize_soft_harmonic:500\n')
        
        #f.write('simulator_run_topoDNA_FENE:100,50,50,append,skip_first\n')
    
        #f.write('sync_simulator_and_system\n')
    
        f.write('write_mono_coords:' + DnaBinFname + ',row\n')

        f.write('write_loops:' + _loops_file(sim_properties, timestep) + '\n')

        f.write('output_state:' + workDir + 'rep_state_{:d}.txt\n'.format(timestep))
        
#         if rep_started:
        
#     sim_properties['last_DNA_step'] = timestep
    
    return None
#########################################################################################
        

#########################################################################################
def writeRiboObstacleFile(time, ribo_site_dict, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    N_edges = sim_properties['lattice_edges']
    
    ribo_center_lattice = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    for ribo_type, type_dict in ribo_site_dict.items():
        
        center_lattice = type_dict['centers']

        ribo_center_lattice = ribo_center_lattice | center_lattice
            
    ribo_center_points = np.argwhere(ribo_center_lattice == True)
    
    ribo_center_points = ribo_center_points - np.array(sim_properties['lattice_center'])
#     print(ribo_center_points)
            
    ribo_center_points = ribo_center_points * (sim_properties['lattice_spacing'] / 1e-9) * 10
    
#     print(ribo_center_points)

#     ribo_cytoplasm_space = int((sim_properties['cyto_radius'])*sim_properties['lattice_spacing']*10/1e-9) - 100

#     ribo_reduced = []

#     for coord in ribo_center_points:
        
#         COMradius = (coord[0]**2 + coord[1]**2 + coord[2]**2)**(1/2)
        
#         if COMradius<ribo_cytoplasm_space:
            
#             ribo_reduced.append(coord)
            
#     ribo_center_points = np.array(ribo_reduced)

    totalRiboCount = len(ribo_center_points)
    
    timestep = int(time/sim_properties['timestep'])
    
    fname = sim_properties['working_directory']+'DNA/' + 'ribo_obstacles_{:d}.bin'.format(timestep)
    
    ribo_center_points = np.reshape(ribo_center_points,(int(len(ribo_center_points)),3),order='F')
    
    with open(fname, 'wb') as f:
        
        ribo_center_points.tofile(f)
        
    return None
#########################################################################################


#########################################################################################
def remapDNA(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """

    genome = sim_properties['genome']
    
    chromosome_length = sim_properties['genome_length']
    
    ori_ter_rotation_factor = int(chromosome_length/2)
    
    workDir = sim_properties['working_directory']+'DNA/'
    RepTopoFname = _ensure_chromo_topo(sim_properties)
    
    repTopoFile = open(RepTopoFname, 'r')
    lines = repTopoFile.readlines()
    
    total_DNA = int(lines[0].split('\n')[0].split('=')[1])
    sim_properties['counts']['chromosome'] = int(total_DNA)
    
    if len(lines)>2:
        
        print('Remapping replicated DNA')
    
        daughterIdxs = lines[2].split('\n')[0].split(',')
        newOriIdx = int(daughterIdxs[3])
        leftArmPoint = int(daughterIdxs[2])
        rightArmPoint = int(daughterIdxs[4])
        
        replicationGoing = True
        
    else:
        
        replicationGoing = False
        
    repTopoFile.close()

    for locusTag, locusDict in genome.items():

        original_start = locusDict['originalStart']
        original_end = locusDict['originalEnd']

        old_starts = locusDict['startIndex']
        old_ends = locusDict['endIndex']

        new_starts = [original_start]
        new_ends = [original_end]
        
        if replicationGoing:

    #         if original_start <= ori_ter_rotation_factor:

            newStart = newOriIdx + original_start - ori_ter_rotation_factor

    #         elif original_start > ori_ter_rotation_factor:

    #             newStart = newOriIdx - original_start

            if leftArmPoint <= newStart <= rightArmPoint:

                new_starts.append(newStart)

    #         if original_end <= ori_ter_rotation_factor:

            newEnd = newOriIdx + original_end - ori_ter_rotation_factor

    #         elif original_end > ori_ter_rotation_factor:

    #             newEnd = newOriIdx - original_end

            if leftArmPoint <= newEnd <= rightArmPoint:

                new_ends.append(newEnd)

        genome[locusTag]['startIndex'] = new_starts
        genome[locusTag]['endIndex'] = new_ends

        genome[locusTag]['prevStartIndex'] = old_starts
        genome[locusTag]['prevEndIndex'] = old_ends
        
#     print(genome)
    
    return genome
#########################################################################################


#########################################################################################
def getReplicatedSegments(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    genomeFile3A = sim_properties['head_directory'] + 'input_data/syn3A.gb'
    genome3A = next(SeqIO.parse(genomeFile3A, "gb"))
    
    genome = sim_properties['genome']
    
    chromosome_length = sim_properties['genome_length']
    
    ori_ter_rotation_factor = int(chromosome_length/2)
    
    workDir = sim_properties['working_directory']+'DNA/'
    RepTopoFname = _ensure_chromo_topo(sim_properties)
    
    repTopoFile = open(RepTopoFname, 'r')
    lines = repTopoFile.readlines()
    
    # If replication is complete, we do not need to advance the forks any further
    if sim_properties['counts']['chromosome'] >= 54338*2:
        return 0,0
    
    # If this is the first instance of replication, the topology will not have a line for the daughter
    # If there is already duaghter information, we use that to determine fork positions
    if len(lines) > 2:
        
        originalIdxs = lines[1].split('\n')[0].split(',')
        daughterIdxs = lines[2].split('\n')[0].split(',')
        originalOri = int(int(originalIdxs[3])*10)
        ccwPoint = int(int(daughterIdxs[1])*10)
        cwPoint = int(int(daughterIdxs[5])*10)
        
        repTopoFile.close()
        
    else:
        
        repTopoFile.close()
        
        repCcw = 40
        repCw = 40
        
        return repCw, repCcw
    
    cwPointIdx = cwPoint - originalOri
    
    print('CW Idx: ', cwPointIdx)
    
    dnasequenceCW = str(genome3A.seq[cwPointIdx+1:cwPointIdx+1+400])
    
    print(dnasequenceCW)
    print(len(dnasequenceCW))
    
    repCw = int(GIP.ReplicationRate(sim_properties, dnasequenceCW))
    
    ccwPointIdx = ccwPoint + originalOri
    
    print('CCW Idx: ', ccwPointIdx)
    
    dnasequenceCCW = str(genome3A.seq[ccwPointIdx-400+1:ccwPointIdx+1])
    
    print(dnasequenceCCW)
    print(len(dnasequenceCCW))
    
    repCcw = int(GIP.ReplicationRate(sim_properties, dnasequenceCCW))
    
    print(repCw, repCcw)
    
    if (repCw == 0) and (repCcw == 0):
        
        return repCw, repCcw
    
    dnaCwReplicated = dnasequenceCW[:(repCw*10)]
    
    dnaCcwReplicated = dnasequenceCCW[:(repCcw*10)]
    
    dNTPcostMap = {'A':['dATP_DNArep_cost', 'dTTP_DNArep_cost'], 'T':['dTTP_DNArep_cost', 'dATP_DNArep_cost'], 
                   'C':['dCTP_DNArep_cost', 'dGTP_DNArep_cost'], 'G':['dGTP_DNArep_cost', 'dCTP_DNArep_cost']}
    
    for base in set(dnaCwReplicated):
        
        baseCount = dnaCwReplicated.count(base)
        
        costIDs = dNTPcostMap[base]
        
        for costID in costIDs:
            
            sim_properties['counts'][costID] = sim_properties['counts'][costID] + baseCount
        
            sim_properties['counts'][costID+'_second'] = sim_properties['counts'][costID+'_second'] + baseCount
        
    sim_properties['counts']['ATP_DNArep_cost'] = sim_properties['counts']['ATP_DNArep_cost'] + int(len(dnaCwReplicated))
    
    sim_properties['counts']['ATP_DNArep_cost_second'] = sim_properties['counts']['ATP_DNArep_cost_second'] + int(len(dnaCwReplicated))
        
    for base in set(dnaCcwReplicated):
        
        baseCount = dnaCcwReplicated.count(base)
        
        costIDs = dNTPcostMap[base]
        
        for costID in costIDs:
            
            sim_properties['counts'][costID] = sim_properties['counts'][costID] + baseCount
        
            sim_properties['counts'][costID+'_second'] = sim_properties['counts'][costID+'_second'] + baseCount
            
    sim_properties['counts']['ATP_DNArep_cost'] = sim_properties['counts']['ATP_DNArep_cost'] + int(len(dnaCcwReplicated))
        
    sim_properties['counts']['ATP_DNArep_cost_second'] = sim_properties['counts']['ATP_DNArep_cost_second'] + int(len(dnaCcwReplicated))
    
    return repCw, repCcw
#########################################################################################


#########################################################################################
def rotateChromosome(time, sim_properties):
    
    workDir = sim_properties['working_directory']+'DNA/'
    
#     timestep = int(time/sim_properties['timestep'])
    
#     DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(timestep)
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
    # Read in monomer coordinates
    with open(DNAfile,'rb') as f:
            
        DNAbin = np.fromfile(f,dtype=np.float64,count=-1)

    DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
    
    chromosome_length = 54338

    # Read replication state topology
    RepTopoFname = _ensure_chromo_topo(sim_properties)

    repTopoFile = open(RepTopoFname, 'r')
    lines = repTopoFile.readlines()

    if len(lines)>2:

        daughterIdxs = lines[2].split('\n')[0].split(',')
        newOriIdx = int(daughterIdxs[3])
        leftArmPoint = int(daughterIdxs[2])
        rightArmPoint = int(daughterIdxs[4])
        
    else:
        return None
        
    # Separate the coordinates for the mother and daughter chromosomes
    mother_monomers = DNAcoords[:54338]
    daughter_monomers = DNAcoords[leftArmPoint-1:rightArmPoint]
    
    # Calculate center of mass of both chromosomes
    mCom = np.average(mother_monomers, axis=0)
    dCom = np.average(daughter_monomers, axis=0)
    
    # Calculate vector between centers of mass
    if mCom[2] > dCom[2]:
        n = mCom - dCom
    else:
        n = dCom - mCom

    # Calculate rotation matrix between center of mass vector and the z axis
    n_z = np.array([0,0,1])

    n_z = n_z/np.linalg.norm(n_z)

    n_norm = n/np.linalg.norm(n)

    axis = np.cross(n_z,n_norm)

    axis = axis/np.linalg.norm(axis)

    theta = np.arccos(np.dot(n_z,n_norm))

    a = np.cos(theta/2)
    b, c, d = axis * np.sin(theta/2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

    rotation_matrix = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    
    sim_properties['rotation_matrix'] = rotation_matrix
    
    # Rotate every monomer for both chromosomes using the rotation matrix
    monos_R = []

    for mono in DNAcoords:

        rotated_mono = np.dot(rotation_matrix,mono)

        monos_R.append(rotated_mono)

    monos_R = np.array(monos_R)
    
    # Update the monomer coordinates file with the rotated chromosomes' state
    DNAfilePre = workDir + 'dna_monomers_{:d}_prerotation.bin'.format(sim_properties['last_DNA_step'])
    DNAfilePost = DNAfile # workDir +  'dna_monomers_{:d}_postrotation.bin'.format(timestep)

    os.system('mv ' + DNAfile + ' ' + DNAfilePre)

    monos_R = np.reshape(monos_R,(int(len(monos_R)),3),order='F')

    with open(DNAfilePost, 'wb') as f:

        monos_R.tofile(f)
        
    sim_properties['rotated_DNA'] = True
        
    return None
        
#########################################################################################


#########################################################################################
def checkDaughtersFullyPartitioned(sim_properties):
    print('Checking DNA Spatial Partitioning Status')

    if sim_properties['last_DNA_step'] is None:
        return False

    workDir = sim_properties['working_directory'] + 'DNA/'
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])

    # Read in monomer coordinates
    with open(DNAfile, 'rb') as f:
        DNAbin = np.fromfile(f, dtype=np.float64, count=-1)

    DNAcoords = DNAbin.reshape((3, DNAbin.shape[0] // 3), order='F').T

    # Read replication state topology
    RepTopoFname = _ensure_chromo_topo(sim_properties)
    with open(RepTopoFname, 'r') as repTopoFile:
        lines = repTopoFile.readlines()

    if len(lines) <= 2:
        return False

    daughterIdxs = lines[2].split('\n')[0].split(',')
    leftArmPoint = int(daughterIdxs[2])
    rightArmPoint = int(daughterIdxs[4])

    # Separate coordinates for the two daughter chromosomes
    daughterA = DNAcoords[:leftArmPoint-1]
    daughterB = DNAcoords[leftArmPoint-1:rightArmPoint]

    # Check z-coordinates for partitioning
    daughterA_z = daughterA[:, 2]
    daughterB_z = daughterB[:, 2]

    condition1 = np.all(daughterA_z < 0) and np.all(daughterB_z > 0)
    condition2 = np.all(daughterA_z > 0) and np.all(daughterB_z < 0)

    return condition1 or condition2
#########################################################################################


#########################################################################################
def rescueDNA(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    print('Rescuing Chromosome from Trapped State')
    print(sim_properties['last_last_DNA_step'])

    workDir = sim_properties['working_directory']+'DNA/'

    last_step = sim_properties['last_DNA_step']
    prev_step = sim_properties['last_last_DNA_step']

    rescueDNAFile = workDir + 'dna_monomers_{:d}.bin'.format(last_step)

    if prev_step is not None:
        # Reuse the previous good DNA step's full state.
        oldDNAFile    = workDir + 'dna_monomers_{:d}.bin'.format(prev_step)
        oldRepState   = workDir + 'rep_state_{:d}.txt'.format(prev_step)
        oldChromoTopo = workDir + 'chromo_topo_{:d}.dat'.format(prev_step)
        oldLoops      = _loops_file(sim_properties, prev_step)
    else:
        # btree failed on the very first DNA hook (no prior step to fall back on):
        # restore the initial pre-replication chromosome. rep_state/loops do not
        # exist yet; chromo_topo is synthesized below for remapDNA.
        oldDNAFile    = workDir + 'x_chain_Syn3A_chromosome_init_rep00001.bin'
        oldRepState   = None
        oldChromoTopo = None
        oldLoops      = None

#     oldQuatFile = workDir + 'dna_quats_{:d}.bin'.format(prev_step)
#     rescueQuatFile = workDir + 'dna_quats_{:d}.bin'.format(last_step)
#     os.system('cp ' + oldQuatFile + ' ' + rescueQuatFile)

    if oldDNAFile and os.path.isfile(oldDNAFile):
        os.system('cp ' + oldDNAFile + ' ' + rescueDNAFile)
    else:
        print('WARNING rescueDNA: no source DNA config found ({}); leaving as-is'.format(oldDNAFile))

    if oldRepState and os.path.isfile(oldRepState):
        os.system('cp ' + oldRepState + ' ' + workDir + 'rep_state_{:d}.txt'.format(last_step))

    if oldChromoTopo and os.path.isfile(oldChromoTopo):
        os.system('cp ' + oldChromoTopo + ' ' + workDir + 'chromo_topo_{:d}.dat'.format(last_step))
    else:
        _write_prereplication_chromo_topo(sim_properties, last_step)

    if oldLoops and os.path.isfile(oldLoops):
        os.system('cp ' + oldLoops + ' ' + _loops_file(sim_properties, last_step))

    return None
#########################################################################################


#########################################################################################
def partitionChromosomes(sim_properties):
#     partitionChromosomes(memcoords, neck_position)
    
    print('PARTITIONING CHROMSOMES')

    workDir = sim_properties['working_directory'] + 'DNA/'
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])

    # Read in monomer coordinates
    with open(DNAfile, 'rb') as f:
        DNAbin = np.fromfile(f, dtype=np.float64, count=-1)

    DNAcoords = DNAbin.reshape((3, DNAbin.shape[0] // 3), order='F').T

    # Read replication state topology
    RepTopoFname = _ensure_chromo_topo(sim_properties)
    with open(RepTopoFname, 'r') as repTopoFile:
        lines = repTopoFile.readlines()

    if len(lines) <= 2:
        return False

    daughterIdxs = lines[2].split('\n')[0].split(',')
    leftArmPoint = int(daughterIdxs[2])
    rightArmPoint = int(daughterIdxs[4])

    # Separate coordinates for the two daughter chromosomes
    C1 = DNAcoords[:leftArmPoint-1]
    C2 = DNAcoords[leftArmPoint-1:rightArmPoint]

#     C1 = DNAcoords[:int(DNAcoords.shape[0]/2)]
#     C2 = DNAcoords[int(DNAcoords.shape[0]/2):]
    print(C1.shape)
    print(C2.shape)
    
    C1_CoM = np.average(C1, axis=0)
    C2_CoM = np.average(C2, axis=0)
    
#     mL_CoM = np.average(memcoords_L, axis=0)
#     mR_CoM = np.average(memcoords_R, axis=0)
    
    mL_CoM = np.array([0,0,-1*sim_properties['divH_Prev']*10])
    mR_CoM = np.array([0,0,sim_properties['divH_Prev']*10])
    
    C1mL = ((C1_CoM[0] - mL_CoM[0])**2 + (C1_CoM[1] - mL_CoM[1])**2 + (C1_CoM[2] - mL_CoM[2])**2)**(0.5)
    C1mR = ((C1_CoM[0] - mR_CoM[0])**2 + (C1_CoM[1] - mR_CoM[1])**2 + (C1_CoM[2] - mR_CoM[2])**2)**(0.5)
    
    if C1mL<=C1mR:
        CL = C1
        CR = C2
        sim_properties['partSides'] = {'m':'L','d':'R'}
    else:
        CL = C2
        CR = C1
        sim_properties['partSides'] = {'m':'R','d':'L'}
        
    max_r_CL = sim_properties['divR']*10 - sim_properties['lattice_spacing']*10/1e-9
    for coord in CL:
        r = ((coord[0] - mL_CoM[0])**2 + (coord[1] - mL_CoM[1])**2 + (coord[2] - mL_CoM[2])**2)**(0.5)
        if r>max_r_CL:
            max_r_CL = float(r)
    print(max_r_CL)
            
    max_r_CR = sim_properties['divR']*10 - sim_properties['lattice_spacing']*10/1e-9
    for coord in CR:
        r = ((coord[0] - mR_CoM[0])**2 + (coord[1] - mR_CoM[1])**2 + (coord[2] - mR_CoM[2])**2)**(0.5)
        if r>max_r_CR:
            max_r_CR = float(r)
            
            
    # Save chromsome coordinates to separate files        
    CL_bin = np.reshape(CL.T,(int(len(CL)*3)),order='F')
    CR_bin = np.reshape(CR.T,(int(len(CR)*3)),order='F')
    
    fnameL = workDir + 'CL.bin'
    with open(fnameL, 'wb') as f:
        CL_bin.tofile(f)
        
    fnameR = workDir + 'CR.bin'
    with open(fnameR, 'wb') as f:
        CR_bin.tofile(f)
    

    # Run C1
    DirectivesFname = writePartitioningChromosomeInputFile(sim_properties, 'L', fnameL, max_r_CL)
    runPartitioningChromosome(DirectivesFname, sim_properties)
#     runChromosome(DirectivesFname, softwareDirectory, working_directory)
    
    # Run C2
    DirectivesFname = writePartitioningChromosomeInputFile(sim_properties, 'R', fnameR, max_r_CR)
    runPartitioningChromosome(DirectivesFname, sim_properties)
#     runChromosome(DirectivesFname, softwareDirectory, working_directory)
    
    # Combine results of running individual chromosomes into a single configuration
    DNAfile = workDir + 'dna_monomers_L.bin'
    with open(DNAfile,'rb') as f:

        DNAbin = np.fromfile(f,dtype=np.float64,count=-1)

        DNAcoordsL = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
        
    DNAfile = workDir + 'dna_monomers_R.bin'
    with open(DNAfile,'rb') as f:

        DNAbin = np.fromfile(f,dtype=np.float64,count=-1)

        DNAcoordsR = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
        
        
    if sim_properties['partSides']['m'] == 'L':
        newDNAcoords = np.concatenate((DNAcoordsL, DNAcoordsR))
    elif sim_properties['partSides']['m'] == 'R':
        newDNAcoords = np.concatenate((DNAcoordsR, DNAcoordsL))
    
    newDNA_bin = np.reshape(newDNAcoords.T,(int(len(newDNAcoords)*3)),order='F')
    
#     DNAfile = working_directory + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
    DNAfileCopy = workDir + 'dna_monomers_{:d}_prepartition.bin'.format(sim_properties['last_DNA_step'])
    
    try:
#         os.popen('mv {} {}'.format(DNAfile,DNAfileCopy))
        os.system('mv {} {}'.format(DNAfile,DNAfileCopy))
    except:
        print('Error moving')
#         print(DNAfile)
    
    print(DNAfile)
    with open(DNAfile, 'wb') as f:
        newDNA_bin.tofile(f)

    return None
#########################################################################################


#########################################################################################
def writePartitioningChromosomeInputFile(sim_properties, chromoID, PrevDnaBinFname, maxR):
    
    headDir = sim_properties['dna_software_directory']
    workDir = sim_properties['working_directory']+'DNA/'
    
    DirectivesFname = workDir + 'partitioning_operations_{:d}_{}.inp'.format(sim_properties['last_DNA_step'], chromoID)
    
    DnaBinFname = workDir + 'dna_monomers_{}.bin'.format(chromoID)
    
#     RiboBinFname = workDir + 'ribosomes_{:d}.bin'.format(timestep)
    
#     PrevDnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(int(timestep-1))
    
#     DnaXyzFname = workDir + 'dna_monomers_{:d}.xyz'.format(timestep)
    
#     DnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(timestep)
    
#     PrevDnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(int(timestep-1))
    
    if chromoID == 'L':
        memZcom = -1*sim_properties['divH_Prev']*10
    elif chromoID == 'R':
        memZcom = sim_properties['divH_Prev']*10

    with open(DirectivesFname, 'w') as f:

        f.write('btree_prng_seed:10\n')
        f.write('replicator_prng_seed:10\n')

        f.write('new_chromo:54338\n')
#         f.write('input_state:' + workDir + 'rep_state.txt\n')

        f.write('load_BD_lengths:' + sim_properties['head_directory'] + 'input_data/in_BD_lengths_LAMMPS_test.txt\n')
    
        f.write('load_mono_coords:' + PrevDnaBinFname + ',row\n')

        f.write('spherical_bdry:{:d},{:d},{:d},{:d}\n'.format(int(maxR), 0, 0, int(memZcom)))
            
        f.write('prepare_simulator:' + workDir + 'log_{}_{:d}.log\n'.format(chromoID,sim_properties['last_DNA_step']))
        f.write('simulator_set_prng_seed:{:d}\n'.format(69))
        f.write('simulator_set_nProc:8\n')
        f.write('simulator_set_DNA_model:' + headDir + 'btree_chromo/' + os.environ.get('DNA_MODEL_DIR_NAME', 'LAMMPS_DNA_model_kk') + '\n')
        f.write('simulator_set_output_details:' + workDir + ',partitioning_{}\n'.format(chromoID))
        f.write('simulator_set_delta_t:1.0E+5\n')
        
#         f.write('simulator_load_loop_params:'+ headDir + 'btree_chromo/test_case/loop_params.txt\n')
        loop_params_path = _write_loop_params_file(sim_properties, _num_smc(sim_properties))
        f.write('simulator_load_loop_params:' + loop_params_path + '\n')

#         f.write('switch_Ori_bdry_attraction:T\n')
#         f.write('switch_Ori_pair_repulsion:T\n')

        # TURN OFF TWISTING #
        f.write('switch_twisting_angles:F\n')
        f.write('switch_ellipsoids:F\n')

        f.write('dump_topology:'+ workDir +'chromo_topo_{}.dat,1\n'.format(chromoID))

        f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{}\n'.format(chromoID))
        
        f.write('simulator_minimize_soft_harmonic:500\n')
        
        f.write('simulator_run_soft_FENE:100,50,50,noappend,first\n')
        
        f.write('sync_simulator_and_system\n')

#         f.write('simulator_relax_progressive:1000,500\n')
    
#         f.write('simulator_relax_progressive:10000,500\n')
        
#         f.write('sync_simulator_and_system\n')

        current_radius = int(maxR)
    
        donePartitioning = False
        
        while not donePartitioning:
            
            f.write('spherical_bdry:{:d},{:d},{:d},{:d}\n'.format(int(current_radius), 0, 0, int(memZcom)))
            
#             f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{}\n'.format(chromoID))

#             f.write('simulator_relax_progressive:1000\n')
            
#             f.write('sync_simulator_and_system\n')

#             f.write('simulator_run_hard_FENE:100,50,50,append,nofirst\n')
            
#             f.write('simulator_expand_bdry_particles:2.0,10,100,10\n')
            
#             f.write('sync_simulator_and_system\n')
            
#             f.write('simulator_relax_progressive:500\n')
            
#             f.write('sync_simulator_and_system\n')
            
#             f.write('simulator_run_hard_FENE:100,50,50,append,nofirst\n')
            
#             f.write('sync_simulator_and_system\n')

            f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{}\n'.format(chromoID))

            f.write('sync_simulator_and_system\n')
        
            f.write('simulator_relax_progressive:1000,500\n')

            f.write('sync_simulator_and_system\n')

            f.write('simulator_run_hard_FENE:10,5,10,append,nofirst\n')
    
            f.write('simulator_expand_bdry_particles:2.0,10,100,10\n')
        
            f.write('sync_simulator_and_system\n')
        
            f.write('simulator_relax_progressive:1000,500\n')

            f.write('sync_simulator_and_system\n')

            f.write('simulator_run_hard_FENE:10,5,10,append,nofirst\n')

            f.write('sync_simulator_and_system\n')
            
            if current_radius<=(sim_properties['divR']*10 - sim_properties['lattice_spacing']*10/1e-9):
                
                donePartitioning = True
                
            else:
                
                current_radius = current_radius - 4*10
            
        
        f.write('simulator_minimize_soft_harmonic:500\n')

        f.write('simulator_run_soft_FENE:100,50,50,append,nofirst\n')

        f.write('sync_simulator_and_system\n')
        
#         f.write('simulator_relax_progressive:1000,500\n')
        
#         f.write('sync_simulator_and_system\n')
        
#         f.write('simulator_run_hard_FENE:10,5,10,append,nofirst\n')
        
#         f.write('sync_simulator_and_system\n')

#         f.write('write_mono_xyz:' + DnaXyzFname + '\n')
        
        f.write('write_mono_coords:' + DnaBinFname + ',row\n')
        
#         f.write('write_mono_quats:' + DnaQuatFname + ',row\n')
        
#         f.write('write_ribo_coords:' + RiboBinFname + ',row\n')
    
    return DirectivesFname
#########################################################################################


#########################################################################################
def runPartitioningChromosome(DirectivesFname, sim_properties):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    headDir = sim_properties['dna_software_directory']
    workDir = sim_properties['working_directory']+'DNA/'
    
#     DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)

#     DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)
    
    DNA_executable = headDir + 'btree_chromo/' + os.environ.get('BTREE_BUILD_DIR', 'build') + '/apps/btree_chromo ' + DirectivesFname
    print(DNA_executable)
    
    os.system(DNA_executable)
    
    return None
#########################################################################################
