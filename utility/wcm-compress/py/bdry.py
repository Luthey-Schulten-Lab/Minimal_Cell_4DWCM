"""Membrane shell + box header of DNA/chromosome.lammpstrj from kept files (see utility/wcm-compress/README.md).
  r_bdry(run_dir)                → BD_l.r_bdry (from the load_BD_lengths file the ops files name; 42.5 fallback)
  directive(run_dir, step)       → the boundary directive of hook <step> (chromosome_operations_<step>.inp; step None → frame 0 init file)
  shell(directive, r_bdry)       → (coords float64 (n,3) in btree order, bbox text lines) via bdry/bdry_gen, cached
  box_lines(run_dir, step, frame0=False) → the three 'ITEM: BOX BOUNDS ff ff ff' value lines of the frame dumped by hook <step>
                                   (btree calc_bbox over shell + the DNA it loaded; see BDRY.md). box_from_bbox(bbox) for a shell alone
  boundary_lines(coords, first_id) → 'id 1 x y z id 1' dump lines ('%g')

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import os, re, subprocess, tempfile, functools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); GEN = os.path.join(HERE, 'bdry', 'bdry_gen')
BDRY_RE = re.compile(r'^((?:spherical|overlapping_spheres|cylindrical|spherocylindrical)_bdry:[^\n]*)$', re.M)

def r_bdry(run_dir):
    run_dir = os.path.abspath(run_dir); D = os.path.join(run_dir, 'DNA')
    ops = sorted(f for f in os.listdir(D) if f.startswith('chromosome_operations_'))
    path = None
    for f in ops[:3]:
        m = re.search(r'^load_BD_lengths:(\S+)$', open(os.path.join(D, f)).read(), re.M)
        if m: path = m.group(1); break
    cands = []
    if path:
        base = os.path.basename(path)
        cands += [os.path.join(run_dir, 'wcm_sidecar', 'input_data', base),                       # lean runs carry it
                  os.path.join(HERE, '..', '..', '..', 'input_data', base), path]                  # this repo's input_data/, the original path
    for c in cands:
        if os.path.isfile(c):
            m = re.search(r'^r_bdry=(\S+)$', open(c).read(), re.M)
            if m: return float(m.group(1))
    return 42.5

def directive(run_dir, step):
    D = os.path.join(run_dir, 'DNA')
    f = os.path.join(D, 'chromosome_operations_%d.inp' % step) if step is not None else next(
        os.path.join(D, x) for x in sorted(os.listdir(D)) if x.endswith('_chromosome_init.inp'))
    m = BDRY_RE.findall(open(f).read())
    return m[-1] if m else None

@functools.lru_cache(maxsize=64)
def shell(directive, r_bdry):
    if not os.path.exists(GEN):
        raise SystemExit('bdry_gen is not built: run utility/wcm-compress/py/bdry/build.sh <btree_chromo source dir> (same compiler as btree_chromo)')
    with tempfile.NamedTemporaryFile(suffix='.bin') as t:
        out = subprocess.run([GEN, directive, repr(r_bdry), t.name], capture_output=True, text=True, check=True).stdout.split('\n')
        co = np.fromfile(t.name, dtype='<f8').reshape(-1, 3)
    assert int(out[0]) == len(co)
    return co, tuple(out[1:4])

def box_from_bbox(bbox):
    return ['%-1.16e %-1.16e' % tuple(float(v) for v in l.split('\t')) for l in bbox]

def _ops(run_dir, step):
    return open(os.path.join(run_dir, 'DNA', 'chromosome_operations_%d.inp' % step)).read()

def loaded_mono(run_dir, step):
    """The float64 bead coordinates btree loaded at the start of hook <step> (its load_mono_coords line, mapped into run_dir)."""
    m = re.search(r'^load_mono_coords:(\S+?),row$', _ops(run_dir, step), re.M)
    return np.fromfile(os.path.join(run_dir, 'DNA', os.path.basename(m.group(1))), dtype='<f8').reshape(-1, 3)

def _g(v): return float('%g' % v)

def box_lines(run_dir, step, frame0=False, rb=None):
    """LAMMPS_sys::calc_bbox (s = 1.2) over the boundary shell and the DNA beads, printed with ostream precision 6 into the data
    file, re-read by LAMMPS, dumped '%-1.16e %-1.16e'. The DNA beads are the coordinates loaded for the hook: at full precision
    for frame 0 (first write, before any sync) and '%g'-rounded for every later frame (the box comes from the write after
    sync_simulator_and_system, which copies LAMMPS's rounded coordinates back). %g rounding is monotonic, so only the extremes
    are rounded."""
    rb = r_bdry(run_dir) if rb is None else rb
    co, _ = shell(directive(run_dir, step), rb)
    m = loaded_mono(run_dir, step)
    dmn, dmx = m.min(0), m.max(0)
    if not frame0: dmn = np.array([_g(v) for v in dmn]); dmx = np.array([_g(v) for v in dmx])
    mn = np.minimum(co.min(0), dmn); mx = np.maximum(co.max(0), dmx)
    mid = 0.5 * (mn + (-1.0) * mx) + mx                    # vqm.v_linterp(0.5, r_min, r_max)
    dr = 0.5 * (mx + (-1.0) * mn)                          # v_xpy(r_max, v_inv(r_min)), v_ax(0.5, .)
    lo = (-1.2) * dr + mid; hi = 1.2 * dr + mid            # v_axpy(-s, dr, r_mid), v_axpy(s, dr, r_mid)
    return ['%-1.16e %-1.16e' % (float('%g' % lo[i]), float('%g' % hi[i])) for i in range(3)]

def boundary_lines(coords, first_id):
    return ['%d 1 %g %g %g %d 1' % (first_id + i, x, y, z, first_id + i) for i, (x, y, z) in enumerate(coords.tolist())]
