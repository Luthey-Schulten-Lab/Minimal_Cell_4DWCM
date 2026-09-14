# Modularization plan

Migration plan for restructuring the 4DWCM host pipeline. Written on branch
`modularize`; no code has been changed yet.

## Goals

1. Remove the `Restart_*` duplication, so the restart path cannot silently
   drift from the fresh-start path.
2. Group the 27 flat modules into packages that match the physical layers of
   the model (lattice utilities, model definition, solvers, coupling,
   chromosome, morphology, I/O).
3. Make the expensive stages individually importable and swappable, so a
   stage can be replaced or skipped without editing the hook.

## Non-goals

- No change to biology. Reaction rates, lattice physics, SMC parameters, and
  RNG seeding are out of scope.
- No performance work. The Stage I–III optimizations stay exactly as they are.
- No renaming of the two entry-point scripts (see constraints).

## Constraints

| Constraint | Consequence |
|---|---|
| `Dockerfile` line 236 does `COPY . /src/4d/` | Any new package directory is copied automatically, but `.dockerignore` must not exclude it. |
| Jobs run `python Whole_Cell_Minimal_Cell.py` | Both entry-point filenames must remain at the repo root and keep their CLI. |
| Referenced in `docker/run_example.sh`, `docker/README.md` (lines 76, 90), `README.md` (lines 90, 107) | Those paths must stay valid or the docs and smoke script break. |
| Job 2741 runs from `/raid/alfiap/.../Optimize_4DWCM_Minimal_Cell` | Separate tree (different filesystem and inode); refactoring this repo cannot disturb it. |
| Paper reproducibility | `VERSIONS.md` pins commits for the published runs; the tags must keep pointing at the pre-refactor tree. |

## Current state

27 modules, 12,156 lines, all at the repo root. The first-party dependency
graph is acyclic and already layered:

```
entry points   Whole_Cell_Minimal_Cell, Restart_Whole_Cell_Minimal_Cell
orchestration  Hook, Restart_Hook, MC_RDME_initialization, Restart_MC_RDME_initialization
services       Communicate, FileSaving, Growth, RibosomesRDME, Division,
               SpatialDnaDynamics, MC_CME, Integrate
model          Rxns_CME, Rxns_RDME, Rxns_ODE, RegionsAndComplexes,
               ImportInitialConditions
leaves         LatticeFunctions, GIP_rates, Diffusion, FreeDTS_functions, InitRdmeDna
```

No import cycles, so modules can be moved bottom-up without breaking
intermediate states.

### Findings that motivate the work

**`Hook.py` and `Restart_Hook.py` are 82% identical.** 72 of ~400 lines differ
(103-line raw diff), and both import the same ten modules. Nearly all of the
diff is docstrings, comments, and `except:` → `except Exception:`. Only six
differences are semantic:

1. `self.restart_time = sim_properties['time']`, absent in the fresh path.
2. The fresh path seeds `fluxes`, `rep_started`, `gamma_V`, `next_gamma_V`,
   and `division_started`; the restart path inherits them from the checkpoint.
3. Schedule offsets: `next_*_time = 1.0` versus `restart_time + 1.0`.
4. `time = t` versus `time = t + self.restart_time`.
5. Save gate: `if time > 0.99` versus `if t > 0.99`.
6. The restart path creates `DNA/loops/`; the fresh path sets
   `last_last_DNA_step = None`.

All six collapse into one class with `restart_time` defaulting to `0.0` and a
single `is_restart` flag.

**`round_sig` is defined four times** — in `LatticeFunctions`, `FileSaving`,
`Hook`, and `Restart_Hook`. `FileSaving` both star-imports `LatticeFunctions`
and redefines the function, so the local definition silently shadows the
imported one.

**Six modules use `from LatticeFunctions import *`** (`Communicate`,
`Division`, `FileSaving`, `Growth`, `RibosomesRDME`, `SpatialDnaDynamics`).
The star surface is only four functions — `round_sig`, `deleteParticle`,
`checkParticle`, `getParticlesInSite` — so these convert to explicit imports
mechanically.

**`constructGIP` is duplicated** between `MC_RDME_initialization` and
`Restart_MC_RDME_initialization`. Several other names (`placeRibosomes`,
`runNewChromosome`, `transcription`, `moveParticles`) are defined in two
modules each, but with different signatures and purposes; those are namespace
collisions rather than duplication, and packages disambiguate them.

**Module-level mutable state** lives in `FileSaving` (save worker pool),
`MC_CME` (`_WORKER_PROC`, `_WORKER_LOCK`, `_WORKER_DEAD`), `RibosomesRDME`
(`_cached_tree`, `_cached_she_coords_size`), and `SpatialDnaDynamics`. These
are process-global singletons; moving the modules is safe, but they must not
end up imported twice under two names, which would duplicate the worker
process and the KDTree cache.

## Target structure

```
Whole_Cell_Minimal_Cell.py           thin shim, CLI unchanged
Restart_Whole_Cell_Minimal_Cell.py   thin shim, CLI unchanged
mincell4d/
    __init__.py
    core/        lattice.py, diffusion.py, gip_rates.py
    model/       rxns_cme.py, rxns_rdme.py, rxns_ode.py, regions.py,
                 initial_conditions.py
    solvers/     cme.py, cme_worker.py, ode.py, rdme_init.py
    coupling/    communicate.py, hook.py
    chromosome/  spatial_dna.py, init_rdme_dna.py
    morphology/  growth.py, division.py, ribosomes.py, freedts.py
    io/          filesaving.py
```

## Phases

Each phase is a separate commit that leaves the tree importable.

### Phase 1 — Unify the hooks

Merge `Restart_Hook` into `Hook` as a single class:

```python
class MyOwnSolver:
    def __init__(self, ..., restart_time=0.0, is_restart=False):
```

Fresh-start state seeding moves behind `if not is_restart:`. The save gate
becomes `if (time - self.restart_time) > 0.99`, which reduces to the existing
behaviour in both paths. `DNA/loops/` creation is made unconditional, since
`exist_ok=True` makes it a no-op on the fresh path.

`Restart_Hook.py` becomes a two-line shim re-exporting `MyOwnSolver` so the
restart entry point keeps working unchanged.

This is the highest-value phase and is independent of the package move.

### Phase 2 — Deduplicate the leaves

Delete the three redundant `round_sig` definitions, keeping the one in
`LatticeFunctions`. Replace all six `from LatticeFunctions import *` with
explicit imports of the four names actually used.

### Phase 3 — Create the package, move leaves upward

Create `mincell4d/` with the subpackage skeleton. Move modules bottom-up
(leaves first, entry points last) using `git mv` so history follows. After
each move, rewrite importers and run the smoke test.

### Phase 4 — Unify the RDME initializers

Same treatment as Phase 1 for `MC_RDME_initialization` and its restart twin,
starting with the duplicated `constructGIP`. These are more divergent (462
differing lines) so this phase is deliberately last and may end up only
sharing `constructGIP` rather than fully merging.

### Phase 5 — Entry points become shims

`Whole_Cell_Minimal_Cell.py` keeps its `argparse` block and delegates to
`mincell4d.run`. The CLI, the Docker invocation, and the SLURM scripts are
untouched.

## Verification

After every phase:

1. **Syntax check** every module: `python -m compileall -q .`
2. **Import smoke test**: import each module in dependency order in a fresh
   interpreter and assert no `ImportError`, no circular import, and that
   `MC_CME._WORKER_PROC` and `RibosomesRDME._cached_tree` each resolve to a
   single object identity regardless of import path.
3. **Import-graph check**: confirm the graph is still acyclic.
4. **Entry-point check**: `python Whole_Cell_Minimal_Cell.py --help` and the
   restart equivalent still parse their arguments.

Note that none of this executes biology. A short run (60 s biological time)
is the only way to confirm the hook still couples correctly, and should be
run before this branch is considered mergeable.

## Risks

| Risk | Mitigation |
|---|---|
| Restart path breaks silently; only shows up on a resumed job | Phase 1 keeps `Restart_Hook.py` as a shim and preserves both call signatures. |
| Doubled CME worker or KDTree cache from dual import paths | Smoke test asserts singleton identity. |
| Docker image stops finding the code | Entry points stay at root; `.dockerignore` reviewed before Phase 3. |
| Published-run reproducibility | Tags in `VERSIONS.md` stay on pre-refactor commits; this branch is not merged until after submission. |
