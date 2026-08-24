# protein_science branch — integration notes

Brief summary of changes on this branch relative to `main` for coupling 4DWCM to **btree_chromo_gpu `protein_science`** (persistent SMC loops).

## 4DWCM (this repo)

- **`SpatialDnaDynamics.py`** — DNA hooks follow `template_replicate.inp`: `load_loops`, `translocate`, `simulator_form_loops:F`, scaled `simulator_run_soft_harmonic` per 4 s hook; `numSmc` from RDME `P_0415` via `_num_smc()`; replication block when `rep_started`.
- **`MC_RDME_initialization.py`** — defaults: `dna_hook_interval_s=4`, `dna_smc_bound_fraction=0.5`, `dna_loop_translocate_bps=500`, BD wall-time scale for hook 2+.
- **`InitRdmeDna.py`** — `gen_sc_chain` obstacle path fix for init DNA.
- **`Restart_MC_RDME_initialization.py`** — restart falls back to `MinCell.lm` if `MinCell_restart_<time>.lm` is missing.
- **`input_data/loop_params.txt`** — protein_science format; legacy copy in `loop_params_legacy.txt`.

## btree_chromo (external build)

Use **`btree_chromo_gpu` `protein_science`** (current tip). WCM-specific SMC
handling is already upstream: `numSmc` from `loop_params.txt` is authoritative
(no replication-length scaling), `sync_M_to_loaded()` + `set_M(numSmc_initial)`
after `read_state`, and `translocate:...,F` passes `tag_extruded=false`.

Point `-dsd` at the directory containing the built `btree_chromo/`.

## Reference trees

[`_reference/`](_reference/) holds local snapshots of **btree_chromo_gpu** and **Minimal_Cell_ChromosomeSegregation** for comparing directives and templates. Nested `.git` dirs are ignored; re-clone upstream if needed.

## Validation (2026-06)

`test_protein_science2`: 900 s bio completed with zero `Create_bonds` errors, accumulating `rep_state` fork lines and growing `data.lammps` bond counts. `numSmc` tracks live `P_0415` (e.g. 50→51 when lattice count rises slowly); not btree-capped at 51.

Simulation output belongs under `Data/` (gitignored).
