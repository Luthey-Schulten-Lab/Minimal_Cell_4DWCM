# btree_chromo patch: SMC count from `loop_params`, not replication length

Stock [btree_chromo_gpu `protein_science`](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu/tree/protein_science) scales the loop-extruder count when loading loops or mapping replication:

```cpp
M = M_init + M_init * (N_mono - N_init) / N_init;
```

That couples SMC number to chromosome length. For whole-cell integration (and other cases where an external model sets the active SMC count), `numSmc` in `loop_params.txt` should be authoritative.

## Apply

```bash
cd /path/to/btree_chromo_gpu   # protein_science branch
patch -p1 < /path/to/Optimize_4DWCM_Minimal_Cell/btree_chromo_wcm/loop_simulator_wcm.patch
patch -p1 < /path/to/Optimize_4DWCM_Minimal_Cell/btree_chromo_wcm/btree_driver_wcm.patch
# rebuild btree_chromo
```

`loop_simulator_wcm.patch` adds `numSmc_initial` / `get_numSmc_initial()` (required by the driver patch) and fixes `read_loop_params` to call `set_M` / `set_N` so `numSmc` from `loop_params.txt` resizes the loop arrays (avoiding a segfault when WCM sets `numSmc` > 50).

After applying both patches, also add in `loop_simulator` (if not already in your tree): **`sync_M_to_loaded()`** sets `M` from the loaded loops file size before **`set_M(numSmc_initial)`** in `load_loops`, so extra SMC slots are created via `birth()` instead of leaving zero indices (fixes LAMMPS `create_bonds` errors when the snapshot has 50 rows but `numSmc` is 51).

The driver patch also fixes `translocate:...,F` to pass `tag_extruded=false` (upstream mistakenly passed `true` in both branches).

## Behavior after patch

| Command | Change |
|---------|--------|
| `load_loops:<file>` | `set_N` → `read_state` → `set_M(numSmc)` from `loop_params` (no replication scaling) |
| `map_replication` | Updates monomer mapping and `set_N` only (no SMC scaling) |

No new directive names — existing `template_replicate.inp` and 4DWCM inputs work unchanged.

The 4DWCM sets `numSmc` each hook from `int((P_0415 / 2) * dna_smc_bound_fraction)` and trims `DNA/loops/loops_<step>.txt` to that count before `load_loops`.
