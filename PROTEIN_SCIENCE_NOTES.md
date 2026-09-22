# protein_science chromosome coupling

Partitioning is driven by **SMC loop extrusion alone** (no fictitious external
force). This branch couples the 4DWCM to
[`btree_chromo_gpu`](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu)
branch **`btree_chromo_v2.0.0`**, which persists SMC loop state across DNA hooks
via `load_loops` / `write_loops` and `translocate`. That lets the model use SMC
dwell times much longer than a single DNA hook interval (~4 s). The loop
persistence originates in Maytin *et al.*'s `protein_science` branch;
`btree_chromo_v2.0.0` carries it plus our boundary and output-publication fixes.

Standalone examples of the same chromosome physics appear in
[Minimal_Cell_ChromosomeSegregation](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_ChromosomeSegregation)
(Maytin *et al.*, *Protein Science*, 2026). Science detail and the default
parameter table are also in the JCP manuscript **SI**.

## Requirements

1. Build/install **`btree_chromo`** from `btree_chromo_gpu` **`btree_chromo_v2.0.0`**
   (current tip). WCM `numSmc` handling is already upstream — do **not** need a
   separate `btree_chromo_wcm` patch. Do **not** use the older
   `simulator_run_loops` API (that simplified SMC number ∝ replicated DNA length
   for the Maytin *et al.* 2026 paper).
2. Use **`input_data/loop_params.txt`** (protein_science format).
   `basal_death_prob` sets dwell time; `numSmc` is the total SMC count across
   chromosomes and is **overwritten each DNA hook** by the 4DWCM Python wrapper
   from RDME Smc counts. Legacy format:
   `input_data/loop_params_legacy.txt` (on `main`).

Point `-dsd` at the parent directory that contains the built `btree_chromo/`.

## Parameters (Maytin *et al.* 2026)

Force-free segregation is controlled by three quantities. Their product divided
by chromosome length is roughly the extruded-loop fraction and should be
**≳ 1** for reliable segregation:

| Parameter | Default in this pipeline | Notes |
|-----------|--------------------------|--------|
| Translocation speed | `dna_loop_translocate_bps = 500` bp/s | Per 4 s hook: both SMC sides advance 100 beads → `translocate:100,T` in `chromosome_operations_*.inp`, then `simulator_form_loops:F` (SMC-head bonds only). |
| Dwell time | `basal_death_prob = 0.0002` in `loop_params.txt` | Unbind probability per 1-bead loop-simulator step → mean 5000 steps ≈ **200 s**. |
| Active SMC number | `int((P_0415 / 2) * dna_smc_bound_fraction)` | Homodimer; default `dna_smc_bound_fraction = 0.5`. With ~200 Smc proteins at *t* = 0 → ~50 active complexes; expression typically increases this over the cycle (replicate-dependent). |

Other knobs set in `MC_RDME_initialization.py`: `dna_hook_interval_s` (default
4.0; override with `DNA_HOOK_INTERVAL_SEC`), BD wall-time scale, soft-harmonic
warmup steps. Protocol implementation: `SpatialDnaDynamics.py`.

## Related code

- `SpatialDnaDynamics.py` — per-hook replicate protocol
- `MC_RDME_initialization.py` — DNA/SMC defaults
- `input_data/loop_params.txt` — btree_chromo loop parameters
