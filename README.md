# Minimal Cell 4DWCM

Four-dimensional whole-cell model (4DWCM) of the genetically minimal bacterium
*JCVI-syn3A*.

| Ref | What |
|-----|------|
| **`main`** | Unoptimized published pipeline (Thornburg *et al.*, *Cell*, 2026) |
| **`1.1.0`** (this branch) | Performance-optimized pipeline with persistent SMC-loop chromosome dynamics |

This branch keeps the mechanistic model and improves wall-clock cost by targeting
cross-scale coupling in the hook (chromosome BD wait, redundant data movement,
and host-side bookkeeping).

Typical full-cycle cost ($7200\,\mathrm{s}$ biological time):

- Unoptimized reference: multi-day (published A100) / ~79 h (NVIDIA B200 baseline)
- This branch (production): approximately **27 hours** on 2 NVIDIA B200 GPUs
- On this branch, chromosome dynamics follow the SMC-mediated segregation framework
of Maytin *et al.* (*Protein Science*, 2026): SMC loops persist across DNA hooks
(dwell time ≫ $4\,\mathrm{s}$ coupling interval). **Partitioning is driven by
SMC alone** (no fictitious external force).

Key defaults (tunable; see [`PROTEIN_SCIENCE_NOTES.md`](PROTEIN_SCIENCE_NOTES.md)
and the manuscript SI):

| Knob | Default | Role |
|------|---------|------|
| Translocation | 500 bp/s (`translocate:100,T` per 4 s hook) | Loop extrusion speed |
| Dwell | `basal_death_prob=0.0002` (~200 s) | How long an SMC stays bound |
| Active SMC | `int((P_0415/2)*0.5)` (~50 at *t*=0) | Bound dimers from RDME Smc count |

Build `btree_chromo` from `btree_chromo_gpu` **`btree_chromo_v2.0.0`** and use
`input_data/loop_params.txt` (not the legacy / `simulator_run_loops` path).

## Dependencies

Install these before running (same conda env for LM + odecell):

| Package | Repository |
|---------|------------|
| Lattice Microbes | https://github.com/Luthey-Schulten-Lab/Lattice_Microbes |
| odecell | https://github.com/Luthey-Schulten-Lab/odecell |
| btree_chromo (Kokkos LAMMPS) | https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu (`btree_chromo_v2.0.0` branch) |
| sc_chain_generation | https://github.com/Luthey-Schulten-Lab/sc_chain_generation |

Use current `btree_chromo_v2.0.0` btree_chromo (SMC count from `loop_params` / WCM
proteome is already upstream). See [`PROTEIN_SCIENCE_NOTES.md`](PROTEIN_SCIENCE_NOTES.md)
for coupling details and [`VERSIONS.md`](VERSIONS.md) for pinned commits.

---

## Docker (recommended for new users)

A public CUDA image builds LM + odecell + sc_chain + Kokkos/LAMMPS +
`btree_chromo` (`btree_chromo_v2.0.0`) + this code. See **[`docker/README.md`](docker/README.md)**.

```bash
# From repo root — Ampere (A100 / many cloud GPUs), default
chmod +x docker/build_*.sh docker/run_example.sh
./docker/build_ampere.sh          # → 4dwcm:ampere  (long first build)

# Other GPU profiles
./docker/build_multi.sh           # → 4dwcm:multi
./docker/build_blackwell.sh       # → 4dwcm:blackwell (B200)

# Smoke run (needs NVIDIA Container Toolkit)
mkdir -p Data
./docker/run_example.sh 4dwcm:ampere docker_smoke 60
```

| Helper | Image | Hardware |
|--------|-------|----------|
| `docker/build_ampere.sh` | `4dwcm:ampere` | `sm_80` (default) |
| `docker/build_multi.sh` | `4dwcm:multi` | multi-arch (slower/larger) |
| `docker/build_blackwell.sh` | `4dwcm:blackwell` | `sm_100` (B200) |

Full-cycle dual-GPU and build-arg details: [`docker/README.md`](docker/README.md).

---

## Quick start

```bash
conda activate <lattice-microbes-env>
# ensure LAMMPS / btree_chromo are on PATH
lammps -h

python Whole_Cell_Minimal_Cell.py \
  -od replicate1 \
  -t 7200 \
  -cd 0 \
  -drs 13 \
  -dsd /path/to/dir/containing/btree_chromo/and/sc_chain_generation/
```

| Flag | Meaning |
|------|---------|
| `-od` / `--outputDir` | Output directory name (no `.`) |
| `-t` / `--simTime` | Biological time (s) |
| `-cd` / `--cudeDevices` | GPU index for RDME |
| `-dsd` / `--dnaSoftwareDirectory` | Parent dir of `btree_chromo/` and `sc_chain_generation/` |
| `-drs` / `--dnaRngSeed` | RNG seed for chromosome programs |
| `-wd` / `--workingDirectory` | Working directory (clusters) |

Restart with `Restart_Whole_Cell_Minimal_Cell.py` using the same `-od` (and
`-t` = **additional** biological seconds to run).

Trajectory and intermediate files are written under `Data/`.

---

## Further wall-clock and output-size reductions

On top of 1.1.0, the hook and its I/O were profiled component by component and trimmed without changing the model
(same random streams, same operations in the same order; each change keeps a switch to the previous code path):

- **Ribosomes:** site scans and placement vectorised; the common no-relocation path runs as compiled Cython
  (`processes/ribo_fast.pyx`); the particle lattice is read on the GPU when the solver has not downloaded it.
- **Hook return codes:** hooks that change only site types tell Lattice Microbes to upload the site types alone.
- **ODE / CME:** input tables, SBML and CME models are parsed or built once per process; odecell's generated flux functor
  reads its parameters through a typed view.
- **Chromosome (DNA) steps:** one resident `btree_chromo --serve` process runs every DNA step (per-step processes as a
  fallback); monomer coordinates are converted and DNA site types updated with array operations.
- **Lean output:** the run no longer writes files it can regenerate byte for byte (~35 GB → ~9 GB per cycle, see below).

With the matching Lattice Microbes and `btree_chromo` speed changes (companion pull requests "speed optimizations for
4DWCM" in those repositories; both are required), a full 7200 s cell cycle took **2.05 h** on 2 NVIDIA B200 GPUs.
The GPU code builds for sm_70 and newer (`docker/build_ampere.sh` for A100, `docker/build_multi.sh` for several generations);
on A100s expect roughly 5–7 h per cycle (estimated from the per-component timings, not measured).

| Switch (environment) | Default | Effect when changed |
|---|---|---|
| `WCM_LEAN_OUTPUT=0` | 1 | write the original layout (`DNA/chromosome.lammpstrj`, `DNA/data.lammps_<step>`, `counts_fluxes_temp/`) |
| `WCM_LEAN_SHADOW=1` | 0 | write the lean and the original layout (to verify the regenerators) |
| `WCM_BTREE_SERVE=0` | 1 | one `btree_chromo` process per DNA step instead of the resident server |
| `WCM_SITE_ONLY_UPLOAD_OFF=1` | unset | hooks always ask for a full lattice upload |
| `WCM_RIBO_COMPILED_OFF=1` / `WCM_RIBO_COMPILED_VERIFY=1` | unset | numpy ribosome path / run both paths and compare |
| `WCM_RIBO_DEVICE_VERIFY=1` | unset | compare the GPU-side ribosome lattice reads with the host computation |
| `WCM_ODE_TYPED_PARAMS_OFF=1`, `WCM_ODE_TYPED_FLUX_OFF=1` | unset | odecell's original generated flux code |

---

## Compressing runs (optional): `utility/wcm-compress/`

By default a run writes a **lean output directory (~9 GB)**: it skips the per-hook LAMMPS data files, the DNA trajectory and
the per-second count backups, which are all regenerated byte for byte from what it keeps (`WCM_LEAN_OUTPUT=0` restores the
original ~35 GB layout). Compression on top of that is **opt-in**: if you do nothing, the run stays its output directory.
To compress a run after it finished:

```bash
export WCM_IMAGE=4dwcm:ampere                                     # or run natively (python3 + numpy + zstd)
utility/wcm-compress/pack_run.sh Data/replicate1                   # → Data/replicate1.wcmpack (~4.3 GB)
utility/wcm-compress/extract_run.sh Data/replicate1.wcmpack Data/replicate1   # full ~35 GB layout, byte for byte
```

`pack_run.sh` unpacks the archive it just wrote and compares every file with the run, bit for bit. Only when that check
passes does it delete the run directory (`VERIFIED.json` records the result); if it fails, the archive is removed and the
run is left untouched. To compress every run automatically, append it to the run command:

```bash
python Whole_Cell_Minimal_Cell.py -od replicate1 -t 7200 -cd 0 -drs 13 -dsd /path/to/software/ \
  && utility/wcm-compress/pack_run.sh Data/replicate1
```

(on Slurm: `&& RUNDIR=$PWD/Data/replicate1 sbatch utility/wcm-compress/pack_job.sbatch`, a CPU-only job). Details, timings
and the format: [`utility/wcm-compress/README.md`](utility/wcm-compress/README.md).

---

## What changed in 1.1.0 (performance)

Relative to `main`, this branch includes profile-guided optimizations of the
hook and chromosome coupling, including:

- Vectorized / cache-aware host-side RDME control paths (ribosomes, I/O, CME counts)
- Overlapped chromosome coupling with a fused GPU conjugate-gradient minimizer
- Reduced loop-synchronization frequency where safe
- In-place chromosome coordinate updates (avoid redundant LAMMPS rebuilds)
- Support for persistent SMC-loop DNA dynamics (`protein_science`)

Science validation against Thornburg *et al.* (2026) is described in the
accompanying manuscript supplementary material.

---

## Citation

- **Model:** Thornburg *et al.*, *Cell* (2026).
- **SMC segregation:** Maytin *et al.*, *Protein Science* (2026).
- **This optimized code:** cite branch/tag **`1.1.0`** (and Zenodo DOI when available).

```text
https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM/tree/1.1.0
```

---

## License / contact

See repository license. Questions: corresponding authors listed in the
performance manuscript (Alfia Parvez, Zaida Luthey-Schulten).
