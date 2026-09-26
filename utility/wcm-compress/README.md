# wcm-compress — lossless packing and unpacking of 4DWCM runs

A full 4DWCM cell cycle in the original layout is about 35 GB. Most of it is redundant: the per-hook LAMMPS data files and
the DNA trajectory are printed copies of `DNA/dna_monomers_<step>.bin`, and `counts_fluxes_temp/` duplicates
`counts_and_fluxes.csv`. This branch therefore writes a lean layout by default (~9 GB, see below), and `wcm-compress`
takes either layout down to ~4–5 GB and back, byte for byte. `wcm-compress` replaces a run with an archive about 7× smaller and restores it byte for byte.

| Script | What it does |
|---|---|
| `pack_run.sh <run_dir> [archive_dir]` | Compresses a run and **proves** the archive: it unpacks it into a scratch directory and checks every file (sha256, size, mode, mtime), directory and link against the run, bit for bit. **Pass:** writes `<archive>/VERIFIED.json` and deletes the run directory. **Fail (or any error):** deletes the archive and leaves the run directory untouched. No reference run is needed: the check compares the archive with the run it was made from. |
| `extract_run.sh <archive_dir> <out_dir>` | Restores the run directory (`CME/`, `counts_fluxes_temp/`, `DNA/`, `fluxes/`, `restart_files/`, the CSV files, `MinCell.lm`, `sim_args.txt`, `sim_properties.pkl`, …). |
| `pack_job.sbatch` | `pack_run.sh` as a CPU-only Slurm job, meant to be submitted at the end of a run job so the GPUs are released right away. |

## Quick start

```bash
# in the Docker image (recommended: same tools as the simulation) ...
export WCM_IMAGE=4dwcm:ampere            # or 4dwcm:blackwell, 4dwcm:multi
# ... or natively, in an environment with python3 + numpy and the zstd CLI (>= 1.4): leave WCM_IMAGE unset

utility/wcm-compress/pack_run.sh Data/replicate1                 # → Data/replicate1.wcmpack (VERIFIED.json); Data/replicate1/ removed
utility/wcm-compress/extract_run.sh Data/replicate1.wcmpack /scratch/replicate1
```

Optional environment: `CORES` (cpuset to pin to, e.g. `0-27`), `JOBS` (worker processes; default all usable cores).

From a Slurm run script, after the simulation finished:

```bash
RUNDIR=$PWD/Data/replicate1 WCM_IMAGE=4dwcm:ampere sbatch utility/wcm-compress/pack_job.sbatch
```

## Compressing automatically after every run (opt-in)

`Whole_Cell_Minimal_Cell.py` does not compress anything by itself: if you do nothing, a run stays its output directory
(~9 GB with the default lean output, ~35 GB with `WCM_LEAN_OUTPUT=0`).
To compress every run as soon as it finishes, add one line after the simulation in your run or job script:

```bash
python Whole_Cell_Minimal_Cell.py -od replicate1 -t 7200 -cd 0 -drs 13 -dsd /path/to/software/ \
  && utility/wcm-compress/pack_run.sh Data/replicate1                   # same node, ~10-15 min, GPU idle meanwhile
# or, on Slurm, hand it to a CPU-only job so the GPUs are released immediately:
#   && RUNDIR=$PWD/Data/replicate1 sbatch utility/wcm-compress/pack_job.sbatch
```

`&&` runs it only when the simulation exited cleanly. If the bit-for-bit check fails, the run directory is kept as it was.

## How much, how fast

Measured on full 7200 s cycles, 28 threads (NVIDIA B200 node):

| Run | Raw | Archive | Ratio | Pack (incl. verification) |
|---|---|---|---|---|
| `modularize` production run (seed 48) | 35.3 GB | 5.08 GB | 6.9× | 15 min |
| this branch with `WCM_LEAN_OUTPUT=0` (seeds 48, 49, 50) | 35.7 GB | 4.56–4.63 GB | 7.8× | 10–15 min |
| lean-output run (default of this branch, see below) | 9.1 GB | 4.27 GB | 2.1× | 7 min |

Unpacking takes 1–3 minutes. For a lean run, `extract_run.sh` then regenerates the full layout: the trajectory in ~1.5 min,
the 1800 data files by `btree_chromo` replay in ~15 min on 3 GPUs. That replay is CPU-heavy (`MAX_GPUS` × `JOBS_PER_GPU`
btree_chromo instances): do not run it on the same CPU socket as a simulation you are timing. The lower bound is the chromosome coordinates themselves:
`dna_monomers_*.bin` are float64 Brownian-dynamics positions whose low mantissa bits are noise, so they stay at
~3.2 GB (1.3×) in any lossless scheme.

## What the archive stores (`py/wcmpack.py`)

Each handler claims the files it can reproduce, **checks its own reconstruction at pack time**, and leaves anything it
cannot reproduce to the generic handler (raw bytes, zstd). An unusual run therefore packs less well, never lossily.

| Files | Stored as |
|---|---|
| `DNA/dna_monomers_<step>.bin` | x/y/z planes; order-preserving integer keys of the float64 bits, second differences along the chain, byte planes, zstd |
| `DNA/data.lammps_<step>`, `DNA/chromosome.lammpstrj` | coordinates dropped (they are `%g` of the monomer files) and rebuilt; the remaining text (ids, types, membrane shell, bonds, angles) de-duplicated |
| `MinCell.lm` (HDF5) | every deflate chunk whose zlib re-compression reproduces its stored bytes is stored decoded — particle-lattice chunks as per-site occupancy plus frequency-ranked species bytes — and re-deflated on unpack; the rest of the file is kept byte for byte |
| `counts_fluxes_temp/*.csv` | rebuilt from the columns of `counts_and_fluxes.csv` |
| everything else | concatenated, zstd `-19 --long=31` |

## Lean-output runs

By default (`WCM_LEAN_OUTPUT=1`) a run of this branch does not write `DNA/chromosome.lammpstrj`,
`DNA/data.lammps_<step>` or `counts_fluxes_temp/*.csv` at all: it produces ~9 GB instead of ~35 GB, plus a small
`wcm_sidecar/` (~20 MB: the frames no kept file describes, the first hook's data file, and the exact `btree_chromo` build).
`WCM_LEAN_OUTPUT=0` restores the original layout; `WCM_LEAN_SHADOW=1` writes both (used to verify the regenerators).
`extract_run.sh` notices the sidecar and regenerates those files after unpacking:

| File | Regenerated by | From |
|---|---|---|
| `DNA/chromosome.lammpstrj` | `py/mono2lammpstrj.py` | monomer coordinates (`%g`), atom types from `chromo_topo_<step>.dat` (btree `prepare_types`) and `loops/loops_<step>.txt`, the membrane shell from btree's own `boundary_surface.cpp` (`py/bdry/bdry_gen`), frame 0 from the sidecar |
| `DNA/data.lammps_<step>` | `py/data_replay.py` | `btree_chromo` replaying each hook's own `chromosome_operations_<step>.inp` up to its data write (needs a GPU; split over up to `MAX_GPUS`) |
| `counts_fluxes_temp/*.csv` | `py/cftemp_regen.py` | `counts_and_fluxes.csv` |

Verification: a full cycle run with `WCM_LEAN_SHADOW=1` (both layouts written) regenerates the trajectory byte for byte
(1801 frames, 8.95 GB), 1799/1800 data files by replay plus the first one from the sidecar, and 7201/7201 counts files.
(`py/compare_trees.py` is the test tool used for that check; packing and extracting never need a reference run —
`pack_run.sh` checks the archive against the run it just packed.)

`mono2lammpstrj.py` is also a visualization tool: `--stride 10 --no-boundary` writes every tenth frame of the chromosome
only (0.7 GB, ~2 s), `--first/--last` a time window; OVITO, VMD and LAMMPS `rerun` read the output.

Build `bdry_gen` once, with the same compiler as `btree_chromo` (inside the image when you use one):

```bash
WCM_IMAGE=4dwcm:ampere utility/wcm-compress/py/bdry/build.sh /path/to/btree_chromo_gpu
```

## Files

| Path | Role |
|---|---|
| `common.sh` | shared settings: `WCM_IMAGE`, `CORES`, `JOBS`, `WCM_GPU_LOCKDIR` |
| `py/wcmpack.py` | the codec (`pack` / `unpack`) |
| `py/manifest.py`, `py/gate.py` | byte manifest of a tree; lossless check against a manifest |
| `py/mono2lammpstrj.py`, `py/bdry.py`, `py/bdry/` | trajectory regenerator and membrane-shell helper |
| `py/data_replay.py` | LAMMPS data file regenerator (btree_chromo replay) |
| `py/cftemp_regen.py` | `counts_fluxes_temp/` regenerator |
| `py/finish_extract.py` | completeness check, modes/mtimes, sidecar removal after extraction |
| `py/compare_trees.py` | test tool only: compare an extracted run with a run that wrote both layouts (not needed to pack or extract) |
