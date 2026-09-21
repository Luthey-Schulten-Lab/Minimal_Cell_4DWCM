# Docker — Minimal Cell 4DWCM

Public CUDA image with Lattice Microbes, odecell, sc_chain_generation,
Kokkos/LAMMPS, and `btree_chromo` (`protein_science`), plus this repository.

**Requirements to build:** Docker, ~50+ GB free disk, multi-hour compile (no GPU needed to *build*).  
**Requirements to run:** [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) + NVIDIA GPU.

## Build profiles

One `Dockerfile`; pick arch with helpers (or pass `--build-arg` yourself).

| Script | Image tag | Target |
|--------|-----------|--------|
| [`build_ampere.sh`](build_ampere.sh) | `4dwcm:ampere` | **Default** — A100 / many cloud GPUs (`sm_80`) |
| [`build_multi.sh`](build_multi.sh) | `4dwcm:multi` | Wider GPUs (`70;75;80;86;89;90`) — slower/larger |
| [`build_blackwell.sh`](build_blackwell.sh) | `4dwcm:blackwell` | NVIDIA B200 (`sm_100`) |

From the **repository root**:

```bash
chmod +x docker/build_*.sh docker/run_example.sh docker/entrypoint.sh

# Recommended for most users
./docker/build_ampere.sh

# Or
./docker/build_multi.sh
./docker/build_blackwell.sh
```

Equivalent Ampere one-liner:

```bash
docker build -f docker/Dockerfile -t 4dwcm:ampere .
```

### Build-args

| Arg | Default | Meaning |
|-----|---------|---------|
| `CUDA_ARCHITECTURES` | `80` | LM + LAMMPS CUDA arches |
| `GPU_ARCH` | `sm_80` | LAMMPS `-DGPU_ARCH` |
| `KOKKOS_ARCH_FLAGS` | `-DKokkos_ARCH_AMPERE80=yes` | Kokkos GPU arch flags |
| `HOST_ARCH` | `ZEN3` | Kokkos host CPU arch |
| `BTREE_NVCC_ARCH` | `sm_80` | `nvcc -arch` for btree fused CG |
| `BTREE_REF` | `d5dcada7…` | btree_chromo commit (see below) |
| `N_PROC_MAKE` | *(nproc)* | Parallel compile jobs |

Example Blackwell override:

```bash
docker build -f docker/Dockerfile -t 4dwcm:blackwell \
  --build-arg CUDA_ARCHITECTURES=100 \
  --build-arg GPU_ARCH=sm_100 \
  --build-arg KOKKOS_ARCH_FLAGS="-DKokkos_ARCH_BLACKWELL100=yes" \
  --build-arg BTREE_NVCC_ARCH=sm_100 \
  .
```

### Reproducing the published runs

The default `BTREE_REF` carries two corrections.

The first is to the topoisomerase pair model: the soft boundary–DNA
cutoff was sized with the ribosome radius (117 Å) rather than the
boundary radius (217 Å), which let DNA pass through the cell envelope
during replication.

The second is to how btree publishes its output. `write_bin` and
`dump_topology` wrote straight to the final path, so the Python driver
polling for those files could open one mid-write and read only part of
the chromosome — rare, but it killed a production run. They now write to
a temporary name and rename it into place, which is atomic.

Every figure and timing number in the paper comes from the *pre-fix*
commit. To rebuild that exact software stack:

```bash
docker build -f docker/Dockerfile -t 4dwcm:paper \
  --build-arg BTREE_REF=1725d1bb0247df1090baf19e31b3124289ff8964 \
  .
```

They differ in that pair cutoff and in output publication. Replication is
unaffected — both reach 108676 chromosome beads — but the pre-fix build
leaks a small fraction of DNA beads outside the membrane during the
replication window
(peak ~0.8% around frame 650), which is visible when rendering.

## Run

### Smoke test (1 GPU, short bio time)

```bash
mkdir -p Data
./docker/run_example.sh 4dwcm:ampere docker_smoke 60
```

Or:

```bash
docker run --rm --gpus all \
  -e DNA_GPU_ID=0 \
  -v "$PWD/Data:/src/4d/Data" \
  4dwcm:ampere \
  python Whole_Cell_Minimal_Cell.py \
    -od docker_smoke -t 60 -cd 0 -drs 1 -dsd /Software/
```

### Full cycle (2 GPUs, production-style)

Expose two GPUs; RDME on device 0, DNA BD on device 1:

```bash
docker run --rm --gpus '"device=0,1"' \
  -e CUDA_VISIBLE_DEVICES=0,1 \
  -e DNA_GPU_ID=1 \
  -v "$PWD/Data:/src/4d/Data" \
  4dwcm:ampere \
  python Whole_Cell_Minimal_Cell.py \
    -od replicate1 -t 7200 -cd 0 -drs 13 -dsd /Software/
```

Optional: pin btree CPU cores with `-e DNA_CPU_CORES=...` (see `PROTEIN_SCIENCE_NOTES.md`).

Interactive shell:

```bash
docker run --rm -it --gpus all -v "$PWD/Data:/src/4d/Data" 4dwcm:ampere bash
```

## Layout inside the image

| Path | Contents |
|------|----------|
| `/src/4d/` | This repository (WORKDIR) |
| `/Software/Lattice_Microbes/` | Public LM (global T/R matrices) |
| `/Software/odecell/` | Metabolic ODE stack |
| `/Software/sc_chain_generation/` | Init chromosome generator |
| `/Software/btree_chromo/` | `protein_science` BD engine |
| `/Software/LAMMPS/` | Kokkos/CUDA LAMMPS build |

Default `-dsd /Software/` matches the baked layout (`btree_chromo/` and `sc_chain_generation/` as siblings).

## Notes

- First build is long (GCC, OpenMPI, LAMMPS, LM, btree). Prefer Ampere unless you need another arch.
- FreeDTS is not baked in (optional morphology path).
- LM is **unmodified** public [Lattice_Microbes](https://github.com/Luthey-Schulten-Lab/Lattice_Microbes) with large species matrices enabled.
- Chromosome engine: [btree_chromo_gpu `protein_science`](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu/tree/protein_science), pinned to a SHA rather than the branch, since the branch moves.
