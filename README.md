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

- Unoptimized reference: multi-day (published A100) / ~79 h (B200 baseline)
- This branch (production): approximately **27 hours** on 2 B200 GPUs
- On this branch, chromosome dynamics follow the SMC-mediated segregation framework
of Maytin *et al.* (*Protein Science*, 2026): SMC loops persist across DNA hooks
(dwell time ≫ $4\,\mathrm{s}$ coupling interval).

---

## Dependencies

Install these before running (same conda env for LM + odecell):

| Package | Repository |
|---------|------------|
| Lattice Microbes | https://github.com/Luthey-Schulten-Lab/Lattice_Microbes |
| odecell | https://github.com/Luthey-Schulten-Lab/odecell |
| btree_chromo (Kokkos LAMMPS) | https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu (`protein_science` branch) |
| sc_chain_generation | https://github.com/Luthey-Schulten-Lab/sc_chain_generation |
| FreeDTS (optional) | https://github.com/weria-pezeshkian/FreeDTS |

Use current `protein_science` btree_chromo (SMC count from `loop_params` / WCM
proteome is already upstream). See [`PROTEIN_SCIENCE_NOTES.md`](PROTEIN_SCIENCE_NOTES.md)
for coupling details and [`VERSIONS.md`](VERSIONS.md) for pinned commits.

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
