# Published simulation versions (JCP performance paper)

This repository supports **two reproducible stacks**. Both build on the
published 4DWCM science model
([Minimal_Cell_4DWCM](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM);
Thornburg et al.). They differ in **chromosome DNA dynamics** and therefore in
wall-clock cost for a full $7200\,\mathrm{s}$ cell cycle.

| Version | Branch (this repo) | Chromosome physics | Typical wall clock (7200 s bio) | Role in the paper |
|---------|--------------------|--------------------|----------------------------------|-------------------|
| **A — 25 h** | `main` | Lighter DNA / FENE-style coupling (no persistent SMC dwell across hooks) | **~25.3 h** | Stage I: coupling-optimized baseline |
| **B — protein_science** | `protein_science` | Persistent SMC loops (dwell ~200 s ≫ 4 s hook); memory-faithful | **~29.4 h (~30 h)** after CPU isolation + single round-trip | Stage III: **final production** config |

Cite **Minimal_Cell_4DWCM** for the model; cite **this repo @ the tag below** for
the performance / architecture release used in the paper.

---

## Version A — 25 hour (coupling-optimized)

**What it is:** Stage I performance engineering (GPU-resident BD / fused
minimizer path, reduced sync, hook-side opts) with the lighter chromosome
formulation. End-to-end ~46 h → **25.3 h**.

### Pin these commits

| Component | Repo | Branch / ref | Commit (validated) |
|-----------|------|--------------|--------------------|
| 4DWCM Python | [Optimize_4DWCM_Minimal_Cell](https://github.com/luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell) | `main` | `5b61880` *(Jan 25 h pack; bump if you re-validate on newer `main`)* |
| Chromosome BD | [btree_chromo-dev](https://github.com/luthey-schulten-chemistry/btree_chromo-dev) | `btree_chromo_2.0` | `257bf35` |
| RDME / CME | [Lattice_Microbes_2.6](https://github.com/luthey-schulten-chemistry/Lattice_Microbes_2.6) | `master` | `75c0205` |
| Container | local Docker | `4dcell-optimize:cuda128` | image tag on the machine |

### Checkout

```bash
git clone git@github.com:luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell.git
cd Optimize_4DWCM_Minimal_Cell
git checkout main   # or: git checkout 5b61880

# btree (25 h engine)
git clone -b btree_chromo_2.0 \
  git@github.com:luthey-schulten-chemistry/btree_chromo-dev.git btree_chromo_2.0
cd btree_chromo_2.0 && git checkout 257bf35 && cd ..
```

### Ready-to-run pack (this lab)

`/raid/alfiap/4dwcm_sim` — self-contained 25 h setup + `run_4dwcm_25h.sh`
(see that folder’s `README.md`).

### Suggested Git tag (when freezing for the paper)

```text
v25h-jcp   → Optimize_4DWCM_Minimal_Cell @ 5b61880 (or revalidated main SHA)
```

---

## Version B — protein_science (~30 h final)

**What it is:** Same coupling optimizations, plus **memory-faithful** DNA
dynamics (persistent SMC loops; Maytin et al. 2026 physics). Without isolation
this rose to ~37.6 h (mid-cycle host contention); with **CPU isolation + single
in-place data round-trip** → **~29.4 h (~30 h)**. This is the **final** pipeline
in the paper.

### Pin these commits

| Component | Repo | Branch / ref | Commit (validated) |
|-----------|------|--------------|--------------------|
| 4DWCM Python | [Optimize_4DWCM_Minimal_Cell](https://github.com/luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell) | `protein_science` | `5fa6919` (CPU isolation + single round-trip) |
| Chromosome BD | [btree_chromo_gpu](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu) | `protein_science` | `ce1d839` *(includes in-place LAMMPS update + fused CG; push if still ahead of origin)* |
| RDME / CME | [Lattice_Microbes_2.6](https://github.com/luthey-schulten-chemistry/Lattice_Microbes_2.6) | `master` | `75c0205` |
| Container | local Docker | `4dcell-optimize:cuda128` | image tag on the machine |

Also apply / use the WCM patches under `btree_chromo_wcm/` as described in the
`protein_science` section of `README.md`.

### Checkout

```bash
git clone git@github.com:luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell.git
cd Optimize_4DWCM_Minimal_Cell
git checkout protein_science   # or: git checkout 5fa6919

# btree (protein_science engine)
git clone -b protein_science \
  https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu.git btree_chromo_gpu
cd btree_chromo_gpu && git checkout ce1d839 && cd ..
# apply btree_chromo_wcm patches per btree_chromo_wcm/README.md
```

### Launchers (this lab)

- A+B validation / full cycle:  
  `btree_chromo_optmiz/shell_scripts/sbatch_ab_validation.sh`
- DNA-hook cadence scan (optional):  
  `btree_chromo_optmiz/shell_scripts/sbatch_dna_hook_scan.sh`

### Suggested Git tag (when freezing for the paper)

```text
v30h-protein-science-jcp   → Optimize_4DWCM_Minimal_Cell @ 5fa6919
```

---

## How to present this on GitHub / Zenodo

1. **Keep `main` = Version A (25 h)** and **`protein_science` = Version B (~30 h)**.  
   Do not merge B into `main` just for the paper; tags make both citable.
2. Create **two release tags** (commands below) after SHAs are final.
3. **Zenodo:** upload one archive that contains *both* tagged trees (or two Zenodo versions) + this `VERSIONS.md` + SLURM scripts. One DOI can list both configurations.
4. **Paper Code Availability** (draft):

> Coupling-optimized 4DWCM sources are available at
> Optimize_4DWCM_Minimal_Cell. The **25 h** configuration is tag `v25h-jcp`
> (branch `main`, lighter DNA dynamics). The **protein_science** configuration
> (persistent SMC memory; final ~30 h pipeline) is tag
> `v30h-protein-science-jcp` (branch `protein_science`). Exact dependency
> commits are listed in `VERSIONS.md`. The underlying whole-cell model is that
> of Thornburg et al. (Minimal_Cell_4DWCM).

### Tag commands (run when ready — do not retag casually)

```bash
cd /path/to/Optimize_4DWCM_Minimal_Cell

git checkout 5b61880   # or revalidated main SHA for 25 h
git tag -a v25h-jcp -m "JCP: 25 h coupling-optimized (lighter DNA) stack"
git push origin v25h-jcp

git checkout 5fa6919   # protein_science final
git tag -a v30h-protein-science-jcp -m "JCP: protein_science ~30 h final (CPU isolation + single round-trip)"
git push origin v30h-protein-science-jcp
```

Mirror tags / SHAs on `btree_chromo_gpu` / `btree_chromo-dev` when you have push access.

---

## Quick decision guide

| You want… | Use |
|-----------|-----|
| Fastest full cycle / coupling-opt paper numbers (25.3 h) | **Version A** (`main` / `v25h-jcp`) |
| Memory-faithful SMC physics (final story, ~30 h) | **Version B** (`protein_science` / `v30h-protein-science-jcp`) |
| Original published science model (no perf paper opts) | [Minimal_Cell_4DWCM](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM) |

---

## Relationship to the paper’s three stages

```
Baseline (~46 h)
    └─ Stage I opts ──────────────► Version A  (~25.3 h)   [main]
           └─ + protein_science DNA ► ~37.6 h (contention)
                  └─ Stage III isolation ► Version B (~29.4 h) [protein_science]
```
