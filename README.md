# MC4D

Ongoing work on optimizing the 4D whole-cell model for the genetically minimal
cell, JCVI-syn3A. This repository is a **performance / architecture** line built
on the published model
([Minimal_Cell_4DWCM](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM)).

## Two published configurations

| Version | Branch | Wall clock (full cycle) | DNA physics |
|---------|--------|-------------------------|-------------|
| **25 h** | [`main`](https://github.com/luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell/tree/main) | ~25.3 h | Lighter / FENE-style (coupling-optimized) |
| **protein_science** (final) | [`protein_science`](https://github.com/luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell/tree/protein_science) | ~29.4 h (~30 h) | Persistent SMC loops (memory-faithful) |

Pinned dependency commits, checkout commands, and suggested release tags are in
**[`VERSIONS.md`](VERSIONS.md)**. Use that file when publishing or citing the
code for the JCP paper.

---

## `protein_science` branch — Partitioning driven by SMC alone (no fictitious force)

The `protein_science` branch couples the 4DWCM to the updated [btree_chromo_gpu `protein_science` branch](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu/tree/protein_science), which persists SMC loop state across DNA hook steps via `load_loops` / `write_loops` and `translocate` (see [Minimal_Cell_ChromosomeSegregation](https://github.com/a2maytin/Minimal_Cell_ChromosomeSegregation) for standalone examples used for the Protein Science paper Maytin et al. 2026). This allows us to simulate SMC looping with the SMC dwell times that are longer than a single DNA hook interval.

**Requirements**

- Build and install `btree_chromo` from the `protein_science` branch of `btree_chromo_gpu`, with the [`btree_chromo_wcm` patch](btree_chromo_wcm/README.md) applied (not the older `simulator_run_loops` API, which had made the simplification for the Maytin et al. 2026 paper that the SMC number is linearly proportional to the DNA replication and doubles by the end of the cell cycle).
- Use the new `input_data/loop_params.txt` format. The `basal_death_prob` tunes the SMC dwell time, and `numSmc` is the total number of SMC across all chromosomes (and is overwritten by the 4DWCM python wrapper each hook step). (The legacy format is kept as `input_data/loop_params_legacy.txt` on `main`.)

**Parameters Relating to SMC Behavior**

Following Maytin et al. 2026, there are three parameters relating to SMC loop extrusion that can be tuned to allow chromosome segregation to occur without the need for an external force - the SMC complex dwell time, SMC complex translocation speed, and number of actively translocating SMC complexes. The product of these three parameters divided by the chromosome(s) length, can be roughly interpreted as the fraction of the chromosome(s) that are extruded into loops, and must be greater than ~1 to ensure reliable segregation.

- **Translocation speed**: This represents how fast the SMC complexes extrude loops, which here we assume to be **500 bp/s**, or 2000 bp every 4 seconds. At the start of each 4 second hook, both sides of the SMC translocate by 100 beads in the `loop_simulator` object. (In the `chromosome_operations` file, this is written as `translocate:100,T`). Then we re-initialize the LAMMPS system, and write the `loop_simulator` state to the LAMMPS system with `simulator_form_loops:F` (the `F` here means to add bonds between SMC heads only and not add any extra bonds within the loop).
- **Dwell time**: This represents how long on average an SMC complex is bound to the chromosome. We implement this by setting the probability that the SMC unbinds each time the `loop_simulator` translocates by 1 step. By default, we set this probability to `0.0002`, which means on average it takes 5000 steps on the `loop_simulator` before unbinding. This corresponds to a dwell time of **200 s**. 
- **Number of active SMC**: This represents how many SMC complexes are bound to the chromosome(s) and are actively extruding loops. We get the total number of SMC complexes from the counts of Smc protein divided by 2 because it is a homodimer, and then multiply by the bound fraction: `int((P_0415 / 2) * dna_smc_bound_fraction)`. The default `dna_smc_bound_fraction=0.5`. Based on proteomics we assume 200 Smc proteins at the beginning of the simulation, or 100 dimers, or **50 active SMC complexes**. This number on average will double to 100 by the end of the cell cycle, but the growth will vary across different replicates. (For example, in the run I tried, we did not get much gene expression for 0415 and the number increased from 50 to 66 after 105 minutes).

**Code Behavior**

Directive layout follows `Minimal_Cell_ChromosomeSegregation/scripts/template_replicate.inp` and `run_btree_chromo_replicate.py`:

- **`load_loops`** — after `simulator_load_loop_params`, load `DNA/loops/loops_<previous_step>.txt` (skipped on the first DNA hook).
- **Initial `simulator_run_soft_harmonic`** — then optional `translocate:360000,F` equilibration on the first hook only.
- **Hook body** — replication (`map_replication`), `write_loops`, `translocate:N,T`, `simulator_form_loops:F`, `simulator_minimize_topoDNA_harmonic`, then **`run_dynamics`** (`simulator_run_soft_harmonic` scaled from the template’s 2 s batch to the 4 s hook).
- Final **`write_loops`** / **`output_state`** for the next hook.

- `numSmc` in `DNA/loop_params.txt` is set each hook to `int((P_0415 / 2) * dna_smc_bound_fraction)` (default `dna_smc_bound_fraction=0.5`). 
- Hook `translocate` steps default to 100 per 4 s interval (`dna_loop_translocate_bps=500`, i.e. 500 bp/s). 
- Apply the small patch in [`btree_chromo_wcm/`](btree_chromo_wcm/README.md) so btree_chromo uses that count directly in `load_loops` and `map_replication` (no replication-length scaling). Loop snapshots are trimmed to match before each `load_loops`.

Reference snapshots of `btree_chromo_gpu` and `Minimal_Cell_ChromosomeSegregation` live in [`_reference/`](_reference/) (see [`PROTEIN_SCIENCE_NOTES.md`](PROTEIN_SCIENCE_NOTES.md)). Run artifacts go in `Data/` (gitignored).

## Required Programs

- [Lattice Microbes](https://github.com/Luthey-Schulten-Lab/Lattice_Microbes) - [https://github.com/Luthey-Schulten-Lab/Lattice_Microbes](https://github.com/Luthey-Schulten-Lab/Lattice_Microbes)
- [odecell](https://github.com/Luthey-Schulten-Lab/odecell) - [https://github.com/Luthey-Schulten-Lab/odecell](https://github.com/Luthey-Schulten-Lab/odecell) (Install this in the same conda environment as Lattice Microbes AFTER building Lattice Microbes)
- [btree_chromo](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu) - [https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu) (Kokkos-enabled LAMMPS; on the `protein_science` git branch here, use the [`protein_science` branch](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu/tree/protein_science) of btree_chromo_gpu)
- [sc_chain_generation](https://github.com/Luthey-Schulten-Lab/sc_chain_generation) - [https://github.com/Luthey-Schulten-Lab/sc_chain_generation](https://github.com/Luthey-Schulten-Lab/sc_chain_generation)
- [FreeDTS](https://github.com/weria-pezeshkian/FreeDTS) - [https://github.com/weria-pezeshkian/FreeDTS](https://github.com/weria-pezeshkian/FreeDTS) (OPTIONAL)

## Running the Model

The model here is runnbale as-is and does not reuire its own installation. Once you have installed the required programs listed above, before running the model, you must activate the conda environment in which you built Lattice Microbes:

```
conda activate envName
```

Then, you must make sure that your LAMMPS installation for ```btree_chromo``` is in your path. You can check this by running:

```
lammps -h
```

Once your environment is ready, you can now run the model.

The python file ```Whole_Cell_Minimal_Cell.py``` is the main executable for the model. The executable has the following user input variables:

|  Variable  | Shorthand | Description |
|------------|-----------|-------------|
| --outputDir | -od | Name of directory that will be created to store trajectory files. DO NOT INCLUDE A '.' |
| --simTime | -t | Amount of biological time to be simulated in units of seconds. |
| --cudeDevices | -cd | Integer index of the GPU to use for the Lattice Microbes RDME solver. |
| --dnaSoftwareDirectory | -dsd | Directory containing the ```btree_chromo/``` and ```sc_chain_generation/``` directories. |
| --dnaRngSeed | -drs | Integer RNG seed that will be used for the chromosome programs. |
| --workingDirectory | -wd | Directory where the simulation is being run in. Defaults to current directory, no input is needed. Only necessary to change on clusters. |

Example executable:

```
python Whole_Cell_Minimal_Cell.py -od replicate1 -t 1200 -cd 1 -drs 13 -dsd /home/zane/Software/
```

## Restarting a simulation

If you want to run more time for a simulation or you need to restart a simulations because of a crash, you can run the available restart python script ```Restart_Whole_Cell_Minimal_Cell.py```

This python executable has the same input variables as the original script. This will only run if ```outputDir``` is a directory that contains all associated files with a previously run simulation from the main script. For example, if you provided the input ```-od replicate1``` previously, giving that same argument to the restart script will take the simulation state from the replicate1 directory.

The ```simTime``` variable is how much more biological time you want to run. If you have simulated 3600 seconds and want to reach a total of 6000 seconds, you would give an input of 2400 for this variable.

## Performance Optimizations

This version includes several performance optimizations compared to the unoptimized version. The optimizations target the main computational bottlenecks: ribosome placement, file I/O, and CME communication.

### 1. Ribosome Placement Optimizations (`RibosomesRDME.py`)

**Optimization 1: KDTree Caching**
- **What**: Caches the KDTree data structure between function calls, only rebuilding when the cytoplasm shell size changes
- **Why**: KDTree construction is expensive and was being rebuilt every hook step
- **Implementation**: Module-level cache (`_cached_tree`, `_cached_she_coords_size`) with size-based invalidation
- **Benefit**: Eliminates redundant KDTree construction, significantly reducing ribosome placement time

**Optimization 2: Vectorized Ribosome Center Updates**
- **What**: Replaces loop-based updates with vectorized NumPy array indexing
- **Original**: Looped through each ribosome coordinate individually
- **Optimized**: Single vectorized operation: `ribo_site_dict['ribos']['centers'][RIBOcoords[:, 1], RIBOcoords[:, 2], RIBOcoords[:, 3]] = True`
- **Benefit**: Faster processing of ribosome center updates

**Optimization 3: Efficient Site Collection**
- **What**: Collects all ribosome sites in a list first, then stacks once with `np.vstack()`
- **Original**: Used `np.concatenate()` in a try/except loop, causing multiple array copies
- **Optimized**: Collects all sites in a list, then performs single `np.vstack()` operation
- **Benefit**: Reduces memory allocations and array copying operations

**Optimization 4: Inline Relocation Checks**
- **What**: Processes sites with inline checks to avoid creating large intermediate boolean arrays
- **Benefit**: Reduces memory usage during ribosome relocation operations

**Optimization 5: Vectorized Cross-Voxel Updates**
- **What**: Vectorized updates for ribosome cross sites with bounds checking using `np.clip()`
- **Original**: Looped through each coordinate individually to update 6 cross sites
- **Optimized**: Uses `np.clip()` for bounds checking and vectorized array indexing for all 6 directions at once
- **Benefit**: Much faster cross-voxel updates, especially for large numbers of ribosomes

**Additional Improvements:**
- Removed timing/debugging code to reduce overhead
- Better code organization with clear optimization comments

### 2. CSV File I/O Optimizations (`FileSaving.py`)

**Optimization: End-Only Concatenation**
- **What**: Individual CSV files for each timepoint are concatenated only at the end of the simulation
- **Original**: Read + Write + Save at each timepoint
- **Optimized**: Saves individual CSV files per timepoint, concatenates only at simulation end
- **Benefits**:
  - Eliminates periodic I/O spikes during long simulations
  - Reduces I/O overhead during simulation execution
  - Better handling of crash scenarios with automatic recovery
  - Preserved individual files for debugging and verification

**Optimization: Batch CSV Concatenation**
- **What**: New `_concatenate_csv_files_optimized()` function processes files in batches
- **Implementation**:
  - Processes files in batches of 20 to reduce memory usage
  - Batch reads multiple files to reduce I/O overhead
  - More efficient merge operations using pandas
- **Benefit**: Faster concatenation when many files need to be merged

**New Functions Added:**
- `_concatenate_csv_files_optimized()`: Batch processing for CSV concatenation
- `recoverFromCrash()`: Handles orphaned CSV files from crashes during restart
- `finalizeCountsAndFluxes()`: Final cleanup function to ensure all temp files are concatenated

**File Structure:**
- `counts_and_fluxes.csv`: Final concatenated file (created at simulation end)
- `counts_fluxes_temp/`: Individual CSV files for each timepoint (preserved after concatenation)

### 3. CME Communication Optimization (`Communicate.py`)

**Optimization: Bulk Species Count Reading**
- **What**: Reads all species counts from CME in one operation instead of one-by-one
- **Original**: Called `PP.getSpecieTrace()` for each species individually in a loop
- **Optimized**: Reads entire `SpeciesCounts` array at once: `species_data = sim[replicate_key]['SpeciesCounts'][:]`
- **Implementation**: Extracts last timepoint for all species at once using array slicing
- **Benefit**: Significantly faster CME communication, especially for simulations with many species

### 4. Code Quality and Type Safety Improvements

**Type Safety Fixes:**
- Added explicit `int()` casts for particle indices when calling `lattice.addParticle()`
- Prevents `TypeError` issues when numpy integers are passed to C++ bindings
- Ensures compatibility with Lattice Microbes C++ API


### Summary of Performance Benefits

1. **Reduced I/O Overhead**: End-only concatenation eliminates periodic I/O spikes
2. **Faster Ribosome Processing**: Vectorized operations and KDTree caching significantly speed up ribosome placement
3. **Lower Memory Usage**: Batch processing and inline checks reduce memory allocations
4. **Better Crash Recovery**: Automatic handling of orphaned files during restart
5. **Faster CME Communication**: Bulk reading instead of per-species calls reduces communication overhead

These optimizations target the main computational bottlenecks identified in profiling, resulting in improved overall simulation performance while maintaining numerical accuracy and correctness.

## Descriptions of Simulation Files

| File Name | Description |
|-----------|-------------|
|```input_data/``` | Directory containing initial conditions, kinetic parameters, and other data used to initialize the whole-cell model. |
| ``` Communicate.py ``` | Procedures used to communicate counts of molecules between different methodologies. For example, extraction of enzyme counts from RDME to use in rates of ODE reactions. |
| ``` Diffusion.py ``` | Sets the diffudion rules for all molecules that are put into the RDME lattice. |
| ``` Division.py ``` | Updates cell morphology and particle positions on the RDME lattice during cell division. |
| ``` FileSaving.py ``` | Records the cell state and generates restart files. |
| ```FreeDTS_functions.py``` | Defines set of functions used to interpret output files of membrane shapes from FreeDTS. |
| ``` GIP_rates.py ``` | Contains functions for rate equations of genetic information processing reactions. |
| ``` Growth.py  ``` | Updates cell morphology and particle positions on the RDME lattice as the cell grows spherically. |
| ``` Hook.py ``` | IMPORTANT Defines main algorithm used when interrupting the RDME solver. Includes executing all other simulations methods and calling functions to communicate betwwen the methods. |
| ``` ImportInitialConditions.py ``` | Loads initial condition data into initialize the system. |
| ``` InitRdmeDNA.py ``` | Creates an initial condition chromosome configuration using ```sc_chain_generation``` and communicates the structure to the RDME lattice. |
| ``` Integrate.py ``` | Runs the ODE integrator for metabolism. |
| ``` LatticeFunctions.py ``` | Defines functions used in communication to manipulate the RDME lattice. |
| ``` MC_CME.py ``` | Creates and runs the global CME simulation for tRNA charging and transcription. |
| ``` MC_RDME_initialization.py ``` | Initializes RDME simulation including site types and reactions. |
| ``` RegionsAndComplexes.py ``` | Initializes shapes of cell regions (e.g. membrane and cytoplasm) onto the RDME lattice. |
| ```Restart_Hook.py``` | Modified version of ```Hook.py``` that is used for simulations that have been restarted. |
| ```Restart_MC_RDME_initialization.py``` | Constructs a whole-cell simulation from the last recorded cell state of a preexisting simulation. |
| ```Restart_Whole_Cell_Minimal_Cell.py``` | Executable for restarting simulations. |
| ``` RibosomesRDME.py ``` | Updates excluded volume of ribosomes so that the ribosomes lattice sites follow the center of mass particle. |
| ``` Run_CME.py ``` | Executable file for global CME. |
| ``` Rxns_CME.py ``` | Defines set of reactions simulated in the global CME. |
| ``` Rxns_ODE.py ``` | Defines set of reactions simulated using deterministic ODEs. |
| ``` Rxns_RDME.py ``` | Defines set of reactions simulated on the RDME lattice. |
| ``` SpatialDnaDynamics.py ``` | Procedures for running ```btree_chromo``` and updating the chromosome state on the RDME lattice. |

## Simulation Files and Input Data by Model Component

| Model Component | Simulation Files | Input Data |
|---------------|-------------|-------------|
|RDME - Reaction-Diffusion Master Equaion| ``` MC_RDME_initialization.py ```, ``` Rxns_RDME.py ```, ``` Diffusion.py ```, ```RegionsAndComplexes.py```, ```RibosomesRDME.py```, ```InitRdmeDna.py``` | ``` oneParamMulder-local_min.json ```, ``` LargeSubunit.xlsx ```|
|CME - Chemical Master Equation|``` MC_CME.py ```, ``` Run_CME.py ```, ``` Rxns_CME.py ``` | ``` kinetic_params.xlsx ```|
|ODE - Ordinary Differential Equations| ``` Rxns_ODE.py ```, ```integrate.py```| ``` kinetic_params.xlsx ```, ``` protein_metabolites.xlsx ```, ``` Syn3A_updated.xml ```|
|BD (DNA) - Brownian Dynamics of DNA| ``` InitRdmeDNA.py ```, ``` SpatialDnaDynamics.py ```| ``` in_BD_lengths_LAMMPS_test.txt ```, ``` loop_params.txt ```, ``` syn3A.gb ```|


## Connectivity of 4DWCM Python Modules

![plot](./4DWCM_Dependency.png)

