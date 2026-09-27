# modelspec

One read of every 4DWCM input: genes, species, reactions, hardcoded roles and the source lines behind them. The model
browser (JCVI-Syn3A 4DWCM), the perturbation checker and the simulator's `--perturb` option all run on it.

## Build (in the container, which has pandas, openpyxl, Biopython and libsbml)

```bash
docker run --rm -u $(id -u):$(id -g) -e HOME=/tmp -v $PWD:$PWD -w $PWD --entrypoint python 4d-cell-numba:latest -m modelspec build
```

This writes `modelspec/out/model_spec.json` and `4DWCM-GUI.html` at the repo root: a self-contained page that opens from
disk on any machine, no network needed.
Mount the repo at its own path, as above, so the page's editor links point at real files. Rebuild after changing
`input_data/`, `roles.py` or the simulator source. The build fails if a function named in `roles.py` has been renamed or
removed.

## Use (standard library only)

```bash
python3 -m modelspec serve                 # browser on localhost:8765; Save writes perturbations/<name>.yaml
python3 -m modelspec check perturbations/ko_0415.yaml
python3 -m modelspec impact 0415 0645      # static impact of knocking these out
```

Opened straight from disk, the browser downloads the YAML file instead of saving it.

## Perturbation file

```yaml
schema: 1
name: "ko_0415"
knockouts:
  - gene: "JCVISYN3A_0415"
    mode: "full"            # full | expression_only | initial_only
knockdowns:
  - gene: "JCVISYN3A_0001"
    promoter_scale: 0.25
initial_protein_counts:
  P_0609: 5
initial_metabolites_mM:
  M_atp_c: 1.0
medium_mM:
  M_glc__D_e: 0.0
reaction_parameters:
  - reaction: "PGI"
    parameter: "kcatF"      # kcatF, kcatR, Km:<species>, k_atp, k_aa, k_tRNA, k_cat, or a Non-Random-Binding name
    scale: 0.1              # or value:
disabled_reactions:
  - "ATPase"
```

What the knockout modes do:

- `full` sets the promoter strength to 0 and starts the protein and mRNA at 0.
- `expression_only` blocks transcription but keeps the protein present at t=0. The model has no protein degradation, so
  that protein only dilutes with growth and division.
- `initial_only` starts the protein and mRNA at 0 but leaves the promoter alone, so the protein is made again.

The promoter strength has a floor of 45, so zeroing an initial count is never enough on its own to knock a gene out.

## Running a perturbation

```bash
python Whole_Cell_Minimal_Cell.py -od run1 -t 60 -cd 0 -drs 1 -dsd /Software/ -p perturbations/ko_0415.yaml
```

Through Slurm, from a private snapshot of this working tree (nothing a run writes lands in the repo):

```bash
slurm/submit_run.sh -b <batch> -n <run name> -t <bio seconds> -g <gpu pair, e.g. 2,3> [-p perturbations/tests/q_ko_rnap.yaml]
```

It copies the sources, the spec and the perturbation files to `/raid/racda/4dwcm_tests/<batch>/_src/<run>/` with a
PROVENANCE.txt and the working-tree diff, writes the run to `<batch>/<run>/` and the logs to `<batch>/logs/`, and uses the
ERA LM and btree builds that the wcm-speed Python needs (the image's public LM makes the run die at the second hook).

`modelspec/apply.py` validates the file against the spec, applies it while the model is built, then reads every edited
value back from the jLM simulation object and writes `perturbation_applied.json` (log + checks) and `perturbation.yaml`
into the run directory. A perturbed run compiles its ODE code in `<run>/ode_build/`: odecell compiles the rate constants
into it, and in a shared working directory runs would load each other's build. Without `-p` nothing of it runs. The restart
script does not take `-p` yet.

`python -m modelspec.tests.test_apply` (container) is the construction-level gate: an empty file builds exactly the
unperturbed model, and each edit type changes exactly what it names. `python3 -m modelspec.compare_runs <runs...>` prints
the key observables of short runs side by side.

## What the browser covers

- **Spreadsheet reactions** (ODE metabolism, CME tRNA charging): parameters, rate law with a symbol legend, enzyme rule.
- **RDME and CME reactions built in Python** (3,556 spatial + 1,012 transcription-elongation reactions): recorded by running
  the simulator's builders on a stand-in Sim (`rdme.py`), grouped into families (one per reaction type, one instance per
  gene). Each family page shows the template, lattice regions, the GIP rate formula with its constants, and the value per
  gene; each gene page lists its own reactions and rates.
- **Constants tab**: jLM rate constants, diffusion coefficients (with the region-transition profile of every species), and
  the constants inside the GIP rate formulas.

Perturbation sections for these: `rdme_rate_constants`, `gene_rate_scales` (gene `"*"` = every gene), `diffusion_scales`,
`gip_constants` (see the docstring of `perturbation.py`). `python -m modelspec build --gui-only` re-inlines the page from
the last spec without the container, for GUI-only changes.

## Files

| file | contents |
|---|---|
| `load.py` | input parsing; each rule names the simulator function it mirrors |
| `roles.py` | roles written in Python rather than read from input_data: RNAP and degradosome subunits, DnaA, replisome, SecY, SMC, ribosome assembly, process dependencies |
| `coderefs.py` | source index: literal species ids, per-gene name patterns, sheet reads, function anchors |
| `impact.py` | static first-order impact: which reactions lose their enzyme, which processes lose a part |
| `rdme.py` | recording stand-in for the jLM Sim; runs the RDME/CME builders and groups the result into reaction families |
| `names.py` | readable names for the BiGG-style reaction ids (the SBML carries none) |
| `perturbation.py` | schema, validation, resolution into explicit edits, YAML read/write |
| `apply.py` | applies a perturbation inside the simulator (`Whole_Cell_Minimal_Cell.py -p`) and reads the values back |
| `compare_runs.py` | side-by-side observables of short runs |
| `gui/` | the browser; `app.js` ports `impact.evaluate` and checks itself against the Python results on load |
| `tests/test_consistency.py` | runs the simulator's own builders on a recording stand-in for the jLM Sim and compares them with the spec |

Run the consistency test in the container:

```bash
docker run --rm -v $PWD:$PWD -w $PWD --entrypoint python 4d-cell-numba:latest -m modelspec.tests.test_consistency
```
