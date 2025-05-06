## README

### Description
Scripts for Task 
1. Print time-dependent concs. and fluxes into PDF files
2. Analyze statistics in Genetic information processes
3. Gibbs free energy change analysis

### Pipeline
1. Use `pkl.sh` to serialize 50 cell replicates' CSV files of multi-omics into a single cell
2. Run `Concs_Fluxes_GIP.ipynb` to do Task 1 and 2
3. For Task 3:
   1. Run `run.sh` to export concs into csv file as import
   2. Run `free_energy_comparison_eQuilibrator.ipynb` to conduct Gibbs free energy change analysis