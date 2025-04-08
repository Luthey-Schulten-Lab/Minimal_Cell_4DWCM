Scripts that can be used to copy files from the DGX to a local workstation and then combine them into one merged directory of simulation files

The .sh scripts require variables to be changed, like the directories to copy from on the DGX and the username and password for sshpass. The password for sshpass is your password on the DGX.

The jupyter notebook combine_two combines the particle lattices from the original and all restart LM files into one lattice trajectory. Currently this only copies the particle lattice, not the site lattice. SInce we are only analyzing particle lattices, I have not yet written a part to copy the site lattice.

The copy to merged jupyter notebooks take the lists of replicates from multiple sets of simulations and copies them reindexed to a merged directory containing the full set of simulations.

Copying the LM files and merging them requires the environment mentioned in lm_analysis
