# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np
import os

import importlib
import sys

import WCM_analysis
importlib.reload(WCM_analysis)

in_dir = './test_data/'
# in_dir = '/home/ben/Data/WCM_runs/current_runs/'
in_label = 'counts_and_fluxes'
min_rep = 1
max_rep = 8

reps = np.arange(min_rep,max_rep+1,dtype=np.int32)
# reps = np.array([1,2,4,5,6,7],dtype=np.int32)

merging_required = False
out_dir = './test_data/'
out_label = 'ensemble_test'


fig_dir = './correlations_fig/'
fig_label = 'test'

###################################
# pre-processing trajectory files #
###################################

if merging_required == True:

    w = WCM_analysis.WCM_ensemble()

    # this assumes trajectories are in the below format
    # (in_dir)/(in_label).(replicate).csv
    w.set_traj_files(in_dir,in_label,reps)

    w.load_trajs()

    w.merge_trajs()

    w.write_merged_ensemble(out_dir,out_label)
    
w = WCM_analysis.WCM_ensemble()

w.read_merged_ensemble(out_dir,out_label)


#####################
# species selection #
#####################

# get the gene information
w.get_filtered_locus_tags('./gene_types.txt')






#########################
# key metabolite levels #
#########################

conc_species = ['M_atp_c','M_pep_c','M_fdp_c','M_pi_c','M_ppi_c']

################
# Key Proteins #
################

group_labels_count = []
grouped_species = dict()

temp_proteins = [445,220,131,727,607,451,606,729,213,221]
count_species = ['P_{:04d}'.format(x) for x in temp_proteins]

for test_specie in count_species:
    temp_label = test_specie
    group_labels_count.append(temp_label)
    grouped_species[temp_label] = dict()
    
    grouped_species[temp_label]['ref_specie'] = None
    grouped_species[temp_label]['species'] = [test_specie]
    grouped_species[temp_label]['weights'] = [1]

print(grouped_species)

#################
# key reactions #
#################

group_labels_rxn = []
grouped_rxns = dict()

temp_label = 'CELL_DEATH'
group_labels_rxn.append(temp_label)
grouped_rxns[temp_label] = dict()
grouped_rxns[temp_label]['rxns'] = ['GLCpts0','PYK','PYK2','PYK3','PYK4','PYK5','PYK6','PYK7','PYK8','PYK9']
grouped_rxns[temp_label]['weights'] = [1,-1,-1,-1,-1,-1,-1,-1,-1,-1]

test_rxns = ['PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK']

for test_rxn in test_rxns:
    temp_label = test_rxn
    group_labels_rxn.append(temp_label)
    grouped_rxns[temp_label] = dict()
    grouped_rxns[temp_label]['rxns'] = [test_rxn]
    grouped_rxns[temp_label]['weights'] = [1]

print(grouped_rxns)

w.plot_time_correlations_vs_final(fig_dir,fig_label,
                                  '.png',
                                  'essentialfactors',
                                  conc_species,
                                  group_labels_count,
                                  grouped_species,
                                  group_labels_rxn,
                                  grouped_rxns)


