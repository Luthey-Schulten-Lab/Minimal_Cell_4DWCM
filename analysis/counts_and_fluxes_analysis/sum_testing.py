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
in_label = 'particle_counts'
min_rep = 1
max_rep = 4

reps = np.arange(min_rep,max_rep+1,dtype=np.int32)

merging_required = False
out_dir = './test_data/'
out_label = 'ensemble_test'


fig_dir = './test_fig/'
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

# test_species = ['M_1ag3p_c','M_12dgr_c','M_10fthfglu3_c']

# regular expression to get all metabolite species
# species_regex = r"(\bM_\w+)"
# test_species = w.get_regex_species_list(species_regex)

# dedicated functions
# test_species = w.get_metabolite_species()
# test_species = w.get_tRNA_species()
# test_species = w.get_single_gene_species('0001')


lt = w.get_all_locus_tags()

group_labels = []
grouped_species = dict()

for i in range(4):

    temp_label = 'group_{:d}'.format(i)
    
    group_labels.append(temp_label)

    grouped_species[temp_label] = dict()

    grouped_species[temp_label]['ref_specie'] = 'R_' + lt[i]
    grouped_species[temp_label]['species'] = []
    grouped_species[temp_label]['weights'] = []
    
    grouped_species[temp_label]['species'].append('RPM_'+lt[i])
    grouped_species[temp_label]['weights'].append(1)
    grouped_species[temp_label]['species'].append('DM_'+lt[i])
    grouped_species[temp_label]['weights'].append(-1)

    grouped_species[temp_label]['weights'] = np.array(grouped_species[temp_label]['weights'],dtype=np.float32)

print(grouped_species)

x = w.get_group_traces(group_labels,
                       grouped_species)

w.plot_groups(fig_dir,fig_label,
              '.png',
              grouped_species,
              group_labels)
