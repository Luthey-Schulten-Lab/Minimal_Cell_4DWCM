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


fig_dir = './RNAP_rRNA_conservation_fig/'
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

filtered_lt = w.get_filtered_locus_tags('./gene_types.txt')

# print(filtered_lt)


####################################################
# transcript conservation for protein-coding genes #
####################################################

lt = filtered_lt['rRNA']
print(lt)

for l in lt:

    print(l)
    print(w.get_single_gene_species(l))

max_RNAP = dict()

max_RNAP['0067'] = 1
max_RNAP['0068'] = 7
max_RNAP['0069'] = 3
max_RNAP['0532'] = 1
max_RNAP['0533'] = 7
max_RNAP['0534'] = 3

group_labels = []
grouped_species = dict()

for i in range(len(lt)):

    temp_label = 'boundRNAP_rRNA_' + lt[i]
    
    group_labels.append(temp_label)

    grouped_species[temp_label] = dict()

    grouped_species[temp_label]['ref_specie'] = None
    grouped_species[temp_label]['species'] = []
    grouped_species[temp_label]['weights'] = []

    if max_RNAP[lt[i]] == 1:
        
        grouped_species[temp_label]['species'].append('RP_'+lt[i])
        grouped_species[temp_label]['weights'].append(1)
        
    else:

        for j_RNAP in range(1,max_RNAP[lt[i]]+1):

            temp_specie = 'RP_' + lt[i] + '_' + str(j_RNAP)
            grouped_species[temp_label]['species'].append(temp_specie)
            grouped_species[temp_label]['weights'].append(j_RNAP)

            temp_specie = 'RP_' + lt[i] + '_t_' + str(j_RNAP)
            grouped_species[temp_label]['species'].append(temp_specie)
            grouped_species[temp_label]['weights'].append(j_RNAP)

            # use below to test for mass conservation
            
            # temp_specie = 'RP_' + lt[i] + '_c1_' + str(j_RNAP)
            # grouped_species[temp_label]['species'].append(temp_specie)
            # grouped_species[temp_label]['weights'].append(-1)

            # temp_specie = 'RP_' + lt[i] + '_c2_' + str(j_RNAP)
            # grouped_species[temp_label]['species'].append(temp_specie)
            # grouped_species[temp_label]['weights'].append(-1)
            

    grouped_species[temp_label]['weights'] = np.array(grouped_species[temp_label]['weights'],dtype=np.float32)


x = w.get_group_traces(group_labels,
                       grouped_species)

print(np.amin(x))

w.plot_groups(fig_dir,fig_label,
          '.png',
          grouped_species,
          group_labels)
