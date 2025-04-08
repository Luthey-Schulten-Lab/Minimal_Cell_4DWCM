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


fig_dir = './transcription_conservation_fig/'
fig_label = 'test'

protein_plotting = True
tRNA_plotting = False
rRNA_plotting = False


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


####################################################
# transcript conservation for protein-coding genes #
####################################################

lt = filtered_lt['protein']
# print(lt)

group_labels = []
grouped_species = dict()

for i in range(len(lt)):

    temp_label = 'transcripts_protein_' + lt[i]
    
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


x = w.get_group_traces(group_labels,
                       grouped_species)

print(np.amin(x))

if protein_plotting == True:
    w.plot_groups(fig_dir,fig_label,
              '.png',
              grouped_species,
              group_labels)


#################################################
# transcript conservation for tRNA-coding genes #
#################################################

lt = filtered_lt['tRNA']
# print(lt)

group_labels = []
grouped_species = dict()

for i in range(len(lt)):

    temp_label = 'transcripts_tRNA_' + lt[i]
    
    group_labels.append(temp_label)

    grouped_species[temp_label] = dict()

    grouped_species[temp_label]['ref_specie'] = 'R_' + lt[i]
    grouped_species[temp_label]['species'] = []
    grouped_species[temp_label]['weights'] = []
    
    grouped_species[temp_label]['species'].append('RPM_'+lt[i])
    grouped_species[temp_label]['weights'].append(1)

    grouped_species[temp_label]['weights'] = np.array(grouped_species[temp_label]['weights'],dtype=np.float32)


x = w.get_group_traces(group_labels,
                       grouped_species)

print(np.amin(x))

if tRNA_plotting == True:
    w.plot_groups(fig_dir,fig_label,
                  '.png',
                  grouped_species,
                  group_labels)


#################################################
# transcript conservation for rRNA-coding genes #
#################################################

# lt = filtered_lt['rRNA']
# # print(lt)

# group_labels = []
# grouped_species = dict()

# for i in range(len(lt)):

#     temp_label = 'transcripts_rRNA_' + lt[i]
    
#     group_labels.append(temp_label)

#     grouped_species[temp_label] = dict()

#     grouped_species[temp_label]['ref_specie'] = 'R_' + lt[i]
#     grouped_species[temp_label]['species'] = []
#     grouped_species[temp_label]['weights'] = []
    
#     grouped_species[temp_label]['species'].append('RPM_'+lt[i])
#     grouped_species[temp_label]['weights'].append(1)

#     grouped_species[temp_label]['weights'] = np.array(grouped_species[temp_label]['weights'],dtype=np.float32)

# x = w.get_group_traces(group_labels,
#                        grouped_species)

# print(np.amin(x))

# if rRNA_plotting == True:
#     w.plot_groups(fig_dir,fig_label,
#                   '.png',
#                   grouped_species,
#                   group_labels)
