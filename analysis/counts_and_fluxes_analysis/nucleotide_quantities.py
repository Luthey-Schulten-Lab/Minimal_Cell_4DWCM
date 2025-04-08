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


fig_dir = './nucleotide_fig/'
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

test_species = w.get_nucleobase_species()

print(test_species)


# test_species = w.get_metabolite_species()


##################
# getting traces #
##################

# time
t = w.get_t()
# counts
# x = w.get_species_traces(test_species)
# concentractions
c = w.get_conc_traces(test_species)
# counts
# x = w.get_avg_species_traces(test_species)
# concentractions
# c = w.get_avg_conc_traces(test_species)


####################
# plotting species #
####################

# plot concentrations in individual plots for each specie
# w.plot_individual_concentrations(fig_dir,
#                                  fig_label,
#                                  '.png',
#                                  test_species)

# plot sample covariances/correlations that calculated using ensemble
# averages at each timepoint before being time-averaged within the
# provided time window
# time_window = np.array([2600,3200],dtype=np.float32)

# w.plot_correlations(fig_dir,
#                     fig_label,
#                     '.png',
#                     time_window,
#                     test_species,
#                     'conc',
#                     'corrX')

test_species = ['M_datp_c','M_dgtp_c','M_dctp_c','M_dttp_c','M_dutp_c']


w.plot_multiple_conc_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_species,
                             'dntp_test')
