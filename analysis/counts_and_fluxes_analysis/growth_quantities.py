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


fig_dir = './growth_fig/'
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


###############
# plot volume #
###############

w.plot_volume(fig_dir,
               fig_label,
               '.png')

w.plot_surfacearea(fig_dir,
                    fig_label,
                    '.png')

w.plot_DNAcontent(fig_dir,
                  fig_label,
                  '.png')
