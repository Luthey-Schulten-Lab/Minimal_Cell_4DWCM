from datetime import datetime
import numpy as np
import os
import importlib
import sys



import WCM_analysis_4DWCM

# Take in parameters
import argparse

ap = argparse.ArgumentParser()


ap.add_argument('-date','--simDate',required=True)
ap.add_argument('-n','--NumberReplicates',required=True)

args = ap.parse_args()

# The date of the data 
simulationdate = args.simDate

# working folder directory
work_dir = f'/data/enguang/4DWCM/{simulationdate}/'

if not os.path.exists(work_dir):
        print('WARNING: {0} does not exist'.format(work_dir))

subsystem = 'pkl'

logfile_dir = work_dir + f'{subsystem}_{simulationdate}.log'

log = open(logfile_dir,'w')
orgin_stdout = sys.stdout
sys.stdout = log

start_time = datetime.now()

print('The merging of CSV files starts at '+ str(start_time))

in_dir = f'/data/enguang/4DWCM/{simulationdate}/'

in_label = 'counts_and_fluxes'
min_rep = 1
max_rep = int(args.NumberReplicates)

reps_total = np.arange(min_rep,max_rep+1,dtype=np.int32)

print('{0} CSV files will be transformed to pkl file'.format(len(reps_total)))

pkl_dir = in_dir

merging_required = True

pkl_label = f'4DWCMensemble_{simulationdate}'


###################################
# pre-processing trajectory files #
###################################


if merging_required == True:

        w = WCM_analysis_4DWCM.WCM_ensemble()

        # this assumes trajectories are in the below format
        # (in_dir)/(in_label).(replicate).csv
        w.set_traj_files(in_dir,in_label,reps_total)

        w.load_trajs()

        w.merge_trajs()

        w.write_merged_ensemble(pkl_dir,pkl_label)

        print('New ensemble of trajectories Created in folder {0}'.format(pkl_dir))


end_time = datetime.now()

print('The merging of CSV files finishes at '+ str(end_time))

print('Time (hour:minutes:seconds) taken to finish the serializing: {0}'.format(end_time - start_time))

# create a WCM_analysis class
w = WCM_analysis_4DWCM.WCM_ensemble()

# Read in the pkl file
w.read_merged_ensemble(pkl_dir,pkl_label)

print('The length of {0} trajs is {1}'.format(w.N_reps, w.t))


sys.stdout = orgin_stdout