"""
This Script is written specifically to output the concentrations for Gibbs Free Energy calculation

Output the time-averaged concentrations of metabolites among different replicates
"""


from datetime import datetime
import numpy as np
import pandas as pd
import os
import importlib
import sys
from Bio import SeqIO
import argparse


import WCM_analysis_4DWCM
import WCM_gene as gene
import WCM_diagnosis as diagnosis
import WCM_mRNA as mRNA


ap = argparse.ArgumentParser()

ap.add_argument('-date','--simDate',required=True)
ap.add_argument('-start_time','--start_time',required=True)
ap.add_argument('-end_time','--end_time',required=True)


args = ap.parse_args()

# The date of the data 
simulation_date = args.simDate
start_time = 60*int(args.start_time)
end_time = 60*int(args.end_time)

pkl_dir = '/data/enguang/4DWCM/{0}/'.format(simulation_date)

pkl_label = '4DWCMensemble_{0}'.format(simulation_date)

# create a WCM_analysis class
w = WCM_analysis_4DWCM.WCM_ensemble()

# Read in the pkl file
w.read_merged_ensemble(pkl_dir,pkl_label)

# Diagnosis
healthy_indices = np.arange(0,w.N_reps)
healthy_list = [index +1 for index in healthy_indices]

# Log file
work_dir = '/data/enguang/4DWCM/{0}/output_concs/'.format(simulation_date)

if not os.path.exists(work_dir):
    os.makedirs(work_dir)

logfile_dir = work_dir + '/output_concs.log'

log = open(logfile_dir,'w')

orgin_stdout = sys.stdout

sys.stdout = log

# old_conc_file = '/data/enguang/CMEODE/Metabolism/parameter_analysis_Zane/average_concentrations.csv'

# old_conc_df = pd.read_csv(old_conc_file,delimiter=',')

# old_meta_list = list(old_conc_df['RXN'])

# old_meta_list = ['M_' + meta for meta in old_meta_list ]

meta_list = w.get_metabolite_species()

# print("In Old WCM but not current one: ", set(old_meta_list) - set(meta_list))

# print("In current WCM but not Old", set(meta_list) - set(old_meta_list))

# meta_concs = w.get_conc_traces(old_meta_list)

# print(np.shape(meta_concs))

print("in Simulation {0}, the number of healthy replicates is {1}".format(simulation_date, w.N_reps))

concs = w.get_conc_traces(meta_list)

mean_concs = np.mean(concs[:,start_time:end_time,:], axis=(1))

# concatenate the ensemble averaged count as the first column
mean_concs = np.concatenate((np.mean(mean_concs, axis=1, keepdims=True), mean_concs), axis=1)

print('The shape of mean_concs is {0}'.format(np.shape(mean_concs)))

meta_list = [meta[2:] for meta in meta_list]

print('The number of exported metabolites is {0}'.format(len(meta_list)))

# Export mean_concs into csv file
first_row = ['rep_' + str(i_rep+1) for i_rep in range(w.N_reps)]

first_row = ['rep_avg'] + first_row

first_row = np.array(first_row)

first_row = first_row.reshape((1, np.shape(first_row)[0]))

print(np.shape(first_row), "first_row")

# First concatenate the header
exported_values_row = np.concatenate((first_row, mean_concs), axis=0)

# exported_values_row.reshape(np.shape(exported_values_row)[0], 1)

print(np.shape(exported_values_row), "exported_values_row")

first_column = np.array(['metabolite'] + meta_list)

first_column = first_column.reshape((np.shape(first_column)[0], 1))

print(np.shape(first_column), "first_column")

exported_values_final = np.concatenate((first_column, exported_values_row), axis=1)

df = pd.DataFrame(exported_values_final)

df.to_csv(work_dir + 'avg_concs_{0}_{1}_{2}.csv'.format(simulation_date, start_time, end_time)
            , index=False, header=False)

sys.stdout = orgin_stdout