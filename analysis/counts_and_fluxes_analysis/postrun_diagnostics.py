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


fig_dir = './diagnostics_fig/'
fig_label = 'test'


############################
# prepare diagnostic tests #
############################

diagnostic_tests = dict()


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


#####################
# growth quantities #
#####################

w.plot_volume(fig_dir,
               fig_label,
               '.png',
              True)

w.plot_surfacearea(fig_dir,
                    fig_label,
                    '.png',
                   True)

w.plot_DNAcontent(fig_dir,
                  fig_label,
                  '.png',
                  True)

lt = w.get_tagged_loci('protein')
l = w.get_loci_lengths(lt)

w.plot_protein_metrics(fig_dir,
                       fig_label,
                       '.png',
                       lt,
                       l,
                       True)

w.plot_protein_transcript_metrics(fig_dir,
                                  fig_label,
                                  '.png',
                                  lt,
                                  l,
                                  True)

test_name = 'volume'
test_description = 'testing the volume'
test_result = True
diagnostic_tests[test_name] = dict()
diagnostic_tests[test_name]['description'] = test_description
diagnostic_tests[test_name]['result'] = test_result


#################################
# RNAP occupancy for rRNA genes #
#################################

# lt = w.get_tagged_loci('rRNA')

# max_RNAP = dict()

# max_RNAP['0067'] = 1
# max_RNAP['0068'] = 7
# max_RNAP['0069'] = 3
# max_RNAP['0532'] = 1
# max_RNAP['0533'] = 7
# max_RNAP['0534'] = 3

# group_labels = []
# grouped_species = dict()

# max_N_chromo = 2

# #for i_chromo in range(1,max_N_chromo+1):

# for i in range(len(lt)):

#     temp_label = 'RDME-gCME_boundRNAP_rRNA_' + lt[i]

#     if max_RNAP[lt[i]] > 1:

#         group_labels.append(temp_label)

#         grouped_species[temp_label] = dict()

#         grouped_species[temp_label]['ref_specie'] = None
#         grouped_species[temp_label]['species'] = []
#         grouped_species[temp_label]['weights'] = []
        
#         for j_RNAP in range(1,max_RNAP[lt[i]]+1):

#             # bound RNAP in RDME
#             temp_specie = 'RP_' + lt[i] + '_' + str(j_RNAP)
#             grouped_species[temp_label]['species'].append(temp_specie)
#             grouped_species[temp_label]['weights'].append(j_RNAP)

#             # bound RNAP in RDME
#             temp_specie = 'RP_' + lt[i] + '_t_' + str(j_RNAP)
#             grouped_species[temp_label]['species'].append(temp_specie)
#             grouped_species[temp_label]['weights'].append(j_RNAP)

#             # use below to test for mass conservation

#             # bound RNAP in gCME
#             temp_specie = 'RP_' + lt[i] + '_c1_' + str(j_RNAP)
#             grouped_species[temp_label]['species'].append(temp_specie)
#             grouped_species[temp_label]['weights'].append(-1)

#             # bound RNAP in gCME
#             temp_specie = 'RP_' + lt[i] + '_c2_' + str(j_RNAP)
#             grouped_species[temp_label]['species'].append(temp_specie)
#             grouped_species[temp_label]['weights'].append(-1)
            

#         grouped_species[temp_label]['weights'] = np.array(grouped_species[temp_label]['weights'],dtype=np.float32)


# x = w.get_group_traces(group_labels,
#                        grouped_species)

# # plot time-dependent RNAP occupancy (RDME) for rRNA genes
# w.plot_groups(fig_dir,fig_label,
#           '.png',
#           grouped_species,
#           group_labels)

# # test conservation between RDME and gCME
# diagnostic_metric = ''
# if (np.sum(x,axis=(0,1,2)) != 0):
#     test_result = False
# else:
#     test_result = True

# test_name = 'RDME_gCME_RNAP_rRNA_conservation'
# test_description = 'RNAP conservation of rRNA genes between RDME and gCME'
# diagnostic_tests[test_name] = dict()
# diagnostic_tests[test_name]['description'] = test_description
# diagnostic_tests[test_name]['result'] = test_result


# for i in range(len(lt)):

#     temp_label = 'RDME_boundRNAP_rRNA_' + lt[i]
    
#     group_labels.append(temp_label)

#     grouped_species[temp_label] = dict()

#     grouped_species[temp_label]['ref_specie'] = None
#     grouped_species[temp_label]['species'] = []
#     grouped_species[temp_label]['weights'] = []

#     # single RNAP occupancy
#     if max_RNAP[lt[i]] == 1:
        
#         grouped_species[temp_label]['species'].append('RP_'+lt[i])
#         grouped_species[temp_label]['weights'].append(1)

#     # multi RNAP occupancy
#     else:

#         for j_RNAP in range(1,max_RNAP[lt[i]]+1):

#             # bound RNAP in RDME
#             temp_specie = 'RP_' + lt[i] + '_' + str(j_RNAP)
#             grouped_species[temp_label]['species'].append(temp_specie)
#             grouped_species[temp_label]['weights'].append(j_RNAP)

#             # bound RNAP in RDME
#             temp_specie = 'RP_' + lt[i] + '_t_' + str(j_RNAP)
#             grouped_species[temp_label]['species'].append(temp_specie)
#             grouped_species[temp_label]['weights'].append(j_RNAP)            

#     grouped_species[temp_label]['weights'] = np.array(grouped_species[temp_label]['weights'],dtype=np.float32)

# # plot time-dependent RNAP occupancy (RDME) for rRNA genes
# w.plot_groups(fig_dir,fig_label,
#           '.png',
#           grouped_species,
#           group_labels)


###############################
# conservation of transcripts #
###############################

# transcript conservation for protein-coding genes
lt = w.get_tagged_loci('protein')
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

# test conservation of transcripts
diagnostic_metric = ''
if (np.amin(x) < 0):
    test_result = False
else:
    test_result = True

test_name = 'transcript_protein_conservation'
test_description = 'transcript conservation for protein-coding genes'
diagnostic_tests[test_name] = dict()
diagnostic_tests[test_name]['description'] = test_description
diagnostic_tests[test_name]['result'] = test_result

# w.plot_groups(fig_dir,fig_label,
#               '.png',
#               grouped_species,
#               group_labels)


# transcript conservation for tRNA-coding genes
lt = w.get_tagged_loci('tRNA')
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

# test conservation of transcripts
diagnostic_metric = ''
if (np.amin(x) < 0):
    test_result = False
else:
    test_result = True

test_name = 'transcript_tRNA_conservation'
test_description = 'transcript conservation for tRNA-coding genes'
diagnostic_tests[test_name] = dict()
diagnostic_tests[test_name]['description'] = test_description
diagnostic_tests[test_name]['result'] = test_result

# w.plot_groups(fig_dir,fig_label,
#               '.png',
#               grouped_species,
#               group_labels)


# transcript conservation for rRNA-coding genes
lt = w.get_tagged_loci('rRNA')
group_labels = []
grouped_species = dict()

for i in range(len(lt)):

    temp_label = 'transcripts_rRNA_' + lt[i]
    
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

# test conservation of transcripts
diagnostic_metric = ''
if (np.amin(x) < 0):
    test_result = False
else:
    test_result = True

test_name = 'transcript_rRNA_conservation'
test_description = 'transcript conservation for rRNA-coding genes'
diagnostic_tests[test_name] = dict()
diagnostic_tests[test_name]['description'] = test_description
diagnostic_tests[test_name]['result'] = test_result

# w.plot_groups(fig_dir,fig_label,
#               '.png',
#               grouped_species,
#               group_labels)


###############
# dNTP levels #
###############

test_species = ['M_datp_c','M_dgtp_c','M_dctp_c','M_dttp_c','M_dutp_c']


w.plot_multiple_conc_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_species,
                             'dntps')

w.plot_individual_concentrations(fig_dir,fig_label,
                                 '.png',
                                 test_species,
                                 True)

#########################
# key metabolite levels #
#########################

test_species = ['M_atp_c','M_pep_c','M_fdp_c','M_pi_c','M_ppi_c']


w.plot_multiple_conc_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_species,
                             'key_metabolites')

w.plot_individual_concentrations(fig_dir,fig_label,
                                 '.png',
                                 test_species,
                                 True)

#######################
# key reaction fluxes #
#######################

# test_rxns = w.get_rxns_list()

# w.plot_individual_rxns(fig_dir,fig_label,
#                        '.png',
#                        test_rxns)

test_rxns = ['PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK']

w.plot_multiple_rxns_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_rxns,
                             'central_metabolism')

w.plot_individual_rxns(fig_dir,
                       fig_label,
                       '.png',
                       test_rxns,
                       True)

test_rxns = ['PUNP1','PUNP2','PUNP3','PUNP4','PUNP5']

w.plot_multiple_rxns_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_rxns,
                             'PUNPs')

test_rxns = ['PYK','PYK2','PYK3','PYK4','PYK5','PYK6','PYK7','PYK8','PYK9']

w.plot_multiple_rxns_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_rxns,
                             'PYKs')

test_rxns = ['PGK','PGK2','PGK3','PGK4']

w.plot_multiple_rxns_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_rxns,
                             'PGKs')

test_rxns = ['ADNabc','DADNabc',\
             'GSNabc','DGSNabc',\
             'DCYTabc',\
             'URIabc',\
             'THMDabc']

w.plot_multiple_rxns_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_rxns,
                             'NucTransport')

test_rxns = ['PIabc','RIBFLVabc',\
             'P5Pabc','5FTHFabc',\
             'NACabc','COAabc',\
             'THMPPabc','SPRMabc']

w.plot_multiple_rxns_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_rxns,
                             'CofactorTransport')

w.plot_individual_rxns(fig_dir,
                       fig_label,
                       '.png',
                       test_rxns,
                       True)

test_rxns = ['ATPase','Kt6','MG2abc','CA2abc']

w.plot_multiple_rxns_per_rep(fig_dir,
                             fig_label,
                             '.png',
                             test_rxns,
                             'IonTransport')

w.plot_individual_rxns(fig_dir,
                       fig_label,
                       '.png',
                       test_rxns,
                       True)

group_labels = []
grouped_rxns = dict()

temp_label = 'CELL_DEATH'

group_labels.append(temp_label)
grouped_rxns[temp_label] = dict()
grouped_rxns[temp_label]['rxns'] = ['GLCpts0','PYK','PYK2','PYK3','PYK4','PYK5','PYK6','PYK7','PYK8','PYK9']
grouped_rxns[temp_label]['weights'] = [1,-1,-1,-1,-1,-1,-1,-1,-1,-1]

w.plot_rxn_groups(fig_dir,fig_label,
                  '.png',
                  grouped_rxns,
                  group_labels,
                  True)

################
# Key Proteins #
################

temp_proteins = [445,220,131,727,607,451,606,729,213,221]
temp_species = ['P_{:04d}'.format(x) for x in temp_proteins]

w.plot_individual_counts(fig_dir,fig_label,
                         '.png',
                         temp_species,
                         True)

w.plot_individual_concentrations(fig_dir,fig_label,
                                 '.png',
                                 temp_species,
                                 True)


###################################
# write out the diagnostic report #
###################################

report_file = fig_dir + fig_label + '_DiagnosticReport.txt'
with open(report_file, 'w') as f:

    test_list = list(diagnostic_tests.keys())
    N_tests = len(test_list)
    count_passed = 0

    for i_test in range(N_tests):
        temp_test = test_list[i_test]

        if diagnostic_tests[temp_test]['result'] == True:
            count_passed += 1

    report_header = '{:d}/{:d} -  diagnostic tests passed'.format(count_passed,
                                                                  N_tests)
    f.write('#'*(len(report_header)+4) + '\n')
    f.write('# ' + report_header + ' #\n')
    f.write('#'*(len(report_header)+4) + '\n\n')

    for i_test in range(N_tests):
        temp_test = test_list[i_test]

        if diagnostic_tests[temp_test]['result'] == True:
            outcome = 'PASS'
        else:
            outcome = 'FAIL'

        f.write(outcome + ' - ' + \
                diagnostic_tests[temp_test]['description'] + '\n')
