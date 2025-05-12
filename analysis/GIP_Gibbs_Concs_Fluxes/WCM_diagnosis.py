"""

Description: 
"""

import numpy as np
from scipy.ndimage import uniform_filter1d 

import WCM_math as math
import WCM_gene as gene

# def atp_shortage(w):
#     """

#     Description: Distinguish healthy/unhealthy based on whether ATP shortage happened or not 
#     """

#     unhealthy_list = []
    
#     Nreps = w.N_reps

#     accumulative_atp = 'M_atp_c_shortage'

#     accumulative_atp_trace = w.get_specie_trace(accumulative_atp)

#     minimum_shortage = np.min(accumulative_atp_trace, axis = 0)

#     unhealthy_boolen_list = minimum_shortage < 0

#     healthy_list = [i_rep+1 for i_rep in range(Nreps) if unhealthy_boolen_list[i_rep] == False]

#     unhealthy_list = [i_rep+1 for i_rep in range(Nreps) if unhealthy_boolen_list[i_rep] == True]

#     print('The number of unhealthy replicates is {0} out of {1} with ratio {2:.2f}'.format(len(unhealthy_list), Nreps, len(unhealthy_list)/Nreps ))

#     return healthy_list, unhealthy_list

def pep_zero(w):
    """
    
    Return the mean concentration of pep between certain time period across different cells
    """
    

    return 


def get_volume_doubling_times(w, reps):
    """
    Get the cell volume doubling times of given rep
    Return a list
    """

    reps = [rep -1 for rep in reps]

    volume = w.volumes[:,reps]

    scaled_volume = volume*np.reciprocal(volume[0,:])

    volume_double_times = []; not_doubled_reps = []

    fold = 2

    for rep in reps:
        y_rep = scaled_volume[:,rep]
        indices = np.where(y_rep > fold)
        first_index = indices[0][0] if indices[0].size > 0 else -1
        if first_index == -1:
            doubling_time = -1
            not_doubled_reps.append(rep+1)
        else: 
            doubling_time = w.t[first_index]
        volume_double_times.append(doubling_time)

    print("{0} replicates {1} in input {2} reps don't double the Volume".format(len(not_doubled_reps), 
                                                                                ','.join(map(str, not_doubled_reps)),
                                                                                len(reps)))

    MMMM = math.get_min_median_mean_max(volume_double_times)/60

    print("Volume Doubling Time: Minimum, Median, Mean, and Maximum is {0} minutes".format(', '.join(f'{t:.2f}' for t in MMMM)))

    return volume_double_times, not_doubled_reps


def get_surface_doubling_times(w, reps, print_flag=True):
    """
    Get the cell surface area doubling time of given reps

    Return a list
    """
    reps = [rep -1 for rep in reps]
    try:
        surface = w.get_specie_trace('SA_nm2')
    except:
        surface = w.get_specie_trace('SA_total')

    scaled_surface = surface*np.reciprocal(surface[0,:])

    surface_doubling_times = []; not_doubled_reps = []

    fold = 2

    for rep in reps:
        y_rep = scaled_surface[:,rep]
        indices = np.where(y_rep > fold)
        first_index = indices[0][0] if indices[0].size > 0 else -1
        if first_index == -1:
            doubling_time = -1
            not_doubled_reps.append(rep+1)
        else:
            doubling_time = w.t[first_index]
        surface_doubling_times.append(doubling_time)

    MMMM = math.get_min_median_mean_max(surface_doubling_times)/60

    if print_flag:
        print("{0} replicates {1} in input {2} reps don't double the Surface Area".format(len(not_doubled_reps), 
                                                                                    ','.join(map(str, not_doubled_reps)),
                                                                                    len(reps)))

        print("Surface Area Doubling Time: Minimum, Median, Mean, and Maximum is {0} minutes".format(', '.join(f'{t:.2f}' for t in MMMM)))

    return surface_doubling_times, not_doubled_reps

def plot_ribosome(w, reps, fig_dir, fig_label):
    """
    plot traces of produced, active, and total ribosome
    plot trace of ratio of active ribosome
    """



    return None

def plot_degradosome(w, reps, fig_dir, fig_label):
    
    TypetoLocusNums = w.TypetoLocusNums

    # free degradosome
    free_degra_count = w.get_specie_trace('Degradosome')

    # Calculate active degradosome count
    ptn_locusNums = TypetoLocusNums['protein']

    degra_bound_ids = ['Degradosome_mRNA_' + locusNum for locusNum in ptn_locusNums]

    degra_bound_counts = w.get_species_traces(degra_bound_ids)

    active_degra_count = np.sum(degra_bound_counts, axis=0)

    total_degra_count = free_degra_count + active_degra_count

    ratio_degra = math.divide_replace(active_degra_count, total_degra_count, 0)

    degra_prefix = ['Free', 'Active', 'Total', 'Active Ratio']
    
    counts = [free_degra_count, active_degra_count, total_degra_count, ratio_degra]

    for i_prefix, prefix in enumerate(degra_prefix):
        ylabel = 'Count [\#]' if not prefix == 'Active Ratio' else 'Ratio [\%]'
        
        title = '{0} degradosome'.format(prefix)
        
        w.plot_in_replicates_single(fig_dir, fig_label, 
                                    '.png', counts[i_prefix], reps,
                                    ylabel, title, True, True)
    
    surface_doubling_times, not_doubled_surface_reps = get_surface_doubling_times(w, reps, False)

    total_degra_count_SA = math.get_doubling_moments_value(surface_doubling_times, total_degra_count)

    xlabel = 'Count [\#]'
    ylabel = 'Frequency'
    title = 'Total Degradosome at SA Doubling time'

    w.plot_hist(fig_dir, fig_label, '.png', 
                total_degra_count_SA, 
                xlabel, ylabel, title, bins=50)
    
    return None

def plot_RNaseY_Dimer(w, reps, fig_dir, fig_label):
    
    TypetoLocusNums = w.TypetoLocusNums

    # free degradosome
    free_count = w.get_specie_trace('RNaseY_Dimer')

    # Calculate active degradosome count
    ptn_locusNums = TypetoLocusNums['protein']

    bound_ids = ['RNaseY_Dimer_mRNA_' + locusNum for locusNum in ptn_locusNums]

    bound_counts = w.get_species_traces(bound_ids)

    active_count = np.sum(bound_counts, axis=0)

    total_count = free_count + active_count

    ratio_degra = active_count/total_count

    degra_prefix = ['Free', 'Active', 'Total', 'Active Ratio']
    
    counts = [free_count, active_count, total_count, ratio_degra]

    for i_prefix, prefix in enumerate(degra_prefix):
        ylabel = 'Count [\#]' if not prefix == 'Active Ratio' else 'Ratio [\%]'
        
        title = '{0} RNaseY_Dimer'.format(prefix)
        
        w.plot_in_replicates_single(fig_dir, fig_label, 
                                    '.png', counts[i_prefix], reps,
                                    ylabel, title, True, True)
    
    return None


def plot_RNAP(w, reps, fig_dir, fig_label):

    TypetoLocusNums = w.TypetoLocusNums

    genes = []

    for geneType in TypetoLocusNums.keys():
        if geneType != 'gene':
            genes.extend(TypetoLocusNums[geneType])

    free_RNAP_count = w.get_specie_trace('RNAP')

    # Active RNAP

    RNAP_bound_ids = ['RNAP_G_' + locusNum for locusNum in genes]

    RNAP_bound_counts = w.get_species_traces(RNAP_bound_ids)

    active_RNAP_count = np.sum(RNAP_bound_counts, axis=0)

    total_RNAP_count = free_RNAP_count + active_RNAP_count

    ratio_RNAP = math.divide_replace(active_RNAP_count, total_RNAP_count, 0)

    RNAP_prefix = ['Free', 'Active', 'Total', 'Active Ratio']

    counts = [free_RNAP_count, active_RNAP_count, total_RNAP_count, ratio_RNAP]

    for i_prefix, prefix in enumerate(RNAP_prefix):
        ylabel = 'Count [\#]' if not prefix == 'Active Ratio' else 'Ratio [\%]'
        
        title = '{0} RNAP'.format(prefix)
        
        w.plot_in_replicates_single(fig_dir, fig_label, 
                                    '.png', counts[i_prefix], reps,
                                    ylabel, title, True, True)
    
    surface_doubling_times, not_doubled_surface_reps = get_surface_doubling_times(w, reps, False)

    total_RNAP_count_SA = math.get_doubling_moments_value(surface_doubling_times, total_RNAP_count)
    
    xlabel = 'Count [\#]'
    ylabel = 'Frequency'
    title = 'Total RNAP at SA Doubling time'

    w.plot_hist(fig_dir, fig_label, '.png', 
                total_RNAP_count_SA, 
                xlabel, ylabel, title, bins=50)
    return None


# def plot_RNAP(w, reps, fig_dir, fig_label):

#     TypetoLocusNums = w.TypetoLocusNums

#     genes = []

#     for geneType in TypetoLocusNums.keys():
#         if geneType != 'gene':
#             genes.extend(TypetoLocusNums[geneType])

#     RNAP_core_count = w.get_specie_trace('RNAP_core')
#     free_RNAP_count = w.get_specie_trace('RNAP') + RNAP_core_count

#     # Active RNAP

#     RNAP_bound_ids = ['RNAP_G_' + locusNum for locusNum in genes]

#     RNAP_bound_counts = w.get_species_traces(RNAP_bound_ids)

#     active_RNAP_count = np.sum(RNAP_bound_counts, axis=0)

#     total_RNAP_count = free_RNAP_count + active_RNAP_count

#     ratio_RNAP = active_RNAP_count/total_RNAP_count

#     RNAP_prefix = ['RNAP_Core','Free', 'Active', 'Total', 'Active Ratio']

#     counts = [RNAP_core_count, free_RNAP_count, active_RNAP_count, total_RNAP_count, ratio_RNAP]

#     for i_prefix, prefix in enumerate(RNAP_prefix):
#         ylabel = 'Count [\#]' if not prefix == 'Active Ratio' else 'Ratio [\%]'
        
#         title = '{0} RNAP'.format(prefix)
        
#         w.plot_in_replicates_single(fig_dir, fig_label, 
#                                     '.png', counts[i_prefix], reps,
#                                     ylabel, title, True, True)
        
#     return None


def plot_cost_smoothed_nucleotide(w, fig_dir, fig_label, cost_names, catagory, nucleotide, windowsize):
        
        ATP_costs = w.get_species_traces(cost_names) # ATP_costs is a 3D array with axies species, time, replicate
        title = '{0} {1} Cost per second'.format(nucleotide, catagory)
        w.plot_ensemble_averaged_multiples(fig_dir, fig_label, '.png',cost_names, 'count', title, True, True )

        ATP_costs_cumsum = np.cumsum(ATP_costs, axis = 1)
        title = '{0} {1} Cost Accumulative'.format(nucleotide, catagory)
        w.plot_ensemble_averaged_multiples_ym(fig_dir, fig_label, '.png', np.mean(ATP_costs_cumsum, axis = 2), cost_names, 'count', title, True, True, linewidth = 1)

        title = '{0} {1} Cost per second Smoothed'.format(nucleotide, catagory)
        averaged_ATP_costs = uniform_filter1d(ATP_costs, size=windowsize, axis=1, mode='reflect')
        w.plot_ensemble_averaged_multiples_ym(fig_dir, fig_label, '.png', np.mean(averaged_ATP_costs, axis = 2), cost_names, 'count', title, True, True )

        return None