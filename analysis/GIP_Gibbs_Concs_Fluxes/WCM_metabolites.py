"""

Description: Functions to analyze the features of metabolites' time traces
"""

import numpy as np

def cal_ribosome(w):

    total_ribosome_count = np.zeros_like(w.get_specie_trace('M_k_c'))

    for locusTag, subdict in w.genomeDict.items():
        locusNum = locusTag.split('_')[1]

        if locusNum in w.TypetoLocusNums['protein']: #mRNAs
            ribo_mRNA = 'RB_' + locusNum
            count = w.get_specie_trace(ribo_mRNA)
            total_ribosome_count += count

    total_ribosome_count += w.get_specie_trace('ribosomeP')


    return total_ribosome_count


def cal_bound_Mg(w, ribo_init):

    Mg_perNTP = 1
    Mg_perribo = 1346

    total_ribo = cal_ribosome(w)

    produced_ribosome_count = total_ribo - ribo_init

    # Extract the bound Mg2+ from total Mg2+

    NTP2chg = {'ATP':4, 'CTP': 4, 'GTP': 4, 'UTP': 4, 
                    'dATP':4, 'dCTP':4, 'dGTP':4, 'dTTP':4}

    bound_Mg_count = np.zeros_like(w.get_specie_trace('M_k_c'))

    bound_Mg_count += produced_ribosome_count*Mg_perribo

    for NTP, chg in NTP2chg.items():
        specie = 'M_' + NTP.lower() + '_c'
        bound_Mg_count += (w.get_specie_trace(specie) - w.get_specie_trace(specie)[0,:])*Mg_perNTP # minus the initial NTPs

    total_Mg_count = w.get_specie_trace('M_mg2_c')

    extracted_Mg_count = total_Mg_count - bound_Mg_count

    free_Mg_count = np.maximum(extracted_Mg_count, 0)

    return bound_Mg_count, extracted_Mg_count, free_Mg_count


# def get_conc_increase(w):
#     """
    
#     """

#     species_map = w.species_map

#     for specie in species_map.keys():

#         trace = w.get_specie_trace(specie)


#     return None

def get_metabolite_rxns(SBML_model, metabolite_id):
    
    """
    Get the list of metabolite involved reactions 
    Return: rxns_sub: dictionary where key is rxn where metabolite as substrate and value is the stoichiometric coefficient
    """

    rxns_sub, rxns_prod = {}, {}

    metabolite = SBML_model.metabolites.get_by_id(metabolite_id)

    for rxn in metabolite.reactions: # metabolite.reactions is a set of reactions
        rxn_id = rxn.id
        stoi = rxn.metabolites[metabolite]
        if stoi < 0:
            rxns_sub[rxn_id] = stoi
        else:
            rxns_prod[rxn_id] = stoi

    return rxns_sub, rxns_prod 

def get_metabolites_rxns(SBML_model, metabolites_id):
    """
    Get the reactions that all reactions contains metabolites in metabolites_id
    """
    rxns = []

    for metabolite in metabolites_id:
        
        rxns_sub, rxns_prod = get_metabolite_rxns(SBML_model, metabolite)

        rxns.extend(list(rxns_sub.keys())); rxns.extend(list(rxns_prod.keys()))

    rxns = list(set(rxns))
    
    return rxns


def get_fluxes_threshold(w, rxns, threshold):
    """
    Get the time spans where the fluxes of rxns exceeds certain threshold
    Return: rxns_time, keys are the rxns where the fluxes exceed threshold

    """

    rxns_times = {}
    
    fluxes = w.get_rxn_traces(rxns)

    for i_rxn, rxn in enumerate(rxns):
        
        flux = fluxes[i_rxn,:,:]

        if np.max(flux) > threshold:
    
            indices = np.asarray(flux>threshold).nonzero()[0] # [0] means at time axis not replicate
    
            rxns_times[rxn] = indices

            print(f"Reaction {rxn} exceeds {threshold} in time {indices}")

    print(f"{', '.join(list(rxns_times.keys()))} exceeds threshold {threshold} out of {len(rxns)} reactions {', '.join(rxns)} ")

    return rxns_times

def get_rxn_metabolites(SBML_model, rxn):
    """
    Return the met id of one particular reaction with prefix 'M_'
    """
    SBML_rxn = SBML_model.reactions.get_by_id(rxn)

    SBML_metabolites = SBML_rxn.metabolites.keys()

    ids = ['M_' + meta.id for meta in SBML_metabolites]

    return ids

def get_rxns_metabolites(SBML_model, rxns):
    
    ids = []
    for rxn in rxns:
        ids.extend(get_rxn_metabolites(SBML_model, rxn))
    
    ids = list(set(ids))

    return ids

def plot_metabolite_fluxes(w, SBML_model, metabolite_id, healthy_list):
    """
    metabolite_id: adp_c, ...
    plot the individual fluxes, summed in and out fluxes, and net flux of one single metabolite
    """

    rxns_sub, rxns_prod = get_metabolite_rxns(SBML_model, metabolite_id)
    fig_dir = w.fig_dir; fig_label = w.fig_label
    all_rxns = list(rxns_sub.keys()) + list(rxns_prod.keys())

    # plot individual fluxes

    # for rxn in all_rxns:
    #   flux = w.get_rxn_trace(rxn)
    #   w.plot_in_replicates_single(fig_dir,fig_label,
    #                            '.png',
    #                           flux, healthy_list, 'Flux', rxn,
    #                            True, True)
    
    # plot multiple rxns
    fluxes = w.get_rxn_traces(all_rxns)

    w.plot_ensemble_averaged_multiples_ym(fig_dir,fig_label,
                               '.png', np.mean(fluxes, axis=2), all_rxns, 'flux', 
                               'rxns' + metabolite_id, True, True, linewidth = 0.5)
    
    # w.plot_multiples_per_rep(fig_dir,fig_label,
    #                                '.png',
    #                                fluxes, all_rxns, healthy_list,
    #                                f'{metabolite_id}_all', 'flux')
    
    # plot in and out
    if metabolite_id != 'pep_c' and metabolite_id != '13dpg_c':
        plot_inout_fluxes(w, rxns_sub, rxns_prod, metabolite_id, healthy_list)

    elif metabolite_id == 'pep_c':
        plot_inout_fluxes_pep(w, healthy_list)

    elif metabolite_id == '13dpg_c':
        plot_inout_fluxes_13dpg(w, healthy_list)
    
    # plot net fluxes
    plot_net_flux(w, rxns_sub, rxns_prod, metabolite_id, healthy_list)

    return None


def plot_inout_fluxes(w, rxns_sub, rxns_prod, metabolite_id, healthy_list):
    """
    Stoichiometric coefficient considered
    """
    
    fig_dir = w.fig_dir; fig_label = w.fig_label

    fluxes_inout = np.zeros((2,w.Nt,len(healthy_list)))

    for rxn, stoi in rxns_sub.items():
        fluxes_inout[0,:,:] += np.abs(stoi) * w.get_rxn_trace(rxn)
    
    for rxn, stoi in rxns_prod.items():
        fluxes_inout[1,:,:] += np.abs(stoi) * w.get_rxn_trace(rxn)

    qauntities = [f'out_{metabolite_id}', f"in_{metabolite_id}"]

    # w.plot_multiples_per_rep(fig_dir,fig_label,
    #                                 '.png',
    #                                 fluxes_inout, qauntities, healthy_list,
    #                                 f'{metabolite_id}_in_out', 'flux')
    
    w.plot_ensemble_averaged_multiples_ym(fig_dir,fig_label,
                               '.png', np.mean(fluxes_inout, axis=2), qauntities, 'flux', 
                               f'{metabolite_id}_in_out', True, True, linewidth = 0.5)
    return None


def plot_inout_fluxes_pep(w, healthy_list):
    """
    in_out fluxes of pep
    """
    fig_dir = w.fig_dir; fig_label = w.fig_label
    
    rxns_pep = ['PYK' + str(i) if i != 1 else 'PYK' for i in range(1, 10)]
    rxns_pep.extend(['GLCpts0', 'ENO'])


    flux_pep = np.zeros((3,w.Nt,len(healthy_list)))

    for rxn in rxns_pep:
        if rxn == 'GLCpts0':
            flux_pep[1,:,:] = w.get_rxn_trace(rxn)
        elif rxn == 'ENO':
            flux_pep[2,:,:] = w.get_rxn_trace(rxn)
        else:
            flux_pep[0,:,:] += w.get_rxn_trace(rxn)

    qauntities_pep = ['Summed_PYKs', 'GLCpts0', 'ENO']

    w.plot_multiples_per_rep(fig_dir,fig_label,
                                    '.png',
                                    flux_pep, qauntities_pep, healthy_list,
                                    'pep_in_out', 'flux')
    return None


def plot_inout_fluxes_13dpg(w, healthy_list):

    fig_dir = w.fig_dir; fig_label = w.fig_label

    rxns_dpg = ['GAPD', 'PGK', 'PGK2', 'PGK3', 'PGK4']

    flux_dpg = np.zeros((2,w.Nt,len(healthy_list)))

    for rxn in rxns_dpg:
        if rxn.startswith('PGK'):
            flux_dpg[0,:,:] += w.get_rxn_trace(rxn)
        else:
            flux_dpg[1,:,:] = w.get_rxn_trace(rxn)

    quantities_dpg = ['Summed_PGK', 'GAPD']

    w.plot_multiples_per_rep(fig_dir,fig_label,
                                    '.png',
                                    flux_dpg, quantities_dpg, healthy_list,
                                    '1,3dpg_in_out', 'flux')
    
    return None



def plot_net_flux(w, rxns_sub, rxns_prod, metabolite_id, healthy_list):
    """
    Stoichiometric coefficient considered
    """
    
    fig_dir = w.fig_dir; fig_label = w.fig_label

    flux_net = np.zeros((1, w.Nt,len(healthy_list)))

    for rxn, stoi in rxns_sub.items():
        flux_net[0,:,:] += stoi * w.get_rxn_trace(rxn)
    
    for rxn, stoi in rxns_prod.items():
        flux_net[0,:,:] += stoi * w.get_rxn_trace(rxn)

    quantity = [f'net_flux_{metabolite_id}']

    # w.plot_multiples_per_rep(fig_dir,fig_label,
    #                                 '.png',
    #                                 flux_net, quantity, healthy_list,
    #                                 f'net_flux_{metabolite_id}', 'flux')
    
    w.plot_ensemble_averaged_multiples_ym(fig_dir,fig_label,
                               '.png', np.mean(flux_net, axis=2), quantity, 'flux', 
                               f'{metabolite_id}_net', True, True, linewidth = 0.5)
    
    return None