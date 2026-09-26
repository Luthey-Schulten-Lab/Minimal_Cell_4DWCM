"""
Authors
-------
Ron Acda — tRNA-charging table parsed once per process
    (using an iterative LLM-guided workflow: https://github.com/quarkron/iterative-hillclimber/tree/main)
Zane Thornburg — original code

General functions to create reactions in a global CME simulation
"""

import numpy as np
import pandas as pd
from processes.Rxns_ODE import _read_excel_cached
from collections import defaultdict, OrderedDict

import utility.GIP_rates as GIP


#########################################################################################
def transcription(csim, sim_properties, locusTag, rnasequence):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    locusNum = locusTag.split('_')[1]
    
    RNAP_gene = 'RP_' + locusNum + '_C1'
    
    new_RNA_ID = 'RP_' + locusNum + '_f' + '_C1'
        
    csim.addReaction(RNAP_gene, new_RNA_ID, GIP.TranscriptionRate(sim_properties, locusTag, rnasequence))
    
    
    RNAP_gene = 'RP_' + locusNum + '_C2'
    
    new_RNA_ID = 'RP_' + locusNum + '_f' + '_C2'
        
    csim.addReaction(RNAP_gene, new_RNA_ID, GIP.TranscriptionRate(sim_properties, locusTag, rnasequence))
    
    return None
#########################################################################################


#########################################################################################
def transcriptionLong(csim, sim_properties, locusTag, rnasequence):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    locusNum = locusTag.split('_')[1]
    
    for chromo in range(1,3):
    
        for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):
            
            off_i = 'RP_' + locusNum + '_c' + str(chromo) + '_open_' + str(i)

            on_i = 'RP_' + locusNum + '_c' + str(chromo)+ '_' + str(i)
            
            if i == sim_properties['long_rna_trsc'][locusTag]['max_rnap']:
            
                off_ip1 = 'RP_' + locusNum  + '_c' + str(chromo) + '_endTrsc'
                
                on_ip1 = 'RP_' + locusNum  + '_c' + str(chromo) + '_done'
                
                trscsequence  = rnasequence[int((i-1)*sim_properties['rnap_spacing']):]
                
            else:
                
                off_ip1 = 'RP_' + locusNum + '_c' + str(chromo) + '_open_' + str(int(i+1))

                on_ip1 = 'RP_' + locusNum + '_c' + str(chromo) + '_' + str(int(i+1))
                
                trscsequence = rnasequence[int((i-1)*sim_properties['rnap_spacing']):int(i*sim_properties['rnap_spacing'])]
                
            substrates = [off_ip1, on_i]
            
            products = [on_ip1, off_i]

            csim.addReaction(tuple(substrates), tuple(products), GIP.TranscriptionRate(sim_properties, locusTag, trscsequence))
    
    return None
#########################################################################################


#########################################################################################
_TRNA_PARAMS = {}


def _trna_params_cached(RXNS_params, key):
    """{reaction name: {parameter type: value}} of the tRNA-charging sheet, as the per-call filters returned them."""
    if key not in _TRNA_PARAMS:
        out = {}
        for rxnID in RXNS_params["Reaction Name"].unique():
            rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]
            d = {}
            for ptype in ('amino acid', 'synthetase', 'k_atp', 'k_aa', 'k_tRNA', 'k_cat'):
                sel = rxn_params.loc[ rxn_params["Parameter Type"] == ptype ]["Value"].values
                if len(sel): d[ptype] = sel[0]
            out[rxnID] = d
        _TRNA_PARAMS[key] = out
    return _TRNA_PARAMS[key]


def tRNAcharging(csim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    # the static tRNA-charging table is parsed from Excel once per process (Rxns_ODE's cache, copy per call)
    RXNS_params = _read_excel_cached(sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx', sheet_name='tRNA Charging')
    
    # the parameter rows of each charging reaction are looked up once per process (the table is static); the values are
    # the same objects the per-call .loc filters returned
    _p = _trna_params_cached(RXNS_params, sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx')

    for tRNA_aa, rnaIDlist in sim_properties['trna_map'].items():
        
#         if tRNA_aa == 'GLN':
            
#             print('Sad')

#         else:
            
        rxnID = tRNA_aa + 'TRS'

        rp = _p[rxnID]

        aaID = rp['amino acid']

        synthetaseID = rp['synthetase']

        synthetaseAtpID = synthetaseID + '_atp'

        csim.addReaction(tuple([synthetaseID, 'M_atp_c']), synthetaseAtpID, rp['k_atp'])

        synthetaseAaID = synthetaseAtpID + '_aa'

        csim.addReaction(tuple([synthetaseAtpID, aaID]), synthetaseAaID, rp['k_aa'])

        for rnaID in rnaIDlist:

            synthetaseTrnaID = synthetaseAaID + '_' + rnaID

            csim.addReaction(tuple([synthetaseAaID, rnaID]), synthetaseTrnaID, rp['k_tRNA'])

            chargedTrnaID = rnaID + '_ch'

            csim.addReaction(synthetaseTrnaID, tuple(['M_amp_c', 'M_ppi_c', synthetaseID, chargedTrnaID]), rp['k_cat'])

            costID = tRNA_aa + '_cost'
            costPaidID = costID + '_paid'

            csim.addReaction(tuple([costID, chargedTrnaID]), tuple([rnaID, costPaidID]), 1e5)
    
    print('tRNA')
    
    return None
#########################################################################################



    
    
    
