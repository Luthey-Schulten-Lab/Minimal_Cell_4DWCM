"""
Authors: Zane Thornburg

General functions to create reactions in a global CME simulation
"""

import numpy as np
import pandas as pd
from collections import defaultdict, OrderedDict

import GIP_rates as GIP


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
def tRNAcharging(csim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    RXNS_params = pd.read_excel(sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx', sheet_name='tRNA Charging')
    
    for tRNA_aa, rnaIDlist in sim_properties['trna_map'].items():
        
#         if tRNA_aa == 'GLN':
            
#             print('Sad')

#         else:
            
        rxnID = tRNA_aa + 'TRS'

        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]

        aaID = rxn_params.loc[ rxn_params["Parameter Type"] == 'amino acid' ]["Value"].values[0]

        synthetaseID = rxn_params.loc[ rxn_params["Parameter Type"] == 'synthetase' ]["Value"].values[0]

        synthetaseAtpID = synthetaseID + '_atp'

        csim.addReaction(tuple([synthetaseID, 'M_atp_c']), synthetaseAtpID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_atp' ]["Value"].values[0])

        synthetaseAaID = synthetaseAtpID + '_aa'

        csim.addReaction(tuple([synthetaseAtpID, aaID]), synthetaseAaID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_aa' ]["Value"].values[0])

        for rnaID in rnaIDlist:

            synthetaseTrnaID = synthetaseAaID + '_' + rnaID

            csim.addReaction(tuple([synthetaseAaID, rnaID]), synthetaseTrnaID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_tRNA' ]["Value"].values[0])

            chargedTrnaID = rnaID + '_ch'

            csim.addReaction(synthetaseTrnaID, tuple(['M_amp_c', 'M_ppi_c', synthetaseID, chargedTrnaID]), rxn_params.loc[ rxn_params["Parameter Type"] == 'k_cat' ]["Value"].values[0])

            costID = tRNA_aa + '_cost'
            costPaidID = costID + '_paid'

            csim.addReaction(tuple([costID, chargedTrnaID]), tuple([rnaID, costPaidID]), 1e5)
    
    print('tRNA')
    
    return None
#########################################################################################



    
    
    
