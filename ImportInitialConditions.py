"""
Authors: Zane Thornburg

Sets the initial conditions for particle counts.  

Does not include initlialization of intermediates, for example bound states of RNAP to genes or Ribosomes to mRNA.
"""

import pandas as pd
import numpy as np
import random

import Diffusion as Diff


#########################################################################################
def initializeParticles(sim, region_dict, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    Diff.defaultDiffusion(sim, region_dict)
    
    Diff.generalDiffusionConstants(sim)
    
    initializeDnaParticles(sim, sim_properties)
    
    initializeProteins(sim, sim_properties)
    
    initializeMRNA(sim, sim_properties)

    initializeTRNA(sim, sim_properties)
    
    initializeRRNA(sim, sim_properties)
    
    initializeRNAP(sim)
    
    initializeRibosomeParticles(sim, region_dict)
    
    initializeDegradosomes(sim, region_dict)
    
    initializePromoterStrengths(sim_properties)
    
    initializeMedium(sim_properties)
    
    initializeMetabolites(sim_properties)
    
    initializeCostCounters(sim_properties)
    
    initializeProteinMetabolites(sim_properties)
    
    initializeTrnaCharging(sim_properties)
    
    initializeLongRnaTracking(sim_properties)
    
    print('Initialized Particle Counts and Diffusion Rules')
    
    return None
#########################################################################################


#########################################################################################
def initializeMedium(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    sim_properties['medium'] = {}
    
    sim_medium = pd.read_excel(sim_properties['head_directory'] + 'input_data/initial_concentrations.xlsx', sheet_name='Simulation Medium')
    
    for row, nutrient in sim_medium.iterrows():
        
        metID = 'M_' + nutrient['Met ID']
        
        sim_properties['medium'][metID] = nutrient['Conc (mM)']
        
#         sim_properties['counts'][metID] = mMtoPart(nutrient['Conc (mM)'])
        
    return None
#########################################################################################


#########################################################################################
def initializeMetabolites(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    metabolite_ic = pd.read_excel(sim_properties['head_directory'] + 'input_data/initial_concentrations.xlsx', sheet_name='Intracellular Metabolites')

    for row, metabolite in metabolite_ic.iterrows():
        
        metID = 'M_' + metabolite['Met ID']
        
        sim_properties['counts'][metID] = mMtoPart(metabolite['Init Conc (mM)'], sim_properties)
        
    return None
#########################################################################################


#########################################################################################
def initializeProteins(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    protein_ic = pd.read_excel(sim_properties['head_directory'] + 'input_data/initial_concentrations.xlsx', sheet_name='Comparative Proteomics')
    
    for row, protein in protein_ic.iterrows():
        
        if row != 0:
        
            locusTag = protein['Locus Tag']

            locusNum = locusTag.split('_')[1]

            proteinID = 'P_'+locusNum
            
            PTN = sim.species(proteinID)

            initial_count = protein['Sim. Initial Ptn Cnt']
            
            sim_properties['counts']['PM_'+locusNum] = 0
            
            if 'S ribosomal' in protein['Gene Product']:
                
                initial_count = max(25, int(initial_count - 500))
            
            if protein['Localization'] == 'trans-membrane':
                
                cytoID = 'C_P_'+locusNum
                
                cPTN = sim.species(cytoID)

                sim.distributeNumber(PTN, sim.region('membrane'), int(initial_count))
                
                Diff.proteinDiffusionCyto(sim, cytoID)

                Diff.proteinDiffusionMem(sim, proteinID)
                
            elif protein['Localization'] == 'peripheral membrane':
                
                sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(initial_count))
                
                Diff.proteinDiffusionPeriph(sim, proteinID)
                
#             elif protein['Gene Product'] 

            else:

                sim.distributeNumber(PTN, sim.region('cytoplasm'), int(initial_count/2))
                sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(initial_count/6))
                sim.distributeNumber(PTN, sim.region('DNA'), int(initial_count/3))

                Diff.proteinDiffusionCyto(sim, proteinID)
        
    return None
#########################################################################################


#########################################################################################
def initializeMRNA(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    mRNA_ic = pd.read_excel(sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx',sheet_name='Diffusion Coefficient')
    
    for row, mRNA in mRNA_ic.iterrows():
        
        locusTag = mRNA['Locus Tag']
        
        if locusTag.startswith('JCVI'):
            
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum
            
            RNA = sim.species(rnaID)
        
            avg_count = mRNA['mRNA Avg Count (#)']

            if avg_count == 0.0:
                avg_count = 0.0001

            init_mRNA_count = np.random.poisson(avg_count*2)

            sim.distributeNumber(RNA, sim.region('cytoplasm'), int(init_mRNA_count))
            
            rnasequence = sim_properties['genome'][locusTag]['RNAsequence']

            Diff.rnaDiffusion(sim, rnaID, rnasequence)
            
            rnaReadID = 'R_'+locusNum+'_d'
            
            rnaRead = sim.species(rnaReadID)
            
            Diff.rnaDiffusion(sim, rnaReadID, rnasequence)
    
#     print('RNA time')
    
    return None
#########################################################################################


#########################################################################################
def initializeTRNA(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    genome = sim_properties['genome']
    
    tRNA_map = {}
    
    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == "tRNA":
            
#             rnaID = locusDict['GeneName']
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum
            
            tRNA = sim.species(rnaID)
            
            sim.distributeNumber(tRNA, sim.region('cytoplasm'), int(200/3))
            sim.distributeNumber(tRNA, sim.region('outer_cytoplasm'), int(200/6))
            sim.distributeNumber(tRNA, sim.region('DNA'), int(200/2))
            
            rnasequence = sim_properties['genome'][locusTag]['RNAsequence']

            Diff.rnaDiffusion(sim, rnaID, rnasequence)
            
            rnaName = locusDict['GeneName'].split('-')[1].upper()
            
            if rnaName not in tRNA_map:
                
                tRNA_map[rnaName] = [rnaID]
                
            else:
                
                tRNA_map[rnaName].append(rnaID)
                
    sim_properties['trna_map'] = tRNA_map
                
#     print('Mapped tRNAs')
#     print(tRNA_map)
            
    return None
#########################################################################################


#########################################################################################
def initializeRRNA(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    genome = sim_properties['genome']
    
    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == "rRNA":
            
#             rnaID = locusDict['GeneName']
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum
            
            rRNA = sim.species(rnaID)
            
#             sim.distributeNumber(tRNA, sim.region('cytoplasm'), int(200/3))
#             sim.distributeNumber(tRNA, sim.region('outer_cytoplasm'), int(200/6))
#             sim.distributeNumber(tRNA, sim.region('DNA'), int(200/2))
            
            rnasequence = sim_properties['genome'][locusTag]['RNAsequence']

            Diff.rnaDiffusion(sim, rnaID, rnasequence)
            
    return None
#########################################################################################


#########################################################################################
def initializeRibosomeParticles(sim, region_dict):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    ribo_part_coords = np.argwhere(region_dict['ribo_centers']['shape']==True)
    
    ribo_part = sim.species('ribosomeP')

    for coord in ribo_part_coords:

        x = coord[0]
        y = coord[1]
        z = coord[2]

        ribo_part.placeParticle(x,y,z,1)
        
    Diff.ribosomeDiffusion(sim, 'ribosomeP')
        
    return None
#########################################################################################


#########################################################################################
def initializeRNAP(sim):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    RNApol = sim.species('RNAP')

    diffRNApol = sim.diffusionConst("RNAP_diff", 0.22e-12)

#     sim.distributeNumber(RNApol, sim.region('cytoplasm'), 187)

    RNApol.diffusionRate(sim.region('DNA'),diffRNApol)
    RNApol.diffusionRate(sim.region('cytoplasm'),diffRNApol)
    RNApol.diffusionRate(sim.region('outer_cytoplasm'),diffRNApol)

    sim.transitionRate(RNApol, sim.region('cytoplasm'), sim.region('DNA'), diffRNApol)
    sim.transitionRate(RNApol, sim.region('DNA'), sim.region('cytoplasm'), diffRNApol)
    sim.transitionRate(RNApol, sim.region('cytoplasm'), sim.region('outer_cytoplasm'), diffRNApol)
    sim.transitionRate(RNApol, sim.region('outer_cytoplasm'), sim.region('cytoplasm'), diffRNApol)
    sim.transitionRate(RNApol, sim.region('outer_cytoplasm'), sim.region('DNA'), diffRNApol)
    sim.transitionRate(RNApol, sim.region('DNA'), sim.region('outer_cytoplasm'), diffRNApol)
    
    print('Added RNAP particles')
    
    return None
#########################################################################################


#########################################################################################
def initializeDegradosomes(sim, region_dict):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    degradosome = sim.species('Degradosome')
    
    peripheral_membrane = np.argwhere(region_dict['outer_cytoplasm']['shape']==True)
    
#     sim.distributeNumber(degradosome, sim.region('outer_cytoplasm'), int(120))

#     for i in range(120):
#         coord = random.choice(peripheral_membrane)

#         x = int(coord[0])
#         y = int(coord[1])
#         z = int(coord[2])

#         degradosome.placeParticle(x,y,z,1)
        
    Diff.degradosomeDiffusion(sim, 'Degradosome')
    
    RnaseY = sim.species('P_0359')
    RnaseJ1 = sim.species('P_0600')
    ptnAssoc = sim.rateConst('degAssoc', 1e7, 2)
    
    sim.region('outer_cytoplasm').addReaction([RnaseY, RnaseJ1], [degradosome], ptnAssoc)
#     Diff.proteinDiffusionPeriph(sim, 'Degradosome')

    print('Degradosomes placed')
    
    return None
#########################################################################################


#########################################################################################
def initializePromoterStrengths(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    proxyPromoterStrengths = {}
    
    genome = sim_properties['genome']
    
    proteomics = pd.read_excel(sim_properties['head_directory'] + 'input_data/initial_concentrations.xlsx', sheet_name='Comparative Proteomics')
    
    for locusTag, locusDict in genome.items():
        
        if locusDict["Type"] == 'protein':
            
            PtnName = proteomics.loc[ proteomics['Locus Tag'] == locusTag ]['Gene Product'].values[0]
            
            PtnCount = proteomics.loc[ proteomics['Locus Tag'] == locusTag ]['Sim. Initial Ptn Cnt'].values[0]
            
            if 'S ribosomal' in PtnName:
                
                proxyPromoterStrengths[locusTag] = min(765, 500+PtnCount)
                
            else:
            
                proxyPromoterStrengths[locusTag] = min(765, max(45, PtnCount))
                
        elif locusDict["Type"] == 'tRNA':
            
#             proxyPromoterStrengths[locusTag] = 180
            proxyPromoterStrengths[locusTag] = 765
            
        elif locusDict["Type"] == 'rRNA':
            
            proxyPromoterStrengths[locusTag] = 765*6.8
            
    sim_properties['promoters'] = proxyPromoterStrengths
    
    return None
#########################################################################################


#########################################################################################
def initializeDnaParticles(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
                        
    genome = sim_properties['genome']
    DNAcoords = sim_properties['DNAcoords']

    for locusTag, locusDict in genome.items():
        
        locusNum = locusTag.split('_')[1]
        
        gene = sim.species('G_' + locusNum + '_C1')
        
        start = int(locusDict['startIndex'][0])
        end = int(locusDict['endIndex'][0])
        
        startXYZ = DNAcoords[start]
        
        gene.placeParticle(int(startXYZ[0]), int(startXYZ[1]), int(startXYZ[2]), 1)
        
    for feature, fdict in sim_properties['chromosome_features'].items():
        
        dnaParticle = sim.species(feature)
        
        partXYZ = DNAcoords[int(fdict['index'])]
        
        dnaParticle.placeParticle(int(partXYZ[0]), int(partXYZ[1]), int(partXYZ[2]), 1)
        
    sim_properties['counts']['chromosome'] = int(len(DNAcoords))
        
    print('Placed chemically relevant DNA particles')
        
    return None
#########################################################################################


#########################################################################################
def initializeCostCounters(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    cost_counters = {'GTP_translat_cost':'M_atp_c', 'ATP_trsc_cost':'M_atp_c', 'ATP_mRNAdeg_cost':'M_atp_c',
                     'ATP_DNArep_cost':'M_atp_c', 'ATP_transloc_cost':'M_atp_c',
                     'ATP_mRNA_cost':'M_atp_c', 'CTP_mRNA_cost':'M_ctp_c', 'UTP_mRNA_cost':'M_utp_c', 'GTP_mRNA_cost':'M_gtp_c',
                     'ATP_tRNA_cost':'M_atp_c', 'CTP_tRNA_cost':'M_ctp_c', 'UTP_tRNA_cost':'M_utp_c', 'GTP_tRNA_cost':'M_gtp_c',
                     'ATP_rRNA_cost':'M_atp_c', 'CTP_rRNA_cost':'M_ctp_c', 'UTP_rRNA_cost':'M_utp_c', 'GTP_rRNA_cost':'M_gtp_c',
                     'dATP_DNArep_cost':'M_datp_c', 'dTTP_DNArep_cost':'M_dttp_c', 'dCTP_DNArep_cost':'M_dctp_c', 'dGTP_DNArep_cost':'M_dgtp_c'}
    
    for cost, metID in cost_counters.items():
        
        sim_properties['counts'][cost] = 0
        
        sim_properties['counts'][cost+'_paid'] = 0
        
        sim_properties['counts'][cost+'_second'] = 0
        
#     sim_properties['cost_counters'] = cost_counters
    
        
    aaCostMap = {"A":"ALA_cost", "R":"ARG_cost", 
        "N":"ASN_cost", "D":"ASP_cost", "C":"CYS_cost", "E":"GLU_cost", "Q":"GLN_cost", "G":"GLY_cost", 
        "H":"HIS_cost", "I":"ILE_cost", "L":"LEU_cost", "K":"LYS_cost", "M":"MET_cost", "F":"PHE_cost", 
        "P":"PRO_cost", "S":"SER_cost", "T":"THR_cost", "W":"TRP_cost", "Y":"TYR_cost", "V":"VAL_cost"}
    
    for aa, aaCost in aaCostMap.items():
            
        sim_properties['counts'][aaCost] = 0
        
        sim_properties['counts'][aaCost+'_paid'] = 0
        
        sim_properties['counts'][aaCost+'_second'] = 0
    
    sim_properties['counts']['FMET_cost'] = 0
    
    sim_properties['counts']['FMET_cost_paid'] = 0
    
    sim_properties['counts']['FMET_cost_second'] = 0
    
    
    NMP_recycle_counters = {'AMP_mRNAdeg_cost':'M_amp_c', 'UMP_mRNAdeg_cost':'M_ump_c', 'CMP_mRNAdeg_cost':'M_cmp_c', 'GMP_mRNAdeg_cost':'M_gmp_c'}

    for recID, metID in NMP_recycle_counters.items():
        
        sim_properties['counts'][recID] = 0
        
        sim_properties['counts'][recID+'_paid'] = 0
        
        sim_properties['counts'][recID+'_second'] = 0
        
#     aaCostMap = {"A":"ALA_cost","R":"ARG_cost","N":"ASN_cost","D":"ASP_cost","C":"CYS_cost",
#                          "E":"GLU_cost","Q":"GLN_cost","G":"GLY_cost","H":"HIS_cost","I":"ILE_cost",
#                          "L":"LEU_cost","K":"LYS_cost","M":"MET_cost","F":"PHE_cost","P":"PRO_cost",
#                          "S":"SER_cost","T":"THR_cost","W":"TRP_cost","Y":"TYR_cost","V":"VAL_cost",
#                          "FM":"FMET_cost"}
    
#     for aa, costID in aaCostMap.items():
        
#         sim_properties['counts'][costID] = 0
    
    return None
#########################################################################################


#########################################################################################
def initializeProteinMetabolites(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """

    data_file = sim_properties['head_directory'] + 'input_data/protein_metabolites.xlsx'
    
    ptnMets = pd.read_excel(data_file, sheet_name='protein metabolites')
    
    for index, row in ptnMets.iterrows():
        
        metabolites = row['Metabolite IDs'].split(',')
        
#         print(metabolites)
        
        for metID in metabolites:
            
            sim_properties['counts'][metID] = 0
        
    return None
#########################################################################################


#########################################################################################
def initializeTrnaCharging(sim_properties):
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

        synthetaseID = rxn_params.loc[ rxn_params["Parameter Type"] == 'synthetase' ]["Value"].values[0]

        synthetaseAtpID = synthetaseID + '_atp'

        sim_properties['counts'][synthetaseAtpID] = 0

        synthetaseAaID = synthetaseAtpID + '_aa'

        sim_properties['counts'][synthetaseAaID] = 0

        for rnaID in rnaIDlist:

            synthetaseTrnaID = synthetaseAaID + '_' + rnaID

            sim_properties['counts'][synthetaseTrnaID] = 0

            chargedTrnaID = rnaID + '_ch'

            sim_properties['counts'][chargedTrnaID] = 0
                
    return None
#########################################################################################


#########################################################################################
def setCmeSpeciesList(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    CmeSpeciesList = []
    
    for specID, count in sim_properties['counts'].items():
        
        if specID.startswith('RP_'):
            
            CmeSpeciesList.append(specID)
            
    CmeSpeciesList.append('M_atp_c')
    CmeSpeciesList.append('M_adp_c')
    CmeSpeciesList.append('M_amp_c')
    CmeSpeciesList.append('M_ppi_c')
    CmeSpeciesList.append('M_pi_c')
       
    trna_RXNS_params = pd.read_excel(sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx', sheet_name='tRNA Charging')
    
    CmeStatesList = []
            
    for tRNA_aa, rnaIDlist in sim_properties['trna_map'].items():
        
#         if tRNA_aa != 'GLN':
        
        rxnID = tRNA_aa + 'TRS'

        rxn_params = trna_RXNS_params.loc[ trna_RXNS_params["Reaction Name"] == rxnID ]

        aaID = rxn_params.loc[ rxn_params["Parameter Type"] == 'amino acid' ]["Value"].values[0]

        CmeSpeciesList.append(aaID)

        synthetaseID = rxn_params.loc[ rxn_params["Parameter Type"] == 'synthetase' ]["Value"].values[0]

        CmeSpeciesList.append(synthetaseID)
        CmeStatesList.append(synthetaseID)

        synthetaseAtpID = synthetaseID + '_atp'

        CmeSpeciesList.append(synthetaseAtpID)

        synthetaseAaID = synthetaseAtpID + '_aa'

        CmeSpeciesList.append(synthetaseAaID)

        for rnaID in rnaIDlist:

            CmeSpeciesList.append(rnaID)
            CmeStatesList.append(rnaID)

            synthetaseTrnaID = synthetaseAaID + '_' + rnaID

            CmeSpeciesList.append(synthetaseTrnaID)

            chargedTrnaID = rnaID + '_ch'

            CmeSpeciesList.append(chargedTrnaID)

            costID = tRNA_aa + '_cost'
            costPaidID = costID + '_paid'

            CmeSpeciesList.append(costID)
            CmeSpeciesList.append(costPaidID)

    CmeSpeciesSet = [*set(CmeSpeciesList)]

    sim_properties['cme_species'] = CmeSpeciesSet
    
    CmeStatesSet = [*set(CmeStatesList)]
    
    CmeStatesTracker = {}
    
    for specID in CmeStatesSet:
        
        CmeStatesTracker[specID] = 0
        
    sim_properties['cme_state_tracker'] = CmeStatesTracker
    
    return None
#########################################################################################


#########################################################################################
def initializeLongRnaTracking(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    long_rna_trsc = {}
    
    long_rna_trsc['genome1'] = {}
    
    long_rna_trsc['genome2'] = {}
    
    genome = sim_properties['genome']
    
    for locusTag, locusDict in genome.items():
        
        locusNum = locusTag.split('_')[1]
        
        rnasequence = locusDict["RNAsequence"]
        
        if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
            
            long_rna_trsc[locusTag] = {}
            
            long_rna_trsc['genome1'][locusTag] = {}
            long_rna_trsc['genome2'][locusTag] = {}
            
            MaxRNAP = int(int(len(rnasequence)/sim_properties['rnap_spacing']))
            
            long_rna_trsc[locusTag]['max_rnap'] = MaxRNAP
            
            long_rna_trsc[locusTag]['rnap_count'] = 0
            
#             long_rna_trsc['genome1'][locusTag]['positions'] = {}
            
#             long_rna_trsc['genome2'][locusTag]['positions'] = {}
            
            long_rna_trsc['genome1'][locusTag]['state'] = 'G_' + locusNum + '_C1'
            
            long_rna_trsc['genome2'][locusTag]['state'] = 'G_' + locusNum + '_C2'
            
#             for i in range(MaxRNAP+1):
                
#                 long_rna_trsc['genome1'][locusTag]['positions'][str(i)] = 0
                
#                 long_rna_trsc['genome2'][locusTag]['positions'][str(i)] = 0
                
    sim_properties['long_rna_trsc'] = long_rna_trsc
    
#     print(sim_properties['long_rna_trsc'])
    
    return None
#########################################################################################


#########################################################################################
def initializeRdmeCounts(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    RDME_counts = sim.particleStatistics(particleLattice=sim.particleLattice, siteLattice=sim.siteLattice)
    
    total_proteins = 0
    
    for name, index in sim_properties['name_to_index'].items():
    
        count = RDME_counts['countBySpecies'][sim.species(name)]
        
        sim_properties['counts'][name] = count
        
        if name.startswith('P_'):
            
            total_proteins = total_proteins + count
            
    sim_properties['counts']['total_proteins_made'] = 0 # int(total_proteins)
    
    return None
#########################################################################################

            
#########################################################################################
def mMtoPart(conc, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's
    
    particles = int(round((conc/1000)*NA*sim_properties['volume_L']))

    return particles
#########################################################################################
