"""
Authors: Zane Thornburg

Functions that communicate particle counts at hook times.
"""

import numpy as np

import pandas as pd

import os

import json
import pickle

import pySTDLM.PostProcessing as PP

import time as timepy

from LatticeFunctions import *

from scipy.optimize import fsolve
from scipy.optimize import least_squares

import MC_CME as MCCME

#########################################################################################
def updateCountsRDME(RDMEsim, sim_properties, lattice):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    Called by:
    Description:
    """
    
    RDME_counts = RDMEsim.particleStatistics(particleLattice=lattice.getParticleLatticeView(), siteLattice=lattice.getSiteLatticeView())
    
    for name, index in sim_properties['name_to_index'].items():
    
        new_count = RDME_counts['countBySpecies'][RDMEsim.species(name)]
        
        sim_properties['counts'][name] = new_count
        
    print('Updated particle counts from RDME')
    
    return None
#########################################################################################


#########################################################################################
def updateCountsCME(sim_properties):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    csimFolder = sim_properties['working_directory']+'CME/'
    CSIMfilename= csimFolder + 'cmeSim.%d.lm'%np.rint(sim_properties['time'])
    
#     print(CMEspecies)

    bad_cme=False
    while not(bad_cme):
        try: #TAB CHANGE
            CMEsim=PP.openLMFile(CSIMfilename)
            CMEspecies = PP.getSpecies(CMEsim)
            count_trace = PP.getSpecieTrace(CMEsim, CMEspecies[0])
            PP.closeLMFile(CMEsim)
            bad_cme=True
        except:
            timepy.sleep(30)
            print('Error in CME run, re-building and re-running simulations')
            MCCME.runGCME(sim_properties)

            
    CMEsim=PP.openLMFile(CSIMfilename)
    
    CMEspecies = PP.getSpecies(CMEsim)
    
    for specie in CMEspecies:
        
        count_trace = PP.getSpecieTrace(CMEsim, specie)

        count = count_trace[-1]
        
        if specie in sim_properties['cme_state_tracker']:
            
#             print('Updating specie: ', specie)
            
            sim_properties['cme_state_tracker'][specie] = sim_properties['counts'][specie] - count

        sim_properties['counts'][specie] = count

    PP.closeLMFile(CMEsim)
    
    print('Updated particle counts from global CME')
    
    try:
        os.remove(CSIMfilename)
        print('Removed gCME File')
    except:
        print('Nothing to delete')
    
    return None
#########################################################################################


#########################################################################################
def updateCountsODE(sim_properties, odeResults, model):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    resFinal = odeResults[-1,:]
    
    # Get the list of metabolites
    mL = model.getMetList()

    # For loop iterating over the species and updating their counts with ODE results
    for ind in range(len(mL)):
        
        specieName = mL[ind].getID()
        
        if not specieName.endswith('_e'):
        
            sim_properties['counts'][mL[ind].getID()] = mMtoPart(resFinal[ind],sim_properties) # Assign updated species counts to particle map using species IDs
    
    return None
#########################################################################################


#########################################################################################
def updateTranscriptionStates(sim_properties, lattice):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    Called by:
    Description:
    """
    
    print('Updating Transcription States and Placing New RNAs')
    
    plattice = lattice.getParticleLatticeView()
    
    genome = sim_properties['genome']
    
#     updateLongGeneStates(sim_properties, plattice)
    
    for locusTag, locusDict in genome.items():
        
        rnasequence = locusDict["RNAsequence"]
        
        if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
            
            updateRNAstateLong(sim_properties, lattice, plattice, locusTag)
            
        else:
        
            updateRNAstateShort(sim_properties, lattice, plattice, locusTag)

    print('Updated Transcription States')
    
    return None
#########################################################################################


#########################################################################################
def updateRNAstateShort(sim_properties, lattice, plattice, locusTag):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    Called by:
    Description:
    """
    
    locusDict = sim_properties['genome'][locusTag]
    
    locusNum = locusTag.split('_')[1]
    
    RnapIdx = sim_properties['name_to_index']['RNAP']
    
    starts = locusDict['startIndex']
    
    ends = locusDict['endIndex']
    
#     print(starts)
    
    for i in range(len(starts)):
        
        chromo = i + 1

        try:
            start = starts[i]
            end = ends[i]
        except:
            continue
        
        finishedRNAP = 'RP_' + locusNum + '_f' + '_C' + str(int(chromo))
        
#         print('RNAP to update: ', finishedRNAP)

        if sim_properties['counts'][finishedRNAP] > 0:
            
            print(finishedRNAP)

            rnasequence = locusDict["RNAsequence"]

    #             if len(rnasequence) > 500:

    #                 RnapID = 'RP_' + locusNum + '_t'

    #             else:

            RnapID = 'RP_' + locusNum + '_C' + str(int(chromo))

            TrscIdx = sim_properties['name_to_index'][RnapID]

            startXYZ = sim_properties['DNAcoords'][start]

            correctChromosome = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), TrscIdx)
            
            if not correctChromosome:
                
                print('Bad check')

            if correctChromosome:

                endXYZ = sim_properties['DNAcoords'][end]

                endSiteOccupancy = lattice.getOccupancy(int(endXYZ[2]), int(endXYZ[1]), int(endXYZ[0]))

                if endSiteOccupancy < 14:

                    deleteParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), TrscIdx)

#                         if len(rnasequence) <= 500:

                    GeneIdx = sim_properties['name_to_index']['G_'+locusNum+'_C'+str(int(chromo))]

                    lattice.addParticle(int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), GeneIdx)

                    RnaIdx = sim_properties['name_to_index']['R_'+locusNum]

                    lattice.addParticle(int(endXYZ[2]), int(endXYZ[1]), int(endXYZ[0]), RnapIdx)

                    lattice.addParticle(int(endXYZ[2]), int(endXYZ[1]), int(endXYZ[0]), RnaIdx)

                    sim_properties['counts'][finishedRNAP] = max(0, sim_properties['counts'][finishedRNAP] - 1)
                    
                    RNAP_made = 'RPM_' + locusNum
    
                    sim_properties['counts'][RNAP_made] = sim_properties['counts'][RNAP_made] + 1

                    calculateTranscriptionCosts(sim_properties, lattice, locusTag)

                    print('Placed New RNA: ', 'R_'+locusNum)
                    
    return None
#########################################################################################


#########################################################################################
def updateRNAstateLong(sim_properties, lattice, plattice, locusTag):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    locusDict = sim_properties['genome'][locusTag]
    
    locusNum = locusTag.split('_')[1]
    
    RnapIdx = sim_properties['name_to_index']['RNAP']
    
    starts = locusDict['startIndex']
    ends = locusDict['endIndex']
    
    for chrom in range(len(starts)):

        try:
            start = starts[chrom]
            end = ends[chrom]
        except:
            continue
            
        chromo = chrom + 1
        
        last_state = sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state']
        
        print(last_state)
        
        if last_state.startswith('G_'):
            
            continue
            
        startXYZ = sim_properties['DNAcoords'][start]
        
        endXYZ = sim_properties['DNAcoords'][end]
        
        RNAP_count = 0
        
        for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):

            on_i = 'RP_' + locusNum + '_c' + str(chromo)+ '_' + str(i)
            
            RNAP_count = RNAP_count + sim_properties['counts'][on_i]
            
            if i == sim_properties['long_rna_trsc'][locusTag]['max_rnap']:
                
                on_ip1 = 'RP_' + locusNum  + '_c' + str(chromo) + '_done'
                
                RNAP_count = RNAP_count + sim_properties['counts'][on_ip1]
                
        RNAP_count = int(RNAP_count)
        
        print(RNAP_count)

#         for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):
            
#             RNAP_count = int(last_state.split('_')[-1])
            
        updateStartState = 'RP_' + locusNum + '_c' + str(chromo) + '_1' # + str(int(1))

        print(updateStartState)

        updateEndState = 'RP_' + locusNum  + '_c' + str(chromo) + '_done'

        if sim_properties['counts'][updateEndState] > 0:

            endSiteOccupancy = lattice.getOccupancy(int(endXYZ[2]), int(endXYZ[1]), int(endXYZ[0]))

            if endSiteOccupancy < 14:

                RNAP_count = RNAP_count - 1
                
                TrscIdx = sim_properties['name_to_index'][last_state]

                deleteParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), TrscIdx)
                
                if RNAP_count == 0:

                    RP_gene = 'G_' + locusNum + '_C' + str(int(chromo))

                elif ('_t_' not in last_state) and (sim_properties['counts'][updateStartState] > 0):

                    RP_gene = 'RP_' + locusNum + '_' + str(int(RNAP_count)) + '_C' + str(int(chromo))

                else:

                    RP_gene = 'RP_' + locusNum + '_t_' + str(int(RNAP_count)) + '_C' + str(int(chromo))

                RP_g_idx = sim_properties['name_to_index'][RP_gene]

                lattice.addParticle(int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), RP_g_idx)

                RnaIdx = sim_properties['name_to_index']['R_'+locusNum]

                lattice.addParticle(int(endXYZ[2]), int(endXYZ[1]), int(endXYZ[0]), RnaIdx)

                lattice.addParticle(int(endXYZ[2]), int(endXYZ[1]), int(endXYZ[0]), RnapIdx)

                sim_properties['counts'][updateEndState] = max(0, sim_properties['counts'][updateEndState] - 1)
                
                sim_properties['counts']['RP_' + locusNum  + '_c' + str(chromo) + '_endTrsc'] = 1
                
                RNAP_made = 'RPM_' + locusNum
    
                sim_properties['counts'][RNAP_made] = sim_properties['counts'][RNAP_made] + 1

                calculateTranscriptionCosts(sim_properties, lattice, locusTag)

                print('Placed New RNA: ', 'R_'+locusNum)

            continue

        if ('_t_' not in last_state) and (sim_properties['counts'][updateStartState] == 0):
            
            print('Updating start ', locusTag)

#             RNAP_gene = 'RP_' + locusNum + '_' + str(RNAP_count)

            TrscIdx = sim_properties['name_to_index'][last_state]
            
#             print(lattice.getOccupancy(int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0])))

            deleteParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), TrscIdx)
            
#             print(lattice.getOccupancy(int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0])))

            updateIdx = sim_properties['name_to_index']['RP_' + locusNum + '_t_' + str(RNAP_count) + '_C' + str(int(chromo))]

            lattice.addParticle(int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), updateIdx)
            
#             print(lattice.getOccupancy(int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0])))

#                 sim_properties['counts'][updateState] = max(0, sim_properties['counts'][updateState] - 1)
                    
    return None
#########################################################################################


#########################################################################################
def updateLongGeneStates(sim_properties, lattice):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    Called by:
    Description:
    """
    
    plattice = lattice.getParticleLatticeView()
    
    genome = sim_properties['genome']
    
    for locusTag, locusDict in genome.items():
        
        locusNum = locusTag.split('_')[1]
        
        rnasequence = locusDict["RNAsequence"]
        
        if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
    
            starts = locusDict['startIndex']

            for chrom in range(len(starts)):
                
                chromo = chrom + 1

                try:
                    start = starts[chrom]
    #                 end = ends[i]
                except:
                    continue

                startXYZ = sim_properties['DNAcoords'][start]

                gene = 'G_' + locusNum + '_C' + str(int(chromo))

                g_idx = sim_properties['name_to_index'][gene]

                stateCheck = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), g_idx)

                if stateCheck:

                    sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = gene

                    continue

                for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):

                    RNAP_gene = 'RP_' + locusNum + '_' + str(i) + '_C' + str(int(chromo))

                    RP_g_idx = sim_properties['name_to_index'][RNAP_gene]

                    stateCheck = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), RP_g_idx)

                    if stateCheck:

                        sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RNAP_gene

                        sim_properties['counts']['RP_' + locusNum + '_c' + str(chromo) + '_' + str(1)] = 1

                        continue

                    RNAP_gene_t = 'RP_' + locusNum + '_t_' + str(i) + '_C' + str(int(chromo))

                    RP_g_t_idx = sim_properties['name_to_index'][RNAP_gene_t]

                    stateCheck = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), RP_g_t_idx)

                    if stateCheck:

                        sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RNAP_gene_t

                        continue
                        
#             parts = getParticlesInSite(plattice, int(z), int(y), int(x))
            
#             for partIdx in parts:
                
#                 deleteParticle(plattice,int(z),int(y),int(x),partIdx)

    return None
#########################################################################################


#########################################################################################
def calculateTranscriptionCosts(sim_properties, lattice, locusTag):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    Called by:
    Description:
    """
    
    plattice = lattice.getParticleLatticeView()
    
    genome = sim_properties['genome']
    
#     for locusTag, locusDict in genome.items():
    locusDict = genome[locusTag]
        
    locusNum = locusTag.split('_')[1]

#     finishedRNAP = 'RP_' + locusNum + '_f'

#     if sim_properties['counts'][finishedRNAP] > 0:

    rnasequence = locusDict['RNAsequence']

    sim_properties['counts']['ATP_trsc_cost'] = sim_properties['counts']['ATP_trsc_cost'] + len(rnasequence)
    
    sim_properties['counts']['ATP_trsc_cost_second'] = sim_properties['counts']['ATP_trsc_cost_second'] + len(rnasequence)

    if locusDict['Type'] == 'protein':

        RNAtype = 'mRNA'

    else:

        RNAtype = locusDict['Type']

    for nucleotide in set(rnasequence):

        nucleotideCount = rnasequence.count(nucleotide)

        sim_properties['counts'][nucleotide + 'TP_' + RNAtype + '_cost'] = sim_properties['counts'][nucleotide + 'TP_' + RNAtype + '_cost'] + nucleotideCount
        
        sim_properties['counts'][nucleotide + 'TP_' + RNAtype + '_cost_second'] = sim_properties['counts'][nucleotide + 'TP_' + RNAtype + '_cost_second'] + nucleotideCount
            
    return None
#########################################################################################            


#########################################################################################
def calculateTranslationCosts(sim_properties, lattice):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    Called by:
    Description:
    """
    
    plattice = lattice.getParticleLatticeView()
    
    genome = sim_properties['genome']
    
    aaCostMap = {"A":"ALA_cost", "R":"ARG_cost", 
    "N":"ASN_cost", "D":"ASP_cost", "C":"CYS_cost", "E":"GLU_cost", "Q":"GLN_cost", "G":"GLY_cost", 
    "H":"HIS_cost", "I":"ILE_cost", "L":"LEU_cost", "K":"LYS_cost", "M":"MET_cost", "F":"PHE_cost", 
    "P":"PRO_cost", "S":"SER_cost", "T":"THR_cost", "W":"TRP_cost", "Y":"TYR_cost", "V":"VAL_cost"}
#     ,"*":"Stop_Codon"}
    
    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == 'protein':

            locusNum = locusTag.split('_')[1]

            translatCostParticle = 'P_' + locusNum + '_TC'

            if sim_properties['counts'][translatCostParticle] > 0:

                TCPidx = sim_properties['name_to_index'][translatCostParticle]

                TCPcoords = np.argwhere(plattice==TCPidx)

                newPtns = len(TCPcoords)

                aasequence = locusDict['AAsequence']
                
                sim_properties['counts']['GTP_translat_cost'] = sim_properties['counts']['GTP_translat_cost'] + newPtns*len(aasequence)*2
                
                sim_properties['counts']['GTP_translat_cost_second'] = sim_properties['counts']['GTP_translat_cost_second'] + newPtns*len(aasequence)*2

                for aa, aaCost in aaCostMap.items():

                    sim_properties['counts'][aaCost] = sim_properties['counts'][aaCost] + newPtns*aasequence.count(aa)
                    
                    sim_properties['counts'][aaCost+'_second'] = sim_properties['counts'][aaCost+'_second'] + newPtns*aasequence.count(aa)

                sim_properties['counts']['FMET_cost'] = sim_properties['counts']['FMET_cost'] + newPtns
                
                sim_properties['counts']['FMET_cost_second'] = sim_properties['counts']['FMET_cost_second'] + newPtns

                for coord in TCPcoords:
                    
#                     print('locusNum')
                    
#                     print(lattice.getOccupancy(int(coord[3]), int(coord[2]), int(coord[1])))

                    deleteParticle(plattice, int(coord[3]), int(coord[2]), int(coord[1]), TCPidx)
                    
#                     print(lattice.getOccupancy(int(coord[3]), int(coord[2]), int(coord[1])))
                    
                    if 'C_P_'+locusNum in sim_properties['name_to_index']:
                        
#                         print('C_P_'+locusNum)
                    
                        lattice.addParticle(int(coord[3]), int(coord[2]), int(coord[1]), int(sim_properties['name_to_index']['C_P_'+locusNum]))
        
                        sim_properties['counts']['PM_'+locusNum] = sim_properties['counts']['PM_'+locusNum] + 1
        
                        sim_properties['counts']['total_proteins_made'] = sim_properties['counts']['total_proteins_made'] + 1
                    else:
                        
#                         print('P_'+locusNum)
                        
                        lattice.addParticle(int(coord[3]), int(coord[2]), int(coord[1]), int(sim_properties['name_to_index']['P_'+locusNum]))
        
                        sim_properties['counts']['PM_'+locusNum] = sim_properties['counts']['PM_'+locusNum] + 1
        
                        sim_properties['counts']['total_proteins_made'] = sim_properties['counts']['total_proteins_made'] + 1
                        
#                     print(lattice.getOccupancy(int(coord[3]), int(coord[2]), int(coord[1])))
                    
#                     deleteParticle(plattice, int(1), int(1), int(1), sim_properties['name_to_index']['P_0001'])

#                     lattice.addParticle(int(1), int(1), int(1), sim_properties['name_to_index']['P_0001'])

    return None
#########################################################################################


#########################################################################################
def calculateDegradationCosts(sim_properties, lattice):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    Called by:
    Description:
    """
    
    plattice = lattice.getParticleLatticeView()
    
    genome = sim_properties['genome']
    
    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == 'protein':

            locusNum = locusTag.split('_')[1]

            degCostParticle = 'DT_' + locusNum

            if sim_properties['counts'][degCostParticle] > 0:

                DTidx = sim_properties['name_to_index'][degCostParticle]

                DTcoords = np.argwhere(plattice==DTidx)

                mRNAdegraded = len(DTcoords)
                
                rnasequence = locusDict['RNAsequence']
            
                sim_properties['counts']['ATP_mRNAdeg_cost'] = sim_properties['counts']['ATP_mRNAdeg_cost'] + mRNAdegraded*len(rnasequence)
                
                sim_properties['counts']['ATP_mRNAdeg_cost_second'] = sim_properties['counts']['ATP_mRNAdeg_cost_second'] + mRNAdegraded*len(rnasequence)

                for nucleotide in set(rnasequence):

                    nucleotideCount = rnasequence.count(nucleotide)

                    sim_properties['counts'][nucleotide + 'MP_mRNAdeg_cost'] = sim_properties['counts'][nucleotide + 'MP_mRNAdeg_cost'] + mRNAdegraded*nucleotideCount
                    
                    sim_properties['counts'][nucleotide + 'MP_mRNAdeg_cost_second'] = sim_properties['counts'][nucleotide + 'MP_mRNAdeg_cost_second'] + mRNAdegraded*nucleotideCount

                for coord in DTcoords:
                    
#                     print('locusNum')
                    
#                     print(lattice.getOccupancy(int(coord[3]), int(coord[2]), int(coord[1])))

                    deleteParticle(plattice, int(coord[3]), int(coord[2]), int(coord[1]), DTidx)
                    
#                     print(lattice.getOccupancy(int(coord[3]), int(coord[2]), int(coord[1])))

                    lattice.addParticle(int(coord[3]), int(coord[2]), int(coord[1]), int(sim_properties['name_to_index']['decay']))
    
                    sim_properties['counts']['DM_' + locusNum] = sim_properties['counts']['DM_' + locusNum] + 1
    
    return None
#########################################################################################


#########################################################################################
def communicateCostsToMetabolism(sim_properties):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    energy_cost_counters = {'GTP_translat_cost':'g', 'ATP_trsc_cost':'a', 'ATP_mRNAdeg_cost':'a', 'ATP_DNArep_cost':'a', 'ATP_transloc_cost':'a'}
    
    nuc_costs = {'ATP_mRNA_cost':'M_atp_c', 'CTP_mRNA_cost':'M_ctp_c', 'UTP_mRNA_cost':'M_utp_c', 'GTP_mRNA_cost':'M_gtp_c',
                     'ATP_tRNA_cost':'M_atp_c', 'CTP_tRNA_cost':'M_ctp_c', 'UTP_tRNA_cost':'M_utp_c', 'GTP_tRNA_cost':'M_gtp_c',
                     'ATP_rRNA_cost':'M_atp_c', 'CTP_rRNA_cost':'M_ctp_c', 'UTP_rRNA_cost':'M_utp_c', 'GTP_rRNA_cost':'M_gtp_c',
                     'dATP_DNArep_cost':'M_datp_c', 'dTTP_DNArep_cost':'M_dttp_c', 'dCTP_DNArep_cost':'M_dctp_c', 'dGTP_DNArep_cost':'M_dgtp_c'}
    
    NMP_recycle_counters = {'AMP_mRNAdeg_cost':'M_amp_c', 'UMP_mRNAdeg_cost':'M_ump_c', 'CMP_mRNAdeg_cost':'M_cmp_c', 'GMP_mRNAdeg_cost':'M_gmp_c'}
    
    
    for cost, nucID in energy_cost_counters.items():
        
        costCount = sim_properties['counts'][cost]
        
#         sim_properties['counts'][cost+'_second'] = int(costCount)
        
        ntpID = 'M_' + nucID + 'tp_c'
        ndpID = 'M_' + nucID + 'dp_c'
        
        if sim_properties['counts'][ntpID] < costCount:
            
            sim_properties['counts'][cost] = costCount - sim_properties['counts'][ntpID] + 1
            
            sim_properties['counts'][cost+'_paid'] = sim_properties['counts'][ntpID] - 1
            
            sim_properties['counts'][ndpID] = sim_properties['counts'][ndpID] + sim_properties['counts'][ntpID] - 1
            sim_properties['counts']['M_pi_c'] = sim_properties['counts']['M_pi_c'] + sim_properties['counts'][ntpID] - 1
            
            sim_properties['counts'][ntpID] = 1
            
        else:
            
            sim_properties['counts'][ntpID] = sim_properties['counts'][ntpID] - costCount
            
            sim_properties['counts'][ndpID] = sim_properties['counts'][ndpID] + costCount
            sim_properties['counts']['M_pi_c'] = sim_properties['counts']['M_pi_c'] + costCount
            
            sim_properties['counts'][cost+'_paid'] = int(costCount)
            
            sim_properties['counts'][cost] = 0
            
            
    for cost, metID in nuc_costs.items():
        
        costCount = sim_properties['counts'][cost]
        
        if sim_properties['counts'][metID] < costCount:
            
            sim_properties['counts'][cost] = costCount - sim_properties['counts'][metID] + 1
            
            sim_properties['counts'][cost+'_paid'] = sim_properties['counts'][metID] - 1
            
            sim_properties['counts']['M_ppi_c'] = sim_properties['counts']['M_ppi_c'] + sim_properties['counts'][metID] - 1
            
            sim_properties['counts'][metID] = 1
            
        else:
            
            sim_properties['counts'][metID] = sim_properties['counts'][metID] - costCount
            
            sim_properties['counts']['M_ppi_c'] = sim_properties['counts']['M_ppi_c'] + costCount
            
            sim_properties['counts'][cost+'_paid'] = int(costCount)
            
            sim_properties['counts'][cost] = 0
            
            
    for recID, metID in NMP_recycle_counters.items():
        
        recycledCount = sim_properties['counts'][recID]
        
        sim_properties['counts'][metID] = sim_properties['counts'][metID] + recycledCount
        
        sim_properties['counts'][recID+'_paid'] = int(recycledCount)
        
        sim_properties['counts'][recID] = 0
        
    
    cost = 'FMET_cost'
    
    costCount = sim_properties['counts'][cost]
    
    metID = 'M_10fthfglu3_c'
    
    if sim_properties['counts'][metID] < costCount:
        
        sim_properties['counts'][cost] = costCount - sim_properties['counts'][metID] + 1
        
        sim_properties['counts'][cost+'_paid'] = sim_properties['counts'][metID] - 1
        
        sim_properties['counts'][metID] = 1
        
        sim_properties['counts']['M_thfglu3_c'] = sim_properties['counts']['M_thfglu3_c'] + sim_properties['counts'][metID] + 1
        
    else:
        
        sim_properties['counts'][metID] = sim_properties['counts'][metID] - costCount
        
        sim_properties['counts']['M_thfglu3_c'] = sim_properties['counts']['M_thfglu3_c'] + costCount

        sim_properties['counts'][cost+'_paid'] = int(costCount)

        sim_properties['counts'][cost] = 0
        
#     resetAminoAcidCostCounters(sim_properties)
    
    return None
#########################################################################################


#########################################################################################
def resetCostCounters(sim_properties):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    cost_counters = {'GTP_translat_cost':'M_gtp_c', 'ATP_trsc_cost':'M_atp_c', 'ATP_mRNAdeg_cost':'M_atp_c',
                     'ATP_DNArep_cost':'M_atp_c', 'ATP_transloc_cost':'M_atp_c',
                     'ATP_mRNA_cost':'M_atp_c', 'CTP_mRNA_cost':'M_ctp_c', 'UTP_mRNA_cost':'M_utp_c', 'GTP_mRNA_cost':'M_gtp_c',
                     'ATP_tRNA_cost':'M_atp_c', 'CTP_tRNA_cost':'M_ctp_c', 'UTP_tRNA_cost':'M_utp_c', 'GTP_tRNA_cost':'M_gtp_c',
                     'ATP_rRNA_cost':'M_atp_c', 'CTP_rRNA_cost':'M_ctp_c', 'UTP_rRNA_cost':'M_utp_c', 'GTP_rRNA_cost':'M_gtp_c',
                     'dATP_DNArep_cost':'M_datp_c', 'dTTP_DNArep_cost':'M_dttp_c', 'dCTP_DNArep_cost':'M_dctp_c', 'dGTP_DNArep_cost':'M_dgtp_c'}
    
    for cost, metID in cost_counters.items():
        
        sim_properties['counts'][cost+'_paid'] = 0
        
        sim_properties['counts'][cost+'_second'] = 0
        
    
    aaCostMap = {"A":"ALA_cost", "R":"ARG_cost", 
        "N":"ASN_cost", "D":"ASP_cost", "C":"CYS_cost", "E":"GLU_cost", "Q":"GLN_cost", "G":"GLY_cost", 
        "H":"HIS_cost", "I":"ILE_cost", "L":"LEU_cost", "K":"LYS_cost", "M":"MET_cost", "F":"PHE_cost", 
        "P":"PRO_cost", "S":"SER_cost", "T":"THR_cost", "W":"TRP_cost", "Y":"TYR_cost", "V":"VAL_cost"}
    
    for aa, aaCost in aaCostMap.items():
        
        sim_properties['counts'][aaCost+'_paid'] = 0
        
        sim_properties['counts'][aaCost+'_second'] = 0
        
    sim_properties['counts']['FMET_cost_paid'] = 0
    
    sim_properties['counts']['FMET_cost_second'] = 0
    
    
    NMP_recycle_counters = {'AMP_mRNAdeg_cost':'M_amp_c', 'UMP_mRNAdeg_cost':'M_ump_c', 'CMP_mRNAdeg_cost':'M_cmp_c', 'GMP_mRNAdeg_cost':'M_gmp_c'}

    for recID, metID in NMP_recycle_counters.items():
        
        sim_properties['counts'][recID+'_paid'] = 0
        
        sim_properties['counts'][recID+'_second'] = 0

    return None
#########################################################################################


#########################################################################################
def checkRepInitState(sim_properties):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    for i in range(20, 31):
        
        repBound = 'ori_rep2_DnaA_'+str(i)
        
        if sim_properties['counts'][repBound] > 0:
            
            repInitState = True
            
            return repInitState
            
        elif sim_properties['counts'][repBound] == 0:
            
            repInitState = False
    
    return repInitState
#########################################################################################


#########################################################################################
def updateSA(sim_properties):
    """
    Inputs:
    
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    genome = sim_properties['genome']
    
    avgProtSA = 28 #24.75 # nm^2, average protein surface area to produce expected 47% coverage for 9.6K membrane proteins
    
#     lipidSizes = {
#     'M_clpn_c':0.4,
#     'M_chsterol_c':0.35, # 0.35, test for Chol. value smaller
#     'M_sm_c':0.45,
#     'M_pc_c':0.55,
#     'M_pg_c':0.6,
#     'M_galfur12dgr_c':0.6,
#     'M_fa_c':0.5, # Should this be here??
#     'M_12dgr_c':0.5, # Scale down, should cdp-dag be added???
#     'M_pa_c':0.5,
#     'M_cdpdag_c':0.5,
#     }
    
    lipidSizes = {
    'M_clpn_c':0.4,
    'M_chsterol_c':0.35, 
    'M_sm_c':0.45,
    'M_pc_c':0.55,
    'M_pg_c':0.6,
    'M_galfur12dgr_c':0.6,
    'M_12dgr_c':0.5, 
    'M_pa_c':0.5,
    'M_cdpdag_c':0.5,
    }
    
    lipidSA = 0
    
    for lipid, size in lipidSizes.items():
        
        lipidSA = lipidSA + sim_properties['counts'][lipid]*size
        
    lipidSA = int(lipidSA*0.513)
    
    sim_properties['counts']['SA_lipid'] = lipidSA
    
    memProtList = ['JCVISYN3A_0005','JCVISYN3A_0008', 'JCVISYN3A_0009', 'JCVISYN3A_0010', 'JCVISYN3A_0011', 'JCVISYN3A_0030', 'JCVISYN3A_0034', 'JCVISYN3A_0060','JCVISYN3A_0095', 'JCVISYN3A_0113','JCVISYN3A_0114','JCVISYN3A_0116','JCVISYN3A_0117','JCVISYN3A_0132', 'JCVISYN3A_0143','JCVISYN3A_0146','JCVISYN3A_0164','JCVISYN3A_0165', 'JCVISYN3A_0166', 'JCVISYN3A_0167', 'JCVISYN3A_0168', 'JCVISYN3A_0169', 'JCVISYN3A_0195', 'JCVISYN3A_0196', 'JCVISYN3A_0197','JCVISYN3A_0235','JCVISYN3A_0239','JCVISYN3A_0248','JCVISYN3A_0249','JCVISYN3A_0296','JCVISYN3A_0304','JCVISYN3A_0314','JCVISYN3A_0317','JCVISYN3A_0326','JCVISYN3A_0332','JCVISYN3A_0338','JCVISYN3A_0345', 'JCVISYN3A_0346','JCVISYN3A_0371','JCVISYN3A_0372','JCVISYN3A_0379','JCVISYN3A_0388','JCVISYN3A_0398','JCVISYN3A_0399','JCVISYN3A_0411','JCVISYN3A_0425', 'JCVISYN3A_0426', 'JCVISYN3A_0427', 'JCVISYN3A_0428','JCVISYN3A_0439','JCVISYN3A_0440','JCVISYN3A_0478','JCVISYN3A_0481','JCVISYN3A_0505','JCVISYN3A_0516','JCVISYN3A_0601','JCVISYN3A_0639', 'JCVISYN3A_0641', 'JCVISYN3A_0642', 'JCVISYN3A_0643', 'JCVISYN3A_0652', 'JCVISYN3A_0685', 'JCVISYN3A_0686', 'JCVISYN3A_0691','JCVISYN3A_0696', 'JCVISYN3A_0706', 'JCVISYN3A_0707', 'JCVISYN3A_0708', 'JCVISYN3A_0774', 'JCVISYN3A_0777','JCVISYN3A_0778','JCVISYN3A_0779', 'JCVISYN3A_0787', 'JCVISYN3A_0789', 'JCVISYN3A_0790', 'JCVISYN3A_0791', 'JCVISYN3A_0792','JCVISYN3A_0795', 'JCVISYN3A_0797', 'JCVISYN3A_0822', 'JCVISYN3A_0827', 'JCVISYN3A_0830', 'JCVISYN3A_0835','JCVISYN3A_0836', 'JCVISYN3A_0839', 'JCVISYN3A_0852','JCVISYN3A_0870', 'JCVISYN3A_0872', 'JCVISYN3A_0876', 'JCVISYN3A_0878', 'JCVISYN3A_0879', 'JCVISYN3A_0881','JCVISYN3A_0908']
    
    memPtnCnt = 0
    
    for locusTag in memProtList:
        
        locusNum = locusTag.split('_')[1]
        
        ptnID = 'P_' + locusNum
        
        memPtnCnt = memPtnCnt + sim_properties['counts'][ptnID]
        
        try:
            
            if sim_properties['secy_states'][ptnID] < sim_properties['counts'][ptnID]:
                
                aasequence = genome[locusTag]['AAsequence']
                
                sim_properties['counts']['ATP_transloc_cost'] = sim_properties['counts']['ATP_transloc_cost'] + int((len(aasequence)/10)*(sim_properties['counts'][ptnID]-sim_properties['secy_states'][ptnID]))
                
                sim_properties['counts']['ATP_transloc_cost_second'] = sim_properties['counts']['ATP_transloc_cost_second'] + int((len(aasequence)/10)*(sim_properties['counts'][ptnID]-sim_properties['secy_states'][ptnID]))
            
        except:
            
#             continue
            
            sim_properties['secy_states'][ptnID] = sim_properties['counts'][ptnID]
            
        sim_properties['secy_states'][ptnID] = int(sim_properties['counts'][ptnID])
        
#     memPtnCnt = memPtnCnt + sim_properties['counts']['Degradosome']
        
    ptnSA = int(memPtnCnt*avgProtSA)
    
    sim_properties['counts']['SA_ptn'] = ptnSA
    
    sim_properties['counts']['SA_total'] = int(lipidSA + ptnSA)
    
    sim_properties['SA'] = int(lipidSA + ptnSA)/1e18
    
    double_V_radius = round_sig(((3*(((4/3)*np.pi*(2.0e-7)**3)*2))/(4*np.pi))**(1/3), sig=4)
    
#     double_V_radius = 2.502e-7
    
#     print('Double radius = ', double_radius)
    
    cyto_radius_nm = min(np.sqrt(sim_properties['SA']/(4*np.pi)), double_V_radius)
    
    cyto_radius_nm_equivalent_sphere = np.sqrt(sim_properties['SA']/(4*np.pi))
    
    ##### Start division if cell has grown to more than double its initial volume #####
    if cyto_radius_nm_equivalent_sphere > double_V_radius:
        
        sim_properties['division_started'] = True
        
        if ('nextDivSA' not in sim_properties.keys()):
            
            sim_properties['nextDivSA'] = sim_properties['SA'] - 1e-16
    
#     cyto_radius = int(round(cyto_radius_nm/sim_properties['lattice_spacing']))
    cyto_radius = int(cyto_radius_nm/sim_properties['lattice_spacing'])
    
#     sim_properties['volume'] = (4/3)*np.pi*(sim_properties['cyto_radius']*sim_properties['lattice_spacing'])**3
    sim_properties['volume'] = (4/3)*np.pi*(cyto_radius_nm)**3
    sim_properties['volume_L'] = sim_properties['volume']*1000
    
    sim_properties['counts']['Volume'] = int(sim_properties['volume_L']*1e21)

    if sim_properties['division_started']:

#         gamma_V = round_sig(sim_properties['volume'] / ((4 / 3) * np.pi * (cyto_radius_nm_equivalent_sphere) ** 3),
#                             sig=4)

#         print('Gamma V: ', gamma_V)

        if (sim_properties['SA'] >= sim_properties['nextDivSA']):
            
            sim_properties['nextDivSA'] = sim_properties['SA'] + 5.0e-15

#             sim_properties['gamma_V'] = round_sig(sim_properties['next_gamma_V'], sig=2)

#             sim_properties['next_gamma_V'] = sim_properties['next_gamma_V'] - 0.01

            # Calculate R and H based on equations for SA and V, of overlapping spheres minus inner caps
            # The equations below are for the SA and V of the *total* cell shape (include both sides of the dividing cell)

            volume_equation = lambda cutoff, radius: ((4 / 3) * np.pi * radius ** 3
                                                          + 2 * np.pi * cutoff * radius ** 2
                                                          - (2 * np.pi / 3) * cutoff ** 3)

            surface_area_equation = lambda cutoff, radius: 4 * np.pi * radius * (cutoff + radius)

            cellV_nm = sim_properties['volume'] * 1e27
            cellSA_nm = sim_properties['SA'] * 1e18
            # Define the system of equations, scaling from meters to nanometers to enhance numerical stability for fsolve
            division_equations = lambda f: [
                volume_equation(f[0], f[1]) - cellV_nm,
                surface_area_equation(f[0], f[1]) - cellSA_nm
            ]

            # Solve for cutoff and radius using initial guesses
            if 'divH' not in sim_properties.keys():
                print('Creating Div')
                result = least_squares(division_equations, [0, sim_properties['cyto_radius_nm']*1e9],  bounds=[0,300])
                sim_properties['divH_Prev'], sim_properties['divR_Prev'] = 0, sim_properties['cyto_radius_nm']*1e9
                
        #         print(result)
            else:
                result = least_squares(division_equations, [sim_properties['divH'], sim_properties['divR']], bounds=[0,300])
                sim_properties['divH_Prev'], sim_properties['divR_Prev'] = float(sim_properties['divH']), float(sim_properties['divR'])
        #         print(result)

            sim_properties['divH'], sim_properties['divR'] = result['x'][0], result['x'][1]  # store the results; units are m
            print(sim_properties['divH'], sim_properties['divR'])

            updateRegions = True

        else:

            updateRegions = False
        
#         sim_properties['gamma_V'] = round_sig(sim_properties['next_gamma_V'], sig=2)
        
#         sim_properties['next_gamma_V'] = sim_properties['next_gamma_V'] - 0.01
        
#         updateRegions = True
        
    else:
        
        if cyto_radius > sim_properties['cyto_radius']:
        
            updateRegions = True

        else:

            updateRegions = False
            
            
    sim_properties['cyto_radius'] = cyto_radius
    sim_properties['cyto_radius_nm'] = cyto_radius_nm
    
    print('SA: ', sim_properties['SA'])
    print('V: ', sim_properties['volume'])
    
    print('cyto radius: ', sim_properties['cyto_radius'])
    print('cyto radius nm: ', sim_properties['cyto_radius_nm'])
    
    return updateRegions
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
    
    particles = max(1, int(round((conc/1000)*NA*sim_properties['volume_L'])))

    return particles
#########################################################################################


#########################################################################################
def SaCtoPart(conc, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    particles = max(1, int(round(conc*sim_properties['SA_total'])))

    return particles
#########################################################################################

# #########################################################################################
# def BADupdateLongGeneStates(sim_properties, lattice):
    
#     plattice = lattice.getParticleLatticeView()
    
#     genome = sim_properties['genome']
    
#     for locusTag, locusDict in genome.items():
        
#         locusNum = locusTag.split('_')[1]
        
#         rnasequence = locusDict["RNAsequence"]
        
#         if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
    
#             starts = locusDict['startIndex']

#             for chrom in range(len(starts)):
                
#                 chromo = chrom + 1

#                 try:
#                     start = starts[chrom]
#     #                 end = ends[i]
#                 except:
#                     continue

#                 startXYZ = sim_properties['DNAcoords'][start]
                
#                 RNAP_count = 0
        
#                 for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):

#                     on_i = 'RP_' + locusNum + '_c' + str(chromo)+ '_' + str(i)

#                     RNAP_count = RNAP_count + sim_properties['counts'][on_i]

#                     if i == sim_properties['long_rna_trsc'][locusTag]['max_rnap']:

#                         on_ip1 = 'RP_' + locusNum  + '_c' + str(chromo) + '_done'

#                         RNAP_count = RNAP_count + sim_properties['counts'][on_ip1]

#                 RNAP_count = int(RNAP_count)

#                 print('Chomosome ' + str(chromo) + ' RNAP count: ', RNAP_count)
                
#                 updateStartState = 'RP_' + locusNum + '_c' + str(chromo) + '_1'
                
#                 if RNAP_count == 0:
                    
#                     gene = 'G_' + locusNum
                    
#                     g_idx = sim_properties['name_to_index'][gene]
                    
#                     RP_gene_1 = 'RP_' + locusNum + '_' + str(int(1))
                    
#                     RP_idx_1 = sim_properties['name_to_index'][RP_gene_1]
                    
#                     geneCheck = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), g_idx)
                    
#                     RpCheck = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), RP_idx_1)
                    
#                     if geneCheck and not RpCheck:
                    
#                         print(gene)

#                         sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = gene
                        
#                         continue
                        
#                     elif RpCheck and not geneCheck:
                        
#                         print(RP_gene_1)
                        
#                         sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene_1
                        
#                         sim_properties['counts']['RP_' + locusNum + '_c' + str(chromo) + '_' + str(1)] = 1
                        
#                         continue
                        
#                     elif geneCheck and RpCheck:
                        
#                         if chromo == 1:
                        
#                             updateStartStateOther = 'RP_' + locusNum + '_c' + str(2) + '_1'
                            
#                         elif chromo == 2:
                            
#                             updateStartStateOther = 'RP_' + locusNum + '_c' + str(1) + '_1'

#                         if sim_properties['counts'][updateStartStateOther] == 0:
                            
#                             print(RP_gene_1)
                        
#                             sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene_1

#                             sim_properties['counts']['RP_' + locusNum + '_c' + str(chromo) + '_' + str(1)] = 1
                            
#                             continue
                            
#                         else:
                            
#                             print(gene)

#                             sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = gene

#                             continue
                        
#                     else:
                        
#                         print('Lost start site')

#                         continue
                    
#                 else:
                    
#                     RP_gene = 'RP_' + locusNum + '_' + str(int(RNAP_count))
                    
#                     RP_idx = sim_properties['name_to_index'][RP_gene]
                    
#                     RP_gene_t = 'RP_' + locusNum + '_t_' + str(int(RNAP_count))
                    
#                     RP_t_idx = sim_properties['name_to_index'][RP_gene_t]

#                     RpCheck = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), RP_idx)
                    
#                     RptCheck = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), RP_t_idx)
                    
#                     Rpp1Check = False
                    
#                     if RNAP_count < sim_properties['long_rna_trsc'][locusTag]['max_rnap']:

#                         RP_gene_p1 = 'RP_' + locusNum + '_' + str(int(RNAP_count+1))
                        
#                         RP_p1_idx = sim_properties['name_to_index'][RP_gene_p1]

#                         Rpp1Check = checkParticle(plattice, int(startXYZ[2]), int(startXYZ[1]), int(startXYZ[0]), RP_p1_idx)
                    
#                     if RptCheck and not (RpCheck or Rpp1Check):
                        
#                         print(RP_gene_t)
                        
#                         sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene_t
                        
#                         continue
                        
#                     elif RpCheck and not (RptCheck or Rpp1Check):
                        
#                         print(RP_gene)
                        
#                         sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene
                        
#                         sim_properties['counts']['RP_' + locusNum + '_c' + str(chromo) + '_' + str(1)] = 1
                        
#                         continue
                        
#                     elif Rpp1Check and not (RpCheck or RptCheck):
                        
#                         print(RP_gene_p1)
                        
#                         sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene_p1
                        
#                         continue
                        
#                     elif (RpCheck and RptCheck):
                        
#                         if chromo == 1:
                        
#                             updateStartStateOther = 'RP_' + locusNum + '_c' + str(2) + '_1'
                            
#                         elif chromo == 2:
                            
#                             updateStartStateOther = 'RP_' + locusNum + '_c' + str(1) + '_1'
                            
#                         if (sim_properties['counts'][updateStartState] == 0) and (sim_properties['counts'][updateStartStateOther] > 0):
                            
#                             print(RP_gene_t)
                        
#                             sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene_t
                            
#                             continue
                            
#                         else:
                            
#                             print(RP_gene)
                        
#                             sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene

#                             sim_properties['counts']['RP_' + locusNum + '_c' + str(chromo) + '_' + str(1)] = 1
                            
#                             continue 
                        
#                     elif (RpCheck and Rpp1Check):
                            
#                         if (sim_properties['counts'][updateStartState] == 0):
                            
#                             print(RP_gene_p1)
                        
#                             sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene_p1

#                             sim_properties['counts']['RP_' + locusNum + '_c' + str(chromo) + '_' + str(1)] = 1
                            
#                             continue
                            
#                         elif (sim_properties['counts'][updateStartState] > 0):
                            
#                             print(RP_gene)
                            
#                             sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene

#                             sim_properties['counts']['RP_' + locusNum + '_c' + str(chromo) + '_' + str(1)] = 1
                            
#                             continue
                        
#                     elif (RptCheck and Rpp1Check):
                        
#                         if chromo == 1:
                        
#                             updateStartStateOther = 'RP_' + locusNum + '_c' + str(2) + '_1'
                            
#                         elif chromo == 2:
                            
#                             updateStartStateOther = 'RP_' + locusNum + '_c' + str(1) + '_1'

#                         if sim_properties['counts'][updateStartStateOther] == 0:
                            
#                             print(RP_gene_p1)
                        
#                             sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene_p1

#                             sim_properties['counts']['RP_' + locusNum + '_c' + str(chromo) + '_' + str(1)] = 1
                            
#                             continue
                            
#                         else:
                            
#                             print(RP_gene_t)

#                             sim_properties['long_rna_trsc']['genome'+str(int(chromo))][locusTag]['state'] = RP_gene_t

#                             continue
                        
#                     else:
                        
#                         print('Lost start site')
                        
#     return None
# #########################################################################################



# #########################################################################################
# def OLDresetAminoAcidCostCounters(sim_properties):
    
#     aaCostMap = {"A":"ALA_cost", "R":"ARG_cost", 
#         "N":"ASN_cost", "D":"ASP_cost", "C":"CYS_cost", "E":"GLU_cost", "Q":"GLN_cost", "G":"GLY_cost", 
#         "H":"HIS_cost", "I":"ILE_cost", "L":"LEU_cost", "K":"LYS_cost", "M":"MET_cost", "F":"PHE_cost", 
#         "P":"PRO_cost", "S":"SER_cost", "T":"THR_cost", "W":"TRP_cost", "Y":"TYR_cost", "V":"VAL_cost"}
    
#     for aa, aaCost in aaCostMap.items():
        
#         sim_properties['counts'][aaCost+'_paid'] = 0
    
#     return None
# #########################################################################################


# #########################################################################################
# def OLDcalculateCosts(sim_properties, lattice):
    
# #     calculateTranscriptionCosts(sim_properties, lattice)
    
#     calculateTranslationCosts(sim_properties, lattice)
    
#     calculateDegradationCosts(sim_properties, lattice)
    
#     return None
# #########################################################################################


# #########################################################################################
# def OLDcalculateDegradationCosts(sim_properties):
    
#     genome = sim_properties['genome']
    
#     for locusTag, locusDict in genome.items():
        
#         if locusDict['Type'] == 'protein':
            
#             locusNum = locusTag.split('_')[1]

#             DegID = 'D_' + locusNum
            
#             if sim_properties['counts'][DegID] < sim_properties['deg_states'][DegID]:
                
#                 mRNAdegraded = sim_properties['deg_states'][DegID] - sim_properties['counts'][DegID]
                
#                 rnasequence = locusDict['RNAsequence']
            
#                 sim_properties['counts']['ATP_mRNAdeg_cost'] = sim_properties['counts']['ATP_mRNAdeg_cost'] + mRNAdegraded*len(rnasequence)

#                 for nucleotide in set(rnasequence):

#                     nucleotideCount = rnasequence.count(nucleotide)

#                     sim_properties['counts'][nucleotide + 'MP_mRNAdeg'] = sim_properties['counts'][nucleotide + 'MP_mRNAdeg'] + mRNAdegraded*nucleotideCount
                    
#             sim_properties['deg_states'][DegID] = sim_properties['counts'][DegID]

#     return None
# #########################################################################################


#########################################################################################
# def OLDcalculateTranscriptionCosts(sim_properties, lattice):
    
#     plattice = lattice.getParticleLatticeView()
    
#     genome = sim_properties['genome']
    
#     for locusTag, locusDict in genome.items():
        
#         locusNum = locusTag.split('_')[1]
        
#         finishedRNAP = 'RP_' + locusNum + '_f'
        
#         if sim_properties['counts'][finishedRNAP] > 0:
            
#             rnasequence = locusDict['RNAsequence']
            
#             sim_properties['counts']['ATP_trsc_cost'] = sim_properties['counts']['ATP_trsc_cost'] + sim_properties['counts'][finishedRNAP]*len(rnasequence)
            
#             if locusDict['Type'] == 'protein':
            
#                 RNAtype = 'mRNA'
                
#             else:
                
#                 RNAtype = locusDict['Type']

#             for nucleotide in set(rnasequence):
        
#                 nucleotideCount = rnasequence.count(nucleotide)
            
#                 sim_properties['counts'][nucleotide + 'TP_' + RNAtype + '_cost'] = sim_properties['counts'][nucleotide + 'TP_' + RNAtype + '_cost'] + sim_properties['counts'][finishedRNAP]*nucleotideCount
            
#     return None
#########################################################################################
