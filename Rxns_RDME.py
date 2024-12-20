"""
Authors: Zane Thornburg

Functions for creating spatially-localized reactions that take place in the RDME simulation
"""

import numpy as np
import pandas as pd
import json
from collections import defaultdict, OrderedDict

import GIP_rates as GIP

import Diffusion as Diff


#########################################################################################
def general_reaction_rates(sim):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    avgdr   = 6.022e23 # molec/mol
    Ecoli_V = 1e-15 #L
    ribo_init = 40*Ecoli_V*avgdr/60/6800
    
    sim.rateConst('RNAP_off', 1e-3, 1)
    sim.rateConst('ribo_bind', ribo_init, 2)
    sim.rateConst('Ribo_off', 1e-3, 1)

    deg_bind_rate = 11*avgdr*Ecoli_V/60/7800 #2400 or 7800 #1/M/s
    sim.rateConst('deg_bind_rate', deg_bind_rate, 2)
    
    secY_init = 1.0e6 #5*Ecoli_V*avgdr/60/6800
    sim.rateConst('secY_on', secY_init, 2)
    sim.rateConst('secY_off', 1e-3, 1)
    
    binding_rate = 4*(180/765)*Ecoli_V*avgdr/11400/60 #/1800/60
    sim.rateConst('RNAP_on', binding_rate, 2)  
    
    sim.rateConst('conversion', 1000000, 1)
    
    decayParticle = sim.species('decay')
    
    Diff.proteinDiffusionCyto(sim, 'decay')
    
    sim.region('cytoplasm').addReaction([decayParticle], [], sim.rc.conversion)
    sim.region('outer_cytoplasm').addReaction([decayParticle], [], sim.rc.conversion)
    sim.region('DNA').addReaction([decayParticle], [], sim.rc.conversion)
    
    return None
#########################################################################################


#########################################################################################
def translation(sim, locusNum, aasequence, membrane=False):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    mRNA = sim.species('R_'+locusNum)
    
    mRNA_read = sim.species('R_'+locusNum+'_d')
    
    sim.region('cytoplasm').addReaction([mRNA_read], [mRNA], sim.rc.conversion)
    sim.region('outer_cytoplasm').addReaction([mRNA_read], [mRNA], sim.rc.conversion)
    sim.region('DNA').addReaction([mRNA_read], [mRNA], sim.rc.conversion)
    
    if membrane:
        ptn = sim.species('C_P_'+locusNum)
    else:
        ptn = sim.species('P_'+locusNum)
        
    ptnCost = sim.species('P_'+locusNum+'_TC')
    
    mRNA_ribo = sim.species('RB_'+locusNum)
    
    Diff.ribosomeDiffusion(sim, 'RB_'+locusNum)
    
    translation_rate = sim.rateConst(locusNum + '_translat', GIP.TranslationRate(aasequence), 1)
    
    ribo_part = sim.species('ribosomeP')
    
    sim.region('ribo_centers').addReaction([mRNA, ribo_part], [mRNA_ribo], sim.rc.ribo_bind)
    sim.region('ribosomes').addReaction([mRNA, ribo_part], [mRNA_ribo], sim.rc.ribo_bind)
    
    sim.region('ribo_centers').addReaction([mRNA_ribo], [mRNA_read, ribo_part, ptnCost], translation_rate)
    sim.region('ribosomes').addReaction([mRNA_ribo], [mRNA_read, ribo_part, ptnCost], translation_rate)
    
    return None
#########################################################################################


#########################################################################################
def transcription(sim, sim_properties, locusTag, rnasequence, restart=False):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    locusNum = locusTag.split('_')[1]
    
    mRNA = sim.species('R_' + locusNum)
    
    gene = sim.species('G_' + locusNum + '_C1')
    
    RNAP_gene = sim.species('RP_' + locusNum + '_C1')
    
    RNAP_finished = 'RP_' + locusNum + '_f' + '_C1'
    
    RNAP_made = 'RPM_' + locusNum
    
    if not restart:
    
        sim_properties['counts'][RNAP_finished] = 0

        sim_properties['counts'][RNAP_made] = 0
    
    RNAP_on = sim.rateConst(locusNum + 'RNAP_on', GIP.RNAP_binding(sim, sim_properties, locusTag), 2)
    
    RNApol = sim.species('RNAP')
        
    sim.region('DNA').addReaction([gene, RNApol], [RNAP_gene], RNAP_on)
    
    gene2 = sim.species('G_' + locusNum + '_C2')
    
    RNAP_gene2 = sim.species('RP_' + locusNum + '_C2')
    
    RNAP_finished = 'RP_' + locusNum + '_f' + '_C2'
    
    if not restart:
    
        sim_properties['counts'][RNAP_finished] = 0
    
    sim.region('DNA').addReaction([gene2, RNApol], [RNAP_gene2], RNAP_on)
    
    return None
#########################################################################################


#########################################################################################
def degradation_mrna(sim, sim_properties, locusNum, rnasequence, restart=False):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    mRNA = sim.species('R_' + locusNum)
    
    Deg_RNA = sim.species('D_' + locusNum)
    
    Deg_track = sim.species('DT_' + locusNum)
    
    if not restart:
    
        sim_properties['counts']['DM_' + locusNum] = 0
    
    Diff.degradosomeDiffusion(sim, 'D_' + locusNum)
    
    Diff.degradosomeDiffusion(sim, 'DT_' + locusNum)
    
    rna_deg_rate = sim.rateConst(locusNum + '_RNAdeg', GIP.mrnaDegradationRate(rnasequence), 1)
    
    degradosome = sim.species('Degradosome')
    
    sim.region('outer_cytoplasm').addReaction([mRNA, degradosome], [Deg_track, Deg_RNA], sim.rc.deg_bind_rate)

    sim.region('outer_cytoplasm').addReaction([Deg_RNA], [degradosome], rna_deg_rate)
    
    return None
#########################################################################################


#########################################################################################
def replicationInitiation(sim, sim_properties, restart=False):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    oriC = sim.species('oriC')
    DnaA = sim.species('P_0001')
    
    oriC_part_Idx = sim_properties['chromosome_features']['oriC']['index']
    
    dsDNA_HA_on = 7800*1000 #M-1 s-1
    dsDNA_LA_on = 35*1000 #M-1 s-1
    
    sim.rateConst('dnaa_ds_ha', dsDNA_HA_on, 2)
    sim.rateConst('dnaa_ds_la', dsDNA_LA_on, 2)
    
    ori_HA = sim.species('ori_HA_DnaA')
    
    ori_LA1 = sim.species('ori_LA1_DnaA')
    
    ori_LA2 = sim.species('ori_LA2_DnaA')
    
    if not restart:
        
        sim_properties['chromosome_features']['ori_HA_DnaA'] = {}
        sim_properties['chromosome_features']['ori_HA_DnaA']['index'] = oriC_part_Idx

        sim_properties['chromosome_features']['ori_LA1_DnaA'] = {}
        sim_properties['chromosome_features']['ori_LA1_DnaA']['index'] = oriC_part_Idx

        sim_properties['chromosome_features']['ori_LA2_DnaA'] = {}
        sim_properties['chromosome_features']['ori_LA2_DnaA']['index'] = oriC_part_Idx
    
    
    sim.region('DNA').addReaction([oriC, DnaA], [ori_HA], sim.rc.dnaa_ds_ha)
    sim.region('DNA').addReaction([ori_HA, DnaA], [ori_LA1], sim.rc.dnaa_ds_la)
    sim.region('DNA').addReaction([ori_LA1, DnaA], [ori_LA2], sim.rc.dnaa_ds_la)
    
    ssDNA_on_rate = 100*1000 #M-1 s-1
    ssDNA_off_rate = 0.55 #s-1

    #ssDNA_on_rate = 140*1000 #M-1 s-1
    #ssDNA_off_rate = 0.42 #s-1
    
    sim.rateConst('dnaa_ss_on', ssDNA_on_rate, 2)
    sim.rateConst('dnaa_ss_off', ssDNA_off_rate, 1)  
    
    replisomeOnRate = 1e6 #M-1 s-1
    sim.rateConst('repOn', replisomeOnRate, 2)
    
    for i in range(1,31):
        
        ssDNAbound = sim.species('ori_ss_DnaA_'+str(i))
        
        if not restart:
        
            sim_properties['chromosome_features']['ori_ss_DnaA_'+str(i)] = {}
            sim_properties['chromosome_features']['ori_ss_DnaA_'+str(i)]['index'] = oriC_part_Idx
        
        if i == 1:
            ssDNAunbound = ori_LA2
        else:
            ssDNAunbound = sim.species('ori_ss_DnaA_'+str(i-1))
            
        sim.region('DNA').addReaction([ssDNAunbound, DnaA], [ssDNAbound], sim.rc.dnaa_ss_on)
        sim.region('DNA').addReaction([ssDNAbound], [ssDNAunbound, DnaA], sim.rc.dnaa_ss_off)
        
    replisome_particle = sim.species('replisome')
    
    #sim.distributeNumber(replisome, sim.region('cytoplasm'), int(16))
    
    #Diff.proteinDiffusionCyto(sim, 'replisome')

    replisome = sim.species('P_0128')
        
    for i in range(20, 31):
        
        ssDNAbound = sim.species('ori_ss_DnaA_'+str(i))
        
        repBound = sim.species('ori_rep_DnaA_'+str(i))
        
        if not restart:
        
            sim_properties['chromosome_features']['ori_rep_DnaA_'+str(i)] = {}
            sim_properties['chromosome_features']['ori_rep_DnaA_'+str(i)]['index'] = oriC_part_Idx
        
        repBound2 = sim.species('ori_rep2_DnaA_'+str(i))
        
        if not restart:
        
            sim_properties['chromosome_features']['ori_rep2_DnaA_'+str(i)] = {}
            sim_properties['chromosome_features']['ori_rep2_DnaA_'+str(i)]['index'] = oriC_part_Idx
        
        sim.region('DNA').addReaction([ssDNAbound, replisome], [repBound], sim.rc.repOn)
        
        sim.region('DNA').addReaction([repBound, replisome], [repBound2], sim.rc.repOn)

    return None
#########################################################################################


#########################################################################################
def translocation_secy(sim, locusNum, aasequence):
    """
    Inputs:
    
    Returns:
    Called by:
    Description:
    """
    
    secy = sim.species('P_0652')
    
    cPTN = sim.species('C_P_'+locusNum)
    PTN = sim.species('P_'+locusNum)
    
    secy_bound = sim.species('S_'+locusNum)

    transloc_rate = sim.rateConst(locusNum + '_insertion', GIP.TranslocationRate(aasequence), 1)

    sim.region('outer_cytoplasm').addReaction([secy,cPTN],[secy_bound],sim.rc.secY_on)
    sim.region('outer_cytoplasm').addReaction([secy_bound],[secy,PTN],transloc_rate)
    
    return None
#########################################################################################


#########################################################################################
def transcriptionLong(sim, sim_properties, locusTag, rnasequence, restart=False):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    locusNum = locusTag.split('_')[1]
    
    mRNA = sim.species('R_' + locusNum)
    
    gene = sim.species('G_' + locusNum+ '_C1')
    
    gene2 = sim.species('G_' + locusNum + '_C2')
    
    RNApol = sim.species('RNAP')
    
    RNAP_made = 'RPM_' + locusNum
    
    if not restart:
    
        sim_properties['counts'][RNAP_made] = 0
    
    RNAP_on = sim.rateConst(locusNum + 'RNAP_on', GIP.RNAP_binding(sim, sim_properties, locusTag), 2)
    
    for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):
    
        RNAP_gene = sim.species('RP_' + locusNum + '_' + str(i) + '_C1')

        RNAP_gene_trsc = sim.species('RP_' + locusNum + '_t_' + str(i) + '_C1')
        
        RNAP_gene2 = sim.species('RP_' + locusNum + '_' + str(i) + '_C2')

        RNAP_gene_trsc2 = sim.species('RP_' + locusNum + '_t_' + str(i) + '_C2')
        
        if i==1:

            sim.region('DNA').addReaction([gene, RNApol], [RNAP_gene], RNAP_on)
            
            sim.region('DNA').addReaction([gene2, RNApol], [RNAP_gene2], RNAP_on)
            
        else:
            
            RNAP_gene_trsc_m1 = sim.species('RP_' + locusNum + '_t_' + str(i-1) + '_C1')
            
            RNAP_gene_trsc_m1_2 = sim.species('RP_' + locusNum + '_t_' + str(i-1) + '_C2')
            
            sim.region('DNA').addReaction([RNAP_gene_trsc_m1, RNApol], [RNAP_gene], RNAP_on)
            
            sim.region('DNA').addReaction([RNAP_gene_trsc_m1_2, RNApol], [RNAP_gene2], RNAP_on)
    
    if not restart:
    
        for chromo in range(1,3):

            for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):

                off_i = 'RP_' + locusNum + '_c' + str(chromo) + '_open_' + str(i)

                on_i = 'RP_' + locusNum + '_c' + str(chromo)+ '_' + str(i)

                sim_properties['counts'][on_i] = 0
                sim_properties['counts'][off_i] = 1

                if i == sim_properties['long_rna_trsc'][locusTag]['max_rnap']:

                    off_ip1 = 'RP_' + locusNum  + '_c' + str(chromo) + '_endTrsc'

                    on_ip1 = 'RP_' + locusNum  + '_c' + str(chromo) + '_done'

                    sim_properties['counts'][on_ip1] = 0
                    sim_properties['counts'][off_ip1] = 1
    
    return None
#########################################################################################


#########################################################################################
def addRibosomeBiogenesis(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    smallAssemblyFile = sim_properties['head_directory'] + "input_data/oneParamMulder-local_min.json"
    
    assemblyData = json.load(open(smallAssemblyFile))
    
    maxImts = 19 #145 # Don't count 30S or 16S
    mses = assemblyData['netmin_rmse']
    spsRemoved = assemblyData['netmin_species']
    spNames = set([ sp['name'] for sp in assemblyData['species'] if sp['name'][0] == 'R' ])
    for err,sps in zip(mses,spsRemoved):

        if len(spNames) <= maxImts:
            break
        spNames.difference_update(set(sps))
        
    assemblyRxns = [r for r in assemblyData['reactions'] if r['intermediate'] in spNames and r['product'] in spNames ]
    
    bindingRates = { p['rate_id']: p for p in assemblyData['parameters'] if p['rate_id'] in set(r['rate_id'] for r in assemblyRxns) }
    
    rnasequence = sim_properties['genome']['JCVISYN3A_0069']['RNAsequence']

    for rxn in assemblyRxns:
        
        riboProteinSubstrate = 'S' + rxn['protein'].split('s')[1]
        
        rxnName = 'Ribo' + riboProteinSubstrate + 'Binding'
        
        bindRate = sim.rateConst(rxnName, float(bindingRates[ rxn['rate_id'] ]['rate'])*1e6, 2)
        
        riboPtnLocusTag = sim_properties['RiboPtnMap'][riboProteinSubstrate]['locusTag']
            
        riboPtnLocusNum = riboPtnLocusTag.split('_')[1]

        riboPtnSubstrate = sim.species('P_'+riboPtnLocusNum)
        
        assemblyProduct = sim.species(rxn['product'])
        
        Diff.rnaDiffusion(sim, rxn['product'], rnasequence)
        
        if rxn['intermediate'] != 'R':
        
            assemblySubstrateNames = [rxn['intermediate']]
            
        elif rxn['intermediate'] == 'R':
            
            assemblySubstrateNames = ['R_0069', 'R_0534']
        
        for substrateName in assemblySubstrateNames:
            
            assemblySubstrate = sim.species(substrateName)
            
            sim.region('cytoplasm').addReaction([riboPtnSubstrate, assemblySubstrate], [assemblyProduct], bindRate)
            sim.region('outer_cytoplasm').addReaction([riboPtnSubstrate, assemblySubstrate], [assemblyProduct], bindRate)
            sim.region('DNA').addReaction([riboPtnSubstrate, assemblySubstrate], [assemblyProduct], bindRate)
            
            
    LargeAssemblyFile = sim_properties['head_directory'] + 'input_data/LargeSubunit.xlsx'
    
    rateConstants = pd.read_excel(LargeAssemblyFile, sheet_name='parameters', index_col=None)
    
    assemblyRxns = pd.read_excel(LargeAssemblyFile, sheet_name='reactions', index_col=None)
    
    rnasequence = sim_properties['genome']['JCVISYN3A_0068']['RNAsequence']
    
    for index, rxn in assemblyRxns.iterrows():
        
        riboPtnID = rxn['substrate']
        
        rate = rateConstants.loc[ rateConstants['Protein'] == riboPtnID ]['Rate'].values[0]
        
        rxnName = riboPtnID + 'binding'
        
        bindRate = sim.rateConst(rxnName, float(rate)*1e6, 2)
        
        if riboPtnID == 'L7':
        
            riboPtnID = 'L7/L12'
            
        if riboPtnID == '5S':
            
            assemblySubstrateNames = ['R_0067', 'R_0532']
            
        else:
            
            riboPtnLocusTag = sim_properties['RiboPtnMap'][riboPtnID]['locusTag']

            riboPtnLocusNum = riboPtnLocusTag.split('_')[1]
            
            assemblySubstrateNames = ['P_'+riboPtnLocusNum]
        
        if rxn['intermediate'] != 'R':
        
            assemblyIntermediateNames = [rxn['intermediate']]
            
        elif rxn['intermediate'] == 'R':
            
            assemblyIntermediateNames = ['R_0068', 'R_0533']
            
        assemblyProduct = sim.species(rxn['product'])
        
        Diff.rnaDiffusion(sim, rxn['product'], rnasequence)
            
        for substrateName in assemblySubstrateNames:
            
            for intermediateName in assemblyIntermediateNames:
            
                assemblySubstrate = sim.species(substrateName)
                
                assemblyIntermediate = sim.species(intermediateName)

                sim.region('cytoplasm').addReaction([assemblySubstrate, assemblyIntermediate], [assemblyProduct], bindRate)
                sim.region('outer_cytoplasm').addReaction([assemblySubstrate, assemblyIntermediate], [assemblyProduct], bindRate)
                sim.region('DNA').addReaction([assemblySubstrate, assemblyIntermediate], [assemblyProduct], bindRate)
    
    LSU_SSU_binding = sim.rateConst('LSU_SSU_bind', 1e7, 2)
    
    LSU = sim.species('R5SL1L2L3L4L5L6L7L9L10L11L13L14L15L16L17L18L19L20L21L22L23L24L27L28L29L31L32L33L34L35L36')
    SSU = sim.species('Rs3s4s5s6s7s8s9s10s11s12s13s14s15s16s17s19s20')
    ribosomeP = sim.species('ribosomeP')
    
    sim.region('cytoplasm').addReaction([SSU, LSU], [ribosomeP], LSU_SSU_binding)
    sim.region('outer_cytoplasm').addReaction([SSU, LSU], [ribosomeP], LSU_SSU_binding)
    sim.region('DNA').addReaction([SSU, LSU], [ribosomeP], LSU_SSU_binding)

    return None
#########################################################################################


#########################################################################################
def addRNAPassembly(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    alpha = sim.species('P_0645')
    beta1 = sim.species('P_0804')
    beta2 = sim.species('P_0803')
    
    a_b1 = sim.species('RNAP_ab1')
    a_b2 = sim.species('RNAP_ab2')
    b1_b2 = sim.species('RNAP_b1b2')
    
    RNAP = sim.species('RNAP')
    
    Diff.proteinDiffusionCyto(sim, 'RNAP_ab1')
    Diff.proteinDiffusionCyto(sim, 'RNAP_ab2')
    Diff.proteinDiffusionCyto(sim, 'RNAP_b1b2')
    
    ptnAssoc = sim.rateConst('ptnAssoc', 1e7, 2)
    
    sim.region('cytoplasm').addReaction([alpha, beta1], [a_b1], ptnAssoc)
    sim.region('outer_cytoplasm').addReaction([alpha, beta1], [a_b1], ptnAssoc)
    sim.region('DNA').addReaction([alpha, beta1], [a_b1], ptnAssoc)
    
    sim.region('cytoplasm').addReaction([alpha, beta2], [a_b2], ptnAssoc)
    sim.region('outer_cytoplasm').addReaction([alpha, beta2], [a_b2], ptnAssoc)
    sim.region('DNA').addReaction([alpha, beta2], [a_b2], ptnAssoc)
    
    sim.region('cytoplasm').addReaction([beta1, beta2], [b1_b2], ptnAssoc)
    sim.region('outer_cytoplasm').addReaction([beta1, beta2], [b1_b2], ptnAssoc)
    sim.region('DNA').addReaction([beta1, beta2], [b1_b2], ptnAssoc)
    
    sim.region('cytoplasm').addReaction([alpha, b1_b2], [RNAP], ptnAssoc)
    sim.region('outer_cytoplasm').addReaction([alpha, b1_b2], [RNAP], ptnAssoc)
    sim.region('DNA').addReaction([alpha, b1_b2], [RNAP], ptnAssoc)
    
    sim.region('cytoplasm').addReaction([beta2, a_b1], [RNAP], ptnAssoc)
    sim.region('outer_cytoplasm').addReaction([beta2, a_b1], [RNAP], ptnAssoc)
    sim.region('DNA').addReaction([beta2, a_b1], [RNAP], ptnAssoc)
    
    sim.region('cytoplasm').addReaction([beta1, a_b2], [RNAP], ptnAssoc)
    sim.region('outer_cytoplasm').addReaction([beta1, a_b2], [RNAP], ptnAssoc)
    sim.region('DNA').addReaction([beta1, a_b2], [RNAP], ptnAssoc)

    return None
#########################################################################################

    
