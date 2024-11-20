# Author: Zane Thornburg

###### Genetic Information Processes Reactions RDME  #######

import numpy as np

from jLM.RegionBuilder import RegionBuilder
from jLM.RDME import Sim as RDMESim
from jLM.RDME import File as RDMEFile
import jLM

import random

import os

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict


import Diffusion as Diff

import ImportInitialConditions as IC

import Rxns_RDME

import pickle


#########################################################################################
def initSimRestart(totalTime, sim_properties_file, workingDirectoryName, headDirectory):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    sim_properties = {}
    
    with open(sim_properties_file, 'rb') as handle:
        sim_properties = pickle.load(handle)
        
#     headDirectory =  os.getcwd() + '/'
    
    sim_properties['head_directory'] = headDirectory
    
    workingDirectory = headDirectory + 'Data/' + workingDirectoryName + '/'

    sim_properties['working_directory'] = workingDirectory
        
    simFolder = sim_properties['working_directory']
    
    filename = simFolder + 'MinCell_restart.lm'

    lattice_spacing = sim_properties['lattice_spacing']
    
    if lattice_spacing == 8e-9:

        N_edges_x = 96 # Number of subvolumes making up and edge of the simulation space N x N x N
        N_edges_y = 96
        N_edges_z = 128
    
    elif lattice_spacing == 10e-9:
        
        N_edges_x = 64 # Number of subvolumes making up and edge of the simulation space N x N x N
        N_edges_y = 64
        N_edges_z = 128
    
    N_edges = [N_edges_x,N_edges_y,N_edges_z]

    sim = RDMESim("JCVI-syn3A",
                  filename,
                  [N_edges_x,N_edges_y,N_edges_z],
                  lattice_spacing,
                  "extracellular",
                  latticeType='Int')
    
    sim_timestep = 50e-6

    sim.timestep = sim_properties['timestep']
    sim.simulationTime=totalTime
    sim.latticeWriteInterval=int(sim_properties['write_interval'])
    sim.speciesWriteInterval=int(sim_properties['write_interval'])
    sim.hookInterval=int(sim_properties['hook_interval'])
    replicates = 1
    sim_properties['total_time'] = sim_properties['total_time'] + totalTime
    
    print('Simulation Initialized')
    
    return sim, sim_properties
#########################################################################################


#########################################################################################
def buildRegionsRestart(sim, sim_properties):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    build = RegionBuilder(sim)
    
    regionIdx = 0
    
    N_edges = sim_properties['lattice_edges']
    
    region_dict = {}
    
#     regions = ["extracellular", "cytoplasm", "outer_cytoplasm", "p_starts", "p_mids", "p_ends", 
#                "ribosomes", "p_start_c", "p_mid_c", "p_end_c", "ribo_centers",
#                "DNA", "membrane"] #, "Z-ring"]

    regions = ["extracellular", "cytoplasm", "outer_cytoplasm",
               "ribosomes", "ribo_centers", "DNA", "membrane"]
    
    for region in regions:
        
        defined_region = sim.region(region)
        region_dict[region] = {}
        region_dict[region]['index'] = regionIdx
        
        region_shape = np.load(sim_properties['working_directory'] + 'restart_files/' + region + '.npy')
        
        region_dict[region]['shape'] = region_shape
        
        regionIdx = regionIdx + 1
        
        
    cytoplasm = build.ellipsoid(radius = sim_properties['cyto_radius'], center = sim_properties['lattice_center'])

    cyto_dilation = build.dilate(cytoplasm, se = build.se26)
    cyto_shell = cyto_dilation & ~cytoplasm
    
    region_dict['outer_cytoplasm']['full_shape'] = cyto_shell
    region_dict['cytoplasm']['full_shape'] = cytoplasm
    
        
    ribo_types = ['ribos'] #,'starts','mids','ends']
    
    ribo_site_dict = {}
    
    for ribo_type in ribo_types:
        ribo_site_dict[ribo_type] = {}
        
    ribo_site_dict['ribos']['center_idx'] = "ribo_centers"
    ribo_site_dict['ribos']['cross_idx'] = "ribosomes"
#     ribo_site_dict['starts']['center_idx'] = "p_start_c"
#     ribo_site_dict['starts']['cross_idx'] = "p_starts"
#     ribo_site_dict['mids']['center_idx'] = "p_mid_c"
#     ribo_site_dict['mids']['cross_idx'] = "p_mids"
#     ribo_site_dict['ends']['center_idx'] = "p_end_c"
#     ribo_site_dict['ends']['cross_idx'] = "p_ends"
    
    region_dict['ribosomes']['shape'] = region_dict['ribosomes']['shape']
    region_dict['ribo_centers']['shape'] = region_dict['ribo_centers']['shape']
    
#     print(ribosome_centers.shape)
#     print(len(np.argwhere(ribosomes==True)))
    
    ribo_site_dict['ribos']['centers'] = region_dict['ribo_centers']['shape']
    ribo_site_dict['ribos']['crosses'] = region_dict['ribosomes']['shape']
#     ribo_site_dict['starts']['centers'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['starts']['crosses'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['mids']['centers'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['mids']['crosses'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['ends']['centers'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['ends']['crosses'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)

#     RiboPtnMap = createRibosomalProteinMap(genome3A)
    
#     N_edges = sim_properties['lattice_edges']
#     DNAsites = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    build.compose(
    (sim.region('extracellular'), region_dict['extracellular']['shape']),
    (sim.region('cytoplasm'), region_dict['cytoplasm']['shape']),
    (sim.region('outer_cytoplasm'), region_dict['outer_cytoplasm']['shape']),
    (sim.region('ribosomes'), region_dict['ribosomes']['shape']),
    (sim.region('DNA'),region_dict['DNA']['shape']),
    (sim.region('ribo_centers'), region_dict['ribo_centers']['shape']),
    (sim.region('membrane'), region_dict['membrane']['shape']))
    
    print('Geometry constructed')
    
    return region_dict, ribo_site_dict
#########################################################################################


#########################################################################################
def initRdmeSpeciesRestart(sim, sim_properties, region_dict):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    print('Loading RDME particle lattice')
    
    plattice = np.load(sim_properties['working_directory'] + 'restart_files/plattice.npy')
    
    for specie, index in sim_properties["name_to_index"].items():
        
        spec = sim.species(specie)
        
        coords = np.argwhere(plattice==int(index))
        
        for coord in coords:
        
            spec.placeParticle(int(coord[1]), int(coord[2]), int(coord[3]), 1)
        
    Diff.defaultDiffusion(sim, region_dict)
    
    Diff.generalDiffusionConstants(sim)
    
    print('RDME particles placed')
    
    return None
#########################################################################################


#########################################################################################
def constructGIP(sim, sim_properties):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
#     createParticleIdxMap(sim, sim_properties)
    
    Rxns_RDME.general_reaction_rates(sim)
    
    addGeneticInformationReactions(sim, sim_properties)
    
    Rxns_RDME.replicationInitiation(sim, sim_properties, restart=True)
    
    print('Added RDME reactions')
    
    return None
#########################################################################################


#########################################################################################
def addGeneticInformationReactions(sim, sim_properties):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    genome = sim_properties['genome']

    for locusTag, locusDict in genome.items():
        
        locusNum = locusTag.split('_')[1]
        
        rnasequence = locusDict["RNAsequence"]
        
        if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
            
            Rxns_RDME.transcriptionLong(sim, sim_properties, locusTag, rnasequence, restart=True)
            
        else:
        
            Rxns_RDME.transcription(sim, sim_properties, locusTag, rnasequence, restart=True)
        
        if locusDict["Type"] == 'protein':
            
            Rxns_RDME.degradation_mrna(sim, sim_properties, locusNum, rnasequence, restart=True)

            aasequence = locusDict["AAsequence"]
            
            membranePtnCheck = 'C_P_'+locusNum in sim_properties['name_to_index']

            if membranePtnCheck:

                Rxns_RDME.translation(sim, locusNum, aasequence, membrane=True)

                Rxns_RDME.translocation_secy(sim, locusNum, aasequence)

            else:

                Rxns_RDME.translation(sim, locusNum, aasequence)

    return None
#########################################################################################


#########################################################################################
def constructAssemblyReactions(sim, sim_properties):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    Rxns_RDME.addRibosomeBiogenesis(sim, sim_properties)
    
    Rxns_RDME.addRNAPassembly(sim, sim_properties)
    
    degradosome = sim.species('Degradosome')
    RnaseY = sim.species('P_0359')
    RnaseJ1 = sim.species('P_0600')
    ptnAssoc = sim.rateConst('degAssoc', 1e7, 2)
    
    sim.region('outer_cytoplasm').addReaction([RnaseY, RnaseJ1], [degradosome], ptnAssoc)
    
    return None
#########################################################################################


#########################################################################################
def initDiffusionRestart(sim, sim_properties):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    genome = sim_properties['genome']
    
    #### Riboosme Diffusion ####
#     Diff.ribosomeDiffusion(sim)
    Diff.ribosomeDiffusion(sim, 'ribosomeP')
    
    #### Degradosome Diffusion ####
    Diff.degradosomeDiffusion(sim, 'Degradosome')
    
    #### RNAP Diffusion ####
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
    
    #### Protein Diffusion ####
    protein_ic = pd.read_excel(sim_properties['head_directory'] + 'input_data/initial_concentrations.xlsx', sheet_name='Comparative Proteomics')
    
    for row, protein in protein_ic.iterrows():
        
        if row != 0:
        
            locusTag = protein['Locus Tag']

            locusNum = locusTag.split('_')[1]

            proteinID = 'P_'+locusNum
            
            PTN = sim.species(proteinID)
            
            if protein['Localization'] == 'trans-membrane':
                
                cytoID = 'C_P_'+locusNum
                
                cPTN = sim.species(cytoID)
                
                Diff.proteinDiffusionCyto(sim, cytoID)

                Diff.proteinDiffusionMem(sim, proteinID)
                
            elif protein['Localization'] == 'peripheral membrane':
                
                Diff.proteinDiffusionPeriph(sim, proteinID)
                
#             elif protein['Gene Product'] 

            else:

                Diff.proteinDiffusionCyto(sim, proteinID)
            
    #### mRNA Diffusion ####        
    mRNA_ic = pd.read_excel(sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx',sheet_name='Diffusion Coefficient')
    
    for row, mRNA in mRNA_ic.iterrows():
        
        locusTag = mRNA['Locus Tag']
        
        if locusTag.startswith('JCVI'):
            
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum
            
            RNA = sim.species(rnaID)
            
            rnasequence = sim_properties['genome'][locusTag]['RNAsequence']

            Diff.rnaDiffusion(sim, rnaID, rnasequence)
            
            rnaReadID = 'R_'+locusNum+'_d'
            
            rnaRead = sim.species(rnaReadID)
            
            Diff.rnaDiffusion(sim, rnaReadID, rnasequence)
            
            
    #### tRNA Diffusion ####   
    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == "tRNA":
            
#             rnaID = locusDict['GeneName']
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum
            
            tRNA = sim.species(rnaID)
            
            rnasequence = sim_properties['genome'][locusTag]['RNAsequence']

            Diff.rnaDiffusion(sim, rnaID, rnasequence)
            
        
    #### rRNA Diffusion ####   
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
            
    print('Initialized diffusion rules for RDME particles')
    
    return None
#########################################################################################


#########################################################################################
def createParticleIdxMap(sim, sim_properties):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    index_to_name = sim.speciesList.__dict__['_id2name']
    
    name_to_index = {}
    
    for index, specieName in index_to_name.items():
        
        name_to_index[specieName] = int(index)
        
    sim_properties["name_to_index"] = name_to_index
        
    return None
#########################################################################################


#########################################################################################
def getMembraneProteinList(name_to_index):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    membrane_proteins = []
    
    for locusTag in genomePtnLocDict:
        
        locusNum = locusTag.split('_')[1]
    
        cytoPtn = 'C_P_'+locumNum
        
        if cytoPtn in name_to_index:
            
            PTN = 'P_'+locusNum
            
            membrane_proteins.append(PTN)
            
    return membrane_proteins
#########################################################################################




