"""
Authors: Zane Thornburg

Initialize main RDME simulation and sim_properties dictionary
Construction of RDME model
"""

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


#########################################################################################
def initSim(hook_step, write_step, totalTime, workingDirectoryName, headDirectory):
    """
    Inputs:
    
    Returns:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Called by:
    Description:
    """
    
    sim_properties = {}
    
#     headDirectory =  os.getcwd() + '/'
    
    sim_properties['head_directory'] = headDirectory

    try:
        os.makedirs(headDirectory + 'Data/')
    except:
        pass

    workingDirectory = headDirectory + 'Data/' + workingDirectoryName + '/'

    simFolder = workingDirectory

    sim_properties['working_directory'] = workingDirectory
    
    DNAfolder = sim_properties['working_directory']+'DNA/'
    
    try:
        os.makedirs(simFolder)
        os.makedirs(DNAfolder)
        print('Created sim directory')
    except:
        print('sim directory already exists')


    filename = simFolder + 'MinCell.lm'

    lattice_spacing = 10e-9
    sim_properties['lattice_spacing'] = lattice_spacing
    
    if lattice_spacing == 8e-9:

        N_edges_x = 96 # Number of subvolumes making up and edge of the simulation space N x N x N
        N_edges_y = 96
        N_edges_z = 128
    
    elif lattice_spacing == 10e-9:
        
        N_edges_x = 64 # Number of subvolumes making up and edge of the simulation space N x N x N
        N_edges_y = 64
        N_edges_z = 128
    
    N_edges = [N_edges_x,N_edges_y,N_edges_z]
    
    sim_properties['lattice_edges'] = N_edges

    N_2_x = N_edges_x/2
    N_2_y = N_edges_y/2
    N_2_z = N_edges_z/2
    
    N_2 = [N_2_x,N_2_y,N_2_z]
    
    sim_properties['lattice_center'] = N_2
    
    sim = RDMESim("JCVI-syn3A",
                  filename,
                  [N_edges_x,N_edges_y,N_edges_z],
                  lattice_spacing,
                  "extracellular",
                  latticeType='Int')

    radius_nm = 2.00e-7
    cyto_radius = int(round(radius_nm/sim.latticeSpacing)) #m converted to lattice sites (8 nm lattice spacing)
    print(cyto_radius)
    
    sim_properties['cyto_radius'] = cyto_radius
    sim_properties['cyto_radius_nm'] = radius_nm
    
    sim_properties['volume'] = (4/3)*np.pi*(sim_properties['cyto_radius']*sim_properties['lattice_spacing'])**3
    sim_properties['volume_L'] = sim_properties['volume']*1000
    
    dna_monomers = 46188

    cyto_vol = (4/3)*np.pi*0.200**3

    cyto_200 = (4/3)*np.pi*0.2**3

    ptn_ratio = (2.3e6*cyto_vol)/(2.3e6*cyto_200)

    sim_center = [N_2_x,N_2_y,N_2_z]
    
    sim_timestep = 50e-6

    sim.timestep = sim_timestep
    sim.simulationTime=totalTime
    sim.latticeWriteInterval=write_step
    sim.speciesWriteInterval=write_step
    sim.hookInterval=hook_step
    replicates = 1
    
    sim_properties['total_time'] = totalTime
    sim_properties['timestep'] = sim_timestep
    sim_properties['write_interval'] = write_step
    sim_properties['hook_interval'] = hook_step
    
    sim_properties['growth_step_counter'] = 11
    
    sim_properties['last_DNA_step'] = None
    
    sim_properties['rep_started'] = False
    
    sim_properties['rnap_spacing'] = int(400)
    
    sim_properties['secy_states'] = {}
#     sim_properties['cost_tracker'] = {}
    
    sim_properties['counts'] = {}
    
    print('Simulation Initialized')
    
    return sim, sim_properties
#########################################################################################


#########################################################################################
def constructGIP(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    createParticleIdxMap(sim, sim_properties)
    
    Rxns_RDME.general_reaction_rates(sim)
    
    addGeneticInformationReactions(sim, sim_properties)
    
    Rxns_RDME.replicationInitiation(sim, sim_properties)
    
    return None
#########################################################################################


#########################################################################################
def addGeneticInformationReactions(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    genome = sim_properties['genome']

    for locusTag, locusDict in genome.items():
        
        locusNum = locusTag.split('_')[1]
        
        rnasequence = locusDict["RNAsequence"]
        
        if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
            
            Rxns_RDME.transcriptionLong(sim, sim_properties, locusTag, rnasequence)
            
        else:
        
            Rxns_RDME.transcription(sim, sim_properties, locusTag, rnasequence)
        
        if locusDict["Type"] == 'protein':
            
            Rxns_RDME.degradation_mrna(sim, sim_properties, locusNum, rnasequence)

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
def createParticleIdxMap(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
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
def mapTranslationStates(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    translat_states = {}
    
    ribosomeParticleIDs = np.full((len(sim_properties['name_to_index'])), False)
    
    ribosomeParticleIDs[sim_properties['name_to_index']['ribosomeP']-1] = True
    
    genome = sim_properties['genome']

    for locusTag, locusDict in genome.items():
        
        if locusDict["Type"] == 'protein':
        
            locusNum = locusTag.split('_')[1]

            mRNA = 'R_'+locusNum
            
#             mRNA_read = 'R_'+locusNum+'_d'

            mRNA_ribo = 'RB_'+locusNum

            ribosomeParticleIDs[sim_properties['name_to_index'][mRNA_ribo]-1] = True

            mRNA_check = 'RB_'+locusNum+'_cp'
            mRNA_poly = 'RB_'+locusNum+'_p'
            mRNA_poly_f = 'RB_'+locusNum+'_pe'

            if mRNA_check in sim_properties['name_to_index']:
                
                translat_states[mRNA] = {}

                translat_states[mRNA]["Idx"] = sim_properties['name_to_index'][mRNA_ribo]
#                 translat_states[mRNA]["ReadRnaIdx"] = sim_properties['name_to_index'][mRNA_read]
                translat_states[mRNA]["CheckIdx"] = sim_properties['name_to_index'][mRNA_check]
                translat_states[mRNA]["PolyIdx"] = sim_properties['name_to_index'][mRNA_poly]
                translat_states[mRNA]["PolyFIdx"] = sim_properties['name_to_index'][mRNA_poly_f]

                ribosomeParticleIDs[sim_properties['name_to_index'][mRNA_check]-1] = True
                ribosomeParticleIDs[sim_properties['name_to_index'][mRNA_poly]-1] = True
                ribosomeParticleIDs[sim_properties['name_to_index'][mRNA_poly_f]-1] = True
            
    sim_properties['riboIDs'] = ribosomeParticleIDs
    sim_properties['translationStates'] = translat_states
    
    print('Mapped Translational States')
            
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


#########################################################################################
def constructAssemblyReactions(sim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    Rxns_RDME.addRibosomeBiogenesis(sim, sim_properties)
    
    Rxns_RDME.addRNAPassembly(sim, sim_properties)
    
    return None
#########################################################################################



