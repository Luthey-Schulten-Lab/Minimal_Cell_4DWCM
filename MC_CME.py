"""
Authors: Zane Thornburg

Create a global CME simulation
"""

##### CME Model #####

import Rxns_CME

from pyLM import CME

import numpy as np

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict

import os

import time as TIME


#########################################################################################
def runGCME(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    csimFolder = sim_properties['working_directory']+'CME/'
    
    # Create simulation
    csim=CME.CMESimulation()
    
    add_CME_species(csim, sim_properties)
    
    constructCME(csim, sim_properties)

    # Set time data
    csim.setWriteInterval(0.1)
    csim.setSimulationTime(1.0)
    
    # Save and run simulation
    CSIMfilename= csimFolder + 'cmeSim.%d.lm'%np.rint(sim_properties['time'])
    
    try:
        os.remove(CSIMfilename)
    except:
        print('Nothing to delete')
    csim.save(CSIMfilename)
    
    pythonExecutable = sim_properties['head_directory'] + 'Run_CME.py' # sim_properties['working_directory'] + 'run_CME.py'
    
    os.system("python %s %s"% (pythonExecutable, CSIMfilename) )

    return None
#########################################################################################


#########################################################################################
def constructCME(csim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    genome = sim_properties['genome']

    for locusTag, locusDict in genome.items():
        
        rnasequence = locusDict["RNAsequence"]
        
        if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
            
            Rxns_CME.transcriptionLong(csim, sim_properties, locusTag, rnasequence)
            
        else:
        
            Rxns_CME.transcription(csim, sim_properties, locusTag, rnasequence)
        
    Rxns_CME.tRNAcharging(csim, sim_properties)
    
    return None
#########################################################################################


#########################################################################################
def add_CME_species(csim, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    start = TIME.time()
    
    specList = sim_properties['cme_species']
    
    for specID in specList:
        
        csim.defineSpecies([specID])
        
        if specID in sim_properties['cme_state_tracker']:
            
            csim.addParticles(specID, count=int(sim_properties['counts'][specID]-sim_properties['cme_state_tracker'][specID]))
            
        else:
        
            csim.addParticles(specID, count=int(sim_properties['counts'][specID]))
        
    print('Time to add species: ', TIME.time()-start)
        
    return None
#########################################################################################
    


        