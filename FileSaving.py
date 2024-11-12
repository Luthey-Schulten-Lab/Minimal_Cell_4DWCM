"""
Authors: Zane Thornburg

Functions that write simulation output to files other than the main .lm simulation file
"""

import numpy as np

import pandas as pd

import os

import json
import pickle

import argparse

import pySTDLM.PostProcessing as PP

from LatticeFunctions import *

#########################################################################################
from math import floor
from math import log10
def round_sig(x, sig=2):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    negative = False
    if x < 0:
        negative = True
    x = abs(x)
    if negative:
        return -1*round(x, sig-int(floor(log10(abs(x))))-1)
    elif x==0.0:
        return 0.0
    else:
        return round(x, sig-int(floor(log10(abs(x))))-1)
#########################################################################################


#########################################################################################
def saveSimArgs(sim_properties, args):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    argsFname = sim_properties['working_directory'] + 'sim_args.txt'
    
    with open(argsFname, 'w') as f:
    
        f.write('~~~~USER INPUT~~~~\n')
        
        for x in vars(args):
            f.write(x+' : ' + str(getattr(args, x)) + '\n')
        
    return None
#########################################################################################


#########################################################################################
def saveParticleCounts(time, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    simFolder = sim_properties['working_directory']
    
    if np.rint(time) == 0:

        particleDF = pd.DataFrame()

        spec_IDs = []

        new_counts = []

        for specID, count in sim_properties['counts'].items():

            spec_IDs.append(specID)

            new_counts.append(count)

        particleDF['Time']= spec_IDs
        particleDF[np.rint(time)] = new_counts

    #             print(pmap)

        particleDF.to_csv(simFolder+'particle_counts.csv',index=False)

    else:

    #             print(pmap)

        new_counts = []

        particleDF = pd.read_csv(simFolder+'particle_counts.csv')

        for specID, count in sim_properties['counts'].items():

            new_counts.append(count)

        particleDF[np.rint(time)] = new_counts

        particleDF.to_csv(simFolder+'particle_counts.csv',index=False)
        
    print('Saved Particle Counts to csv')
        
    return None      
#########################################################################################


#########################################################################################
def saveODEfluxes(time, sim_properties, odeResults, model, solver):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    simFolder = sim_properties['working_directory']
    
    if np.rint(time) == 1:
        
        resStart = odeResults[0,:]
        
        currentFluxes = solver.calcFlux(0, resStart)
#                         print("Is INF:", np.where( np.isinf(finalFluxes) ) )

        # Create list of reactions and fluxes
        fluxList = []
        fluxIDs = []
        for indx,rxn in enumerate(model.getRxnList()):
            fluxIDs.append(rxn.getID())
            fluxList.append( (round_sig(currentFluxes[indx], sig=3)) )

        fluxDF = pd.DataFrame()
        
        fluxDF['Time'] = fluxIDs
        fluxDF[np.rint(0)] = fluxList

        fluxDF.to_csv(simFolder+'metabolic_fluxes.csv',index=False)

    
    resFinal = odeResults[-1,:]
    
    currentFluxes = solver.calcFlux(0, resFinal)
    
    fluxList = []

    for indx,rxn in enumerate(model.getRxnList()):

        fluxList.append( (round_sig(currentFluxes[indx], sig=3)) )

    fluxDF = pd.read_csv(simFolder+'metabolic_fluxes.csv')

    fluxDF[np.rint(time)] = fluxList

    fluxDF.to_csv(simFolder+'metabolic_fluxes.csv',index=False)
        
    return None
#########################################################################################


#########################################################################################
def saveCountsAndFluxes(time, sim_properties, odeResults, model, solver):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    simFolder = sim_properties['working_directory']
    
    if np.rint(time) == 0:

        particleDF = pd.DataFrame()

        spec_IDs = []

        new_counts = []

        for specID, count in sim_properties['counts'].items():

            spec_IDs.append(specID)

            new_counts.append(count)

        particleDF['Time']= spec_IDs
        particleDF[np.rint(time)] = new_counts

    #             print(pmap)

        particleDF.to_csv(simFolder+'init_particle_counts.csv',index=False)
        
        print('Saved Particle Counts to csv')

    else:

    #             print(pmap)
    
        if np.rint(time) == 1:

            resStart = odeResults[0,:]

            currentFluxes = solver.calcFlux(0, resStart)
    #                         print("Is INF:", np.where( np.isinf(finalFluxes) ) )

            # Create list of reactions and fluxes
            fluxList = []
            fluxIDs = []
            for indx,rxn in enumerate(model.getRxnList()):
                fluxIDs.append('F_' + rxn.getID())
                fluxList.append( (round_sig(currentFluxes[indx], sig=3)) )

            fluxDF = pd.DataFrame()

            fluxDF['Time'] = fluxIDs
            fluxDF['0.0'] = fluxList
            
            particleDF = pd.read_csv(simFolder+'init_particle_counts.csv')
            
            resultsDF = pd.concat([particleDF,fluxDF],ignore_index=True)
            resultsDF.reset_index()

            resultsDF.to_csv(simFolder+'counts_and_fluxes.csv',index=False)

        new_results = []

        resultsDF = pd.read_csv(simFolder+'counts_and_fluxes.csv')

        for specID, count in sim_properties['counts'].items():

            new_results.append(count)

        resFinal = odeResults[-1,:]

        currentFluxes = solver.calcFlux(0, resFinal)

        for indx,rxn in enumerate(model.getRxnList()):

            new_results.append( (round_sig(currentFluxes[indx], sig=3)) )

        resultsDF[np.rint(time)] = new_results

        resultsDF.to_csv(simFolder+'counts_and_fluxes.csv',index=False)
        
        print('Saved Counts and Fluxes')
        
    return None
#########################################################################################


#########################################################################################
def writeSimProp(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    sim_properties_file = sim_properties['working_directory'] + 'sim_properties.pkl'
    
    try:
        os.remove(sim_properties_file)
    except:
        print('Nothing to delete')

    with open(sim_properties_file, 'wb') as SPfile: 
        
        pickle.dump(sim_properties, SPfile)
            
    print('Saved Sim Properties')
            
    return None
#########################################################################################


#########################################################################################
def writeLatticeFiles(sim_properties, lattice, region_dict):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    plattice = lattice.getParticleLatticeView()
    
    plattice_file = sim_properties['working_directory'] + 'restart_files/plattice.npy'
    
    try:
        os.remove(plattice_file)
    except:
        print('Nothing to delete')

    np.save(plattice_file, plattice)

    for region, type_dict in region_dict.items():
    
        region_file = sim_properties['working_directory'] + 'restart_files/' + region + '.npy'
    
        try:
            os.remove(region_file)
        except:
            print('Nothing to delete')
            
        np.save(region_file, type_dict['shape'])
        
        if 'cytoplasm' in region:
            
            region_file = sim_properties['working_directory'] + 'restart_files/' + region + '_full.npy'
            
            np.save(region_file, type_dict['full_shape'])
            
    print('Saved Lattice Files')
            
    return None
#########################################################################################


#########################################################################################
# def OLDsaveODEfluxes(sim_properties, odeResults, model, solver):
    
#     resFinal = odeResults[-1,:]

#     resStart = odeResults[0,:]
    
#     simFolder = sim_properties['working_directory']

# #     if (np.rint(time)/60).is_integer():
#     minute = np.rint(sim_properties['time'])
#     currentFluxes = solver.calcFlux(0, resStart)
# #                         print("Is INF:", np.where( np.isinf(finalFluxes) ) )

#     # Create list of reactions and fluxes
#     fluxList = []
#     for indx,rxn in enumerate(model.getRxnList()):
#         fluxList.append( (rxn.getID(), currentFluxes[indx]) )

#     fluxDF = pd.DataFrame(fluxList)

#     fluxFileName = simFolder+'fluxes/fluxDF_'+str(minute)+'_start.csv'

#     fluxDF.to_csv(fluxFileName,header=False)

# #             minute = int(int(time)/60)
#     currentFluxes = solver.calcFlux(0, resFinal )
# #                         print("Is INF:", np.where( np.isinf(finalFluxes) ) )

#     # Create list of reactions and fluxes
#     fluxList = []
#     for indx,rxn in enumerate(model.getRxnList()):
#         fluxList.append( (rxn.getID(), currentFluxes[indx]) )

#     fluxDF = pd.DataFrame(fluxList)

#     fluxFileName = simFolder+'fluxes/fluxDF_'+str(minute)+'_end.csv'

#     fluxDF.to_csv(fluxFileName,header=False)

#     fluxFileName = simFolder+'fluxes/fluxDF_'+str(minute)+'.csv'

#     fluxDF.to_csv(fluxFileName,header=False)

#     print('Saved fluxes at ' + str(minute) + ' minutes.')

# #     if time > totalTime-delt:
# #         print(time)
# #         finalFluxes = solver.calcFlux(0, resFinal )
# # #                         print("Is INF:", np.where( np.isinf(finalFluxes) ) )

# #         # Create list of reactions and fluxes
# #         fluxList = []
# #         for indx,rxn in enumerate(model.getRxnList()):
# #             fluxList.append( (rxn.getID(), finalFluxes[indx]) )

# #         fluxDF = pd.DataFrame(fluxList)

# #         fluxDF.to_csv(simFolder+'fluxDF_final.csv',header=False)

# #         print('Saved final fluxes.')
        
#     return None
#########################################################################################