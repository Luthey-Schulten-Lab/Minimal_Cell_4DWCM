"""
Authors: Zane Thornburg
        alfiap (2026-05) -- async np.save for restart files (Patch 3)

Functions that write simulation output to files other than the main .lm simulation file
"""

import numpy as np

import pandas as pd

import os

import json
import pickle

import argparse
import glob
import atexit

from concurrent.futures import ThreadPoolExecutor

import pySTDLM.PostProcessing as PP

from LatticeFunctions import *


#########################################################################################
# Patch 3 (async np.save): the per-bucket lattice/region writes were ~0.78 s per
# call x 14400 calls = ~3 h on the master thread in job 1638.  np.save on
# multi-MB int8 arrays is dominated by pwrite() syscalls; offloading them to a
# small worker pool lets the next hook start while the previous save flushes.
#
# A copy() is taken before submitting so the array we hand to the worker is
# decoupled from the live LM lattice (which the next hook will mutate).  This
# costs a small amount of extra RAM on the way out but preserves the exact same
# byte-for-byte output as the synchronous path.
#
# At the start of every writeLatticeFiles call we wait for the previous batch
# to drain, so we never accumulate more than one batch of pending writes and
# the on-disk state at any "complete_steps>=80" tick is consistent.
#########################################################################################

_SAVE_EXECUTOR = None
_PENDING_SAVES = []


def _get_save_executor():
    global _SAVE_EXECUTOR
    if _SAVE_EXECUTOR is None:
        _SAVE_EXECUTOR = ThreadPoolExecutor(max_workers=2, thread_name_prefix='lattice_save')
        atexit.register(_drain_pending_saves)
    return _SAVE_EXECUTOR


def _atomic_np_save(path, arr):
    tmp = path + '.tmp'
    try:
        os.remove(tmp)
    except FileNotFoundError:
        pass
    np.save(tmp, arr)
    os.replace(tmp + '.npy' if not tmp.endswith('.npy') else tmp, path)


def _save_worker(path, arr):
    """Run inside the worker thread.  Just writes; errors propagate through the future."""
    np.save(path, arr)


def _drain_pending_saves():
    global _PENDING_SAVES
    if not _PENDING_SAVES:
        return
    for fut in _PENDING_SAVES:
        try:
            fut.result(timeout=300)
        except Exception as exc:
            print(f'WARNING: async lattice save failed: {exc}')
    _PENDING_SAVES = []


def _submit_save(path, arr):
    executor = _get_save_executor()
    fut = executor.submit(_save_worker, path, arr.copy())
    _PENDING_SAVES.append(fut)
#########################################################################################

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

        particleDF.to_csv(simFolder+'particle_counts.csv',index=False)

    else:

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
        
    fluxes = np.empty((0, len(model.getRxnList())), dtype=float)

    for i in range(odeResults.shape[0]):
        
        flux_transient = solver.calcFlux(0, odeResults[i,:]) # List - FIXED BUG

        fluxes = np.append(fluxes, [np.asarray(flux_transient)], axis=0)
        
    resFinal = odeResults[-1,:]
    
    currentFluxes = solver.calcFlux(0, resFinal)
    print('End Fluxes: ', currentFluxes.shape)
    
    averageFluxes = np.mean(fluxes, axis=0)
    print('Avg Fluxes: ', averageFluxes.shape)
    
    fluxList = []

    for indx,rxn in enumerate(model.getRxnList()):

        fluxList.append( (round_sig(averageFluxes[indx], sig=3)) )

    fluxDF = pd.read_csv(simFolder+'metabolic_fluxes.csv')

    fluxDF[np.rint(time)] = fluxList

    fluxDF.to_csv(simFolder+'metabolic_fluxes.csv',index=False)
        
    return None
#########################################################################################


#########################################################################################
def _concatenate_csv_files_optimized(temp_dir, final_csv, sim_properties):
    """
    OPTIMIZED: Concatenate CSV files more efficiently using batch processing
    
    Improvements:
    1. Batch reads multiple files at once (reduces I/O overhead)
    2. Processes in chunks of 20 files to reduce memory spikes
    3. Only reads final CSV if necessary
    4. Uses efficient merge operations
    5. Keeps individual files in temp folder after concatenation (preserved for backup)
    
    Inputs:
    temp_dir - Directory containing individual CSV files
    final_csv - Path to final concatenated CSV file
    sim_properties - Dictionary for tracking
    
    Returns:
    None
    """
    
    # Get all individual CSV files, sorted by time
    csv_files = sorted(glob.glob(temp_dir + 'counts_fluxes_*.csv'), 
                      key=lambda x: float(x.split('_')[-1].replace('.csv', '')))
    
    if not csv_files:
        return
    
    # OPTIMIZATION: Process in smaller batches if there are many files
    BATCH_SIZE = 20  # Process 20 files at a time to reduce memory usage
    
    # Read existing final CSV if it exists
    if os.path.exists(final_csv):
        try:
            combined_df = pd.read_csv(final_csv)
            existing_cols = set(combined_df.columns) - {'Time'}
        except Exception as e:
            print(f'WARNING: Could not read existing final CSV: {e}. Starting fresh.')
            combined_df = None
            existing_cols = set()
    else:
        combined_df = None
        existing_cols = set()
    
    # Batch read and process files
    new_dataframes = []
    
    for i in range(0, len(csv_files), BATCH_SIZE):
        batch_files = csv_files[i:i+BATCH_SIZE]
        batch_dfs = []
        
        # Read batch of files
        for csv_file in batch_files:
            try:
                temp_df = pd.read_csv(csv_file)
                time_col = [c for c in temp_df.columns if c != 'Time'][0]
                
                # Skip if already exists (restart scenario)
                if time_col in existing_cols:
                    print(f'SKIP: Time column {time_col} already exists. Skipping {csv_file}')
                    continue
                
                # Keep only Time and the data column
                batch_dfs.append(temp_df[['Time', time_col]])
                existing_cols.add(time_col)
            except Exception as e:
                print(f'WARNING: Failed to read {csv_file}: {e}')
                continue
        
        if batch_dfs:
            # Merge batch horizontally first
            if len(batch_dfs) > 1:
                batch_merged = batch_dfs[0]
                for df in batch_dfs[1:]:
                    batch_merged = pd.merge(batch_merged, df, on='Time', how='outer')
                new_dataframes.append(batch_merged)
            else:
                new_dataframes.append(batch_dfs[0])
    
    # Single merge operation for all new data
    if new_dataframes:
        if combined_df is not None:
            # Merge all new data with existing
            for new_df in new_dataframes:
                combined_df = pd.merge(combined_df, new_df, on='Time', how='outer')
        else:
            # Start with first dataframe
            combined_df = new_dataframes[0]
            for new_df in new_dataframes[1:]:
                combined_df = pd.merge(combined_df, new_df, on='Time', how='outer')
    
    if combined_df is None:
        print('WARNING: No data to concatenate.')
        return
    
    # Sort columns: Time first, then numeric time columns in order
    time_cols = [c for c in combined_df.columns if c != 'Time']
    time_cols_sorted = sorted(time_cols, key=lambda x: float(x))
    combined_df = combined_df[['Time'] + time_cols_sorted]
    
    # Write final concatenated file
    combined_df.to_csv(final_csv, index=False)
    
    # Keep individual files in temp folder (do not delete or move them)
    print(f'Concatenated {len(csv_files)} files into {final_csv}. Individual files remain in {temp_dir}')
    
    return None
#########################################################################################


#########################################################################################
def recoverFromCrash(sim_properties):
    """
    Recovery function for crash scenarios: Restore row_order and concat_counter if needed.
    
    Call this function during restart initialization to ensure data consistency.
    
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    None
    
    Description:
    - Restores row_order from final CSV if not in sim_properties
    - Resets concat_counter based on last time in final CSV
    """
    
    simFolder = sim_properties['working_directory']
    final_csv = simFolder + 'counts_and_fluxes.csv'
    temp_dir = simFolder + 'counts_fluxes_temp/'
    
    # Create temp directory if it doesn't exist
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    # Restore row_order from final CSV if not in sim_properties
    if 'row_order' not in sim_properties and os.path.exists(final_csv):
        try:
            final_df = pd.read_csv(final_csv)
            sim_properties['row_order'] = list(final_df['Time'])
            print('CRASH RECOVERY: Restored row_order from final CSV file.')
        except Exception as e:
            print(f'WARNING: Could not restore row_order from final CSV: {e}')
    
    # Reset concat_counter based on last time in final CSV
    if os.path.exists(final_csv):
        try:
            final_df = pd.read_csv(final_csv)
            time_cols = [c for c in final_df.columns if c != 'Time']
            if time_cols:
                # Get the last time column (highest numeric value)
                time_cols_sorted = sorted(time_cols, key=lambda x: float(x))
                last_time = float(time_cols_sorted[-1])
                # Estimate concat_counter (approximate, but safe)
                sim_properties['concat_counter'] = int(last_time) + 1
                print(f'CRASH RECOVERY: Reset concat_counter to {sim_properties["concat_counter"]} based on last time {last_time}.')
        except Exception as e:
            print(f'WARNING: Could not reset concat_counter: {e}')
            sim_properties['concat_counter'] = 0
    
    return None
#########################################################################################


#########################################################################################
def saveCountsAndFluxes(time, sim_properties, odeResults, model, solver):
    """
    Optimized version: Write individual CSV files per timepoint, concatenate at end only
    
    Inputs:
    time - Current biological time
    sim_properties - Dictionary of simulation variables and state trackers
    odeResults - ODE integration results
    model - Metabolic model object
    solver - ODE solver object
    
    Returns:
    None
    
    Description:
    Writes individual CSV files for each timepoint (fast, no read needed), then
    concatenates them ONLY at the end of the simulation. This avoids periodic
    I/O overhead during simulation.
    
    Note: recoverFromCrash() should be called at restart initialization to handle
    orphaned files from previous crashes.
    """
    
    simFolder = sim_properties['working_directory']
    final_csv = simFolder + 'counts_and_fluxes.csv'
    temp_dir = simFolder + 'counts_fluxes_temp/'
    
    # Create temp directory if it doesn't exist
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    # Initialize concat_counter if not set (from recoverFromCrash or previous run)
    if 'concat_counter' not in sim_properties:
        sim_properties['concat_counter'] = 0
    
    if np.rint(time) == 0:
        # RESTART CHECK: If init_particle_counts.csv exists and row_order is set, skip initialization
        init_file = simFolder+'init_particle_counts.csv'
        if os.path.exists(init_file) and 'row_order' in sim_properties:
            print('SKIP: init_particle_counts.csv already exists and row_order is set (restart scenario). Skipping time==0 initialization.')
            return None
        
        # Initialize: Create particle counts DataFrame structure
        particleDF = pd.DataFrame()
        spec_IDs = []
        new_counts = []
        
        for specID, count in sim_properties['counts'].items():
            spec_IDs.append(specID)
            new_counts.append(count)
        
        particleDF['Time'] = spec_IDs
        particleDF[np.rint(time)] = new_counts
        particleDF.to_csv(init_file, index=False)
        print('Saved Particle Counts to csv')
        
        # Store the row order (species IDs) - will be extended with flux IDs at time=1
        sim_properties['row_order'] = spec_IDs
        
    else:
        if np.rint(time) == 1:
            # RESTART CHECK: If final_csv exists and has 0.0 and 1.0, skip initialization
            if os.path.exists(final_csv):
                try:
                    final_df = pd.read_csv(final_csv)
                    if '0.0' in final_df.columns and '1.0' in final_df.columns:
                        print('SKIP: Timepoints 0.0 and 1.0 already exist in final CSV (restart scenario). Skipping initialization.')
                        # Restore row_order from final CSV
                        sim_properties['row_order'] = list(final_df['Time'])
                        # Set concat_counter appropriately
                        time_cols = [c for c in final_df.columns if c != 'Time']
                        if time_cols:
                            time_cols_sorted = sorted(time_cols, key=lambda x: float(x))
                            last_time = float(time_cols_sorted[-1])
                            sim_properties['concat_counter'] = int(last_time) + 1
                        return None
                except Exception as e:
                    print(f'WARNING: Could not check final CSV at time==1: {e}')
            
            # First time: Get flux IDs and create complete row order
            resStart = odeResults[0,:]
            currentFluxes = solver.calcFlux(0, resStart)
            
            # Create flux IDs
            fluxIDs = []
            for indx, rxn in enumerate(model.getRxnList()):
                fluxIDs.append('F_' + rxn.getID())
            
             # Create complete row order: species + fluxes
            # RESTART CHECK: If init_particle_counts.csv doesn't exist but final_csv does, restore from final_csv
            init_file = simFolder+'init_particle_counts.csv'
            if not os.path.exists(init_file) and os.path.exists(final_csv):
                try:
                    final_df = pd.read_csv(final_csv)
                    # Extract species IDs (everything before flux IDs which start with 'F_')
                    all_ids = list(final_df['Time'])
                    species_ids = [s for s in all_ids if not s.startswith('F_')]
                    print(f'RESTART: Restored species IDs from final CSV (found {len(species_ids)} species).')
                except Exception as e:
                    raise ValueError(f"init_particle_counts.csv missing and could not restore from final CSV: {e}")
            else:
                particleDF = pd.read_csv(init_file)
                species_ids = list(particleDF['Time'])
            
            sim_properties['row_order'] = species_ids + fluxIDs
            
            # Write time 0.0 data (from init_particle_counts.csv + initial fluxes)
            fluxList = [round_sig(currentFluxes[indx], sig=3) for indx in range(len(model.getRxnList()))]
            
            # Get the data column name dynamically (the column that's not 'Time')
            data_cols = [c for c in particleDF.columns if c != 'Time']
            if len(data_cols) == 0:
                raise ValueError("No data column found in init_particle_counts.csv")
            data_col = data_cols[0]  # Get first data column (should be '0' or '0.0')
            time0_data = list(particleDF[data_col]) + fluxList
            
            time0_df = pd.DataFrame({
                'Time': sim_properties['row_order'],
                '0.0': time0_data
            })
            time0_df.to_csv(temp_dir + 'counts_fluxes_0.0.csv', index=False)
            
            # Write time 1.0 data (current time)
            new_results = []
            for specID in species_ids:
                new_results.append(sim_properties['counts'][specID])
            
            # Calculate average fluxes for time 1.0
            fluxes = np.empty((0, len(model.getRxnList())), dtype=float)
            for i in range(odeResults.shape[0]):
                flux_transient = solver.calcFlux(0, odeResults[i,:])
                fluxes = np.append(fluxes, [np.asarray(flux_transient)], axis=0)
            averageFluxes = np.mean(fluxes, axis=0)
            
            for indx in range(len(model.getRxnList())):
                new_results.append(round_sig(averageFluxes[indx], sig=3))
            
            time1_df = pd.DataFrame({
                'Time': sim_properties['row_order'],
                '1.0': new_results
            })
            time1_df.to_csv(temp_dir + 'counts_fluxes_1.0.csv', index=False)
            
            # No concatenation at time==1 anymore - will happen at end
            sim_properties['concat_counter'] = 2
            print('Initialized individual CSV files (concatenation will occur at simulation end)')
        
        elif np.rint(time) > 1:  # For time > 1: Write individual file
            # Safety check
            if 'row_order' not in sim_properties:
                # Try to restore from final CSV if it exists
                if os.path.exists(final_csv):
                    try:
                        final_df = pd.read_csv(final_csv)
                        sim_properties['row_order'] = list(final_df['Time'])
                        print('RESTART: Restored row_order from final CSV file.')
                    except Exception as e:
                        raise ValueError(f"row_order not initialized and could not restore from final CSV: {e}")
                else:
                    raise ValueError("row_order not initialized. This should happen at time==1.")
            
            time_int = np.rint(time)
            # Use integer format for consistency (time_int is already rounded)
            time_str = str(int(time_int))
            
            # RESTART CHECK: Skip if this timepoint already exists in final CSV (restart scenario)
            if os.path.exists(final_csv):
                try:
                    final_df = pd.read_csv(final_csv)
                    # Check both integer and float string formats for consistency
                    time_str_float = str(float(time_int))
                    if time_str in final_df.columns or time_str_float in final_df.columns:
                        print(f'SKIP: Timepoint {time_int} already exists in final CSV (restart scenario). Skipping write.')
                        # Update concat_counter if needed
                        if sim_properties.get('concat_counter', 0) <= int(time_int):
                            sim_properties['concat_counter'] = int(time_int) + 1
                        return None
                except Exception as e:
                    print(f'WARNING: Could not check final CSV for existing timepoint: {e}')
            
            # Collect new results in the correct order
            new_results = []
            species_ids = [s for s in sim_properties['row_order'] if s in sim_properties['counts']]
            
            for specID in species_ids:
                new_results.append(sim_properties['counts'][specID])
            
            # Calculate average fluxes
            fluxes = np.empty((0, len(model.getRxnList())), dtype=float)
            for i in range(odeResults.shape[0]):
                flux_transient = solver.calcFlux(0, odeResults[i,:])
                fluxes = np.append(fluxes, [np.asarray(flux_transient)], axis=0)
            
            averageFluxes = np.mean(fluxes, axis=0)
            
            # Add flux values (they come after species counts in row_order)
            for indx in range(len(model.getRxnList())):
                new_results.append(round_sig(averageFluxes[indx], sig=3))
            
            # Write individual CSV file for this timepoint (FAST - no read needed!)
            temp_file = temp_dir + f'counts_fluxes_{time_str}.csv'
            
            # Skip if temp file already exists (shouldn't happen, but safety check)
            if os.path.exists(temp_file):
                print(f'WARNING: Temp file {temp_file} already exists. Skipping write to avoid overwrite.')
                return None
            
            timepoint_df = pd.DataFrame({
                'Time': sim_properties['row_order'],
                time_str: new_results
            })
            timepoint_df.to_csv(temp_file, index=False)
            
            sim_properties['concat_counter'] = sim_properties.get('concat_counter', 0) + 1
            
            # ONLY concatenate at the end of simulation
            time_int_val = int(time_int)
            write_counter = sim_properties['concat_counter']
            total_time = sim_properties.get('total_time', None)
            
            is_final = False
            if total_time is not None:
                # Allow small tolerance for floating point comparison
                if abs(time - total_time) < 0.1 or time >= total_time:
                    is_final = True
            
            if is_final:
                # Final concatenation: Get ALL remaining temp files (handles any that weren't caught earlier)
                remaining_files = sorted(glob.glob(temp_dir + 'counts_fluxes_*.csv'), 
                                       key=lambda x: float(x.split('_')[-1].replace('.csv', '')))
                if remaining_files:
                    _concatenate_csv_files_optimized(temp_dir, final_csv, sim_properties)
                    print(f'FINAL: Concatenated all remaining CSV files at simulation end (time={time_int_val}, write #{write_counter}, {len(remaining_files)} files)')
                else:
                    print(f'FINAL: No remaining temp files to concatenate at end (time={time_int_val})')
            else:
                print(f'Saved individual CSV file (time={time_int_val}, write #{write_counter}, will concatenate at end)')
        
    return None
#########################################################################################


#########################################################################################
def finalizeCountsAndFluxes(sim_properties):
    """
    Final cleanup function: Concatenate all remaining individual CSV files into final CSV.
    
    Call this function at the end of simulation to ensure all temp files are concatenated.
    This should be called even if saveCountsAndFluxes doesn't hit the exact total_time.
    This is especially important for restart scenarios where the simulation might end
    without hitting the exact total_time condition.
    
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    None
    """
    simFolder = sim_properties['working_directory']
    final_csv = simFolder + 'counts_and_fluxes.csv'
    temp_dir = simFolder + 'counts_fluxes_temp/'
    
    if not os.path.exists(temp_dir):
        print('FINALIZE: No temp directory found. Nothing to finalize.')
        return None
    
    # Check for any remaining temp files
    csv_files = sorted(glob.glob(temp_dir + 'counts_fluxes_*.csv'), 
                      key=lambda x: float(x.split('_')[-1].replace('.csv', '')))
    
    if csv_files:
        print(f'FINALIZE: Found {len(csv_files)} remaining temp CSV files. Concatenating...')
        _concatenate_csv_files_optimized(temp_dir, final_csv, sim_properties)
        print('FINALIZE: Completed final concatenation.')
    else:
        print('FINALIZE: No remaining temp files to concatenate.')
    
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
    
    # Patch 3: drain any in-flight saves from the previous bucket so the on-disk
    # state is consistent before we kick off a new batch.
    _drain_pending_saves()

    plattice = lattice.getParticleLatticeView()
    plattice_file = sim_properties['working_directory'] + 'restart_files/plattice.npy'
    _submit_save(plattice_file, np.asarray(plattice))

    for region, type_dict in region_dict.items():
        region_file = sim_properties['working_directory'] + 'restart_files/' + region + '.npy'
        _submit_save(region_file, np.asarray(type_dict['shape']))

        if 'cytoplasm' in region:
            region_file = sim_properties['working_directory'] + 'restart_files/' + region + '_full.npy'
            _submit_save(region_file, np.asarray(type_dict['full_shape']))

    print('Saved Lattice Files (async, %d in flight)' % len(_PENDING_SAVES))

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