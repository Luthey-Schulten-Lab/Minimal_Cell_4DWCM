"""
Authors: Zane Thornburg
        alfiap (2026-05) -- persistent gCME worker (kills `os.system` overhead)

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
import sys
import atexit
import subprocess
import threading

import time as TIME


#########################################################################################
# Persistent CME solver worker
#
# Replaces `os.system("python Run_CME.py <fname>")` per gCME call (1 Hz) with a
# single long-lived Python subprocess that already has `lm` and `GillespieDSolver`
# imported.  Saves the lm-import + CUDA-init overhead on every call.
#
# Per-call breakdown from job 1638:
#   * Python setup:      ~0.08 s  (csim build)        -- not affected
#   * .lm file save:     ~0.90 s  (pyLM HDF5 writer)  -- not affected
#   * `os.system` solve: ~0.65 s                      -- TARGETED
# Of that ~0.65 s, ~0.5 s is python startup + lm import; the worker eliminates
# that on every call after the first.
#
# Falls back gracefully to the old `os.system` path if the worker fails to
# start, dies, or returns an error -- this preserves the original semantics
# in any failure mode.
#########################################################################################

_WORKER_PROC = None
_WORKER_LOCK = threading.Lock()
_WORKER_DEAD = False  # latched if worker dies; we revert to os.system from then on


def _worker_script_path(sim_properties):
    return sim_properties['head_directory'] + 'Run_CME_Worker.py'


def _ensure_worker(sim_properties):
    """Start the persistent CME worker once. Idempotent. Thread-safe."""
    global _WORKER_PROC, _WORKER_DEAD

    if _WORKER_DEAD:
        return None
    if _WORKER_PROC is not None and _WORKER_PROC.poll() is None:
        return _WORKER_PROC

    with _WORKER_LOCK:
        if _WORKER_PROC is not None and _WORKER_PROC.poll() is None:
            return _WORKER_PROC

        worker_script = _worker_script_path(sim_properties)
        if not os.path.isfile(worker_script):
            print(f"CME_WORKER: script not found at {worker_script} -- falling back to os.system", flush=True)
            _WORKER_DEAD = True
            return None

        try:
            _WORKER_PROC = subprocess.Popen(
                [sys.executable, '-u', worker_script],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
                universal_newlines=True,
            )
        except Exception as exc:
            print(f"CME_WORKER: failed to start ({exc}) -- falling back to os.system", flush=True)
            _WORKER_DEAD = True
            return None

        # Drain the "started" banner so subsequent reads only see DONE/ERR lines.
        ready = False
        for _ in range(50):
            line = _WORKER_PROC.stdout.readline()
            if not line:
                break
            print('  ' + line.rstrip())
            if 'CME_WORKER: started' in line:
                ready = True
                break
        if not ready:
            print("CME_WORKER: did not announce ready -- falling back to os.system", flush=True)
            try:
                _WORKER_PROC.kill()
            except Exception:
                pass
            _WORKER_PROC = None
            _WORKER_DEAD = True
            return None

        atexit.register(_shutdown_worker)
        print("CME_WORKER: persistent worker is ready", flush=True)
        return _WORKER_PROC


def _shutdown_worker():
    global _WORKER_PROC
    proc = _WORKER_PROC
    if proc is None:
        return
    try:
        if proc.poll() is None:
            proc.stdin.write("EXIT\n")
            proc.stdin.flush()
            proc.wait(timeout=10)
    except Exception:
        try:
            proc.kill()
        except Exception:
            pass
    finally:
        _WORKER_PROC = None


def _run_via_worker(sim_properties, csim_filename):
    """Send filename to worker; return True on success, False to trigger fallback."""
    global _WORKER_PROC, _WORKER_DEAD

    proc = _ensure_worker(sim_properties)
    if proc is None:
        return False

    try:
        proc.stdin.write(csim_filename + "\n")
        proc.stdin.flush()
    except Exception as exc:
        print(f"CME_WORKER: stdin write failed ({exc}) -- falling back to os.system", flush=True)
        _WORKER_DEAD = True
        return False

    while True:
        line = proc.stdout.readline()
        if not line:
            print("CME_WORKER: stdout EOF (worker died) -- falling back to os.system", flush=True)
            _WORKER_DEAD = True
            return False
        line = line.rstrip()
        if line.startswith("DONE "):
            return True
        if line.startswith("ERR "):
            print(f"CME_WORKER: solver error -> {line}", flush=True)
            return False
        # otherwise: a stray print from the worker -- echo it and keep reading
        print('  ' + line, flush=True)
#########################################################################################


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
    
    # Time Python setup
    setupstart = TIME.time()
    csim=CME.CMESimulation()
    add_CME_species(csim, sim_properties)
    constructCME(csim, sim_properties)
   
    print('='*60)
    print('CME Model Statistics:')
    print(f'  Number of species: {len(csim.species_id)}')
    print(f'  Number of reactions: {len(csim.reactions)}')
    print('='*60)

    csim.setWriteInterval(0.1)
    csim.setSimulationTime(1.0)
    setuptime = TIME.time() - setupstart
    
    # Save and run simulation
    CSIMfilename= csimFolder + 'cmeSim.%d.lm'%np.rint(sim_properties['time'])
    
    try:
        os.remove(CSIMfilename)
    except:
        print('Nothing to delete')
    
    # Time file save
    savestart = TIME.time()
    csim.save(CSIMfilename)
    savetime = TIME.time() - savestart
    
    # Time C++ solver execution
    solverstart = TIME.time()
    used_worker = _run_via_worker(sim_properties, CSIMfilename)
    if not used_worker:
        # Fallback path -- preserves identical semantics to the legacy code.
        pythonExecutable = sim_properties['head_directory'] + 'Run_CME.py'
        os.system("python %s %s" % (pythonExecutable, CSIMfilename))
    solvertime = TIME.time() - solverstart
    if used_worker:
        print('  CME solver path: persistent worker')
    else:
        print('  CME solver path: os.system fallback')
    
    print('  CME solver breakdown:')
    print('    Python setup:  {:.4f}s'.format(setuptime))
    print('    File save:     {:.4f}s'.format(savetime))
    print('    C++ solver:    {:.4f}s ({:.1f}% of CME solver time)'.format(solvertime, 100*solvertime/(setuptime+savetime+solvertime)))

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
    


        