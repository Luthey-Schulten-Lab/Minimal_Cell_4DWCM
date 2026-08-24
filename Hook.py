"""
RDME hook solver for the 4DWCM.

Authors
-------
Alfia Parvez
Zane Thornburg

Coordinates growth/division, chromosome BD, ribosome placement, global CME,
and ODE metabolism at each RDME interrupt.
"""

import Growth as growth
import Division as division
import RibosomesRDME as ribosomesRDME
import SpatialDnaDynamics as DNA
import MC_CME as MCCME
import Rxns_ODE as ODE
import ImportInitialConditions as IC
import Integrate as integrate
import Communicate as communicate
import FileSaving as save

from jLM.RegionBuilder import RegionBuilder
from jLM.RDME import Sim as RDMESim
from jLM.RDME import File as RDMEFile
import jLM

import time as TIME

import numpy as np
import os
from math import floor, log10

from datetime import datetime


class MyOwnSolver:
    
    def __init__(self, lmFile, sim_properties, region_dict, ribo_site_dict, termination_time=None):
        
        super(MyOwnSolver, self).__init__()
    
        if isinstance(lmFile, (RDMEFile, RDMESim)):
            self.rdme = lmFile
        else:
            self.rdme = RDME.File(lmFile)
            
        self.sim_properties = sim_properties

        # The time a which hook solver has been stepped into, initial value = 0
        self.hookStartTime = 60
        self.nextWriteStep = int(0)
        self.writeInterval = int(sim_properties['write_interval'])
        
        self.complete_steps = 0
        
        self.termination_time = termination_time
        
        if termination_time is not None:
            start_dt = datetime.now()
            start_timestamp = round(start_dt.timestamp())
            self.start_datetime = datetime.fromtimestamp(int(start_timestamp))
            self.final_simulation_time = None

        self.sim_properties['fluxes'] = None

        self.sim_properties['rep_started'] = False
        
        self.sim_properties['next_gamma_V'] = 0.99
        
        self.sim_properties['gamma_V'] = 1.0
        
        self.sim_properties['division_started'] = False
        
        self.next_growth_time = 1.0
        
        self.next_gCME_time = 1.0
        
        self.next_DNA_time = 0.0
        
        self.next_metabolism_time = 1.0
        
        self.translation_update_step = 8
        self.ribo_step = 0

        # Prefer Cython ODE solver with a build-once cache (see Integrate.setSolverCached).
        self._ode_use_cython = True
        self._ode_cython_failed = False
        self._ode_solver_cache = {}

        # Ribosome EV placement every hook; translation_update_step gates polysome updates.
        self.ribo_place_every_n_hooks = 1
        self._ribo_place_counter = 0

        self.endLastHook = TIME.time()
        
        self.region_dict = region_dict
        self.ribo_site_dict = ribo_site_dict
        
        self.N_edges = sim_properties['lattice_edges']
        
        self.ribo_IDs = sim_properties['riboIDs'] #ribo_IDs
        
        self.sim_properties['last_last_DNA_step'] = None
        
        try:
            csim_folder=sim_properties['working_directory']+'CME/'
            os.mkdir(csim_folder)
            print('Created global CME directory')
        except:
            print('CME sim directory already exists')
            
        try:
            flux_folder=sim_properties['working_directory']+'fluxes/'
            os.mkdir(flux_folder)
            print('Created metabolic fluxes directory')
        except:
            print('Metabolic fluxes directory already exists')
            
        try:
            flux_folder=sim_properties['working_directory']+'restart_files/'
            os.mkdir(flux_folder)
            print('Created restart files directory')
        except:
            print('Restart files directory already exists')
            

        print('Initialized Solver')

        return None
        
    def hookSimulation(self, t, lattice):
        """Called each RDME interrupt. Return 1 or 2 to push lattice changes to the GPU."""

        
        time = t
        
        print('Current biological time: ', time)
        
        self.sim_properties['time'] = time
        
        print('Time between hook steps: ', TIME.time()-self.endLastHook)
        
        if time == 0:
            
            communicate.updateCountsRDME(self.rdme, self.sim_properties, lattice)
            
            IC.setCmeSpeciesList(self.sim_properties)

        if (time > 0):

            # Optional wall-clock cutoff (maximumHours) — terminate cleanly before queue kill.
            
            if self.termination_time is not None:
                
                current_datetime = datetime.now()
                
                elapsed = current_datetime-self.start_datetime
                elapsed_hours = elapsed.total_seconds()/3600
                
                if elapsed_hours>=self.termination_time:
                    
                    if self.final_simulation_time is None:
                        
                        self.final_simulation_time = self.next_DNA_time + 7.5
                
                    if (time >= self.final_simulation_time):
                        
                        raise Exception(f"Reached end of allowed simulation time after {elapsed_hours} hours. Solver will terminate without saving.")
            
            # Chromosome BD + morphology
            if (time >= self.next_DNA_time):
                
                print("Updating SA and Volume")
                
                dnastart = TIME.time()
                
                updateRegions = communicate.updateSA(self.sim_properties)

                print("Initializing Chromosome Update")

                if not self.sim_properties['rep_started']:

                    repInitState = communicate.checkRepInitState(self.sim_properties)

                    if repInitState:

                        self.sim_properties['rep_started'] = True

                        print('REPLICATION STARTED')
                        
                if self.sim_properties['division_started'] and updateRegions:
                    
                    region_dict, DNA_lattice_coords, genome = DNA.updateChromosomeDivision(time, lattice, self.sim_properties, self.region_dict, self.ribo_site_dict)
                    
                else:

                    region_dict, DNA_lattice_coords, genome = DNA.updateChromosome(time, lattice, self.sim_properties, self.region_dict, self.ribo_site_dict, updateRegions)

                self.region_dict = region_dict

                self.sim_properties['DNAcoords'] = DNA_lattice_coords

                self.sim_properties['genome'] = genome
                
                print("Chromosome Update Complete")

                if updateRegions:
                    
                    print('Updating morphology')
                    
                    if self.sim_properties['division_started']:
                        
                        region_dict, ribo_site_dict = division.divide_cell(self.rdme, lattice, self.sim_properties, self.region_dict, self.ribo_site_dict)
                        
                        self.ribo_site_dict = ribo_site_dict
                        
                    else:

                        region_dict = growth.grow_cell(self.rdme, lattice, self.sim_properties, self.region_dict, self.ribo_site_dict, self.sim_properties['cyto_radius'])

                    self.region_dict = region_dict

                _dna_hook_s = float(self.sim_properties.get('dna_hook_interval_s', 4.0))
                self.next_DNA_time = self.next_DNA_time + _dna_hook_s
                
                print('DNA time: ', TIME.time()-dnastart,
                      '(next DNA @ {:.1f}s, interval={:.1f}s)'.format(
                          self.next_DNA_time, _dna_hook_s))

                print('Updated cell architecture')


            # Ribosome excluded volume / polysomes
            ribostart = TIME.time()

            self.ribo_step = self.ribo_step + 1

            if self.ribo_step >= self.translation_update_step:

                updateTranslat = True

                self.ribo_step = 0

            else:

                updateTranslat = False

            # Place EV every N hooks; always place when translation state updates.
            self._ribo_place_counter += 1
            should_place = (
                updateTranslat
                or (self._ribo_place_counter >= self.ribo_place_every_n_hooks)
            )

            if should_place:
                print("Moving ribosomes")
                self._ribo_place_counter = 0

                ribo_site_dict = ribosomesRDME.placeRibosomes(
                    lattice, self.sim_properties, self.region_dict, self.ribo_site_dict,
                    updateTranslat=updateTranslat,
                )

                region_dict = ribosomesRDME.updateRiboSites(
                    lattice, ribo_site_dict, self.region_dict, self.sim_properties,
                )

                self.region_dict = region_dict
                self.ribo_site_dict = ribo_site_dict

                print('Ribo time: ', TIME.time()-ribostart)
                print("Moved ribosomes")
            else:
                print('Ribo time (skipped): ', TIME.time()-ribostart)


            # Global CME
            if (time >= self.next_gCME_time):

                cmestart = TIME.time()

                prestart = TIME.time()
                communicate.updateCountsRDME(self.rdme, self.sim_properties, lattice)
                communicate.calculateTranslationCosts(self.sim_properties, lattice)
                communicate.updateLongGeneStates(self.sim_properties, lattice)
                pretime = TIME.time() - prestart

                cmesolverstart = TIME.time()
                MCCME.runGCME(self.sim_properties)
                cmesolvertime = TIME.time() - cmesolverstart

                poststart = TIME.time()

                countsCMEstart = TIME.time()
                communicate.updateCountsCME(self.sim_properties)
                countsCMEtime = TIME.time() - countsCMEstart
                
                transStatesstart = TIME.time()
                communicate.updateTranscriptionStates(self.sim_properties, lattice)
                transStatestime = TIME.time() - transStatesstart

                # Skip redundant post-CME updateCountsRDME (next pre-CME rescans).
                countsRDMetime = 0.0

                posttime = TIME.time() - poststart

                self.next_gCME_time = self.next_gCME_time + 1.0
                
                totaltime = TIME.time() - cmestart
                print('CME time breakdown:')
                print('  Pre-CME Python:     {:.4f}s ({:.1f}%)'.format(pretime, 100*pretime/totaltime))
                print('  CME solver (C++):    {:.4f}s ({:.1f}%)'.format(cmesolvertime, 100*cmesolvertime/totaltime))
                print('  Post-CME Python:     {:.4f}s ({:.1f}%)'.format(posttime, 100*posttime/totaltime))
                print('    updateCountsCME:      {:.4f}s ({:.1f}% of post-CME)'.format(countsCMEtime, 100*countsCMEtime/posttime if posttime > 0 else 0))
                print('    updateTranscriptionStates: {:.4f}s ({:.1f}% of post-CME)'.format(transStatestime, 100*transStatestime/posttime if posttime > 0 else 0))
                print('    updateCountsRDME:     {:.4f}s ({:.1f}% of post-CME)'.format(countsRDMetime, 100*countsRDMetime/posttime if posttime > 0 else 0))
                print('  CME time:      {:.4f}s'.format(totaltime))



            # Metabolism (ODE)
            if (time >= self.next_metabolism_time):

                odestart = TIME.time()

                communicate.calculateDegradationCosts(self.sim_properties, lattice)

                communicate.communicateCostsToMetabolism(self.sim_properties)

                print('Initializing ODE simulation')
                model = ODE.initModel(self.sim_properties)
                print('Initialized ODE simulation')

                initVals = integrate.getInitVals(model)

                # Cython ODE path with noCython fallback for the rest of the run on failure.
                if self._ode_use_cython and not self._ode_cython_failed:
                    try:
                        solver = integrate.setSolverCached(model, self._ode_solver_cache)
                    except Exception as exc:
                        print(f"ODE: Cython solver build failed ({exc}); falling back to noCython for the remainder of the run.")
                        self._ode_cython_failed = True
                        self._ode_solver_cache = {}
                        solver = integrate.noCythonSetSolver(model)
                else:
                    solver = integrate.noCythonSetSolver(model)

                odeResults = integrate.runODE(initVals, solver, model)

                communicate.updateCountsODE(self.sim_properties, odeResults, model)

                self.next_metabolism_time = self.next_metabolism_time + 1.0

                print('ODE time: ', TIME.time()-odestart)

        
        if (self.complete_steps>=80) or (self.complete_steps==0):
            
            if time>0.99:
                
                filestart = TIME.time()

                # Skip extra updateCountsRDME here; counts were refreshed earlier in this hook.
                save.saveCountsAndFluxes(time, self.sim_properties, odeResults, model, solver)

                communicate.resetCostCounters(self.sim_properties)

                save.writeSimProp(self.sim_properties)

                save.writeLatticeFiles(self.sim_properties, lattice, self.region_dict)
                print('Communication and file write time: ', TIME.time()-filestart)
                
            print('Return 2 time: ', time)
            
            self.complete_steps = 1
            
            self.endLastHook = TIME.time()
            
            return 2
        
        else:
            print('Return 1 time: ', time)
            self.complete_steps = self.complete_steps + 1
            self.endLastHook = TIME.time()
            return 1


def round_sig(x, sig=2):
    """Round ``x`` to ``sig`` significant figures."""
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
