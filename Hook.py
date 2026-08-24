"""
Authors: Zane Thornburg

Hook algorithm
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

        # Patch 2 (ODE Cython): prefer the Cython-compiled solver (~2-3x faster
        # per step than the no-Cython lsoda path).
        #
        # The historical IndexError (`stepClass.calcFlux_c` indexing
        # self.params[N] in a size-N array, verified 2026-05 in smoke 1641) was
        # NOT an odecell bug: Integrate.setSolver() was constructing the functor
        # with model.getInitVals() (the metabolite initial-concentration vector)
        # instead of the optimization-space parameter values. That sized
        # self.params to the metabolite count, so calcFlux_c ran off the end at
        # the first parameter index beyond it. Fixed in Integrate.setSolver()
        # (now uses getOptSpace() values captured before prepareFunctor()).
        # The try/except fallback below is retained as a safety net.
        #
        # Patch 2b (build-once ODE functor): odecell's buildCall(cythonBuild=True)
        # regenerates and *recompiles* the flux module from scratch on every call
        # (~95 s of Cython->C->gcc per metabolism step, dwarfing the millisecond
        # integration). It recompiles because the per-step enzyme concentrations
        # were baked into the generated code as literals. Now that "Enzyme" is in
        # the optimization space (Rxns_ODE.addParameter("Enzyme", ..., lb!=ub)),
        # every per-step-varying quantity reaches the functor through self.params,
        # so the compiled module is identical across steps. Integrate.setSolverCached
        # compiles ONCE into this cache and afterwards only re-instantiates the
        # compiled stepClass with the current parameter values (milliseconds),
        # recompiling solely when the opt-space / metabolite structure changes
        # (e.g. after replication). This turns ~95 s/step into ~95 s once.
        self._ode_use_cython = True
        self._ode_cython_failed = False
        self._ode_solver_cache = {}

        # Ribosome excluded-volume updates: every hook (12.5 ms bio = 80/s),
        # matching Thornburg et al. Stokes–Einstein cadence (Methods). Patch 6
        # had set this to 2 (every other hook / 25 ms); that exceeds the ~17 ms
        # diffusion bound and is disabled here. `translation_update_step = 8`
        # still gates polysome/translation-state updates inside placeRibosomes.
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
#             print(csim_folder)
            os.mkdir(csim_folder)
            print('Created global CME directory')
        except:
            print('CME sim directory already exists')
            
        try:
            flux_folder=sim_properties['working_directory']+'fluxes/'
#             print(csim_folder)
            os.mkdir(flux_folder)
            print('Created metabolic fluxes directory')
        except:
            print('Metabolic fluxes directory already exists')
            
        try:
            flux_folder=sim_properties['working_directory']+'restart_files/'
#             print(csim_folder)
            os.mkdir(flux_folder)
            print('Created restart files directory')
        except:
            print('Restart Files Directory directory already exists')
            

        print('Initialized Solver')

        return None
        
    # The hookSimulation method defined here will be called at every frame 
    # write time.  The return value is either 0, 1, or 2, which will indicate 
    # if we changed the state or not and need the lattice to be copied back 
    # to the GPU before continuing.  If you do not return 1 or 2, your changes 
    # will not be reflected.
    def hookSimulation(self, t, lattice):
        
        time = t
        
        print('Current biological time: ', time)
        
        self.sim_properties['time'] = time
        
        print('Time between hook steps: ', TIME.time()-self.endLastHook)
        
        ##### Execution of Growth #####
        
        if time == 0:
            
            communicate.updateCountsRDME(self.rdme, self.sim_properties, lattice)
            
            IC.setCmeSpeciesList(self.sim_properties)
        
#         if (time >= self.hookStartTime):
        if (time > 0):
        
            ##### For jobs running on some machines, there might be a maximum allowed time. #####
            ##### We terminate the simulation early to prevent recording errors. #####
            
            if self.termination_time is not None:
                
                current_datetime = datetime.now()
                
                elapsed = current_datetime-self.start_datetime
                elapsed_hours = elapsed.total_seconds()/3600
                
                if elapsed_hours>=self.termination_time:
                    
                    if self.final_simulation_time is None:
                        
                        self.final_simulation_time = self.next_DNA_time + 7.5
                
                    if (time >= self.final_simulation_time):
                        
                        raise Exception(f"Reached end of allowed simulation time after {elapsed_hours} hours. Solver will terminate without saving.")
            
            ##### Update Chromosome Configuration and Cell Architecture #####

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

                # Advance by dna_hook_interval_s (default 4.0; override via
                # DNA_HOOK_INTERVAL_SEC at init — see MC_RDME_initialization).
                _dna_hook_s = float(self.sim_properties.get('dna_hook_interval_s', 4.0))
                self.next_DNA_time = self.next_DNA_time + _dna_hook_s
                
                print('DNA time: ', TIME.time()-dnastart,
                      '(next DNA @ {:.1f}s, interval={:.1f}s)'.format(
                          self.next_DNA_time, _dna_hook_s))

                print('Updated cell architecture')


            ##### RDME Modification to Update Ribosome Positions and Polysome Structures #####

            ribostart = TIME.time()

            self.ribo_step = self.ribo_step + 1

            if self.ribo_step >= self.translation_update_step:

                updateTranslat = True

                self.ribo_step = 0

            else:

                updateTranslat = False

            # Place EV every ribo_place_every_n_hooks (default 1 = every hook).
            # Always place when updateTranslat fires so new polysomes are tracked.
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


            ##### Cell-wide (Global) CME for Well-Stirred Stochastic Reactions #####

            if (time >= self.next_gCME_time):
                
                cmestart = TIME.time()
                
                # Pre-CME Python operations
                prestart = TIME.time()
                communicate.updateCountsRDME(self.rdme, self.sim_properties, lattice)
                communicate.calculateTranslationCosts(self.sim_properties, lattice)
                communicate.updateLongGeneStates(self.sim_properties, lattice)
                pretime = TIME.time() - prestart
                
                # Actual CME simulation (Python setup + C++ solver)
                cmesolverstart = TIME.time()
                MCCME.runGCME(self.sim_properties)
                cmesolvertime = TIME.time() - cmesolverstart
                
                               # Post-CME Python operations
                poststart = TIME.time()
                
                # Time each post-CME operation individually
                countsCMEstart = TIME.time()
                communicate.updateCountsCME(self.sim_properties)
                countsCMEtime = TIME.time() - countsCMEstart
                
                transStatesstart = TIME.time()
                communicate.updateTranscriptionStates(self.sim_properties, lattice)
                transStatestime = TIME.time() - transStatesstart

                # Patch 5 (skip redundant 2nd updateCountsRDME): job 1638 spent
                # ~3.0 h re-scanning the full particle lattice immediately after
                # updateTranscriptionStates.  The very next pre-CME hook (1 s of
                # biological time later) re-scans anyway, so the only window of
                # potential staleness is < 1 s.  The lone consumer of
                # sim_properties['counts'] in that window is a possible DNA
                # update (every 4 s of bio); sim_properties['counts'] for
                # non-CME species is already correct from the pre-CME scan, and
                # for CME species it has just been refreshed by updateCountsCME
                # plus the on-lattice particle adjustments inside
                # updateTranscriptionStates.  Net: zero biological effect, ~3 h
                # saved.  Keep the timing variable for compatibility with the
                # printed breakdown below.
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



            ##### ODE Integrator for Metabolism #####

            if (time >= self.next_metabolism_time):

                odestart = TIME.time()

                communicate.calculateDegradationCosts(self.sim_properties, lattice)

                communicate.communicateCostsToMetabolism(self.sim_properties)

                print('Initializing ODE simulation')
                model = ODE.initModel(self.sim_properties)
                print('Initialized ODE simulation')

                ### Want to get the current values, not necessarily the initial values
                initVals=integrate.getInitVals(model)

                ### Patch 2: prefer the Cython-compiled solver (~2-3x faster
                ### per step than the no-Cython lsoda path).  On any failure
                ### we latch to the noCython path for the rest of the run.
                if self._ode_use_cython and not self._ode_cython_failed:
                    try:
                        # Build/compile once, reuse the compiled module thereafter.
                        solver = integrate.setSolverCached(model, self._ode_solver_cache)
                    except Exception as exc:
                        print(f"ODE: Cython solver build failed ({exc}); falling back to noCython for the remainder of the run.")
                        self._ode_cython_failed = True
                        self._ode_solver_cache = {}
                        solver = integrate.noCythonSetSolver(model)
                else:
                    solver = integrate.noCythonSetSolver(model)

                ### Run the integrator
                odeResults = integrate.runODE(initVals, solver, model)

                communicate.updateCountsODE(self.sim_properties, odeResults, model)

#                         save.saveODEfluxes(time, self.sim_properties, odeResults, model, solver)

                self.next_metabolism_time = self.next_metabolism_time + 1.0

                print('ODE time: ', TIME.time()-odestart)

        
        if (self.complete_steps>=80) or (self.complete_steps==0):
            
            if time>0.99:
                
                filestart = TIME.time()

                # Patch 8 (skip file-save updateCountsRDME): the file-save
                # block fires every 80 hooks, which lines up with gCME's 1 Hz
                # cadence in the same hook.  Pre-CME already scanned the lattice
                # ~50 ms earlier, updateCountsCME just refreshed CME species,
                # and updateCountsODE just refreshed metabolites.  The only
                # lattice mutation between the pre-CME scan and here is
                # updateTranscriptionStates adding/removing mRNA particles --
                # but mRNA species are CME species and were already refreshed
                # by updateCountsCME from the gCME .lm file.  Skipping this
                # third per-second lattice scan saves ~1.5 s x 7200 = ~3 h.
                # communicate.updateCountsRDME(self.rdme, self.sim_properties, lattice)
                
#                 communicate.calculateCosts(self.sim_properties, lattice)
                
#                 save.saveParticleCounts(time, self.sim_properties)
        
                save.saveCountsAndFluxes(time, self.sim_properties, odeResults, model, solver)
        
                communicate.resetCostCounters(self.sim_properties)
                
                save.writeSimProp(self.sim_properties)
                
                save.writeLatticeFiles(self.sim_properties, lattice, self.region_dict)
                print('Communication and file write time: ', TIME.time()-filestart)
                
#                 communicate.communicateCostsToMetabolism(self.sim_properties)
                
            print('Return 2 time: ', time)
            
            self.complete_steps = 1
            
            self.endLastHook = TIME.time()
            
            return 2
        
        else:
            print('Return 1 time: ', time)
            self.complete_steps = self.complete_steps + 1
            self.endLastHook = TIME.time()
            return 1
#             return 1
 
    
# Just a convenient function for rounding a value "x" 
# to a given number of significant figures "sigs"
def round_sig(x, sig=2):
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
