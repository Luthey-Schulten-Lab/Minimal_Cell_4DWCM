
import argparse
import os

#########################################################################################
ap = argparse.ArgumentParser()

ap.add_argument('-od', '--outputDir', required = True) # output directory name
ap.add_argument('-t', '--simTime', type=int, required= True) # Simulation Time: in seconds, e.g. 120

ap.add_argument('-cd', '--cudaDevices', type=int, default=0) # cuda device for simulation

ap.add_argument('-drs', '--dnaRngSeed', type=int, default=1) # RNG seed modifier for btree chromo

ap.add_argument("-dsd", "--dnaSoftwareDirectory", default='/home/zane/Software/BRGDNA/') # Directory containing btree_chromo/ used to generate chromosome configurations 

ap.add_argument("-m", "--membrane", type=int, default=1)

args = ap.parse_args()
#########################################################################################


#########################################################################################
hook_step = int(250) #steps
write_step = int(20000) #steps
totalTime = args.simTime #s
workingDirectoryName = args.outputDir
#########################################################################################


#########################################################################################
import jLM

from jLM.Solvers import makeSolver

import lm

from lm import IntMpdRdmeSolver
# from lm import MGPUIntMpdRdmeSolver

# import MC_RDME_initialization as initialization

import RegionsAndComplexes as InitGeom

import Restart_MC_RDME_initialization as MCRDME

import ImportInitialConditions as IC

import Communicate as communicate

import FileSaving as save
#########################################################################################


#########################################################################################
headDirectory =  os.getcwd() + '/'
workingDirectory = headDirectory + 'Data/' + workingDirectoryName + '/'
sim_properties_file = workingDirectory + 'sim_properties.pkl'
backup_sim_properties = workingDirectory + 'sim_properties_at_restart.pkl'
os.system("cp %s %s"% (sim_properties_file, backup_sim_properties))
#########################################################################################


# #########################################################################################
# sim, sim_properties = initialization.initSim(hook_step, write_step, totalTime, workingDirectoryName)

# save.saveSimArgs(sim_properties, args)

# sim_properties['dna_rng_seed'] = int(args.dnaRngSeed)

# sim_properties['dna_software_directory'] = str(args.dnaSoftwareDirectory)

# sim_properties['membrane_directory'] =sim_properties['head_directory'] + 'input_data/membrane/cell{:d}/'.format(args.membrane)
    
# sim_properties['division_started'] = False

# region_dict, ribo_site_dict = InitGeom.buildRegions(sim, sim_properties)

# IC.initializeParticles(sim, region_dict, sim_properties)

# MCRDME.createParticleIdxMap(sim, sim_properties)

# MCRDME.constructGIP(sim, sim_properties)

# MCRDME.createParticleIdxMap(sim, sim_properties)

# MCRDME.constructAssemblyReactions(sim, sim_properties)

# MCRDME.createParticleIdxMap(sim, sim_properties)

# MCRDME.mapTranslationStates(sim_properties)

# IC.initializeRdmeCounts(sim, sim_properties)

# communicate.updateSA(sim_properties)

# # save.saveParticleCounts(0, sim_properties)
# save.saveCountsAndFluxes(0, sim_properties, None, None, None)
# #########################################################################################


#########################################################################################
sim, sim_properties = MCRDME.initSimRestart(totalTime, sim_properties_file, workingDirectoryName)

sim_properties['dna_software_directory'] = str(args.dnaSoftwareDirectory)

sim_properties['membrane_directory'] =sim_properties['head_directory'] + 'input_data/membrane/cell{:d}/'.format(args.membrane)
    
region_dict, ribo_site_dict = MCRDME.buildRegionsRestart(sim, sim_properties)
    
MCRDME.initRdmeSpeciesRestart(sim, sim_properties, region_dict)

MCRDME.initDiffusionRestart(sim, sim_properties)

MCRDME.constructGIP(sim, sim_properties)

MCRDME.constructAssemblyReactions(sim, sim_properties)
#########################################################################################


#########################################################################################
import Restart_Hook as Hook

mc4dSolver = Hook.MyOwnSolver

Solver = makeSolver(IntMpdRdmeSolver, mc4dSolver)
solver = Solver(sim, sim_properties, region_dict, ribo_site_dict)

sim.finalize()

sim.run(solver=solver, cudaDevices=[int(args.cudaDevices)])
#########################################################################################

