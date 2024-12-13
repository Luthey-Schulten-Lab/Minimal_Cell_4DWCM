"""
Authors: Zane Thornburg

Updates chromosome configuration from interaction with btree_chromo by Benjamin Gilbert
Includes DNA replication
"""

import numpy as np

import os

import subprocess

from subprocess import Popen, PIPE

import time as timepy

from Bio import SeqIO
from Bio.Seq import Seq

from LatticeFunctions import *

import GIP_rates as GIP

import FreeDTS_functions as fdf


#########################################################################################
def updateChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict, updateRegions):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    if sim_properties['last_DNA_step'] is None:
        
        DNA_lattice_coords = sim_properties['DNAcoords']
        
        genome = sim_properties['genome']
        
    elif sim_properties['last_DNA_step'] is not None:
    
        checkLastChromosome(sim_properties)
        
        if ('rotated_DNA' not in sim_properties.keys()) and sim_properties['division_started']:
            
            rotateChromosome(time, sim_properties)

#         deleteOldChromosome(lattice)
    
        region_dict, DNA_lattice_coords = placeNewChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict)

        genome = remapDNA(sim_properties)

        sim_properties['genome'] = genome
        
        moveDnaParticles(sim_properties, lattice, DNA_lattice_coords)
        
    if np.rint(time) < sim_properties['total_time']:

        writeRiboObstacleFile(time, ribo_site_dict, sim_properties)

    #     writeMembraneBoundaryFile(time, region_dict, sim_properties)

        writeChromosomeInputFile(time, sim_properties, updateRegions)

        runNewChromosome(time, sim_properties)
        
        print('Started running new chromosome in background')
        
    return region_dict, DNA_lattice_coords, genome
#########################################################################################


#########################################################################################
def updateChromosomeDivision(time, lattice, sim_properties, region_dict, ribo_site_dict):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    print('DNA Division Step')
    
    if sim_properties['last_DNA_step'] is None:
        
        DNA_lattice_coords = sim_properties['DNAcoords']
        
        genome = sim_properties['genome']
        
    elif sim_properties['last_DNA_step'] is not None:
    
        checkLastChromosome(sim_properties)
        
        if 'rotated_DNA' not in sim_properties.keys():
            
            rotateChromosome(time, sim_properties)

#         deleteOldChromosome(lattice)
    
        region_dict, DNA_lattice_coords = placeNewChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict)

        genome = remapDNA(sim_properties)

        sim_properties['genome'] = genome
        
        moveDnaParticles(sim_properties, lattice, DNA_lattice_coords)

        sim_properties['DNAcoords'] = DNA_lattice_coords
        
    if np.rint(time) < sim_properties['total_time']:

#         writeRiboObstacleFile(time, ribo_site_dict, sim_properties)

    #     writeMembraneBoundaryFile(time, region_dict, sim_properties)
    
#         bp_grow_steps = makeBdryCompressionFiles(time, sim_properties, slice_size=5)
        
#         DirectivesFname = writeBrdyGrowChromosomeInputFile(time, sim_properties, bp_grow_steps)

        writeDivisionChromosomeInputFile(time, sim_properties)

        runDivChromosome(time, sim_properties)
    
        timestep = int(time/sim_properties['timestep'])
    
        lastStep = sim_properties['last_DNA_step']

        sim_properties['last_last_DNA_step'] = lastStep

        sim_properties['last_DNA_step'] = timestep
    
        region_dict, DNA_lattice_coords = placeNewChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict)

        genome = remapDNA(sim_properties)

        sim_properties['genome'] = genome
        
        moveDnaParticles(sim_properties, lattice, DNA_lattice_coords)
        
        print('Updated chromosomes with dividing membrane')

    return region_dict, DNA_lattice_coords, genome
#########################################################################################


#########################################################################################
def checkLastChromosome(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    print("Checking for configuration from previous DNA step")
    
    workDir = sim_properties['working_directory']+'DNA/'
    
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
    print(DNAfile)
    
    last_DNA_complete = os.path.isfile(DNAfile)
    
    if not last_DNA_complete:
        
        print("Waiting on BRGDNA to complete configuration")
        
    DNA_wait = 0
    
    while not last_DNA_complete:
        
        if DNA_wait>=300:
            
            rescueDNA(sim_properties)
            
            print('Created Rescue Files')
            
            return None
        
        last_DNA_complete = os.path.isfile(DNAfile)
        
        timepy.sleep(10)
        
        DNA_wait = DNA_wait + 10
        
    print("Waited seconds: "+str(DNA_wait))
        
    print("Chromosome configuration ready to load")
        
    return None
#########################################################################################


#########################################################################################
def placeNewChromosome(time, lattice, sim_properties, region_dict, ribo_site_dict):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    timestep = int(time/sim_properties['timestep'])
    
    workDir = sim_properties['working_directory']+'DNA/'
    
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
    fileType = DNAfile.split('.')[1] #'bin'
    
    print(fileType)
    
    N_edges = sim_properties['lattice_edges']
    
    N_edges_x = N_edges[0] # Number of subvolumes making up and edge of the simulation space N x N x N
    N_edges_y = N_edges[1]
    N_edges_z = N_edges[2]

    N_2_x = N_edges_x/2
    N_2_y = N_edges_y/2
    N_2_z = N_edges_z/2
    
    DNAsites = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    DNA_lattice_coords = {}
    
    if fileType=='bin':
        
        with open(DNAfile,'rb') as f:
            
            DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
            
        DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
        
        print(DNAcoords.shape)
        
        for i in range(len(DNAcoords)):
#         for DNAparticle in DNAcoords:
            DNAparticle = DNAcoords[i]
            
            x = DNAparticle[0]
            y = DNAparticle[1]
            z = DNAparticle[2]
            
            x_lattice = (x*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_x
            y_lattice = (y*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_y
            z_lattice = (z*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_z
            
            DNAsites[int(x_lattice),int(y_lattice),int(z_lattice)] = True
            
            DNA_lattice_coords[i+1] = [int(x_lattice),int(y_lattice),int(z_lattice)]
            
#             lattice.addParticle(int(z_lattice),int(y_lattice),int(x_lattice),3)
            
        
    elif fileType=='xyz':
        
        with open(DNAfile,'rb') as f:
            
            for line_number,line in enumerate(f):
                
                if line_number == 0:
                    
                    DNAparticleCount=int(line)
                    
                elif line_number == 1:
                    
                    continue
                    
                else:
                    
                    particleIdx = int(line_number-2)
                    
                    atomic_symbol, x, y, z = line.split()
                    
                    x_lattice = (float(x)*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_x
                    y_lattice = (float(y)*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_y
                    z_lattice = (float(z)*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_z

                    DNAsites[int(x_lattice),int(y_lattice),int(z_lattice)] = True
                    
                    DNA_lattice_coords[particleIdx+1] = [int(x_lattice),int(y_lattice),int(z_lattice)]
        
#     sim_properties['DNAcoords'] = DNA_lattice_coords

    oldDNA = np.argwhere(region_dict['DNA']['shape']==True)
    
    for site in oldDNA:

        if region_dict['cytoplasm']['full_shape'][site[0], site[1], site[2]] == True:

            lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['cytoplasm']['index'])

            continue

        if region_dict['outer_cytoplasm']['full_shape'][site[0], site[1], site[2]] == True:

            lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['outer_cytoplasm']['index'])

            continue
            
        if region_dict['membrane']['shape'][site[0], site[1], site[2]] == True:

            lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['membrane']['index'])

            continue

        for ribo_type, type_dict in ribo_site_dict.items():
            
            if region_dict[type_dict['cross_idx']]['shape'][site[0], site[1], site[2]] == True:
                
                crossID = type_dict['cross_idx']
                
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict[crossID]['index'])
        
            if region_dict[type_dict['center_idx']]['shape'][site[0], site[1], site[2]] == True:
                
                centerID = type_dict['center_idx']

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict[centerID]['index'])
            
    
    region_dict['DNA']['shape'] = DNAsites
    
    region_dict['cytoplasm']['shape'] = region_dict['cytoplasm']['full_shape'] & ~region_dict["DNA"]["shape"]
    region_dict['outer_cytoplasm']['shape'] = region_dict['outer_cytoplasm']['full_shape'] & ~region_dict["DNA"]["shape"]
    
    newDNAsites = np.argwhere(region_dict['DNA']['shape']==True)
    
    for site in newDNAsites:
        
        lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['DNA']['index'])
    
    print("Updated DNA Sites")
    
    return region_dict, DNA_lattice_coords
#########################################################################################
    
    
#########################################################################################
def moveDnaParticles(sim_properties, lattice, DNA_lattice_coords):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
                        
    genome = sim_properties['genome']
    OldDnaParticleCoords = sim_properties['DNAcoords']
    
    plattice = lattice.getParticleLatticeView()

    for locusTag, locusDict in genome.items():
        
        locusNum = locusTag.split('_')[1]
        
#         trsc2ID = sim_properties['name_to_index']['RP_' + locusNum + '_t']
        
        starts = locusDict['startIndex']
        ends = locusDict['endIndex']
        
        old_starts = locusDict['prevStartIndex']
        old_ends = locusDict['prevEndIndex']
        
#         print(OldDnaParticleCoords)

        for i in range(len(starts)):
        
            chromo = i + 1
            
            geneID = sim_properties['name_to_index']['G_' + locusNum + '_C' + str(int(chromo))]
        
            start = starts[i]
            
            if (i+1)>len(old_starts):
                
                newStartXYZ = DNA_lattice_coords[start]

                lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), geneID)

                print('Placed new gene for ', locusTag)
                    
            else:
                
                old_start = old_starts[i]

                oldStartXYZ = OldDnaParticleCoords[old_start]

                newStartXYZ = DNA_lattice_coords[start]
                
                rnasequence = locusDict["RNAsequence"]
                
                if (locusDict['Type'] == 'rRNA') and (len(rnasequence) > 2*sim_properties['rnap_spacing']):
                    
                    if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID):

                        deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID)

                        lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), geneID)

                        continue

                    for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):

                        trscID = sim_properties['name_to_index']['RP_' + locusNum + '_' + str(i) + '_C' + str(int(chromo))]

                        if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID):

                            deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID)

                            lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), trscID)

                            break

                        trsc2ID = sim_properties['name_to_index']['RP_' + locusNum + '_t_' + str(i) + '_C' + str(int(chromo))]

                        if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID):

                            deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID)

                            lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), trsc2ID)

                            break
                    
                else:
                    
                    trscID = sim_properties['name_to_index']['RP_' + locusNum + '_C' + str(int(chromo))]

                    if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID):

                        deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID)

                        lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), geneID)
                        
                        continue

                    if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID):

                        deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID)

                        lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), trscID)
                        
                        continue
                    
#                 if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID):

#                     deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID)

#                     lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), trsc2ID)

#                 else:

#                     print('Something is wrong with gene: ', locusTag)
        
#     sim_properties['DNAcoords'] = DNA_lattice_coords

    for feature, fdict in sim_properties['chromosome_features'].items():
        
        oldPartXYZ = OldDnaParticleCoords[int(fdict['index'])]
        
        partIdx = sim_properties['name_to_index'][feature]
        
        if checkParticle(plattice, int(oldPartXYZ[2]), int(oldPartXYZ[1]), int(oldPartXYZ[0]), partIdx):
        
            newPartXYZ = DNA_lattice_coords[int(fdict['index'])]

            deleteParticle(plattice, int(oldPartXYZ[2]), int(oldPartXYZ[1]), int(oldPartXYZ[0]), partIdx)

            lattice.addParticle(int(newPartXYZ[2]), int(newPartXYZ[1]), int(newPartXYZ[0]), partIdx)
        
    print('Placed chemically relevant DNA particles')
        
    return None
#########################################################################################
    
    
#########################################################################################   
def runNewChromosome(time, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    timestep = int(time/sim_properties['timestep'])
    
#     os.system('LAMMPS_NEW_LOAD twistable_BD_OMP')
    
#     headDir = sim_properties['head_directory'] + 'external_programs/'
    
    headDir = sim_properties['dna_software_directory']
    workDir = sim_properties['working_directory']+'DNA/'
    
    DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)
    
#     DNA_executable = headDir + 'btree_chromo/build/apps/program ' + DirectivesFname
#     print(DNA_executable)
    
#     os.system(DNA_executable)

#     if sim_properties['last_DNA_step'] is not None:

#         LAMMPStraj= workDir + 'chromosome_{:d}.lammpstrj'.format(sim_properties['last_DNA_step'])

#         os.system('rm ' + LAMMPStraj)
    
    DNAargs = [headDir + 'btree_chromo/build/apps/btree_chromo', DirectivesFname]
    
#     Popen(DNAargs, stdin=None, stdout=None, stderr=None)

    Popen(DNAargs, stdin=None, stdout=subprocess.DEVNULL, stderr=None)
    
    lastStep = sim_properties['last_DNA_step']
    
    sim_properties['last_last_DNA_step'] = lastStep
    
    sim_properties['last_DNA_step'] = timestep
    
    print("Set next DNA configuration to run in background")
    
    return
#########################################################################################

    
#########################################################################################
def writeChromosomeInputFile(time, sim_properties, updateRegions):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
#     headDir = sim_properties['head_directory'] + 'external_programs/'
    headDir = sim_properties['dna_software_directory']
    workDir = sim_properties['working_directory']+'DNA/'
    
    rep_started = sim_properties['rep_started']
    
    cyto_radius_angstroms = int((sim_properties['cyto_radius'])*sim_properties['lattice_spacing']*10/1e-9)
    timestep = int(time/sim_properties['timestep'])
    
    rng_number = int(timestep/10000+sim_properties['dna_rng_seed'])
    
    processor_number = 8 #25
    
    loop_number = int(sim_properties['counts']['P_0415']/2)
    
    DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)
    
    DnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(timestep)
    
    try:
        PrevDnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    except:
        PrevDnaBinFname = workDir + 'x_chain_Syn3A_chromosome_init_rep00001.bin'
    
    RiboFname = workDir + 'ribo_obstacles_{:d}.bin'.format(timestep)
    
    MemBoundaryFname = workDir + 'mem_boundary_{:d}.bin'.format(timestep)
    
    DnaXyzFname = workDir + 'dna_monomers_{:d}.xyz'.format(timestep)
    
    DnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(timestep)
    
    try:
        PrevDnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(sim_properties['last_DNA_step'])
    except:
        PrevDnaQuatFname = None
    
    with open(DirectivesFname, 'w') as f:

        f.write('btree_prng_seed:10\n')
        f.write('replicator_prng_seed:10\n')

#         f.write('new_chromo:54338\n')
        
        if not rep_started:
            f.write('new_chromo:54338\n')
#             sim_properties['rep_started'] = True
        else:
            f.write('input_state:' + workDir + 'rep_state_{:d}.txt\n'.format(sim_properties['last_DNA_step']))

#         f.write('load_BD_lengths:' + headDir + 'btree_chromo/test_case/in_BD_lengths_LAMMPS_test.txt\n')
        f.write('load_BD_lengths:' + sim_properties['head_directory'] + 'input_data/in_BD_lengths_LAMMPS_test.txt\n')
        
        f.write('load_mono_coords:' + PrevDnaBinFname + ',row\n')
        
#         if PrevDnaQuatFname is not None:
#             f.write('load_mono_quats:' + PrevDnaQuatFname + ',row\n')
        
#         if not sim_properties['division_started']:
            
        if not sim_properties['division_started']:
            f.write('spherical_bdry:{:d},0,0,0\n'.format(cyto_radius_angstroms))
            
        if sim_properties['division_started']:
            cyto_radius_angstroms = int((sim_properties['divR']) * 10)
            cyto_height_angstroms = int((sim_properties['divH']) * 10)
            f.write('overlapping_spheres_bdry:{:d},{:d},0,0,0,0,0,1\n'.format(cyto_height_angstroms, cyto_radius_angstroms))

#         if (sim_properties['gamma_V'] == 1.0) and sim_properties['division_started']:
#             f.write('spherical_bdry:{:d},0,0,0\n'.format(cyto_radius_angstroms))
            
#         if (sim_properties['gamma_V'] < 1.0) and sim_properties['division_started']:
#             memFname = sim_properties['DNA_membrane_file']
#             f.write('load_bdry_coords:' + memFname + ',row\n')
        
#         f.write('load_bdry_coords:' + MemBoundaryFname + ',row\n')

#         f.write('write_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
            
        f.write('prepare_simulator:' + workDir + 'log_{:d}.log\n'.format(timestep))
        f.write('simulator_set_prng_seed:{:d}\n'.format(rng_number))
        f.write('simulator_set_nProc:{:d}\n'.format(processor_number))
        f.write('simulator_set_DNA_model:' + headDir + 'btree_chromo/LAMMPS_DNA_model_kk\n')
#         f.write('simulator_set_output_details:' + workDir + ',chromosome_{:d}\n'.format(timestep))
        f.write('simulator_set_output_details:' + workDir + ',chromosome\n'.format(timestep))

        f.write('simulator_set_delta_t:1.0E+5\n')
        
        # TURN OFF TWISTING #
        f.write('switch_twisting_angles:F\n')
        f.write('switch_ellipsoids:F\n')
        
#         f.write('simulator_load_loop_params:'+ headDir + 'btree_chromo/test_case/loop_params.txt\n')
#         f.write('simulator_load_loop_params:'+ headDir + 'loop_params.txt\n')
        f.write('simulator_load_loop_params:'+ sim_properties['head_directory'] + 'input_data/loop_params.txt\n')
        
#         f.write('switch_Ori_bdry_attraction:T\n')
#         f.write('switch_Ori_pair_repulsion:T\n')
        
        if rep_started:
            
            if sim_properties['last_DNA_step'] is None:
                repCw = 40
                repCcw = 40
            else:
                repCw, repCcw = getReplicatedSegments(sim_properties)
                
            f.write('set_initial_state\n')
            f.write('transform:m_cw' + str(repCw) + '_ccw' + str(repCcw) + '\n')
            f.write('set_final_state\n')
            f.write('map_replication\n')
            
        f.write('dump_topology:'+ workDir +'chromo_topo_{:d}.dat,1\n'.format(timestep))

        f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
        
        f.write('simulator_relax_progressive:1000,500\n')
        
        f.write('sync_simulator_and_system\n')
        
#         f.write('switch_fork_partition_repulsion:T\n')
        
#         f.write('sync_simulator_and_system\n')
        
        f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
        
        if sim_properties['last_DNA_step'] is None:
            f.write('switch_fork_partition_repulsion:T\n')
            f.write('simulator_run_loops:{:d},10000,10000,10000,noappend,first\n'.format(loop_number))
        else:
            if not checkDaughtersFullyPartitioned(sim_properties):
                print('DNA not fully partitioned yet; running loops with fork partition force ON')
                f.write('switch_fork_partition_repulsion:T\n')
            else:
                print('DNA is fully partitioned; running loops with fork partition force OFF')
            f.write('simulator_run_loops:{:d},10000,10000,10000,append,skip_first\n'.format(loop_number))
        
#         f.write('repeat:2\n')
#         f.write('simulator_run_loops:{:d},1000,100,100,append,nofirst\n'.format(loop_number))
#         f.write('end_repeat\n')
#         f.write('simulator_run_soft_FENE:20000,10000,10000,append,skip_first\n')
        f.write('sync_simulator_and_system\n')
    
        #if not updateRegions:
        #    f.write('load_ribo_coords:' + RiboFname + ',row\n')
        
        #f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
        
        #f.write('simulator_minimize_soft_harmonic:500\n')
        
        #f.write('simulator_run_topoDNA_FENE:100,50,50,append,skip_first\n')
    
        #f.write('sync_simulator_and_system\n')
    
        f.write('output_state:' + workDir + 'rep_state_{:d}.txt\n'.format(timestep))
        
#         f.write('write_mono_quats:' + DnaQuatFname + ',row\n')
        
#         f.write('write_mono_xyz:' + DnaXyzFname + '\n')
        
        f.write('write_mono_coords:' + DnaBinFname + ',row\n')
        
#         if rep_started:
        
#     sim_properties['last_DNA_step'] = timestep
    
    return None
#########################################################################################
        

#########################################################################################
def writeRiboObstacleFile(time, ribo_site_dict, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    N_edges = sim_properties['lattice_edges']
    
    ribo_center_lattice = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    for ribo_type, type_dict in ribo_site_dict.items():
        
        center_lattice = type_dict['centers']

        ribo_center_lattice = ribo_center_lattice | center_lattice
            
    ribo_center_points = np.argwhere(ribo_center_lattice == True)
    
    ribo_center_points = ribo_center_points - np.array(sim_properties['lattice_center'])
#     print(ribo_center_points)
            
    ribo_center_points = ribo_center_points * (sim_properties['lattice_spacing'] / 1e-9) * 10
    
#     print(ribo_center_points)

#     ribo_cytoplasm_space = int((sim_properties['cyto_radius'])*sim_properties['lattice_spacing']*10/1e-9) - 100

#     ribo_reduced = []

#     for coord in ribo_center_points:
        
#         COMradius = (coord[0]**2 + coord[1]**2 + coord[2]**2)**(1/2)
        
#         if COMradius<ribo_cytoplasm_space:
            
#             ribo_reduced.append(coord)
            
#     ribo_center_points = np.array(ribo_reduced)

    totalRiboCount = len(ribo_center_points)
    
    timestep = int(time/sim_properties['timestep'])
    
    fname = sim_properties['working_directory']+'DNA/' + 'ribo_obstacles_{:d}.bin'.format(timestep)
    
    ribo_center_points = np.reshape(ribo_center_points,(int(len(ribo_center_points)),3),order='F')
    
    with open(fname, 'wb') as f:
        
        ribo_center_points.tofile(f)
        
    return None
#########################################################################################


#########################################################################################
def remapDNA(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """

    genome = sim_properties['genome']
    
    chromosome_length = sim_properties['genome_length']
    
    ori_ter_rotation_factor = int(chromosome_length/2)
    
    workDir = sim_properties['working_directory']+'DNA/'
    RepTopoFname = workDir +'chromo_topo_{:d}.dat'.format(sim_properties['last_DNA_step'])
    
    repTopoFile = open(RepTopoFname, 'r')
    lines = repTopoFile.readlines()
    
    total_DNA = int(lines[0].split('\n')[0].split('=')[1])
    sim_properties['counts']['chromosome'] = int(total_DNA)
    
    if len(lines)>2:
        
        print('Remapping replicated DNA')
    
        daughterIdxs = lines[2].split('\n')[0].split(',')
        newOriIdx = int(daughterIdxs[3])
        leftArmPoint = int(daughterIdxs[2])
        rightArmPoint = int(daughterIdxs[4])
        
        replicationGoing = True
        
    else:
        
        replicationGoing = False
        
    repTopoFile.close()

    for locusTag, locusDict in genome.items():

        original_start = locusDict['originalStart']
        original_end = locusDict['originalEnd']

        old_starts = locusDict['startIndex']
        old_ends = locusDict['endIndex']

        new_starts = [original_start]
        new_ends = [original_end]
        
        if replicationGoing:

    #         if original_start <= ori_ter_rotation_factor:

            newStart = newOriIdx + original_start - ori_ter_rotation_factor

    #         elif original_start > ori_ter_rotation_factor:

    #             newStart = newOriIdx - original_start

            if leftArmPoint <= newStart <= rightArmPoint:

                new_starts.append(newStart)

    #         if original_end <= ori_ter_rotation_factor:

            newEnd = newOriIdx + original_end - ori_ter_rotation_factor

    #         elif original_end > ori_ter_rotation_factor:

    #             newEnd = newOriIdx - original_end

            if leftArmPoint <= newEnd <= rightArmPoint:

                new_ends.append(newEnd)

        genome[locusTag]['startIndex'] = new_starts
        genome[locusTag]['endIndex'] = new_ends

        genome[locusTag]['prevStartIndex'] = old_starts
        genome[locusTag]['prevEndIndex'] = old_ends
        
#     print(genome)
    
    return genome
#########################################################################################


#########################################################################################
def getReplicatedSegments(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    genomeFile3A = sim_properties['head_directory'] + 'input_data/syn3A.gb'
    genome3A = next(SeqIO.parse(genomeFile3A, "gb"))
    
    genome = sim_properties['genome']
    
    chromosome_length = sim_properties['genome_length']
    
    ori_ter_rotation_factor = int(chromosome_length/2)
    
    workDir = sim_properties['working_directory']+'DNA/'
    RepTopoFname = workDir +'chromo_topo_{:d}.dat'.format(sim_properties['last_DNA_step'])
    
    repTopoFile = open(RepTopoFname, 'r')
    lines = repTopoFile.readlines()
    
    if len(lines) > 2:
        
        originalIdxs = lines[1].split('\n')[0].split(',')
        daughterIdxs = lines[2].split('\n')[0].split(',')
        originalOri = int(int(originalIdxs[3])*10)
        ccwPoint = int(int(daughterIdxs[1])*10)
        cwPoint = int(int(daughterIdxs[5])*10)
        
        repTopoFile.close()
        
    else:
        
        repTopoFile.close()
        
        repCcw = 40
        repCw = 40
        
        return repCw, repCcw
    
    cwPointIdx = cwPoint - originalOri
    
#     print('CW Idx: ', cwPointIdx)
    
    dnasequenceCW = str(genome3A.seq[cwPointIdx+1:cwPointIdx+1+400])
    
#     print(dnasequenceCW)
#     print(len(dnasequenceCW))
    
    repCw = int(GIP.ReplicationRate(sim_properties, dnasequenceCW))
    
    ccwPointIdx = ccwPoint + originalOri
    
#     print('CCW Idx: ', ccwPointIdx)
    
    dnasequenceCCW = str(genome3A.seq[ccwPointIdx-400+1:ccwPointIdx+1])
    
#     print(dnasequenceCCW)
#     print(len(dnasequenceCCW))
    
    repCcw = int(GIP.ReplicationRate(sim_properties, dnasequenceCCW))
    
    if (repCw == 0) and (repCcw == 0):
        
        return repCw, repCcw
    
    dnaCwReplicated = dnasequenceCW[:(repCw*10)]
    
    dnaCcwReplicated = dnasequenceCCW[:(repCcw*10)]
    
    dNTPcostMap = {'A':['dATP_DNArep_cost', 'dTTP_DNArep_cost'], 'T':['dTTP_DNArep_cost', 'dATP_DNArep_cost'], 
                   'C':['dCTP_DNArep_cost', 'dGTP_DNArep_cost'], 'G':['dGTP_DNArep_cost', 'dCTP_DNArep_cost']}
    
    for base in set(dnaCwReplicated):
        
        baseCount = dnaCwReplicated.count(base)
        
        costIDs = dNTPcostMap[base]
        
        for costID in costIDs:
            
            sim_properties['counts'][costID] = sim_properties['counts'][costID] + baseCount
        
            sim_properties['counts'][costID+'_second'] = sim_properties['counts'][costID+'_second'] + baseCount
        
    sim_properties['counts']['ATP_DNArep_cost'] = sim_properties['counts']['ATP_DNArep_cost'] + int(len(dnaCwReplicated))
    
    sim_properties['counts']['ATP_DNArep_cost_second'] = sim_properties['counts']['ATP_DNArep_cost_second'] + int(len(dnaCwReplicated))
        
    for base in set(dnaCcwReplicated):
        
        baseCount = dnaCcwReplicated.count(base)
        
        costIDs = dNTPcostMap[base]
        
        for costID in costIDs:
            
            sim_properties['counts'][costID] = sim_properties['counts'][costID] + baseCount
        
            sim_properties['counts'][costID+'_second'] = sim_properties['counts'][costID+'_second'] + baseCount
        
    sim_properties['counts']['ATP_DNArep_cost_second'] = sim_properties['counts']['ATP_DNArep_cost_second'] + int(len(dnaCcwReplicated))
    
    return repCw, repCcw
#########################################################################################


#########################################################################################
def writeDivisionChromosomeInputFile(time, sim_properties):

    headDir = sim_properties['dna_software_directory']
    workDir = sim_properties['working_directory']+'DNA/'
    
    rep_started = sim_properties['rep_started']
    
    cyto_radius_angstroms = int((sim_properties['cyto_radius'])*sim_properties['lattice_spacing']*10/1e-9)
    timestep = int(time/sim_properties['timestep'])
    
    rng_number = int(timestep/10000+sim_properties['dna_rng_seed'])
    
    processor_number = 8 #25
    
    loop_number = int(sim_properties['counts']['P_0415']/2)
    
    DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)
    
    DnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(timestep)
    
    try:
        PrevDnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    except:
        PrevDnaBinFname = workDir + 'x_chain_Syn3A_chromosome_init_rep00001.bin'
    
    RiboFname = workDir + 'ribo_obstacles_{:d}.bin'.format(timestep)
    
    MemBoundaryFname = workDir + 'mem_boundary_{:d}.bin'.format(timestep)
    
    DnaXyzFname = workDir + 'dna_monomers_{:d}.xyz'.format(timestep)
    
    DnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(timestep)
    
    try:
        PrevDnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(sim_properties['last_DNA_step'])
    except:
        PrevDnaQuatFname = None
    
    with open(DirectivesFname, 'w') as f:

        f.write('btree_prng_seed:10\n')
        f.write('replicator_prng_seed:10\n')
        
        if not rep_started:
            f.write('new_chromo:54338\n')
        else:
            f.write('input_state:' + workDir + 'rep_state_{:d}.txt\n'.format(sim_properties['last_DNA_step']))

        f.write('load_BD_lengths:' + sim_properties['head_directory'] + 'input_data/in_BD_lengths_LAMMPS_test.txt\n')
        
        f.write('load_mono_coords:' + PrevDnaBinFname + ',row\n')
            
        # Create the division cell shape
        f.write('overlapping_spheres_bdry:{:d},{:d},0,0,0,0,0,1\n'.format(int(sim_properties['divH_Prev']*10), int(sim_properties['divR_Prev']*10)))
            
        # Set simulator parameters and paths
        f.write('prepare_simulator:' + workDir + 'log_{:d}.log\n'.format(timestep))
        f.write('simulator_set_prng_seed:{:d}\n'.format(rng_number))
        f.write('simulator_set_nProc:{:d}\n'.format(processor_number))
        f.write('simulator_set_DNA_model:' + headDir + 'btree_chromo/LAMMPS_DNA_model_kk\n')

        f.write('simulator_set_output_details:' + workDir + ',chromosome\n'.format(timestep))

        f.write('simulator_set_delta_t:1.0E+5\n')
        
        # TURN OFF TWISTING #
        f.write('switch_twisting_angles:F\n')
        f.write('switch_ellipsoids:F\n')
        
        f.write('simulator_load_loop_params:'+ sim_properties['head_directory'] + 'input_data/loop_params.txt\n')
            
        f.write('dump_topology:'+ workDir +'chromo_topo_{:d}.dat,1\n'.format(timestep))

        f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
        
        # Minimize in previous division shape
        f.write('simulator_relax_progressive:1000,500\n')
        
        f.write('simulator_run_soft_FENE:100,50,50,append,skip_first\n')
        
        f.write('sync_simulator_and_system\n')
        
        # Calculate change in membrane shape
        bead_size = 3.4
        
        divHdiff = sim_properties['divH'] - sim_properties['divH_Prev']
        divRdiff = sim_properties['divR'] - sim_properties['divR_Prev']
        
        #Progressively change membrane shape to new geometry
        if divHdiff>bead_size or divRdiff>bead_size:
            
            Hsteps = int(divHdiff/bead_size+1)
            Rsteps = int(divRdiff/bead_size+1)
            
            div_steps = max(Hsteps,Rsteps)
            
            Hdelt = divHdiff/div_steps
            Rdelt = divRdiff/div_steps
            
            for i in range(int(div_steps)):
            
                divHA = int((sim_properties['divH_Prev'] + (i+1)*Hdelt)*10)
                divRA = int((sim_properties['divR_Prev'] + (i+1)*Rdelt)*10)

                f.write('overlapping_spheres_bdry:{:d},{:d},0,0,0,0,0,1\n'.format(divHA,divRA))
                f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
                f.write('simulator_relax_progressive:1000,500\n')
                f.write('simulator_run_soft_FENE:100,50,50,append,skip_first\n')
                f.write('sync_simulator_and_system\n')
            
        
        # Minimize in the new division shape
        divHA = int(sim_properties['divH']*10)
        divRA = int(sim_properties['divR']*10)

        f.write('overlapping_spheres_bdry:{:d},{:d},0,0,0,0,0,1\n'.format(divHA,divRA))
        f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
        f.write('simulator_relax_progressive:1000,500\n')
        f.write('simulator_run_soft_FENE:100,50,50,append,skip_first\n')
        f.write('sync_simulator_and_system\n')
        
        # Run looping
        f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))

        if not checkDaughtersFullyPartitioned(sim_properties):
            f.write('switch_fork_partition_repulsion:T\n')
        f.write('simulator_run_loops:{:d},10000,10000,10000,append,skip_first\n'.format(loop_number))
        
        f.write('sync_simulator_and_system\n')
        
        # Run some BD steps so that the minimized chromosome state is recorded
        f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
        
        f.write('simulator_minimize_soft_harmonic:500\n')
        f.write('simulator_run_topoDNA_FENE:1000,500,500,append,skip_first\n')
    
        f.write('sync_simulator_and_system\n')
    
        # Write out the monomer coordinates of the configuration in the new membrane shape
        f.write('output_state:' + workDir + 'rep_state_{:d}.txt\n'.format(timestep))
        
        f.write('write_mono_coords:' + DnaBinFname + ',row\n')
    
    return None
#########################################################################################


#########################################################################################
def runDivChromosome(time, sim_properties):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    timestep = int(time/sim_properties['timestep'])
    
#     timestep = sim_properties['timestep']
    
#     os.system('LAMMPS_NEW_LOAD twistable_BD_OMP')
    
    headDir = sim_properties['dna_software_directory']
    workDir = sim_properties['working_directory']+'DNA/'
    
#     DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)

    DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)
    
    DNA_executable = headDir + 'btree_chromo/build/apps/btree_chromo ' + DirectivesFname
    print(DNA_executable)
    
    os.system(DNA_executable)
    
    return None
#########################################################################################


#########################################################################################
def rotateChromosome(time, sim_properties):
    
    workDir = sim_properties['working_directory']+'DNA/'
    
#     timestep = int(time/sim_properties['timestep'])
    
#     DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(timestep)
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
    # Read in monomer coordinates
    with open(DNAfile,'rb') as f:
            
        DNAbin = np.fromfile(f,dtype=np.float64,count=-1)

    DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
    
    chromosome_length = 54338

    # Read replication state topology
    RepTopoFname = workDir +'chromo_topo_{:d}.dat'.format(sim_properties['last_DNA_step'])

    repTopoFile = open(RepTopoFname, 'r')
    lines = repTopoFile.readlines()

    if len(lines)>2:

        daughterIdxs = lines[2].split('\n')[0].split(',')
        newOriIdx = int(daughterIdxs[3])
        leftArmPoint = int(daughterIdxs[2])
        rightArmPoint = int(daughterIdxs[4])
        
    # Separate the coordinates for the mother and daughter chromosomes
    mother_monomers = DNAcoords[:54338]
    daughter_monomers = DNAcoords[leftArmPoint-1:rightArmPoint]
    
    # Calculate center of mass of both chromosomes
    mCom = np.average(mother_monomers, axis=0)
    dCom = np.average(daughter_monomers, axis=0)
    
    # Calculate vector between centers of mass
    if mCom[2] > dCom[2]:
        n = mCom - dCom
    else:
        n = dCom - mCom

    # Calculate rotation matrix between center of mass vector and the z axis
    n_z = np.array([0,0,1])

    n_z = n_z/np.linalg.norm(n_z)

    n_norm = n/np.linalg.norm(n)

    axis = np.cross(n_z,n_norm)

    axis = axis/np.linalg.norm(axis)

    theta = np.arccos(np.dot(n_z,n_norm))

    a = np.cos(theta/2)
    b, c, d = axis * np.sin(theta/2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

    rotation_matrix = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    
    sim_properties['rotation_matrix'] = rotation_matrix
    
    # Rotate every monomer for both chromosomes using the rotation matrix
    monos_R = []

    for mono in DNAcoords:

        rotated_mono = np.dot(rotation_matrix,mono)

        monos_R.append(rotated_mono)

    monos_R = np.array(monos_R)
    
    # Update the monomer coordinates file with the rotated chromosomes' state
    DNAfilePre = workDir + 'dna_monomers_{:d}_prerotation.bin'.format(sim_properties['last_DNA_step'])
    DNAfilePost = DNAfile # workDir +  'dna_monomers_{:d}_postrotation.bin'.format(timestep)

    os.system('mv ' + DNAfile + ' ' + DNAfilePre)

    monos_R = np.reshape(monos_R,(int(len(monos_R)),3),order='F')

    with open(DNAfilePost, 'wb') as f:

        monos_R.tofile(f)
        
    sim_properties['rotated_DNA'] = True
        
    return None
        
#########################################################################################


#########################################################################################
def checkDaughtersFullyPartitioned(sim_properties):
    print('Checking DNA Spatial Partitioning Status')

    workDir = sim_properties['working_directory'] + 'DNA/'
    DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])

    # Read in monomer coordinates
    with open(DNAfile, 'rb') as f:
        DNAbin = np.fromfile(f, dtype=np.float64, count=-1)

    DNAcoords = DNAbin.reshape((3, DNAbin.shape[0] // 3), order='F').T

    # Read replication state topology
    RepTopoFname = workDir + 'chromo_topo_{:d}.dat'.format(sim_properties['last_DNA_step'])
    with open(RepTopoFname, 'r') as repTopoFile:
        lines = repTopoFile.readlines()

    if len(lines) <= 2:
        return False

    daughterIdxs = lines[2].split('\n')[0].split(',')
    leftArmPoint = int(daughterIdxs[2])
    rightArmPoint = int(daughterIdxs[4])

    # Separate coordinates for the two daughter chromosomes
    daughterA = DNAcoords[:leftArmPoint-1]
    daughterB = DNAcoords[leftArmPoint-1:rightArmPoint]

    # Check z-coordinates for partitioning
    daughterA_z = daughterA[:, 2]
    daughterB_z = daughterB[:, 2]

    condition1 = np.all(daughterA_z < 0) and np.all(daughterB_z > 0)
    condition2 = np.all(daughterA_z > 0) and np.all(daughterB_z < 0)

    return condition1 or condition2
#########################################################################################


#########################################################################################
def rescueDNA(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    print('Rescuing Chromosome from Trapped State')
    
    workDir = sim_properties['working_directory']+'DNA/'
    
    oldDNAFile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_last_DNA_step'])
    
    rescueDNAFile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
    os.system('cp ' + oldDNAFile + ' ' + rescueDNAFile)
    
    oldQuatFile = workDir + 'dna_quats_{:d}.bin'.format(sim_properties['last_last_DNA_step'])
    
    rescueQuatFile = workDir + 'dna_quats_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
    os.system('cp ' + oldQuatFile + ' ' + rescueQuatFile)
    
    oldRepState = workDir + 'rep_state_{:d}.txt'.format(sim_properties['last_last_DNA_step'])
                                                          
    rescueRepState = workDir + 'rep_state_{:d}.txt'.format(sim_properties['last_DNA_step'])
                                                          
    os.system('cp ' + oldRepState + ' ' + rescueRepState)
    
    oldChromoTopo = workDir + 'chromo_topo_{:d}.dat'.format(sim_properties['last_last_DNA_step'])
    
    rescueChromoTopo = workDir + 'chromo_topo_{:d}.dat'.format(sim_properties['last_DNA_step'])
    
    os.system('cp ' + oldChromoTopo + ' ' + rescueChromoTopo)
    
    return None
#########################################################################################


#########################################################################################
# def NEWcheckLastChromosome(sim_properties):
    
#     print("Checking for configuration from previous DNA step")
    
#     workDir = sim_properties['working_directory']+'DNA/'
#     DnaLogFile = workDir + 'log_{:d}.log'.format(sim_properties['last_DNA_step'])
    
#     DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
#     last_DNA_complete = os.path.isfile(DNAfile)
    
#     if not last_DNA_complete:
        
#         print("Waiting on BRGDNA to complete configuration")
        
#     DNA_wait = 0
    
#     while not last_DNA_complete:
        
#         if DNA_wait>=1800:
            
#             rescueDNA(sim_properties)
            
#             print('Created Rescue Files')
            
#             return None
        
#         last_DNA_complete = os.path.isfile(DNAfile)
        
#         timepy.sleep(10)
        
#         DNA_wait = DNA_wait + 10
        
#     print("Waited seconds: "+str(DNA_wait))
        
#     print("Chromosome configuration ready to load")
        
#     return None
#########################################################################################

# #########################################################################################
# def writeBrdyGrowChromosomeInputFile(time, sim_properties, bp_grow_steps):
#     """
#     Inputs:
#     sim_properties - Dictionary of simulation variables and state trackers
    
#     Returns:
#     Called by:
#     Description:
#     """
    
# #     timestep = int(time/sim_properties['timestep'])
    
#     headDir = sim_properties['dna_software_directory']
#     workDir = sim_properties['working_directory']+'DNA/'
    
#     timestep = int(time/sim_properties['timestep'])
    
#     rng_number = int(timestep/10000+sim_properties['dna_rng_seed'])
    
#     processor_number = 8 #25
    
#     loop_number = int(sim_properties['counts']['P_0415']/2)
    
#     DirectivesFname = workDir + 'chromosome_operations_{:d}.inp'.format(timestep)
    
#     DnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(timestep)
    
#     try:
#         PrevDnaBinFname = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
#     except:
#         PrevDnaBinFname = workDir + 'x_chain_Syn3A_chromosome_init_rep00001.bin'
        
#     try:
#         PrevDnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(sim_properties['last_DNA_step'])
#     except:
#         PrevDnaQuatFname = None
    
# #     RiboFname = workDir + 'ribo_obstacles_{:d}.bin'.format(timestep)
    
#     MemBoundaryFname = workDir + 'mem_boundary_{:d}.bin'.format(timestep)
    
#     DnaXyzFname = workDir + 'dna_monomers_{:d}.xyz'.format(timestep)
    
#     DnaQuatFname = workDir + 'dna_quats_{:d}.bin'.format(timestep)
    
#     with open(DirectivesFname, 'w') as f:

#         f.write('btree_prng_seed:10\n')
#         f.write('replicator_prng_seed:10\n')

# #         f.write('new_chromo:5000\n')
#         f.write('input_state:' + workDir + 'rep_state_{:d}.txt\n'.format(sim_properties['last_DNA_step']))

#         f.write('load_BD_lengths:' + sim_properties['head_directory'] + 'input_data/in_BD_lengths_LAMMPS_test.txt\n')
# #         f.write('load_BD_lengths:/home/zane/in_BD_lengths_LAMMPS_test.txt\n')
        
#         f.write('load_mono_coords:' + PrevDnaBinFname + ',row\n')
        
#         f.write('load_mono_quats:' + PrevDnaQuatFname + ',row\n')
        
#         memFname = workDir + 'mem_boundary_{:d}_1.bin'.format(timestep)
            
#         f.write('load_bdry_coords:' + memFname + ',row\n')

# #         f.write('spherical_bdry:{:d},0,0,0\n'.format(cyto_radius_angstroms))
            
#         f.write('prepare_simulator:' + workDir + 'log_{:d}.log\n'.format(timestep))
#         f.write('simulator_set_prng_seed:{:d}\n'.format(rng_number))
#         f.write('simulator_set_nProc:{:d}\n'.format(processor_number))
#         f.write('simulator_set_DNA_model:' + headDir + 'btree_chromo/LAMMPS_DNA_model\n')
# #         f.write('simulator_set_output_details:' + workDir + ',chromosome_{:d}\n'.format(timestep))
#         f.write('simulator_set_output_details:' + workDir + ',chromosome\n'.format(timestep))

#         f.write('simulator_set_delta_t:1.0E+5\n')
        
# #         f.write('simulator_load_loop_params:'+ headDir + 'btree_chromo/test_case/loop_params.txt\n')
# #         f.write('simulator_load_loop_params:'+ headDir + 'loop_params.txt\n')
#         f.write('simulator_load_loop_params:'+ sim_properties['head_directory'] + 'input_data/loop_params.txt\n')
        
# #         f.write('switch_Ori_bdry_attraction:T\n')
# #         f.write('switch_Ori_pair_repulsion:T\n')

#         f.write('dump_topology:'+ workDir +'chromo_topo_{:d}.dat,1\n'.format(timestep))

#         f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
        
#         f.write('simulator_relax_progressive:1000,500\n')
        
#         f.write('sync_simulator_and_system\n')
        
#         for i in range(1,bp_grow_steps+1):
            
#             memFname = workDir + 'mem_boundary_{:d}_{:d}.bin'.format(timestep,i)
            
#             f.write('load_bdry_coords:' + memFname + ',row\n')
            
#             f.write('sys_write_sim_read_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))
# #             sys_write_sim_read_LAMMPS_data_at_timestep:/home/ben/Data/btree_chromo/bdry_particle_expansion/data.5000mono_50ribo_bdry900A
#             f.write('simulator_expand_bdry_particles:2.0,10,100,10\n')
#             f.write('sync_simulator_and_system\n')
#             f.write('simulator_run_hard_FENE:10,5,10,append,skip_first\n')
#             f.write('sync_simulator_and_system\n')
        
# #         f.write('write_LAMMPS_data:' + workDir + 'data.lammps_{:d}\n'.format(timestep))

# #         f.write('switch_Ori_pair_repulsion:T\n')

# #         loop_number=20
        
# #         f.write('simulator_run_loops:{:d},1000,100,100,append,nofirst\n'.format(loop_number))
        
#         f.write('sync_simulator_and_system\n')

# #         f.write('write_mono_xyz:' + DnaXyzFname + '\n')
        
#         f.write('write_mono_coords:' + DnaBinFname + ',row\n')
        
#         f.write('write_mono_quats:' + DnaQuatFname + ',row\n')
        
#         f.write('output_state:' + workDir + 'rep_state_{:d}.txt\n'.format(timestep))
        
# #         f.write('write_ribo_coords:' + RiboBinFname + ',row\n')
    
#     return DirectivesFname
# #########################################################################################


# #########################################################################################
# def makeBdryCompressionFiles(time, sim_properties, slice_size=10):
#     """
#     Inputs:
#     sim_properties - Dictionary of simulation variables and state trackers
    
#     Returns:
#     Called by:
#     Description:
#     """
    
#     print('NEW MEMBRANE SHAPE')
    
#     timestep = int(time/sim_properties['timestep'])
    
#     workDir = sim_properties['working_directory']+'DNA/'
    
#     DNAfile = workDir + 'dna_monomers_{:d}.bin'.format(sim_properties['last_DNA_step'])
    
#     with open(DNAfile,'rb') as f:

#         DNAbin = np.fromfile(f,dtype=np.float64,count=-1)

#         DNAcoords = np.divide(DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T,10)

#     zmax = 0
#     zmin = 0

#     for coord in DNAcoords:

#         if coord[2]>zmax:
#             zmax = float(coord[2])

#         if coord[2]<zmin:
#             zmin = float(coord[2])

#         xy_radius = (coord[0]**2 + coord[1]**2)**(1/2)

# #         if xy_radius > DNAxyMax:

# #             DNAxyMax = float(xy_radius)

#     print(zmin,zmax)
    
#     gamma_V = str(round_sig(sim_properties['gamma_V'], sig=2))

#     memFname = sim_properties['membrane_directory'] + 'output{:s}_mash4/InnerBM.dat'.format(gamma_V)
#     memcoords = fdf.readTS2CG(memFname, rescaleFactor=1)
#     print(len(memcoords))

#     zlow = float(zmin)
#     zhigh = zlow + slice_size
#     initial_scaling = 0
    
#     while (zhigh<(zmax+slice_size)):
        
#         DNAxyMax = 0
        
#         for coord in DNAcoords:
            
#             if zlow<=coord[2]<=zhigh:

#                 xy_radius = (coord[0]**2 + coord[1]**2)**(0.5)

#                 if xy_radius > DNAxyMax:

#                     DNAxyMax = float(xy_radius)
                    
                    
#         DNAxyMax = DNAxyMax+5+1.7 #+4.25 #Account for size of boundary and DNA particles

#         min_rxy = float(DNAxyMax)
        
#         for coord in memcoords:
            
#             if zlow<=coord[2]<=zhigh:
                
#                 xy_radius = (coord[0]**2 + coord[1]**2)**(0.5)

#                 if xy_radius<min_rxy:

#                     min_rxy = float(xy_radius)
                    
#         layer_scaling_factor = DNAxyMax/min_rxy
                    
#         initial_scaling = max(initial_scaling, layer_scaling_factor)
        
#         zlow = float(zhigh)
#         zhigh = zhigh + slice_size

                    
#     max_r = 0            
#     for coord in memcoords:
# #         if zlow<coord[2]<zhigh:
#         xyz_radius = (coord[0]**2 + coord[1]**2 + coord[2]**2)**(0.5)
#         if xyz_radius > max_r:
#             max_r = float(xyz_radius)
                    
#     print(max_r)

# #     initial_scaling = DNAxyMax/min_rxy
#     print(initial_scaling)
#     bp_increase = 5 #4.25 #nm
#     max_r_scaled = initial_scaling*max_r
#     max_r_diff = max_r_scaled - max_r
#     print(max_r_diff)
#     bp_grow_steps = int(max_r_diff/bp_increase)
#     scaling = float(initial_scaling)
#     curr_max_r = float(max_r_scaled)
#     for i in range(1,bp_grow_steps+1):
#     #     makeMembrane()
#         writeMembraneBoundaryFileBdryGrow(memcoords, workDir, timestep, rescale_factor=scaling, rescale_count=i)
#         if i < bp_grow_steps:
#             curr_max_r = curr_max_r - bp_increase
#             scaling = curr_max_r/max_r
        
#     memFname = sim_properties['membrane_directory'] + 'output{:s}_mash3/InnerBM.dat'.format(gamma_V)
#     memcoords = fdf.readTS2CG(memFname, rescaleFactor=1)
#     print(len(memcoords))
#     lowResFname = writeMembraneBoundaryFileBdryGrow(memcoords, workDir, timestep, rescale_factor=scaling, rescale_count=0)
    
#     sim_properties['DNA_membrane_file'] = lowResFname
        
# #     print(curr_max_r)
# #     print(scaling)

#     return bp_grow_steps
# #########################################################################################


# #########################################################################################
# def writeMembraneBoundaryFileBdryGrow(memBoundaryCoords, working_directory, timestep, rescale_factor=1, rescale_count=1):
#     """
#     Inputs:
    
#     Returns:
#     Called by:
#     Description:
#     """
    
#     workDir = working_directory

#     # print(memBoundaryCoords)
#     memBoundaryCoords = np.array(memBoundaryCoords)*10*rescale_factor
# #     print(memBoundaryCoords.shape)

#     memBoundaryCoords = np.reshape(memBoundaryCoords,(int(len(memBoundaryCoords)),3),order='F')
# #     print(memBoundaryCoords.shape)
    
# #     print(memBoundaryCoords)

#     fname = workDir + 'mem_boundary_{:d}_{:d}.bin'.format(timestep,rescale_count)

#     with open(fname, 'wb') as f:

#         memBoundaryCoords.tofile(f)
        
#     return fname
# #########################################################################################


# OldMoveDnaParticlesrRNA
# if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID):

#     deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), geneID)

#     lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), geneID)

#     continue

# for i in range(1, sim_properties['long_rna_trsc'][locusTag]['max_rnap']+1):

#     trscID = sim_properties['name_to_index']['RP_' + locusNum + '_' + str(i)]

#     if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID):

#         deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trscID)

#         lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), trscID)

#         break

#     trsc2ID = sim_properties['name_to_index']['RP_' + locusNum + '_t_' + str(i)]

#     if checkParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID):

#         deleteParticle(plattice, int(oldStartXYZ[2]), int(oldStartXYZ[1]), int(oldStartXYZ[0]), trsc2ID)

#         lattice.addParticle(int(newStartXYZ[2]), int(newStartXYZ[1]), int(newStartXYZ[0]), trsc2ID)

#         break


# def writeMembraneBoundaryFile(time, region_dict, sim_properties):
    
#     timestep = int(time/sim_properties['timestep'])
    
#     memBoundaryCoords = np.argwhere(region_dict['membrane']['shape']==True)
    
#     memBoundaryCoords = memBoundaryCoords - np.array(sim_properties['lattice_center'])

#     memBoundaryCoords = memBoundaryCoords * (sim_properties['lattice_spacing'] / 1e-9) * 10
    
#     fname = sim_properties['working_directory'] + 'mem_boundary_{:d}.bin'.format(timestep)
    
#     memBoundaryCoords = np.reshape(memBoundaryCoords,(int(len(memBoundaryCoords)),3),order='F')
    
#     with open(fname, 'wb') as f:
        
#         memBoundaryCoords.tofile(f)
        
#     return
