"""
Authors: Zane Thornburg

Create an initial configuration for a single chromsome using sc_chain_generation program by Benjamin Gilbert
"""

import numpy as np

import os

from scipy import spatial

from scipy.optimize import fsolve


#########################################################################################
def InitDnaSites(sim_properties, ribo_site_dict):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    DNAsites - 3D binary array defining RDME lattice sites occupied by DNA
    
    Called by:
    regions_and_complexes.py - buildRegions
    
    Description:
    Use sc_chain_generation by Benjamin Gilbert to create an initial chromosome configuration at 10 bp resolution
    """
    
    writeChromosomeInputFile(sim_properties['cyto_radius'], sim_properties)
    
    writeRiboObstacleFile(ribo_site_dict, sim_properties)
    
    runNewChromosome(sim_properties)
    
    DNAsites = assignDnaSites(sim_properties)
        
    return DNAsites
#########################################################################################


#########################################################################################
def runNewChromosome(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    None
    
    Called by:
    InitDnaSites
    
    Description:
    """
    
    input_fname = sim_properties['working_directory']+'DNA/' + 'Syn3A_chromosome_init.inp'
    
    log_fname = sim_properties['working_directory'] + 'log_init.log'
    
    DNA_executable = sim_properties['dna_software_directory'] + 'sc_chain_generation/fortran/gen_sc_chain --i_f=%s --o_d=%s --o_l=%s --s=10 --l=%s --n_t=8 --bin --xyz' % (input_fname, sim_properties['working_directory']+'DNA/', 'Syn3A_chromosome_init', log_fname)
    print(DNA_executable)
    
    os.system(DNA_executable)
    
    return None
#########################################################################################

    
#########################################################################################
def writeChromosomeInputFile(cyto_radius, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    None
    
    Called by:
    InitDnaSites
    
    Description:
    """
    
    cyto_radius_angstroms = cyto_radius*sim_properties['lattice_spacing']*10/1e-9 - sim_properties['lattice_spacing']*10/1e-9
    
    fname = sim_properties['working_directory']+'DNA/' + 'Syn3A_chromosome_init.inp'
    
    with open(fname, 'w') as f:
        
        f.write('run = 2\n')
        f.write('min_rep = 1\n')
        f.write('max_rep = 1\n')
        f.write('\n')
        f.write('seed = {:d}\n'.format(int(sim_properties['dna_rng_seed'])))
        f.write('r = 1.7E+1\n')
        f.write('r_sc = 2.4E+1\n')
        f.write('R_o = 1.0E+2\n')
        f.write('R_b = {:.2e}\n'.format(cyto_radius_angstroms))
        f.write('\n')
        f.write('cyclic_shift = 0\n')
        f.write('knot_prevent = 1\n')
        f.write('\n')
        f.write('N_o = -1\n')
        f.write('obstacle_file = ' + sim_properties['working_directory']+'DNA/' + 'ribo_obstacles_init.inp\n')
        f.write('\n')
        f.write('N_chains = 1\n')
        f.write('\n')
        f.write('N_stages = 4\n')
        f.write('12,2000\n')
        f.write('6,8000\n')
        f.write('3,18000\n')
        f.write('1,54338\n')
                
    return None
#########################################################################################

        
#########################################################################################
def writeRiboObstacleFile(ribo_site_dict, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    None
    
    Called by:
    InitDnaSites
    
    Description:
    """
    
    N_edges = sim_properties['lattice_edges']
    
    ribo_center_lattice = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    for ribo_type, type_dict in ribo_site_dict.items():
        
        try:
        
            center_lattice = type_dict['centers']

            ribo_center_lattice = ribo_center_lattice | center_lattice
            
        except:
            
            continue
            
            
    ribo_center_points = np.argwhere(ribo_center_lattice == True)
#     print(ribo_center_points)
            
    ribo_center_points = ribo_center_points * (sim_properties['lattice_spacing'] / 1e-9) * 10
    
    print(len(ribo_center_points))

    totalRiboCount = len(ribo_center_points)
    
    fname = sim_properties['working_directory']+'DNA/' + 'ribo_obstacles_init.inp'
    
    with open(fname, 'w') as f:
        
        f.write('N = {:d}\n'.format(totalRiboCount))
        
        riboIdx = 1
        
        for coord in ribo_center_points:
#             print(riboIdx,coord[0],coord[1],coord[2])
            f.write('{:d},{:.2e},{:.2e},{:.2e}\n'.format(riboIdx,coord[0],coord[1],coord[2]))
            riboIdx = riboIdx + 1
        
    return None
#########################################################################################


#########################################################################################
def assignDnaSites(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    DNAsites - 3D binary array of lattice sites occupied by DNA
    
    Called by:
    InitDnaSites
    
    Description:
    """
    
#     timestep = int(time/sim_properties['timestep'])
    
    DNAfile = sim_properties['working_directory']+'DNA/' + 'x_chain_Syn3A_chromosome_init_rep00001.xyz'
    
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
        
        for DNAparticle in DNAcoords:
            
            x = DNAparticle[0]
            y = DNAparticle[1]
            z = DNAparticle[2]
            
            x_lattice = (x*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_x
            y_lattice = (y*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_y
            z_lattice = (z*1e-9)//(10*sim_properties['lattice_spacing'])+N_2_z
            
            DNAsites[int(x_lattice),int(y_lattice),int(z_lattice)] = True
            
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
        
#                     lattice.addParticle(int(z_lattice),int(y_lattice),int(x_lattice),int(3))
            
            print('Initialized DNA lattice sites')
        
    sim_properties['DNAcoords'] = DNA_lattice_coords
            
#     region_dict['DNA']['shape'] = DNAsites
    
    return DNAsites
#########################################################################################


# def placeDnaParticles(sim, sim_properties):
                
#     DNAfile = sim_properties['working_directory']+'DNA/' + 'x_chain_Syn3A_chromosome_init_rep00001.xyz'
    
#     fileType = DNAfile.split('.')[1] #'bin'
    
#     print(fileType)
    
#     dna_part = sim.species('DNAparticle')
                
#     N_edges = sim_properties['lattice_edges']
    
#     N_edges_x = N_edges[0] # Number of subvolumes making up and edge of the simulation space N x N x N
#     N_edges_y = N_edges[1]
#     N_edges_z = N_edges[2]

#     N_2_x = N_edges_x/2
#     N_2_y = N_edges_y/2
#     N_2_z = N_edges_z/2
    
#     if fileType=='bin':
        
#         with open(DNAfile,'rb') as f:
            
#             DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
            
#         DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
        
#         print(DNAcoords.shape)
        
# #         DNAcoordsLattice = []
        
#         for DNAparticle in DNAcoords:
            
#             x = DNAparticle[0]
#             y = DNAparticle[1]
#             z = DNAparticle[2]
            
#             x_lattice = (x*1e-9)//(10*sim.latticeSpacing)+N_2_x
#             y_lattice = (y*1e-9)//(10*sim.latticeSpacing)+N_2_y
#             z_lattice = (z*1e-9)//(10*sim.latticeSpacing)+N_2_z
            
#             dna_part.placeParticle(int(x_lattice),int(y_lattice),int(z_lattice),1)
            
# #         print(DNAcoords.shape)
# #         print(DNAcoords)
        
# #         print(DNAcoordsLattice)
        
# #         print(np.min(DNAcoords))
# #         print(np.max(DNAcoords))
        
# #         print(np.min(DNAcoordsLattice))
# #         print(np.max(DNAcoordsLattice))
        
#     elif fileType=='xyz':
        
#         dna_part_added = 0
        
#         with open(DNAfile,'rb') as f:
            
#             for line_number,line in enumerate(f):
                
#                 if line_number == 0:
                    
#                     DNAparticleCount=int(line)
                    
#                 elif line_number == 1:
                    
#                     continue
                    
#                 else:
                    
#                     atomic_symbol, x, y, z = line.split()
                    
#                     x_lattice = (float(x)*1e-9)//(10*sim.latticeSpacing)+N_2_x
#                     y_lattice = (float(y)*1e-9)//(10*sim.latticeSpacing)+N_2_y
#                     z_lattice = (float(z)*1e-9)//(10*sim.latticeSpacing)+N_2_z
# #                     print(int(x_lattice),int(y_lattice),int(z_lattice))

#         #             DNAcoordsLattice.append([x_lattice,y_lattice,z_lattice])
#                     dna_part.placeParticle(int(x_lattice),int(y_lattice),int(z_lattice),1)
        
#                     dna_part_added = dna_part_added + 1
            
#         print(dna_part_added)
        
#     return
                