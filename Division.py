"""
Authors: Zane Thornburg

Functions to update cell morphology during division in 3D
"""

from jLM.RegionBuilder import RegionBuilder
from jLM.RDME import Sim as RDMESim
from jLM.RDME import File as RDMEFile
import jLM

import numpy as np

from scipy import spatial

from scipy.optimize import fsolve

import time as TIME

import RibosomesRDME as ribosomesRDME

import FreeDTS_functions as fdf

from LatticeFunctions import *


#########################################################################################
def divide_cell(RDMEsim, lattice, sim_properties, region_dict, ribo_site_dict):
    """
    Inputs:
    RDMEsim - jLM RDME simulation object
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """

    print('DIVIDING CELL')
#     print('Cell Surface Area: ', cell_SA)
#     print('Cell Volume: ', cell_V)

    index_to_name = RDMEsim.speciesList.__dict__['_id2name']
    
    old_region_shapes = {}
    
    for region, type_dict in region_dict.items():
        
        if 'cytoplasm' in region:
            
            old_region_shapes[region] = type_dict['full_shape']
            
#         elif region == 'membrane':
            
#             try:
#                 print('Old Membrane Shape No DNA')
#                 old_region_shapes[region] = type_dict['shape_noDNA']
#             except:
#                 old_region_shapes[region] = type_dict['shape']
            
        else:
        
            old_region_shapes[region] = type_dict['shape']
    
    sim_properties, region_dict, ribo_site_dict = buildNewDivRegions(RDMEsim, sim_properties, lattice, sim_properties['lattice_center'], region_dict, ribo_site_dict)

    for region, type_dict in region_dict.items():

        regionIdx = type_dict['index']
        
        region_sites = np.argwhere(type_dict['shape']==True)
        
        for site in region_sites:
            
            lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), regionIdx)
            
            
    old_plattice = np.array(lattice.getParticleLatticeView())
        
    
    moveParticles('cytoplasm', lattice, sim_properties, region_dict, old_region_shapes, old_plattice, index_to_name)
    
    moveParticles('outer_cytoplasm', lattice, sim_properties, region_dict, old_region_shapes, old_plattice, index_to_name)
    
    moveParticles('membrane', lattice, sim_properties, region_dict, old_region_shapes, old_plattice, index_to_name)
    
    
    
#     moveParticles('DNA', lattice, region_dict, old_region_shapes)
    
#     moveParticles('ribosomes', lattice, region_dict, old_region_shapes)
    
#     moveParticles('ribo_centers', lattice, region_dict, old_region_shapes)

    print('DIVISION STEP COMPLETE')
    
#     if cutoff < 10:
        
#         cutoff = 10
#         radius = 200

    return region_dict, ribo_site_dict
#########################################################################################


#########################################################################################
def moveParticles(region, lattice, sim_properties, region_dict, old_region_shapes, old_plattice, index_to_name):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    None
    
    Called by:
    Description:
    """
    
    start = TIME.time()
    
    ribo_IDs = sim_properties['riboIDs']
    
    print('Moving Particles in Region: ', region)
    
    siteIdx = region_dict[region]['index']
    
    if 'cytoplasm' in region:
        
        region_coords = np.argwhere(region_dict[region]['full_shape']==True)
        
#     elif region == 'membrane':
        
#         print('Membrane No DNA')
    
#         region_coords = np.argwhere(region_dict[region]['shape_noDNA']==True)
#         print(len(region_coords))
        
    else:

        region_coords = np.argwhere(region_dict[region]['shape']==True)
        print(len(region_coords))

    region_tree = spatial.KDTree(region_coords)
    
    old_region_coords = np.argwhere(old_region_shapes[region]==True)
    print(len(old_region_coords))
    
    plattice = lattice.getParticleLatticeView()
    
    true_count = 0
    
    if 'cytoplasm' in region:
        
        region_shape = region_dict[region]['full_shape']
        
    else:
        
        region_shape = region_dict[region]['shape']
    
    for siteCoord in old_region_coords:
        
        x = siteCoord[0]
        y = siteCoord[1]
        z = siteCoord[2]
        
#         if region_dict[region]['shape'][int(x), int(y), int(z)] == True:
#             true_count = true_count + 1
        
        if region_shape[int(x), int(y), int(z)] != True:
            
#             print('Moving particles at site: ', siteCoord)
            
            parts = getParticlesInSite(old_plattice, int(z), int(y), int(x))
            
            dist, regionSiteIdx = region_tree.query(siteCoord)
            
            newSiteCoord = region_coords[regionSiteIdx]
            
            for partIdx in parts:
                
#                 if ('G_' in index_to_name[partIdx]) or ('RP_' in index_to_name[partIdx]):
                    
#                     print(index_to_name[partIdx])
#                     print(sim_properties['name_to_index'][index_to_name[partIdx]])
#                     print(partIdx)
                
                if (ribo_IDs[partIdx-1] != True) and ('G_' not in index_to_name[partIdx]) and ('RP_' not in index_to_name[partIdx]):
                
                    deleteParticle(plattice,int(z),int(y),int(x),partIdx)

                    Occ = lattice.getOccupancy(int(newSiteCoord[2]),int(newSiteCoord[1]),int(newSiteCoord[0]))

                    if Occ < 11:

                        lattice.addParticle(int(newSiteCoord[2]),int(newSiteCoord[1]),int(newSiteCoord[0]),int(partIdx))

                    else:

#                         print("Lattice Site Full: ", newSiteCoord)
#                         print("Searching for Available Alternative Desitnation Site")

                        placed = False

                        place_counter = 1

                        while not placed:

                            siteDists, siteIdxs = region_tree.query(siteCoord, k=int(10*place_counter))

                            for regionCoordIdx in siteIdxs:

                                test_position = region_coords[regionCoordIdx]

                                NewOcc = lattice.getOccupancy(int(test_position[2]),int(test_position[1]),int(test_position[0]))

                                if NewOcc < 11:

#                                     print("New Destionation Site Found: ", test_position)

                                    lattice.addParticle(int(test_position[2]),int(test_position[1]),int(test_position[0]),int(partIdx))

                                    placed = True

                                    break

                            place_counter = place_counter + 1
                        
    print('Total time to move particles: ', TIME.time()-start)
    print('Moved Particles for Region: ', region)
    
    return None
#########################################################################################


#########################################################################################
def buildNewDivRegions(RDMEsim, sim_properties, lattice, sim_center, region_dict, ribo_site_dict):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Called by:
    Description:
    """

    build = RegionBuilder(RDMEsim)
    
    N_edges = sim_properties['lattice_edges']

    cellV = sim_properties['volume']/((1e-9)**3)/2

    cellSA = sim_properties['SA']/((1e-9)**2)/2
    
    def division_equations(f):

        cutoff, radius = f

        return ((2/3)*np.pi*radius**3 + np.pi*cutoff*radius**2 - np.pi/3*cutoff**3 - cellV, 2*np.pi*radius*(cutoff+radius) - cellSA)

    cutoff, radius = fsolve(division_equations, (c2c,div_radius))

    print(fsolve(division_equations, (c2c,div_radius)))

    print(cutoff, radius)

    new_radius = radius*(1e-9)/sim_properties['lattice_spacing']
    print(new_radius)

    c2c_dist = (cutoff*2*(1e-9))/sim_properties['lattice_spacing']
    print(c2c_dist)

    cell1_center = [sim_center[0], sim_center[1], sim_center[2]-np.rint(c2c_dist/2)]

    cytoplasm1 = build.ellipsoid(radius = new_radius, center = cell1_center)

    cell2_center = [sim_center[0], sim_center[1], sim_center[2]+np.rint(c2c_dist/2)]

    cytoplasm2 = build.ellipsoid(radius = new_radius, center = cell2_center)

    cytoplasm = cytoplasm1 | cytoplasm2

    cyto_dilation = build.dilate(cytoplasm, se = build.se26)
    cyto_shell = cyto_dilation & ~cytoplasm
    cyto_dilation = build.dilate(cyto_dilation, se = build.se26)
    membrane = cyto_dilation & ~cyto_shell & ~cytoplasm
    extracellular = ~cyto_dilation
    
#     cytoplasm = cytoplasm & ~region_dict["DNA"]["shape"]
#     cyto_shell = cyto_shell & ~region_dict["DNA"]["shape"]

#     ribo_site_dict = ribosomesRDME.placeRibosomes(region_dict, ribo_site_dict, N_edges, ribo_IDs, lattice, growth_step=True)
    
#     region_dict['extracellular']['shape'] = extracellular
#     region_dict['membrane']['shape'] = membrane
#     region_dict['outer_cytoplasm']['shape'] = cyto_shell
#     region_dict['cytoplasm']['shape'] = cytoplasm

#     extracellular = ~membrane & ~cyto_shell & ~cytoplasm

    region_dict['extracellular']['shape'] = extracellular
    region_dict['membrane']['shape'] = membrane
    region_dict['outer_cytoplasm']['full_shape'] = cyto_shell
    region_dict['cytoplasm']['full_shape'] = cytoplasm
    
    cytoplasm_noDNA = cytoplasm & ~region_dict["DNA"]["shape"]
    cyto_shell_noDNA = cyto_shell & ~region_dict["DNA"]["shape"]
    membrane_noDNA = membrane & ~region_dict["DNA"]["shape"]
    
    region_dict['outer_cytoplasm']['shape'] = cyto_shell_noDNA
    region_dict['cytoplasm']['shape'] = cytoplasm_noDNA
    region_dict['membrane']['shape_noDNA'] = membrane_noDNA
    
    ribo_site_dict = ribosomesRDME.placeRibosomes(lattice, sim_properties, region_dict, ribo_site_dict, updateTranslat=False)
    
    region_dict = ribosomesRDME.updateRiboSites(lattice, ribo_site_dict, region_dict)
    
    for ribo_type, type_dict in ribo_site_dict.items():
        
        region_dict[type_dict['center_idx']]['shape'] = type_dict['centers']
        region_dict[type_dict['cross_idx']]['shape'] = type_dict['crosses']

    print('Geometry constructed')

    return sim_properties, region_dict, ribo_site_dict
#########################################################################################
