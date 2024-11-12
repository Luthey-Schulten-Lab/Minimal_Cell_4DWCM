"""
Authors: Zane Thornburg

Functions to update cell morphology before division has started
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

from LatticeFunctions import *


#########################################################################################
def grow_cell(RDMEsim, lattice, sim_properties, region_dict, ribo_site_dict, cyto_radius):
    """
    Inputs:
    
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """

    print('GROWING CELL')

    #update the site types
    print("start to add layer",cyto_radius)

    old_region_shapes = {}
    
    for region, type_dict in region_dict.items():
        
        old_region_shapes[region] = type_dict['shape']
    
    region_dict = buildNewRegions(RDMEsim, sim_properties['lattice_center'], cyto_radius, region_dict, ribo_site_dict)

    for region, type_dict in region_dict.items():

        regionIdx = type_dict['index']
        
        region_sites = np.argwhere(type_dict['shape']==True)
        
        for site in region_sites:
            
            lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), regionIdx)

#         region_count=region_count+1
        
    moveParticles('membrane', lattice, region_dict, old_region_shapes)
    
    moveParticles('outer_cytoplasm', lattice, region_dict, old_region_shapes)

    print('GROWTH STEP COMPLETE')

    return region_dict
#########################################################################################


#########################################################################################
def buildNewRegions(RDMEsim, sim_center, cyto_radius, region_dict, ribo_site_dict):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """

    build = RegionBuilder(RDMEsim)

    cytoplasm = build.ellipsoid(radius = cyto_radius, center = sim_center)

    cyto_dilation = build.dilate(cytoplasm, se = build.se26)
    cyto_shell = cyto_dilation & ~cytoplasm
    cyto_dilation = build.dilate(cyto_dilation, se = build.se26)
    membrane = cyto_dilation & ~cyto_shell & ~cytoplasm
    extracellular = ~cyto_dilation
    
    region_dict['outer_cytoplasm']['full_shape'] = cyto_shell
    region_dict['cytoplasm']['full_shape'] = cytoplasm
    
    cytoplasm = cytoplasm & ~region_dict["DNA"]["shape"]
    cyto_shell = cyto_shell & ~region_dict["DNA"]["shape"]
    
    region_dict['extracellular']['shape'] = extracellular
    region_dict['membrane']['shape'] = membrane
    region_dict['outer_cytoplasm']['shape'] = cyto_shell
    region_dict['cytoplasm']['shape'] = cytoplasm

#     for ribo_type, type_dict in ribo_site_dict.items():
        
#         region_dict[type_dict['center_idx']]['shape'] = type_dict['centers']
#         region_dict[type_dict['cross_idx']]['shape'] = type_dict['crosses']

    print('Geometry constructed')

    return region_dict
#########################################################################################


#########################################################################################
def moveParticles(region, lattice, region_dict, old_region_shapes):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    
    Returns:
    Called by:
    Description:
    """
    
    start = TIME.time()
    
    print('Moving Particles in Region: ', region)
    
    siteIdx = region_dict[region]['index']

    region_coords = np.argwhere(region_dict[region]['shape']==True)
    print(len(region_coords))

    region_tree = spatial.KDTree(region_coords)
    
    old_region_coords = np.argwhere(old_region_shapes[region]==True)
    print(len(old_region_coords))
    
    plattice = lattice.getParticleLatticeView()
    
    true_count = 0
    
    for siteCoord in old_region_coords:
        
        x = siteCoord[0]
        y = siteCoord[1]
        z = siteCoord[2]
        
#         if region_dict[region]['shape'][int(x), int(y), int(z)] == True:
#             true_count = true_count + 1
        
        if region_dict[region]['shape'][int(x), int(y), int(z)] != True:
            
#             print('Moving particles at site: ', siteCoord)
            
            parts = getParticlesInSite(plattice, int(z), int(y), int(x))
            
            dist, regionSiteIdx = region_tree.query(siteCoord)
            
            newSiteCoord = region_coords[regionSiteIdx]
            
            for partIdx in parts:
                
                deleteParticle(plattice,int(z),int(y),int(x),partIdx)
            
                Occ = lattice.getOccupancy(int(newSiteCoord[2]),int(newSiteCoord[1]),int(newSiteCoord[0]))

                if Occ < 15:

                    lattice.addParticle(int(newSiteCoord[2]),int(newSiteCoord[1]),int(newSiteCoord[0]),int(partIdx))

                else:

                    print("Lattice Site Full: ", newSiteCoord)
                    print("Searching for Available Alternative Desitnation Site")

                    placed = False

                    place_counter = 1

                    while not placed:

                        siteDists, siteIdxs = region_tree.query(siteCoord, k=int(10*place_counter))

                        for regionCoordIdx in siteIdxs:

                            test_position = region_coords[regionCoordIdx]

                            NewOcc = lattice.getOccupancy(int(test_position[2]),int(test_position[1]),int(test_position[0]))

                            if NewOcc < 15:

                                print("New Destionation Site Found: ", test_position)

                                lattice.addParticle(int(test_position[2]),int(test_position[1]),int(test_position[0]),int(partIdx))

                                placed = True

                                break

                        place_counter = place_counter + 1
                        
    print('Total time to move particles: ', TIME.time()-start)
    print('Moved Particles for Region: ', region)
    
    return None
#########################################################################################

