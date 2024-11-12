"""
Authors: Zane Thornburg

Updating ribosome positions on the RDME lattice
"""

import numpy as np

from scipy import spatial

from scipy.optimize import fsolve

import time as TIME

from LatticeFunctions import *

#########################################################################################
def placeRibosomes(lattice, sim_properties, region_dict, ribo_site_dict, updateTranslat=True, growth_step=False):
    """
    Inputs:
    lattice - LM lattice object including particle and site lattice
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    N_edges = sim_properties['lattice_edges']
    ribo_IDs = sim_properties['riboIDs']
    
    membrane = region_dict['membrane']['shape']
    extracellular = region_dict['extracellular']['shape']
    cyto_shell = region_dict['outer_cytoplasm']['shape']
    
    plattice = lattice.getParticleLatticeView()

    mvPartIdx = 2
    
    ribosome_centers = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    ribosomes = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    she_coords = np.argwhere(cyto_shell==True)
    
    tree = spatial.KDTree(she_coords)
    
    ribo_center_points = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    RIBOidx = sim_properties['name_to_index']['ribosomeP']

    RIBOcoords = np.argwhere(plattice==RIBOidx)
    
    for coord in RIBOcoords:
        
        ribo_site_dict['ribos']['centers'][int(coord[1]), int(coord[2]), int(coord[3])] = True
    
    for ribo_type, type_dict in ribo_site_dict.items():

        centerID = type_dict['center_idx']

        ribo_center_sites = np.argwhere(type_dict['centers']==True)
        
        try:
            old_ribo_sites = np.concatenate((old_ribo_sites, ribo_center_sites))
        except:
            old_ribo_sites = ribo_center_sites
            
        crossID = type_dict['cross_idx']

        ribo_cross_sites = np.argwhere(type_dict['crosses']==True)
            
        old_ribo_sites = np.concatenate((old_ribo_sites, ribo_cross_sites))

#     print(len(old_ribo_sites))
    start = TIME.time()
    ribo_count = 0
    for riboSite in old_ribo_sites:
        
        x, y, z = riboSite[:]

        if lattice.getOccupancy( int(z), int(y), int(x))>0:
                
            parts = getParticlesInSite(plattice, int(z), int(y), int(x))

            if len(parts)>0:

                for particleIdx in parts:

                    if ribo_IDs[particleIdx-1] == True:

                        x_int, y_int, z_int = int(x), int(y), int(z)

                        if (membrane[x_int,y_int,z_int]==True) or (extracellular[x_int,y_int,z_int]==True):

                            dist, sheCoordIdx = tree.query([x_int,y_int,z_int])

                            new_position = she_coords[sheCoordIdx]
                            
                            deleteParticle(plattice,int(z_int),int(y_int),int(x_int),int(particleIdx))
                            
                            Occ = lattice.getOccupancy(int(new_position[2]),int(new_position[1]),int(new_position[0]))

                            if Occ < 15:

                                lattice.addParticle(int(new_position[2]), int(new_position[1]), int(new_position[0]), int(particleIdx))

                            else:
                                
                                print("Lattice Site Full: ", new_position)
                                print("Searching for Available Alternative Ribosome Desitnation Site")

                                placed = False

                                place_counter = 1

                                while not placed:

                                    siteDists, siteIdxs = tree.query([x_int,y_int,z_int], k=int(10*place_counter))

                                    for regionCoordIdx in siteIdxs:

                                        test_position = she_coords[regionCoordIdx]

                                        NewOcc = lattice.getOccupancy(int(test_position[2]), int(test_position[1]), int(test_position[0]))

                                        if NewOcc < 15:

                                            print("New Destionation Site Found: ", test_position)

                                            lattice.addParticle(int(test_position[2]), int(test_position[1]), int(test_position[0]), int(particleIdx))
                                            
                                            new_position = test_position

                                            placed = True

                                            break

                                    place_counter = place_counter + 1


                            ribo_center_points[int(new_position[0]),int(new_position[1]),int(new_position[2])] = True
                            
                            ribo_count = ribo_count + 1

                        else:

                            ribo_center_points[x_int,y_int,z_int] = True
                            
                            ribo_count = ribo_count + 1
          
    ribo_center_points = np.argwhere(ribo_center_points==True)
        
    ribo_types = ['ribos']
    
    for ribo_type in ribo_types:

        ribo_site_dict[ribo_type]['centers'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
        ribo_site_dict[ribo_type]['crosses'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
        

        
    mono_ribos = []

    for coord in ribo_center_points:

        x = int(coord[0])
        y = int(coord[1])
        z = int(coord[2])

        xyz = [x,y,z]
            
        mono_ribos.append([x,y,z])

        
    for ribo_type, type_dict in ribo_site_dict.items():

        if ribo_type == 'ribos':

            for center_point in mono_ribos:

                x_int = int(center_point[0])
                y_int = int(center_point[1])
                z_int = int(center_point[2])

                type_dict['centers'][x_int,y_int,z_int] = True

                type_dict['crosses'][x_int+1,y_int,z_int] = True
                type_dict['crosses'][x_int-1,y_int,z_int] = True
                type_dict['crosses'][x_int,y_int+1,z_int] = True
                type_dict['crosses'][x_int,y_int-1,z_int] = True
                type_dict['crosses'][x_int,y_int,z_int+1] = True
                type_dict['crosses'][x_int,y_int,z_int-1] = True


    return ribo_site_dict
#########################################################################################


#########################################################################################
def updateRiboSites(lattice, ribo_site_dict, region_dict):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    for ribo_type, type_dict in ribo_site_dict.items():
                        
        crossID = type_dict['cross_idx']
        centerID = type_dict['center_idx']

        old_centers = np.argwhere(region_dict[centerID]['shape']==True)

        for site in old_centers:

            if region_dict['cytoplasm']['shape'][site[0], site[1], site[2]] == True:

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['cytoplasm']['index'])

                continue

            if region_dict['outer_cytoplasm']['shape'][site[0], site[1], site[2]] == True:

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['outer_cytoplasm']['index'])

                continue
                
            if region_dict['DNA']['shape'][site[0], site[1], site[2]] == True:

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['DNA']['index'])

                continue

        old_cross = np.argwhere(region_dict[crossID]['shape']==True)

        for site in old_cross:

            if region_dict['cytoplasm']['shape'][site[0], site[1], site[2]] == True:

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['cytoplasm']['index'])

                continue

            if region_dict['outer_cytoplasm']['shape'][site[0], site[1], site[2]] == True:

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['outer_cytoplasm']['index'])

                continue
                
            if region_dict['DNA']['shape'][site[0], site[1], site[2]] == True:

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['DNA']['index'])

                continue

            if region_dict['extracellular']['shape'][site[0], site[1], site[2]] == True:

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['extracellular']['index'])

                continue

        region_dict[centerID]['shape'] = type_dict['centers']
        region_dict[crossID]['shape'] = type_dict['crosses']


    for ribo_type, type_dict in ribo_site_dict.items():

        crossID = type_dict['cross_idx']

        ribo_sites = np.argwhere(type_dict['crosses']==True)

        for site in ribo_sites:

            if (region_dict['membrane']['shape'][site[0], site[1], site[2]] == False) and (region_dict['DNA']['shape'][site[0], site[1], site[2]] == False):

                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict[crossID]['index'])


    for ribo_type, type_dict in ribo_site_dict.items():

        centerID = type_dict['center_idx']

        ribo_center_sites = np.argwhere(type_dict['centers']==True)

        for site in ribo_center_sites:

            lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict[centerID]['index'])
            
    
#     print("Updated Regions with New Ribosome Positions")
    
    return region_dict
#########################################################################################


    