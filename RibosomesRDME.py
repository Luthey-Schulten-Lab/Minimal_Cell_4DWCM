"""
Authors: Zane Thornburg

Updating ribosome positions on the RDME lattice
"""

import numpy as np
from scipy import spatial
from LatticeFunctions import *

# Module-level cache for KDTree (Optimization 1)
# Note: Only caching tree and size to avoid expensive array comparisons
_cached_tree = None
_cached_she_coords_size = None

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
    global _cached_tree, _cached_she_coords_size
    
    N_edges = sim_properties['lattice_edges']
    ribo_IDs = sim_properties['riboIDs']
    
    membrane = region_dict['membrane']['shape']
    extracellular = region_dict['extracellular']['shape']
    cyto_shell = region_dict['outer_cytoplasm']['shape']
    
    plattice = lattice.getParticleLatticeView()
    
    ribo_center_points = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
    
    # Optimization 1: KDTree caching - only rebuild when cyto_shell size changes
    # Note: We use size-based invalidation to avoid expensive array comparisons
    she_coords = np.argwhere(cyto_shell==True)
    current_size = len(she_coords)
    
    if (_cached_tree is None or _cached_she_coords_size != current_size):
        if current_size > 0:
            tree = spatial.KDTree(she_coords)
            _cached_tree = tree
            _cached_she_coords_size = current_size
        else:
            # No valid shell coordinates - reset and return
            ribo_site_dict['ribos']['centers'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
            ribo_site_dict['ribos']['crosses'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
            return ribo_site_dict
    else:
        tree = _cached_tree
    
    RIBOidx = sim_properties['name_to_index']['ribosomeP']
    RIBOcoords = np.argwhere(plattice==RIBOidx)
    
    # Optimization 2: Vectorized ribosome center updates
    if len(RIBOcoords) > 0:
        ribo_site_dict['ribos']['centers'][RIBOcoords[:, 1], RIBOcoords[:, 2], RIBOcoords[:, 3]] = True
    
    # Optimization 3: Efficient collection of old ribosome sites
    # Collect in a list first, then vstack once (more efficient than multiple vstacks)
    old_ribo_sites_list = []
    for ribo_type, type_dict in ribo_site_dict.items():
        ribo_center_sites = np.argwhere(type_dict['centers']==True)
        if len(ribo_center_sites) > 0:
            old_ribo_sites_list.append(ribo_center_sites)
        
        ribo_cross_sites = np.argwhere(type_dict['crosses']==True)
        if len(ribo_cross_sites) > 0:
            old_ribo_sites_list.append(ribo_cross_sites)
    
    if not old_ribo_sites_list:
        # No ribosomes to process
        ribo_site_dict['ribos']['centers'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
        ribo_site_dict['ribos']['crosses'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
        return ribo_site_dict
    
    # Stack all sites at once (more efficient than multiple vstacks)
    old_ribo_sites = np.vstack(old_ribo_sites_list)
    
    # Optimization 4: Process all sites with inline relocation checks
    # (Inline checks avoid creating large intermediate boolean arrays)
    for riboSite in old_ribo_sites:
        x, y, z = riboSite[:]
        x_int, y_int, z_int = int(x), int(y), int(z)
        
        if lattice.getOccupancy(z_int, y_int, x_int) > 0:
            parts = getParticlesInSite(plattice, z_int, y_int, x_int)
            
            if len(parts) > 0:
                for particleIdx in parts:
                    if ribo_IDs[particleIdx-1]:
                        # Check if relocation needed (inline check like original)
                        if (membrane[x_int, y_int, z_int] == True) or (extracellular[x_int, y_int, z_int] == True):
                            # Need relocation - use KDTree
                            dist, sheCoordIdx = tree.query([x_int, y_int, z_int])
                            new_position = she_coords[sheCoordIdx]
                            
                            deleteParticle(plattice, z_int, y_int, x_int, particleIdx)
                            
                            new_z, new_y, new_x = int(new_position[2]), int(new_position[1]), int(new_position[0])
                            Occ = lattice.getOccupancy(new_z, new_y, new_x)
                            
                            if Occ < 15:
                                lattice.addParticle(new_z, new_y, new_x, int(particleIdx))
                            else:
                                print("Lattice Site Full: ", new_position)
                                print("Searching for Available Alternative Ribosome Destination Site")
                                
                                placed = False
                                place_counter = 1
                                
                                while not placed:
                                    k_neighbors = min(10 * place_counter, len(she_coords))
                                    siteDists, siteIdxs = tree.query([x_int, y_int, z_int], k=k_neighbors)
                                    
                                    # Handle single result case
                                    if k_neighbors == 1:
                                        siteIdxs = [siteIdxs]
                                    
                                    for regionCoordIdx in siteIdxs:
                                        test_position = she_coords[regionCoordIdx]
                                        test_z, test_y, test_x = int(test_position[2]), int(test_position[1]), int(test_position[0])
                                        
                                        NewOcc = lattice.getOccupancy(test_z, test_y, test_x)
                                        
                                        if NewOcc < 15:
                                            print("New Destination Site Found: ", test_position)
                                            lattice.addParticle(test_z, test_y, test_x, int(particleIdx))
                                            new_position = test_position
                                            placed = True
                                            break
                                    
                                    place_counter += 1
                                
                                if not placed:
                                    print("Warning: Could not place ribosome particle", particleIdx)
                                    continue
                            
                            ribo_center_points[int(new_position[0]), int(new_position[1]), int(new_position[2])] = True
                        else:
                            # No relocation needed
                            ribo_center_points[x_int, y_int, z_int] = True
    
    # Reset ribosome site dictionaries
    ribo_site_dict['ribos']['centers'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
    ribo_site_dict['ribos']['crosses'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
    
    # Get ribosome center points
    ribo_center_points_coords = np.argwhere(ribo_center_points)
    
    if len(ribo_center_points_coords) == 0:
        return ribo_site_dict
    
    # Optimization 5: Vectorized cross-voxel updates with bounds checking
    x = ribo_center_points_coords[:, 0]
    y = ribo_center_points_coords[:, 1]
    z = ribo_center_points_coords[:, 2]
    
    # Set centers
    ribo_site_dict['ribos']['centers'][x, y, z] = True
    
    # Set crosses with np.clip for bounds checking
    x_plus = np.clip(x + 1, 0, N_edges[0] - 1)
    x_minus = np.clip(x - 1, 0, N_edges[0] - 1)
    y_plus = np.clip(y + 1, 0, N_edges[1] - 1)
    y_minus = np.clip(y - 1, 0, N_edges[1] - 1)
    z_plus = np.clip(z + 1, 0, N_edges[2] - 1)
    z_minus = np.clip(z - 1, 0, N_edges[2] - 1)
    
    ribo_site_dict['ribos']['crosses'][x_plus, y, z] = True
    ribo_site_dict['ribos']['crosses'][x_minus, y, z] = True
    ribo_site_dict['ribos']['crosses'][x, y_plus, z] = True
    ribo_site_dict['ribos']['crosses'][x, y_minus, z] = True
    ribo_site_dict['ribos']['crosses'][x, y, z_plus] = True
    ribo_site_dict['ribos']['crosses'][x, y, z_minus] = True

    return ribo_site_dict
#########################################################################################


#########################################################################################
def updateRiboSites(lattice, ribo_site_dict, region_dict, sim_properties=None):
    """
    Inputs:
    lattice - LM lattice object
    ribo_site_dict - Dictionary of ribosome site information
    region_dict - Dictionary of region information
    sim_properties - Dictionary of simulation properties (optional, for compatibility)
    
    Returns:
    region_dict - Updated region dictionary
    
    Called by:
    Hook.py, Division.py, Restart_Hook.py
    
    Description:
    Updates ribosome sites on the lattice based on current ribosome positions
    """
    
    for ribo_type, type_dict in ribo_site_dict.items():
        crossID = type_dict['cross_idx']
        centerID = type_dict['center_idx']

        old_centers = np.argwhere(region_dict[centerID]['shape']==True)
        
        if len(old_centers) > 0:
            # Vectorized region checking
            cyto_mask = region_dict['cytoplasm']['shape'][old_centers[:, 0], old_centers[:, 1], old_centers[:, 2]]
            outer_cyto_mask = region_dict['outer_cytoplasm']['shape'][old_centers[:, 0], old_centers[:, 1], old_centers[:, 2]]
            dna_mask = region_dict['DNA']['shape'][old_centers[:, 0], old_centers[:, 1], old_centers[:, 2]]
            
            # Process each region type (priority: cytoplasm > outer_cytoplasm > DNA)
            cyto_sites = old_centers[cyto_mask]
            for site in cyto_sites:
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['cytoplasm']['index'])
            
            outer_cyto_sites = old_centers[outer_cyto_mask & ~cyto_mask]
            for site in outer_cyto_sites:
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['outer_cytoplasm']['index'])
            
            dna_sites = old_centers[dna_mask & ~cyto_mask & ~outer_cyto_mask]
            for site in dna_sites:
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['DNA']['index'])

        old_cross = np.argwhere(region_dict[crossID]['shape']==True)
        
        if len(old_cross) > 0:
            # Vectorized region checking for crosses
            cyto_mask = region_dict['cytoplasm']['shape'][old_cross[:, 0], old_cross[:, 1], old_cross[:, 2]]
            outer_cyto_mask = region_dict['outer_cytoplasm']['shape'][old_cross[:, 0], old_cross[:, 1], old_cross[:, 2]]
            dna_mask = region_dict['DNA']['shape'][old_cross[:, 0], old_cross[:, 1], old_cross[:, 2]]
            extra_mask = region_dict['extracellular']['shape'][old_cross[:, 0], old_cross[:, 1], old_cross[:, 2]]
            
            # Process each region type (priority: cytoplasm > outer_cytoplasm > DNA > extracellular)
            cyto_sites = old_cross[cyto_mask]
            for site in cyto_sites:
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['cytoplasm']['index'])
            
            outer_cyto_sites = old_cross[outer_cyto_mask & ~cyto_mask]
            for site in outer_cyto_sites:
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['outer_cytoplasm']['index'])
            
            dna_sites = old_cross[dna_mask & ~cyto_mask & ~outer_cyto_mask]
            for site in dna_sites:
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['DNA']['index'])
            
            extra_sites = old_cross[extra_mask & ~cyto_mask & ~outer_cyto_mask & ~dna_mask]
            for site in extra_sites:
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict['extracellular']['index'])

        region_dict[centerID]['shape'] = type_dict['centers']
        region_dict[crossID]['shape'] = type_dict['crosses']

    # Set cross sites
    for ribo_type, type_dict in ribo_site_dict.items():
        crossID = type_dict['cross_idx']
        ribo_sites = np.argwhere(type_dict['crosses']==True)
        
        if len(ribo_sites) > 0:
            # Vectorized checking
            membrane_mask = region_dict['membrane']['shape'][ribo_sites[:, 0], ribo_sites[:, 1], ribo_sites[:, 2]]
            dna_mask = region_dict['DNA']['shape'][ribo_sites[:, 0], ribo_sites[:, 1], ribo_sites[:, 2]]
            valid_mask = ~membrane_mask & ~dna_mask
            
            valid_sites = ribo_sites[valid_mask]
            for site in valid_sites:
                lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict[crossID]['index'])

    # Set center sites
    for ribo_type, type_dict in ribo_site_dict.items():
        centerID = type_dict['center_idx']
        ribo_center_sites = np.argwhere(type_dict['centers']==True)
        
        for site in ribo_center_sites:
            lattice.setSiteType(int(site[2]), int(site[1]), int(site[0]), region_dict[centerID]['index'])
    
    return region_dict
#########################################################################################