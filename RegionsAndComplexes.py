"""
Authors: Zane Thornburg

Initialize cell architecture and ribosome positions in the RDME simulation
"""

from jLM.RegionBuilder import RegionBuilder
import jLM

import InitRdmeDna

import pandas as pd
import numpy as np
import random
from collections import defaultdict, OrderedDict

from Bio import SeqIO
from Bio.Seq import Seq


#########################################################################################
def buildRegions(sim, sim_properties):
    """
    Inputs:
    sim - jLM RDME simulation object
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    region_dict - dictionary containing 3D binary arrays that define the shapes of regions on the RDME lattice
    ribo_site_dict - dictionary containing ribosome position information and 3D binary arrays for lattice shapes
    
    Called by:
    Whole_Cell_Minimal_Cell.py
    
    Description:
    Construct the initial region shapes for the RDME simulation
    """
    
    regionIdx = 0
    
    N_edges = sim_properties['lattice_edges']
    
    region_dict = {}
    
#     regions = ["extracellular", "cytoplasm", "outer_cytoplasm", "p_starts", "p_mids", "p_ends", 
#                "ribosomes", "p_start_c", "p_mid_c", "p_end_c", "ribo_centers",
#                "DNA", "membrane"] #, "Z-ring"]

    regions = ["extracellular", "cytoplasm", "outer_cytoplasm",
               "ribosomes", "ribo_centers", "DNA", "membrane"]
    
    for region in regions:
        
        defined_region = sim.region(region)
        region_dict[region] = {}
        region_dict[region]['index'] = regionIdx
        region_dict[region]['shape'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
        
        regionIdx = regionIdx + 1
        
#     index_to_name = sim.regionList.__dict__['_id2name']

#     for index, region in index_to_name:

    genomeFile3A = sim_properties['head_directory'] + 'input_data/syn3A.gb'
    genome3A = next(SeqIO.parse(genomeFile3A, "gb"))
#     genome = mapDNA(genome3A)
    
    sim_properties['genome'], sim_properties['genome_length'] = mapDNA(genome3A)
    
    addChromosomeFeatures(sim_properties)
        
    ribo_types = ['ribos'] #,'starts','mids','ends']
    
    ribo_site_dict = {}
    
    for ribo_type in ribo_types:
        ribo_site_dict[ribo_type] = {}
        
    ribo_site_dict['ribos']['center_idx'] = "ribo_centers"
    ribo_site_dict['ribos']['cross_idx'] = "ribosomes"
#     ribo_site_dict['starts']['center_idx'] = "p_start_c"
#     ribo_site_dict['starts']['cross_idx'] = "p_starts"
#     ribo_site_dict['mids']['center_idx'] = "p_mid_c"
#     ribo_site_dict['mids']['cross_idx'] = "p_mids"
#     ribo_site_dict['ends']['center_idx'] = "p_end_c"
#     ribo_site_dict['ends']['cross_idx'] = "p_ends"

    build = RegionBuilder(sim)

    cytoplasm = build.ellipsoid(radius = sim_properties['cyto_radius'], center = sim_properties['lattice_center'])
    
    cyto_dilation = build.dilate(cytoplasm, se = build.se26)
    cyto_shell = cyto_dilation & ~cytoplasm
    cyto_dilation = build.dilate(cyto_dilation, se = build.se26)
    membrane = cyto_dilation & ~cyto_shell & ~cytoplasm
    extracellular = ~cyto_dilation
    
    region_dict['outer_cytoplasm']['full_shape'] = cyto_shell
    region_dict['cytoplasm']['full_shape'] = cytoplasm
    
    ribosome_centers, ribosomes = placeRibosomes(sim, cytoplasm, sim_properties['lattice_edges'])
    
    region_dict['ribosomes']['shape'] = ribosomes
    region_dict['ribo_centers']['shape'] = ribosome_centers
    
    print(ribosome_centers.shape)
    print(len(np.argwhere(ribosomes==True)))
    
    ribo_site_dict['ribos']['centers'] = ribosome_centers
    ribo_site_dict['ribos']['crosses'] = ribosomes
#     ribo_site_dict['starts']['centers'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['starts']['crosses'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['mids']['centers'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['mids']['crosses'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['ends']['centers'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
#     ribo_site_dict['ends']['crosses'] = np.full((N_edges[0], N_edges[1], N_edges[2]), False)

#     RiboPtnMap = createRibosomalProteinMap(genome3A)
    
    sim_properties['RiboPtnMap'] = createRibosomalProteinMap(genome3A)
    
    DNAsites = InitRdmeDna.InitDnaSites(sim_properties, ribo_site_dict)
    
#     N_edges = sim_properties['lattice_edges']
#     DNAsites = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    
    region_dict['DNA']['shape'] = DNAsites
    
    cytoplasm = cytoplasm & ~DNAsites
    cyto_shell = cyto_shell & ~DNAsites
    
    region_dict['extracellular']['shape'] = extracellular
    region_dict['membrane']['shape'] = membrane
    region_dict['outer_cytoplasm']['shape'] = cyto_shell
    region_dict['cytoplasm']['shape'] = cytoplasm
    
    build.compose(
    (sim.region('extracellular'), extracellular),
    (sim.region('cytoplasm'), cytoplasm),
    (sim.region('outer_cytoplasm'), cyto_shell),
    (sim.region('ribosomes'), ribosomes),
    (sim.region('DNA'),DNAsites),
    (sim.region('ribo_centers'), ribosome_centers),
    (sim.region('membrane'), membrane))
    
    print('Geometry constructed')
    
    return region_dict, ribo_site_dict
#########################################################################################
        

#########################################################################################
def placeRibosomes(sim, cytoplasm, N_edges):
    """
    Inputs:
    sim - jLM RDME simulation object
    cytoplasm - 3D binary array defining cytoplasm region shape on RDME lattice
    
    Returns:
    ribosome_centers - 3D binary array defining positions of ribosome center lattice sites
    ribosomes - 3D binary array defining positions of cross sites of ribosome shapes on RDME lattice
    
    Called by:
    buildRegions
    
    Description:
    Randomly determines positions for 500 ribosomes within the cytoplasm
    """
    
    ribosome_centers = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    ribosomes = np.full((N_edges[0], N_edges[1], N_edges[2]), False)

    riboNum = 500
    
    cyto_coords = np.argwhere(cytoplasm==True)

    ribo_centers = []

    print(len(cyto_coords))
    print(cyto_coords)
    print(ribo_centers)

    for i in range(riboNum):
        rand_Idx = np.random.randint(len(cyto_coords), size=1)[0]
        ribo_centers.append(cyto_coords[rand_Idx])
        cyto_coords = np.delete(cyto_coords, rand_Idx, 0)

    print(len(cyto_coords))
    print(cyto_coords)
    print(len(ribo_centers))

    for center_point in ribo_centers:
        
        x_int = center_point[0]
        y_int = center_point[1]
        z_int = center_point[2]
        
        ribosome_centers[x_int,y_int,z_int] = True

        ribosomes[x_int+1,y_int,z_int] = True
        ribosomes[x_int-1,y_int,z_int] = True
        ribosomes[x_int,y_int+1,z_int] = True
        ribosomes[x_int,y_int-1,z_int] = True
        ribosomes[x_int,y_int,z_int+1] = True
        ribosomes[x_int,y_int,z_int-1] = True
    
    return ribosome_centers, ribosomes
#########################################################################################


#########################################################################################
def mapDNA(genome):
    """
    Inputs:
    genome - genome dictionary from sim_properties, contains all gene names, positions, and sequences
    
    Returns:
    DNAmap - Dictionary that defines all gene names, positions, and sequences
    chromosome_length - length of chromosome in bp
    
    Called by:
    buildRegions
    
    Description:
    Maps gene positions to 10 bp resolution for communication with the 10 bp resolution chromosome configurations
    Records the RNA and AA sequences for all genes for quick referencing when calculating costs
    """
    
    DNAmap = {}
    
    chromosome_length = int((genome.features[0].location.end+1)/10)
    
    ori_ter_rotation_factor = int(chromosome_length/2)
    
    print(ori_ter_rotation_factor)
    
    for feature in genome.features:
        
        strand = feature.strand
        
        if strand == 1:
        
            start = int(feature.location.start/10)
            
            if start<=ori_ter_rotation_factor:
                
                start = start + ori_ter_rotation_factor
                
            elif start>ori_ter_rotation_factor:
                
                start = start - ori_ter_rotation_factor
             
            end = int(feature.location.end/10)
            
            if end<=ori_ter_rotation_factor:
                
                end = end + ori_ter_rotation_factor
            
            elif end>ori_ter_rotation_factor:
                
                end = end - ori_ter_rotation_factor
            
        elif strand == -1:
            
            start = int(feature.location.end/10)
            
            if start<=ori_ter_rotation_factor:
                
                start = start + ori_ter_rotation_factor
                
            elif start>ori_ter_rotation_factor:
                
                start = start - ori_ter_rotation_factor
            
            end = int(feature.location.start/10)
            
            if end<=ori_ter_rotation_factor:
                
                end = end + ori_ter_rotation_factor
            
            elif end>ori_ter_rotation_factor:
                
                end = end - ori_ter_rotation_factor
        
        if feature.type == 'CDS':
            
            if('protein_id' in feature.qualifiers.keys()):
            
                locusTag = feature.qualifiers['locus_tag'][0]

                DNAmap[locusTag] = {}

                DNAmap[locusTag]['Type'] = 'protein'

                DNAmap[locusTag]['startIndex'] = [int(start)]
                
                DNAmap[locusTag]['originalStart'] = int(start)

                DNAmap[locusTag]['endIndex'] = [int(end)]
                
                DNAmap[locusTag]['originalEnd'] = int(end)

                DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())

                DNAmap[locusTag]['AAsequence'] = str(feature.location.extract(genome.seq).transcribe().translate(table=4))
                
                DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            
        elif feature.type == 'tRNA':
            
            locusTag = feature.qualifiers['locus_tag'][0]
            
            DNAmap[locusTag] = {}
            
            DNAmap[locusTag]['Type'] = 'tRNA'
            
            DNAmap[locusTag]['startIndex'] = [int(start)]
                
            DNAmap[locusTag]['originalStart'] = int(start)

            DNAmap[locusTag]['endIndex'] = [int(end)]

            DNAmap[locusTag]['originalEnd'] = int(end)
            
            DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
            
            DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            
        elif feature.type == 'rRNA':
            
            locusTag = feature.qualifiers['locus_tag'][0]
            
            DNAmap[locusTag] = {}
            
            DNAmap[locusTag]['Type'] = 'rRNA'
            
            DNAmap[locusTag]['startIndex'] = [int(start)]
                
            DNAmap[locusTag]['originalStart'] = int(start)

            DNAmap[locusTag]['endIndex'] = [int(end)]

            DNAmap[locusTag]['originalEnd'] = int(end)
            
            DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
            
            DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            
    return DNAmap, chromosome_length
#########################################################################################


#########################################################################################
def addChromosomeFeatures(sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    None
    
    Called by:
    buildRegions

    Description:
    Defines a dictionary that will be used to track chromosome features other than gene start and end sites 
    (e.g. the origin particle)
    """
    
    chromoFeatures = {}
    
    chromoFeatures['oriC'] = {}
    
    chromoFeatures['oriC']['index'] = int(sim_properties['genome_length']/2 + 2)
    
    sim_properties['chromosome_features'] = chromoFeatures
    
    return None
#########################################################################################


#########################################################################################
def createRibosomalProteinMap(genome):
    """
    Inputs:
    genome - genome dictionary from sim_properties, contains all gene names, positions, and sequences
    
    Returns:
    RiboPtnMap - dictionary defining attributes of ribosomal protein genes
    
    Called by:
    buildRegions
    
    Description:
    Writes dictionary that enables quick calling of ribosomal protein gene locus tags and names
    """

    RiboPtnMap = {}
    
    for feature in genome.features:
        
        if('protein_id' in feature.qualifiers.keys()):
            
            if 'S ribosomal' in feature.qualifiers['product'][0]:
                
                locusTag = feature.qualifiers['locus_tag'][0]
                
                RibosomalProteinID = feature.qualifiers['product'][0].split(' ')[-1]

                RiboPtnMap[RibosomalProteinID] = {}
                
                RiboPtnMap[RibosomalProteinID]['locusTag'] = locusTag
                
    return RiboPtnMap
#########################################################################################
