"""
Authors: Zane Thornburg

General functions for defining the diffusion rules for all macromolecules.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict

import numpy as np


#########################################################################################
def defaultDiffusion(sim, region_dict):
    """
    Inputs:
    sim - jLM RDME simulation object
    region_dict - Dictionary of regions on the RDME site lattice and their shapes
    
    Returns:
    None
    
    Called by:
    importInitialConditions.py - initializeParticles
    
    Description:
    Sets the default diffusion to zero among all regions for all particles until otherwise defined later
    """
    
    for region1, rDict1 in region_dict.items():
        
        for region2, rDict2 in region_dict.items():
        
            sim.transitionRate(None, sim.region(region1), sim.region(region2), sim.diffusionZero)
    
    return None
#########################################################################################


#########################################################################################
def generalDiffusionConstants(sim):
    """
    Inputs:
    sim - jLM RDME simulation object
    
    Returns:
    None
    
    Called by:
    
    Description:
    Defines general diffusion coefficient values that are used throughout the RDME simulation
    """
    
#     diffPtnCyt = sim.diffusionConst("diffPtnCyt", 1e-12)
# diffPtnDNA = sim.diffusionConst("diffPtnDna", 0.5e-12)
    
    sim.diffusionConst("diffPtn", 1e-12)
    sim.diffusionConst("diffPtnDna", 0.5e-12)
    sim.diffusionConst("diffMemPtn", 1e-13)
    sim.diffusionConst("diffDeg", 0.031e-12)
    
    return None
#########################################################################################


#########################################################################################
def proteinDiffusionCyto(sim, ptnID):
    """
    Inputs:
    sim - jLM RDME simulation object
    ptnID - (string) Name of protein particle in the RDME simulation
    
    Returns:
    None
    
    Called by:
    
    Description:
    Sets diffusion coefficients for cytoplasmic proteins
    """
    
    # cytoplasmic protein diffusion
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('ribo_centers'), sim.dc.diffPtn)

    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('ribosomes'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('cytoplasm'), sim.dc.diffPtn)

    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)

    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('DNA'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('DNA'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('DNA'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('DNA'), sim.dc.diffPtnDna)
    
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('ribosomes'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('ribosomes'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('ribosomes'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('ribosomes'), sim.dc.diffPtnDna)
    
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('DNA'), sim.dc.diffPtnDna)
    
    return None
#########################################################################################


#########################################################################################
def proteinDiffusionMem(sim, ptnID):
    """
    Inputs:
    sim - jLM RDME simulation object
    ptnID - (string) Name of protein particle in the RDME simulation
    
    Returns:
    None
    
    Called by:
    
    Description:
    Sets diffusion coefficients for transmembrane proteins
    """
    
    # membrane protein diffusion
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffMemPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('membrane'), sim.dc.diffMemPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('membrane'), sim.region('membrane'), sim.dc.diffMemPtn)
    
    # To make sure membrane proteins make it back to membrane if they leave
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('membrane'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('membrane'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('membrane'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('membrane'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('DNA'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('DNA'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('DNA'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('ribosomes'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('cytoplasm'), sim.dc.diffPtn)    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('DNA'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('DNA'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('ribo_centers'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('ribosomes'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)

    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)

    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('DNA'), sim.dc.diffPtn)
    
    return None
#########################################################################################


#########################################################################################
def proteinDiffusionPeriph(sim, ptnID):
    """
    Inputs:
    sim - jLM RDME simulation object
    ptnID - (string) Name of protein particle in the RDME simulation
    
    Returns:
    None
    
    Called by:
    
    Description:
    Sets diffusion coefficients for proteins that are restricted to the peripheral membrane region (outer_cytoplasm)
    """
    
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffMemPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('DNA'), sim.dc.diffPtnDna)
    
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('DNA'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('DNA'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)

    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('outer_cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('cytoplasm'), sim.dc.diffPtn)

    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('DNA'), sim.dc.diffPtnDna)

    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('DNA'), sim.dc.diffPtnDna)
    
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('ribo_centers'), sim.dc.diffPtn)

    sim.transitionRate(sim.species(ptnID), sim.region('ribo_centers'), sim.region('ribosomes'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(ptnID), sim.region('cytoplasm'), sim.region('ribosomes'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('outer_cytoplasm'), sim.region('ribosomes'), sim.dc.diffPtnDna)
    sim.transitionRate(sim.species(ptnID), sim.region('ribosomes'), sim.region('ribosomes'), sim.dc.diffPtnDna)

    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('ribo_centers'), sim.dc.diffPtnDna)

    sim.transitionRate(sim.species(ptnID), sim.region('DNA'), sim.region('ribosomes'), sim.dc.diffPtnDna)
    
    return None
#########################################################################################


#########################################################################################
def degradosomeDiffusion(sim, degID):
    """
    Inputs:
    sim - jLM RDME simulation object
    degID - (string) Name of degradosome particle (or bound state) in the RDME simulation
    
    Returns:
    None
    
    Called by:
    
    Description:
    Sets diffusion coefficients for degradosome particle or mRNA bound state
    """
    
    sim.transitionRate(sim.species(degID), sim.region('outer_cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffDeg)
    
    sim.transitionRate(sim.species(degID), sim.region('cytoplasm'), sim.region('outer_cytoplasm'), sim.dc.diffDeg)
    sim.transitionRate(sim.species(degID), sim.region('DNA'), sim.region('outer_cytoplasm'), sim.dc.diffDeg)
    
    sim.transitionRate(sim.species(degID), sim.region('ribosomes'), sim.region('outer_cytoplasm'), sim.dc.diffDeg)

    sim.transitionRate(sim.species(degID), sim.region('ribo_centers'), sim.region('outer_cytoplasm'), sim.dc.diffDeg)

    sim.transitionRate(sim.species(degID), sim.region('outer_cytoplasm'), sim.region('DNA'), sim.dc.diffDeg)
    
    sim.transitionRate(sim.species(degID), sim.region('cytoplasm'), sim.region('cytoplasm'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(degID), sim.region('DNA'), sim.region('DNA'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(degID), sim.region('cytoplasm'), sim.region('DNA'), sim.dc.diffPtn)
    sim.transitionRate(sim.species(degID), sim.region('DNA'), sim.region('cytoplasm'), sim.dc.diffPtn)
    
    sim.transitionRate(sim.species(degID), sim.region('ribosomes'), sim.region('ribosomes'), sim.dc.diffPtnDna)
    
    sim.transitionRate(sim.species(degID), sim.region('ribosomes'), sim.region('cytoplasm'), sim.dc.diffDeg)
    sim.transitionRate(sim.species(degID), sim.region('ribo_centers'), sim.region('cytoplasm'), sim.dc.diffDeg)    
    sim.transitionRate(sim.species(degID), sim.region('ribosomes'), sim.region('DNA'), sim.dc.diffDeg)
    sim.transitionRate(sim.species(degID), sim.region('ribo_centers'), sim.region('DNA'), sim.dc.diffDeg)
    sim.transitionRate(sim.species(degID), sim.region('ribosomes'), sim.region('ribo_centers'), sim.dc.diffDeg)
    sim.transitionRate(sim.species(degID), sim.region('ribo_centers'), sim.region('ribosomes'), sim.dc.diffDeg)
    
    return None
#########################################################################################

    
#########################################################################################
def rnaDiffusion(sim, RNAID, RNAseq):
    """
    Inputs:
    sim - jLM RDME simulation object
    RNAID - (string) Name of RNA particle in the RDME simulation (e.g. R_0001)
    RNAseq - RNA sequence corresponding to the particle with name RNAID, capitalized (AUCG)
    
    Returns:
    None
    
    Called by:
    
    Description:
    Sets the diffusion coefficients of an RNA particle with name RNAID
    """
    
    diffRNA = sim.diffusionConst(RNAID + '_diff', RNA_diff_coeff(RNAseq))
    diffRNADNA = sim.diffusionConst(RNAID + '_diffDna', RNA_diff_coeff(RNAseq)/2)
    
        # mRNA diffusion
    sim.transitionRate(sim.species(RNAID), sim.region('ribosomes'), sim.region('ribo_centers'), diffRNA)

    sim.transitionRate(sim.species(RNAID), sim.region('ribo_centers'), sim.region('ribosomes'), diffRNA)

    sim.transitionRate(sim.species(RNAID), sim.region('ribosomes'), sim.region('ribosomes'), diffRNA)
    
    sim.transitionRate(sim.species(RNAID), sim.region('ribosomes'), sim.region('cytoplasm'), diffRNA)

    sim.transitionRate(sim.species(RNAID), sim.region('ribo_centers'), sim.region('cytoplasm'), diffRNA)
    
    sim.transitionRate(sim.species(RNAID), sim.region('cytoplasm'), sim.region('ribo_centers'), diffRNA)

    sim.transitionRate(sim.species(RNAID), sim.region('cytoplasm'), sim.region('ribosomes'), diffRNA)
    
    sim.transitionRate(sim.species(RNAID), sim.region('ribosomes'), sim.region('outer_cytoplasm'), diffRNA)

    sim.transitionRate(sim.species(RNAID), sim.region('ribo_centers'), sim.region('outer_cytoplasm'), diffRNA)
    
    sim.transitionRate(sim.species(RNAID), sim.region('outer_cytoplasm'), sim.region('ribo_centers'), diffRNA)

    sim.transitionRate(sim.species(RNAID), sim.region('outer_cytoplasm'), sim.region('ribosomes'), diffRNA)
    
    sim.transitionRate(sim.species(RNAID), sim.region('cytoplasm'), sim.region('cytoplasm'), diffRNA)
    sim.transitionRate(sim.species(RNAID), sim.region('outer_cytoplasm'), sim.region('outer_cytoplasm'), diffRNA)
    sim.transitionRate(sim.species(RNAID), sim.region('cytoplasm'), sim.region('outer_cytoplasm'), diffRNA)
    sim.transitionRate(sim.species(RNAID), sim.region('outer_cytoplasm'), sim.region('cytoplasm'), diffRNA)
    
    sim.transitionRate(sim.species(RNAID), sim.region('cytoplasm'), sim.region('DNA'), diffRNADNA)
    sim.transitionRate(sim.species(RNAID), sim.region('outer_cytoplasm'), sim.region('DNA'), diffRNADNA)
    sim.transitionRate(sim.species(RNAID), sim.region('ribo_centers'), sim.region('DNA'), diffRNADNA)
    sim.transitionRate(sim.species(RNAID), sim.region('ribosomes'), sim.region('DNA'), diffRNADNA)
    
    sim.transitionRate(sim.species(RNAID), sim.region('DNA'), sim.region('cytoplasm'), diffRNA)
    sim.transitionRate(sim.species(RNAID), sim.region('DNA'), sim.region('outer_cytoplasm'), diffRNA)
    sim.transitionRate(sim.species(RNAID), sim.region('DNA'), sim.region('ribo_centers'), diffRNA)

    sim.transitionRate(sim.species(RNAID), sim.region('DNA'), sim.region('ribosomes'), diffRNA)
    
    sim.transitionRate(sim.species(RNAID), sim.region('DNA'), sim.region('DNA'), diffRNADNA)
    
    return None
#########################################################################################


#########################################################################################  
def ribosomeDiffusion(sim, riboID):
    """
    Inputs:
    sim - jLM RDME simulation object
    riboID - (string) Name of ribosome particle in the RDME simulation (e.g. ribosomeP, RB_0001)
    
    Returns:
    None
    
    Called by:
    
    Description:
    Sets diffusion coefficients for ribosome particles
    """
    
#     diffPtn = sim.diffusionConst("ptn_diff", 1e-12)
#     diffRNA = sim.diffusionConst("mRNA_diff",0.01e-12)
    diffRibo = sim.diffusionConst("ribo_diff", 0.001e-12)
    diffTRibo = sim.diffusionConst("tribo_diff", 0.0001e-12)
    
    sim.transitionRate(sim.species(riboID), sim.region('ribosomes'), sim.region('ribo_centers'), diffRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_starts'), sim.region('p_start_c'), diffRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_mids'), sim.region('p_mid_c'), diffRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_ends'), sim.region('p_end_c'), diffRibo)

    sim.transitionRate(sim.species(riboID), sim.region('ribo_centers'), sim.region('ribosomes'), diffRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_start_c'), sim.region('p_starts'), diffRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_mid_c'), sim.region('p_mids'), diffRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_end_c'), sim.region('p_ends'), diffRibo)
    
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('ribo_centers'), sim.region('DNA'), diffTRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_start_c'), sim.region('DNA'), diffTRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_mid_c'), sim.region('DNA'), diffTRibo)
#     sim.transitionRate(sim.species('ribosomeP'), sim.region('p_end_c'), sim.region('DNA'), diffTRibo)


# # sim.transitionRate(sim.species('ribosomeP'), sim.region('ribosomes'), sim.region('ribo_centers'), diffRibo)

#     sim.transitionRate(sim.species('ribosomeP'), sim.region('ribo_centers'), sim.region('ribosomes'), diffRibo)
    
    return None
#########################################################################################

    
#########################################################################################
def RNA_diff_coeff(rnasequence):
    """
    Inputs:
    rnasequence - (string) nucleotide sequence of an RNA, capitalized (AUCG)
    
    Returns:
    mrna_diff_coeff - estimated diffusion coeficient of RNA with input sequence rnasequence (micron^2/s)
    
    Called by:
    
    Description:
    Calculates estimated diffusion coefficient of input rnasequence by estimating its hydrodynamic radius
    """
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
        
    n_tot = sum(list(baseCount.values()))
    
    N_A = baseCount["A"]

    N_U = baseCount["U"]

    N_C = baseCount["C"]

    N_G = baseCount["G"]
    
    molec_mass = 337 #309 #g/mol/nucleotide
    density = 1.75*1000000 #g/m^3
    N_A = 6.023e23 #mol^-1
    
    R_H = ((3*molec_mass*n_tot)/(4*np.pi*N_A*density))**(1/3)

    visc = 1.17 #0.15 #7.1 #17.5 #0.05 #0.001 #Pa*s
    kB = 1.380e-23
    Temp = 310
    
    mrna_diff_coeff = kB*Temp/(6*np.pi*visc*R_H)
#     print(mrna_diff_coeff)
    
    return mrna_diff_coeff
#########################################################################################

