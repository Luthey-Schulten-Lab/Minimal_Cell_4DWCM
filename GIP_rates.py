"""
Authors: Zane Thornburg

Rate constant caluculations for genetic information processing reactions
"""

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np

#################################################################################################

# FIX PART TO MM#

#################################################################################################


##################
# Define how to calculate transcription rate constants for transcription reactions.
# Uses mature transcript length and proteomics for promoter strength.


#########################################################################################
def TranscriptionRate(sim_properties, locusTag, rnasequence):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
        
    proxyPromoterStrength = sim_properties['promoters'][locusTag]
    
    kcat_mod = min(rnaPolKcat*(proxyPromoterStrength/180),85)

    kcat_mod = max(10,kcat_mod)

    # Add total number of monomers to parameter dict
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    pmap = sim_properties['counts']
    
    atp = partTomM(max(1,pmap['M_atp_c']), sim_properties)
    ctp = partTomM(max(1,pmap['M_ctp_c']), sim_properties)
    gtp = partTomM(max(1,pmap['M_gtp_c']), sim_properties)
    utp = partTomM(max(1,pmap['M_utp_c']), sim_properties)
    
    NMonoSum = NMono_A*rnaPolKd/atp + NMono_C*rnaPolKd/ctp + NMono_U*rnaPolKd/utp + NMono_G*rnaPolKd/gtp
    
    k_transcription = kcat_mod / ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
#########################################################################################


#########################################################################################
def RNAP_binding(sim, sim_properties, locusTag):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    proxyPromoterStrength = sim_properties['promoters'][locusTag]
    
    binding_rate = 10*(180/765)*(proxyPromoterStrength/180)*Ecoli_V*avgdr/11400/60
    #4
    
    return binding_rate
#########################################################################################


#########################################################################################
def TranslationRate(aasequence):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    # Add translation reaction
    
    # Considers amino acids up to the first stop codon.
#     aasequence = aasequence[0:aasequence.find("*")]
    
    # Check that we know all residues used in the sequence
#     if ( set(aasequence) - set(aaMap.keys()) ):
#         raise Exception("Unknown residue(s) in Protein sequence {}".format(set(aasequence) - set(aaMap.keys())) )
    
    # Count how many times each residue is used
    aaCount = defaultdict(int)
    for aa in set(aasequence):
        aaCount[aa] = aasequence.count(aa)
    
    NMono_A = aaCount["A"]
    NMono_R = aaCount["R"]
    NMono_N = aaCount["N"]
    NMono_D = aaCount["D"]
    NMono_C = aaCount["C"]
    NMono_E = aaCount["E"]
    NMono_Q = aaCount["Q"]
    NMono_H = aaCount["H"]
    NMono_I = aaCount["I"]
    NMono_L = aaCount["L"]
    NMono_K = aaCount["K"]
    NMono_M = aaCount["M"]
    NMono_P = aaCount["P"]
    NMono_S = aaCount["S"]
    NMono_T = aaCount["T"]
    NMono_W = aaCount["W"]
    NMono_Y = aaCount["Y"]
    NMono_G = aaCount["G"]
    NMono_F = aaCount["F"]
    NMono_V = aaCount["V"]
    
    NStop = aaCount["*"]
    
    if NStop > 1:
        print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")
    
    NMonoDict = [NMono_A,NMono_R,NMono_N,NMono_D,NMono_C,NMono_E,NMono_Q,NMono_H,
                 NMono_I,NMono_L,NMono_K,NMono_M,NMono_P,NMono_S,NMono_T,NMono_W,
                 NMono_Y,NMono_G,NMono_F,NMono_V]
    
    NMonoSum = 0
    
    for nmono in range(0,len(NMonoDict)):
        NMonoSum = NMonoSum + NMonoDict[nmono]*riboKd/ctRNAconc
        
    n_tot = sum(list(aaCount.values()))
    
    kcat_mod = riboKcat #*ribo_num #*0.4
    
    k_translation = kcat_mod / ((riboKd**2)/(ctRNAconc**2) + NMonoSum + n_tot - 1)
    
    return k_translation
#########################################################################################


#########################################################################################
def mrnaDegradationRate(rnasequence):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    kcat = 88 #1/s

    n_tot = sum(list(baseCount.values()))

    k_deg = kcat / n_tot 
    
    return k_deg
#########################################################################################


#########################################################################################
def TranslocationRate(aasequence):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    # Count how many times each residue is used
    aaCount = defaultdict(int)
    for aa in set(aasequence):
        aaCount[aa] = aasequence.count(aa)
    
    ptnLen = sum(list(aaCount.values()))
    
    k_transloc = 50/ptnLen #secyNum*
    
    return k_transloc
#########################################################################################


#########################################################################################
def ReplicationRate(sim_properties, dnasequence):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(dnasequence):
        baseCount[base] = dnasequence.count(base)
        
#     proxyPromoterStrength = sim_properties['promoters'][locusTag]
    
    kcat = 100 * 4 # bp/4*s

    # Add total number of monomers to parameter dict
    
#     CMono1 = baseMap[ rnasequence[0] ]
    
#     CMono2 = baseMap[ rnasequence[1] ]

    n_tot = sum(list(baseCount.values()))
    
    if n_tot == 0:
        
        return 0

    NMono_A = baseCount["A"]
    
    NMono_T = baseCount["T"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    pmap = sim_properties['counts']
    
    datp = partTomM(max(1,pmap['M_datp_c']), sim_properties)
    dctp = partTomM(max(1,pmap['M_dctp_c']), sim_properties)
    dgtp = partTomM(max(1,pmap['M_dgtp_c']), sim_properties)
    dttp = partTomM(max(1,pmap['M_dttp_c']), sim_properties)
    
    NMonoSum = NMono_A*dnaPolKd/datp + NMono_C*dnaPolKd/dctp + NMono_T*dnaPolKd/dttp + NMono_G*dnaPolKd/dgtp
    
    base_rep_rate = kcat / n_tot
    
    k_replication = kcat / (NMonoSum + n_tot - 1)
    
    rep_rate = int((k_replication / base_rep_rate) * (kcat / 10))
    
    return rep_rate
#########################################################################################



baseMap = OrderedDict({ "A":"M_atp_c", "U":"M_utp_c", "G":"M_gtp_c", "C":"M_ctp_c" })
# baseMapToMonoP = OrderedDict({ "A":"M_amp_c", "U":"M_ump_c", "G":"M_gmp_c", "C":"M_cmp_c" })

dnaPolKd = 0.01 #mM

# Global parameters for transcription
rnaPolKcat = 20 # nt/s
rnaPolK0 = 1e-4 #mM
rnaPolKd = 0.1 #mM

rrnaPolKcat = 85 # nt/s

krnadeg = 0.00578/2 # 1/s
# rna_deg_rate = sim.rateConst('RNAdeg', krnadeg, 2)

ptnDegRate = 7.70e-06 # 1/s

ATPconc = 4 #mM
UTPconc = 3 #mM
CTPconc = 1 #mM
GTPconc = 2 #mM

# Cell radius (meters):
# r_cell = 2.5*(10**-7)
r_cell = 2.0*(10**-7) # m

CytoVolume = (4*np.pi/3)*1000*r_cell**3 # L
cellVolume = CytoVolume

subvolume_vol = 1000*(8e-9)**3 # L

# print(cellVolume)

# Avogadro:
avgdr   = 6.022e23 # molec/mol
Avognum = avgdr

NaV = Avognum * subvolume_vol

countToMiliMol = 1000/(avgdr*cellVolume)

RnaPconc = 187*countToMiliMol # mM

# Global parameter for degradation of mRNAs
rnaDegRate = 0.00578/2 # 1/s

# degrad_bind_rate = 11/60/countToMiliMol*1000 #1/M/s
Ecoli_V = 1e-15 #L
degrad_bind_rate = 11*avgdr*Ecoli_V/60/2400 #1/M/s 7800

# Create a map for rna sequence to NTP concentration.
baseMap = OrderedDict({ "A":ATPconc, "U":UTPconc, "G":GTPconc, "C":CTPconc })

# Create Dictionaries to map tRNAs to associated aa abbreviations in protein sequences.
aaMap = OrderedDict({"A":"M_ala__L_c", "R":"M_arg__L_c", 
    "N":"M_asn__L_c", "D":"M_asp__L_c", "C":"M_cys__L_c", "E":"M_glu__L_c", "Q":"M_gln__L_c", "G":"M_gly_c", 
    "H":"M_his__L_c", "I":"M_ile__L_c", "L":"M_leu__L_c", "K":"M_lys__L_c", "M":"M_met__L_c", "F":"M_phe__L_c", 
    "P":"M_pro__L_c", "S":"M_ser__L_c", "T":"M_thr__L_c", "W":"M_trp__L_c", "Y":"M_tyr__L_c", "V":"M_val__L_c",
    "*":"Stop_Codon"})

aaTRNAMap = OrderedDict({"A":"M_alatrna_c", "R":"M_argtrna_c", 
    "N":"M_asntrna_c", "D":"M_asptrna_c", "C":"M_cystrna_c", "E":"M_glutrna_c", "Q":"M_glntrna_c", "G":"M_glytrna_c", 
    "H":"M_histrna_c", "I":"M_iletrna_c", "L":"M_leutrna_c", "K":"M_lystrna_c", "M":"M_mettrna_c", "F":"M_phetrna_c", 
    "P":"M_protrna_c", "S":"M_sertrna_c", "T":"M_thrtrna_c", "W":"M_trptrna_c", "Y":"M_tyrtrna_c", "V":"M_valtrna_c"})

aaTRNAFreeMap = OrderedDict({"A":"M_trnaala_c", "R":"M_trnaarg_c", 
    "N":"M_trnaasn_c", "D":"M_trnaasp_c", "C":"M_trnacys_c", "E":"M_trnaglu_c", "Q":"M_trnagln_c", "G":"M_trnagly_c", 
    "H":"M_trnahis_c", "I":"M_trnaile_c", "L":"M_trnaleu_c", "K":"M_trnalys_c", "M":"M_trnamet_c", "F":"M_trnaphe_c", 
    "P":"M_trnapro_c", "S":"M_trnaser_c", "T":"M_trnathr_c", "W":"M_trnatrp_c", "Y":"M_trnatyr_c", "V":"M_trnaval_c"})


defaultPtnCount = 1

# Global parameters for translation
riboKcat = 12 # 1/s
riboK0 = 4*25e-6 # mM
riboKd = 0.001 # mM

# ribo_init = 30*Ecoli_V*avgdr/60/6800

ribosomeConc = 500*countToMiliMol # mM

# Concentration of charged tRNA
ctRNAconc = 200*countToMiliMol # mM


#########################################################################################
def partTomM(particles, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's

    conc = (particles*1000.0)/(NA*sim_properties['volume_L'])

    return conc
#########################################################################################
# Inputs:
# Returns:
# Called by:
# Description:
#########################################################################################

