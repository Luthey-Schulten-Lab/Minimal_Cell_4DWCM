"""

Description: functions to analyze gene for CMEODE WCM
"""

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import pandas as pd 
from Bio.Seq import Seq

def mapDNA(gbfile_path):
    """
    Return: Multi-layer dictionary containing all the genome information
    """

    genome = next(SeqIO.parse(gbfile_path, "gb"))
    DNAmap = {}
    chromosome_length = int((genome.features[0].location.end+1)/10)
    ori_ter_rotation_factor = int(chromosome_length/2)
    # print(ori_ter_rotation_factor)   
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
                # Excluding 3 pseudo genes
                # Totoal coding genes are 455
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
            
        else:
            # 2 misc binding features, 2 ncRNA genes, 1 tm RNA gene, 
            try:
                locusTag = feature.qualifiers['locus_tag'][0]           
                DNAmap[locusTag] = {}          
                DNAmap[locusTag]['Type'] = feature.type
                DNAmap[locusTag]['startIndex'] = [int(start)]
                DNAmap[locusTag]['originalStart'] = int(start)
                DNAmap[locusTag]['endIndex'] = [int(end)]
                DNAmap[locusTag]['originalEnd'] = int(end)          
                DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
                DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            except:
                continue
            
    return DNAmap

def categorizeGenes(gbfile_path):
    """
    Input:
        gbfile_path: the path to the gene bank file
    Output:
        TypetoLocusNums: Python Dictionary where key is gene type, value is the list of locusNums
        geneTypes: six gene types
    Description: Categorize genes into protein, ncRNA, gene, rRNA, tRNA, tmRNA
    """
    
    genomeDict = mapDNA(gbfile_path)

    geneTypes = []
    geneNumbers = []
    TypetoLocusNums = {}

    for locusTag, locusDict in genomeDict.items():
        locusNum = locusTag.split('_')[1]

        if locusDict['Type'] not in geneTypes:
            geneTypes.append(locusDict['Type'])
            TypetoLocusNums[locusDict['Type']] = []

        TypetoLocusNums[locusDict['Type']].append(locusNum)
    
    for type, locusNums in TypetoLocusNums.items():
        geneNumbers.append(len(locusNums))

    # map(str, list) will convert the list's elements into strings.
    print('Six types of genes in Syn3A are {0} with respective numbers {1}.'.format(', '.join(geneTypes), ', '.join(map(str, geneNumbers))))

    return TypetoLocusNums, geneTypes


# Map locusNum to Gene sequence with intergeneic region
def initializeLocusNumtoGeneSeq(gbfile_path):
    """
    Called at the very beginning of the simulation once

    Description: Set up sub dictionary sim_properties['LocusNumtoIndex'] where the values are the genes' locusNums and keys are the start and end position of each gene
    """

    genome3A = next(SeqIO.parse(gbfile_path, "gb"))

    dna = genome3A
    gene_secondhalf = []
    gene_firsthalf = []
    LocusNumtoIndex = {}
    LocusNumtoGeneSeq = {}

    endposition = 543086
    # the index of the end position of 0910 is 543086
    # The length of the whole chromosome is len(dna.seq) is 543379
    position = 0

    gene_list = []
    for i in range(len(dna.features)):
        if ('product' in dna.features[i].qualifiers.keys()):
            #print(i) # This first statement works
            #print(dna.features[i].qualifiers['product'])
            if dna.features[i].qualifiers['product'][0]:# Figure out how to sort out for ribosomal operons?
                #print(dna.features[i].qualifiers['product'])
                gene_list.append(i)

    for gene in gene_list:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real
        if start < len(dna.seq)/2:
            gene_firsthalf.append(gene)
            if start == 0:
                LocusNumtoIndex[locusNum] = [position, end]
                position = end
            
            else:
                LocusNumtoIndex[locusNum] = [position, end]
                position = end
        else:
            gene_secondhalf.append(gene)

    position = endposition
    gene_secondhalf.reverse()

    for gene in gene_secondhalf:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real
        if end == endposition:

            LocusNumtoIndex[locusNum] = [start, position]
            position = start
            
        else:
            LocusNumtoIndex[locusNum] = [start, position]
            position = start

    for LocusNum, Index in LocusNumtoIndex.items():
        start, end = Index[0], Index[1]
        gene_seq = Seq(str(dna.seq[start:end]))
        LocusNumtoGeneSeq[LocusNum] = gene_seq
    

    return LocusNumtoIndex, LocusNumtoGeneSeq
 
def riboTolocusMap(genome):
    """
    Input; genome of syn3A
    
    Description: Map ribosome protein names to LocusTags e.g. 'L27' to 'JCVISYN3A_0499'
        Change r-protein name L7/L12 to L7L12
    """

    ribotoLocus = {}

    for locusTag, subdic in genome.items():
        if subdic['Type'] != 'gene':
            GeneName = subdic['GeneName']
            if GeneName.find('ribosomal protein') > 0:
                ribo = GeneName.split(' ')[-1]
                ribo = ribo.replace('/','')
                ribotoLocus[ribo] = []

    for locusTag, subdic in genome.items():
        if subdic['Type'] != 'gene':
            GeneName = subdic['GeneName']
            if GeneName.find('ribosomal protein') > 0:
                ribo = GeneName.split(' ')[-1]
                ribo = ribo.replace('/','')

                ribotoLocus[ribo].append(locusTag)

    ribotoLocus['L27'] = ['JCVISYN3A_0499']

    rPtn_locusTags = [locusTag for locusTags in ribotoLocus.values() for locusTag in locusTags]

    return ribotoLocus, rPtn_locusTags

def locusToriboMap(ribotoLocus):

    locusTagToribo = {}

    for ribo, locusTags in ribotoLocus.items():
        for locusTag in locusTags:
            locusTagToribo[locusTag] = ribo
        

    return locusTagToribo

def get_tRNAMap(genome):
    """
    Set up the dictionary aa_tRNA with amino acids and correspoding tRNAs {'LEU': ['R_0070', 'R_0423', 'R_0506'],...}
    tRNA_aa: {'R_0070':LEU, }
    """
    aa_tRNA = {}

    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == "tRNA":
            
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum

            rnaName = locusDict['GeneName'].split('-')[1].upper()
            
            if rnaName not in aa_tRNA:
                
                aa_tRNA[rnaName] = [rnaID]
                
            else:
                
                aa_tRNA[rnaName].append(rnaID)
    
    tRNA_aa = {}

    for aa, tRNAs in aa_tRNA.items():
        for tRNA in tRNAs:
            tRNA_aa[tRNA] = aa


    return aa_tRNA, tRNA_aa

def get_aa_synthetase(kinetic_params_path):
    """
    Get the dictionary mapping amino acids to protein synthetase
    """
    tRNA_params = pd.read_excel(kinetic_params_path, sheet_name='tRNA Charging')

    aa_synthetase = {}

    for index, row in tRNA_params.iterrows():
        if row['Parameter Type'] == 'synthetase':
            AA = row['Reaction Name'][:3]
            aa_synthetase[AA] = row['Value']

    return aa_synthetase

def get_aa_tRNAcharging_intermediates(AA, gbfile_path, kinetic_params_path):
       """
       AA: Upper case of AA, e.g. ASP, ILE
       Return: list of synthetase, atp_synthetase, atp_synthetase_aa, atp_synthetase_tRNA, free_tRNA, charged_tRNA
       """

       genome = mapDNA(gbfile_path)

       aa_tRNA , tRNA_aa = get_tRNAMap(genome)

       aa_synthetase = get_aa_synthetase(kinetic_params_path)

       synthetase = aa_synthetase[AA]; synthetase_atp = synthetase + '_atp'; synthetase_atp_aa = synthetase_atp + '_aa'

       aa_intermediates = [synthetase, synthetase_atp, synthetase_atp_aa]

       for tRNA in aa_tRNA[AA]:
           synthetaseTrnaID = synthetase_atp_aa + '_' + tRNA

           charged_tRNA = tRNA + '_ch'

           aa_intermediates.extend([tRNA, synthetaseTrnaID, charged_tRNA ])

       return aa_intermediates

def get_aa_transport(gbfile_path, kinetic_params_path):
    """
    return: Map the aa names to transport reaction names 
    """

    genome = mapDNA(gbfile_path)
    aa_tRNA, tRNA_aa = get_tRNAMap(genome)

    transport_params = pd.read_excel(kinetic_params_path, sheet_name='Transport')
    aa_names = [aa.lower() for aa in aa_tRNA.keys()]

    free_aas = ['M_{0}__L_c'.format(aa) if aa != 'gly' else 'M_gly_c' for aa in aa_names]

    aa_transport = {}; 

    for free_aa in free_aas:
        aa_transport[free_aa] = []
        for index, row in transport_params.iterrows():
            if row['Related Species'] == free_aa:
                aa_transport[free_aa].append(row['Reaction Name'])

    return aa_transport


def get_filamentlength_new(w):
    """
    Input: 
    
    Description: Get the trajectories of filament length of dnaA on mother and two daughter chromosomes for all replicates
    """

    time = w.get_t()
    
    reps = np.arange(1, w.N_reps+1)
    
    filas_mother = np.zeros((len(time), len(reps)))
    filas_d1 = np.zeros((len(time), len(reps)))
    filas_d2 = np.zeros((len(time), len(reps)))


    for i_rep in range(len(reps)):
        fila_mother = np.zeros(len(time))
        fila_d1 = np.zeros(len(time))
        fila_d2 = np.zeros(len(time))
        
        rep = reps[i_rep]

        for i in range(1,31):
            unbound_specie = 'ssdnaAFila_' + str(i)
            unbound2_specie = 'ssdnaAFila2_' + str(i)
            unbound3_specie = 'ssdnaAFila3_' + str(i)
            fila_mother += w.get_specie_trace_single_rep(unbound_specie, rep)*i
            fila_d1 += w.get_specie_trace_single_rep(unbound2_specie, rep)*i
            fila_d2 += w.get_specie_trace_single_rep(unbound3_specie, rep)*i
        
        filas_mother[:,i_rep] = fila_mother
        filas_d1[:,i_rep] = fila_d1
        filas_d2[:,i_rep] = fila_d2

    return filas_mother, filas_d1, filas_d2


def get_filamentlength(w):
    """
    Input: 

    Description: Get the trajectories of filament length of dnaA on mother and two daughter chromosomes for all replicates
    """

    time = w.get_t()
    
    reps = np.arange(1, w.N_reps+1)
    
    filas_mother = np.zeros((len(time), len(reps)))
    filas_d1 = np.zeros((len(time), len(reps)))
    filas_d2 = np.zeros((len(time), len(reps)))


    for i_rep in range(len(reps)):
        fila_mother = np.zeros(len(time))
        fila_d1 = np.zeros(len(time))
        fila_d2 = np.zeros(len(time))
        
        rep = reps[i_rep]

        for i in range(1,31):
            unbound_specie = 'ssDNAunboundSite_' + str(i)
            unbound2_specie = 'ssDNAunboundSite2_' + str(i)
            unbound3_specie = 'ssDNAunboundSite3_' + str(i)
            fila_mother += w.get_specie_trace_single_rep(unbound_specie, rep)*i
            fila_d1 += w.get_specie_trace_single_rep(unbound2_specie, rep)*i
            fila_d2 += w.get_specie_trace_single_rep(unbound3_specie, rep)*i
        
        filas_mother[:,i_rep] = fila_mother
        filas_d1[:,i_rep] = fila_d1
        filas_d2[:,i_rep] = fila_d2

    return filas_mother, filas_d1, filas_d2

def get_filamentlength_mother(w):
    """
    Input: 

    Description: Get the trajectories of filament length of dnaA on mother and two daughter chromosomes for all replicates
    """

    time = w.get_t()
    
    reps = np.arange(1, w.N_reps+1)
    
    filas_mother = np.zeros((len(time), len(reps)))


    for i_rep in range(len(reps)):
        fila_mother = np.zeros(len(time))
        
        rep = reps[i_rep]

        for i in range(1,31):
            unbound_specie = 'ssDNAunboundSite_' + str(i)
            unbound2_specie = 'ssDNAunboundSite2_' + str(i)
            unbound3_specie = 'ssDNAunboundSite3_' + str(i)
            fila_mother += w.get_specie_trace_single_rep(unbound_specie, rep)*i
        
        filas_mother[:,i_rep] = fila_mother

    return filas_mother


def analyze_initiation(w, reps):
    """
    Input: w
            reps: serial numbers of replicates, e.g. [1,4,9,..,12]
    Return:
            ini_rounds: [1,0,2,..,3]
            ini_times: [[800], [], [700, 5000], ..., [500, 4000, 8000]]
    Description:
    """
    
    reps = [rep - 1 for rep in reps]

    rep_num = len(reps)

    rep_str = ', '.join(str(rep) for rep in reps)

    print('Following is the initiation analysis of {0} replicates {1}.'.format(rep_num, rep_str))

    ini_check = w.get_specie_trace('RepInitCheck')
    
    ini_rounds = []
    
    ini_times = []

    for rep in reps:

        ini_check_rep = ini_check[:,rep]

        diff = ini_check_rep[1:] - ini_check_rep[:-1]
        
        locs = np.where(diff == 1)[0]
        
        ini_rounds.append(len(locs))

        ini_times.append(locs)

        time_str = ''

        for loc in locs:
            time_str = time_str + str(loc/60) + ','

        
        print('Replicate {0} finished {1} rounds of replication initiation at {2} minutes'.format(rep, len(locs), time_str))
    
    print('Average rounds of initiation is {0:.2f} over the whole cell cycle'.format(sum(ini_rounds)/len(ini_rounds)))
    
    print('Minimum rounds of replication initiation is {0}, maximum rounds is {1}'.format(min(ini_rounds), max(ini_rounds)))

    print('Distribution of initiation: 0 round {0:.2f}, 1 round {1:.2f}, 2 rounds {2:.2f}, 3 rounds {3:.2f}, 4 rounds {4:.2f}, 5 rounds {5:.2f}'
          .format(ini_rounds.count(0)/len(ini_rounds), ini_rounds.count(1)/len(ini_rounds), 
                  ini_rounds.count(2)/len(ini_rounds), ini_rounds.count(3)/len(ini_rounds),
                  ini_rounds.count(4)/len(ini_rounds),ini_rounds.count(5)/len(ini_rounds)))
    
    ini_mother_times, ini_daughter1_times, ini_daughter2_times = convert_ini_times(ini_rounds, ini_times)

    # print('Average time to finish first, second, and third round initiation are {0}, {1} and {2} seconds'.format(sum(ini_mother_times)/len(ini_mother_times), 
    #                                             sum(ini_daughter1_times)/len(ini_daughter1_times), sum(ini_daughter2_times)/len(ini_daughter2_times)))

    return ini_rounds, ini_times, ini_mother_times, ini_daughter1_times, ini_daughter2_times


def convert_ini_times(ini_rounds, ini_times):
    """
    
    Description: Convert replication initiation information into a dictionary
    """
    
    max_ini = max(ini_rounds)

    ini_mother_times = []
    ini_daughter1_times = []
    ini_daughter2_times = []

    for ini_time in ini_times:
        try:
            ini_mother_times.append(ini_time[0])
            ini_daughter1_times.append(ini_time[1])
            ini_daughter2_times.append(ini_time[2])
        except:
            None
            # print('Not enough initiation')
    
    return ini_mother_times, ini_daughter1_times, ini_daughter2_times



def analyze_filament(w, reps):
    
    filas_mother, fila_d1, filas_d2 = get_filamentlength(w)



    return None


