"""

Description: functions to analyze ptns for CMEODE WCM
"""
import pandas as pd
import numpy as np
import WCM_gene as gene
import WCM_math as math

def gettRNAmap(genome):
    """
    Input: genome

    Description: 
    Set up the sub dictionary with amino acids and correspoding tRNAs {'LEU': ['R_0070', 'R_0423', 'R_0506'],...}
    """

    tRNA_map = {}

    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == "tRNA":
            
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum

            rnaName = locusDict['GeneName'].split('-')[1].upper()
            
            if rnaName not in tRNA_map:
                
                tRNA_map[rnaName] = [rnaID]
                
            else:
                
                tRNA_map[rnaName].append(rnaID)
                

    return tRNA_map


def getSynthetaseToAA(kineticfile_dir):
    synthetaseToAA = {}

    tRNA_params = pd.read_excel(kineticfile_dir, sheet_name='tRNA Charging')

    for index, row in tRNA_params.iterrows():
        if row['Parameter Type'] == 'synthetase':
            synthetase = row['Value']
            AA = row['Reaction Name'][0:3]
            synthetaseToAA[synthetase] = AA

    return synthetaseToAA


def gettRNALocusNums(tRNA_map):
    """
    
    Description: return a list of tRNA locusNums
    """
    tRNALocusNums = []


    return tRNALocusNums



def get_categories_ptn(init_conc_path):
    """
    Based on the proteomics in Excel
    Return: the list of locusNums that with experimental protemomic counts larger than 10 (non-ribosomal ptns)
        ptn_list contains 338 ptns, ribosomal_ptn_list contains 52, ten_list 65
    """


    ptn_list = []
    
    ribosomal_ptn_list = []
    
    ten_ptn_list = []

    proteomics = pd.read_excel(init_conc_path, sheet_name='Comparative Proteomics')
    
    for row, protein in proteomics.iterrows():
        if row != 0:
            locusTag = protein['Locus Tag']
            locusNum = locusTag.split('_')[1]
            PtnName = protein['Gene Product']
    #         print(PtnName)
            if 'S ribosomal' not in PtnName: #exclude the ribosomal proteins
                initial_count = protein['Sim. Initial Ptn Cnt']
                if initial_count > 10:
                    ptn_list.append(locusNum)
                else:
                    ten_ptn_list.append(locusNum)
            else:
                ribosomal_ptn_list.append(locusNum)

    # print(f"ribosome Ptn: {ribosomal_ptn_list}")

    return ptn_list, ribosomal_ptn_list, ten_ptn_list

def getPtnLocations(init_conc_path):
    """
    
    Description: get the proteins at different locations
    """
    proteomics = pd.read_excel(init_conc_path, sheet_name='Comparative Proteomics', skiprows=[1])
    locations = ['trans-membrane', 'peripheral membrane', 'cytoplasm', 'lipoprotein', 'unidentified' ]

    locations_ptns = {}

    for location in locations:
        locations_ptns[location] = []
        for index, row in proteomics.iterrows():
            if row['Localization'] == location:
                locusTag = row['Locus Tag']
                locations_ptns[location].append(locusTag)

    for location in locations:
        print(f"{len(locations_ptns[location])} {location} proteins in Syn3A")

    return locations_ptns

def get_cplx_dict(cplx_path, locations_ptns):
    """
    Return a dictionary of complexes involved in the simulation

    Return: cplx_dict, two layer dictionary: {complex name, subdict}, {'init_count':0, 'Stoi': {'JCVISYN3A_0001':1, ...}, 'transMP_count': 2}
    """

    transmemPtnList = locations_ptns['trans-membrane']

    cplx_df = pd.read_excel(cplx_path, sheet_name='Complexes')
    abc_df = pd.read_excel(cplx_path, sheet_name='ABC Transporter')
    ATP_df = pd.read_excel(cplx_path, sheet_name='ATPSynthase', dtype=str)

    cplx_dict = {}

    def getABCTransporter(abc_complex, abc_df, subunit):

        locusNums = []
        
        for index, row in abc_df.iterrows():
            if row['Name'] == abc_complex:
                if row['Subunit'].startswith(subunit):
                    locusNums.append(row['Gene'])
        
        if locusNums != ['None']:
            locusTags = ['JCVISYN3A_' + locusNum for locusNum in locusNums]
        else:
            locusTags = []

        return locusTags



    for index, row in cplx_df.iterrows():
        name = row['Name']
        cplx_dict[name] = {}
        # init_count
        cplx_dict[name]['init_count'] = row['Count']
        # Stoi
        type = row['Type']; cplx_dict[name]['Stoi'] = {}
        locusNums = row['Genes'].split(';'); locusTags = ['JCVISYN3A_' +locusNum for locusNum in locusNums]

        if type == 'dimer':
            for locusTag in locusTags:
                cplx_dict[name]['Stoi'][locusTag] = locusTags.count(locusTag)

        elif type == 'dimer_activate':
            for locusTag in locusTags:
                cplx_dict[name]['Stoi'][locusTag] = 1

        elif type == 'ECFModule':
            for locusTag in locusTags:
                cplx_dict[name]['Stoi'][locusTag] = 1

        elif type == 'ECF':
            cplx_dict[name]['Stoi'] = dict(cplx_dict['ECF']['Stoi'])
            cplx_dict[name]['Stoi'][locusTags[0]] = 1

        elif type == 'ABCTransporter':
            subunits = ['TMD', 'NBD', 'SBP']
            for subunit in subunits:
                locusTags = getABCTransporter(name, abc_df, subunit)
                for locusTag in set(locusTags):
                    cplx_dict[name]['Stoi'][locusTag] = locusTags.count(locusTag)

        elif type == 'ATPSynthase':
            for locusTag in locusTags:
                stoi = ATP_df[ATP_df['Gene'] == locusTag.split('_')[1]]['Stoi']
                cplx_dict[name]['Stoi'][locusTag] = int(stoi)
        elif type == 'RNAP':
            for locusTag in locusTags:
                cplx_dict[name]['Stoi'][locusTag] = locusTags.count(locusTag) 
            cplx_dict[name]['Stoi']['JCVISYN3A_0645'] = 2
        elif type == 'SecYEG':
            for locusTag in locusTags:
                cplx_dict[name]['Stoi'][locusTag] = locusTags.count(locusTag)
        else:
            print(f'WARNING: Unknown Assembly Pathway {type}')

    for name, subdict in cplx_dict.items():
        Stoi = subdict['Stoi']
        transMP_count = 0
        for locusTag, stoi in Stoi.items():
            if locusTag in transmemPtnList:
                transMP_count += stoi
        subdict['transMP_count'] = transMP_count


    return cplx_dict



def getPtnInitCount(init_conc_path):
    """
    Get the initial count of ptn from Excel sheet
    """

    locusTagtoPtnInitCount = {}

    proteomics = pd.read_excel(init_conc_path, sheet_name='Comparative Proteomics', skiprows=[1])

    for index, row in proteomics.iterrows():
        locusTag = row['Locus Tag']
        init_count = row['Sim. Initial Ptn Cnt']
        locusTagtoPtnInitCount[locusTag] = init_count

    print(f"{len(locusTagtoPtnInitCount)} Proteins in proteomics Excel Sheet")

    return locusTagtoPtnInitCount



def correctInitPtnCount(PtnIniCount, locusTag, rPtn_locusTags, print_flag=False):
    """
    Input: 
        PtnIniCount: Count in proteomics Excel sheet
        locusTag
        cplx_dict: Python dictionary of complexes
        rPtn_locusTags

    Correct the initial protein counts for protein in complexes

    Return: 
        corrected_PtnIniCount: count of free protein after correction
        total_PtnIniCount: count of total ptn in the initial state after correction
        complex: complex that contains this protein 
    """


    complex = []

    total_PtnIniCount = 0
    
    old_PtnIniCount = PtnIniCount

    corrected_PtnIniCount = PtnIniCount

    # Assumption about the counts of free ribosomal protein pool (both SSU and LSU) is 5% of the complete assembled ribosome
    if locusTag in rPtn_locusTags:
        corrected_PtnIniCount = max(25, int(PtnIniCount - 500))
        complex.append('ribosomeP')
        total_PtnIniCount = corrected_PtnIniCount + 500
        if print_flag:
            print(f"P_{locusTag.split('_')[1]} in Ribosome Biogenesis Initial Count Corrected from {PtnIniCount} to {corrected_PtnIniCount}")
    else:
        total_PtnIniCount += corrected_PtnIniCount
        
    if print_flag:
        print(f"{locusTag.replace('JCVISYN3A', 'P')} in {', '.join(complex)} Initial Count Corrected from {old_PtnIniCount} to {corrected_PtnIniCount}")
        
        print(f"{locusTag.replace('JCVISYN3A', 'P')} Total Initial Count is {total_PtnIniCount}")

    return corrected_PtnIniCount, total_PtnIniCount, complex


def getScaledPtnRatio(w, init_conc_path, category):
    """
    Input:
        init_conc_path: path to the proteomics Excel sheet
        category: 
        ptn_list: list of locusNums of proteins that fall into category
        ptns_ratio: 3D numpy array and the values are the ratio of total protein counts over each initial total protein count     
        surface_doubling_times: list of SA doubling time

    Output:

    Description:
        Check the scaled protein abundance at the SA doubling times
    """

    normal_ptn_list, ribosomal_ptn_list, ten_ptn_list = get_categories_ptn(init_conc_path)

    rPtn_locusTags = ['JCVISYN3A_' + locusNum for locusNum in ribosomal_ptn_list]

    locusTagtoPtnInitCount = getPtnInitCount(init_conc_path)
    
    # ptns_lists = [normal_ptn_list, ribosomal_ptn_list, ten_ptn_list]

    if category == 'normal':
        ptns_list = normal_ptn_list
    elif category == 'rPtn':
        ptns_list = ribosomal_ptn_list
    elif category == 'ten':
        ptns_list = ten_ptn_list
    elif category == 'trans-membrane':
        TM_locusTags = getPtnLocations(init_conc_path)['trans-membrane']
        ptns_list = [locusTag.split('_')[1] for locusTag in TM_locusTags]
    elif category == 'entire':
        entire_locusTags = locusTagtoPtnInitCount.keys()
        ptns_list = [locusTag.split('_')[1] for locusTag in entire_locusTags]

    locusTags_ptns = ['JCVISYN3A_' + locusNum for locusNum in ptns_list]

    Produced_ptns = ['Produced_P_{0}'.format(locusNum) for locusNum in ptns_list ]

    Produced_ptns_counts = w.get_species_traces(Produced_ptns)

    ptns_ratio = np.zeros_like(Produced_ptns_counts)

    for i_ptn, locusTag in enumerate(locusTags_ptns):
        Produced_ptn_count = Produced_ptns_counts[i_ptn,:,:]

        PtnInitCount = locusTagtoPtnInitCount[locusTag]

        corrected_PtnIniCount, total_PtnIniCount, complex = correctInitPtnCount(PtnInitCount, locusTag, rPtn_locusTags)

        scaled_produced_ptn_count = Produced_ptn_count / total_PtnIniCount

        scaled_total_ptn_count = scaled_produced_ptn_count + 1

        ptns_ratio[i_ptn,:,:] = scaled_total_ptn_count

    return ptns_list, ptns_ratio


def getAbnormalPtn(init_conc_path, category, ptns_list, ptns_ratio, surface_doubling_times, genomeDict, threshold=[1.5,3]):
    abnormal_ptns = []; abnormal_ratios = []
    
    abnormal_ptns_dicts = []
    # get rPtn_locusTags
    normal_ptn_list, ribosomal_ptn_list, ten_ptn_list = get_categories_ptn(init_conc_path)
    rPtn_locusTags = ['JCVISYN3A_' + locusNum for locusNum in ribosomal_ptn_list]

    proteomics = pd.read_excel(init_conc_path, sheet_name='Comparative Proteomics', skiprows=[1])

    for i_ptn, locusNum in enumerate(ptns_list):
        ptn_ratio = ptns_ratio[i_ptn, :, :] # 2D array

        scaledPtn_SA_doubling = np.mean(math.get_doubling_moments_value(surface_doubling_times,ptn_ratio), axis=0) # Scalar value

        if scaledPtn_SA_doubling < threshold[0] or scaledPtn_SA_doubling > threshold[1]:
            abnormal_ptns.append(locusNum); abnormal_ratios.append(scaledPtn_SA_doubling)

            abnormal_ptn_dict = {}
            locusTag = 'JCVISYN3A_' + locusNum
            row = proteomics[proteomics['Locus Tag'] == locusTag]
            PtnIniCount = row['Sim. Initial Ptn Cnt'].values[0]

            corrected_PtnIniCount, total_PtnIniCount, complex = correctInitPtnCount(PtnIniCount, locusTag, rPtn_locusTags, print_flag=False)
            
            abnormal_ptn_dict['Locus Tag'] = locusTag
            abnormal_ptn_dict['Gene Name'] = row['Gene Name'].values[0]
            abnormal_ptn_dict['Gene Product'] = row['Gene Product'].values[0]
            abnormal_ptn_dict['Protein Length'] = len(genomeDict[locusTag]['AAsequence'])

            abnormal_ptn_dict['Exp. Ptn Cnt'] = row['Exp. Ptn Cnt'].values[0]
            # ptn['Sim. Initial Ptn Cnt'] = PtnIniCount
            abnormal_ptn_dict['Sim. Free Initial Ptn Cnt'] = corrected_PtnIniCount
            abnormal_ptn_dict['Sim. Total Initial Ptn Cnt'] = total_PtnIniCount
            
            abnormal_ptn_dict['Scaled Abundance'] = f"{scaledPtn_SA_doubling:.3f}"


            abnormal_ptns_dicts.append(abnormal_ptn_dict)
    
    abnornal_ptn_df = pd.DataFrame(abnormal_ptns_dicts)
    # print(abnornal_ptn_df)
    if abnornal_ptn_df.empty:
        sorted_abnornal_ptn_df = pd.DataFrame()
    else:
        sorted_abnornal_ptn_df = abnornal_ptn_df.sort_values(by='Scaled Abundance')

    return abnormal_ptns, abnormal_ratios, sorted_abnornal_ptn_df
    
def getRxnsDict():
    """
    Description: for the old kinetic_params Excel sheet
    {'BPNT': {'locusNums': ['0139']},
    'FAKr': {'locusNums': ['0420', '0616', '0617'], 'GPR': 'and'},...}

    """
    sheet_names = ['Central', 'Nucleotide', 'Lipid', 
                    'Cofactor', 'Transport']
    kinetic_params = '/data/enguang/CMEODE/input_data/kinetic_params.xlsx'
    rxns_dict = {}
    for sheet_name in sheet_names:
        params = pd.read_excel(kinetic_params, sheet_name=sheet_name)
        for index, row in params.iterrows():
            if row['Parameter Type'] == 'Eff Enzyme Count':
                rxn_id = row['Reaction Name']
                ptns = row['Value']
                if ptns != 'default':
                    ptns = ptns.split('-')
                    locusNums = [ptn.split('_')[1] for ptn in ptns]
                    rxns_dict[rxn_id] = {}
                    rxns_dict[rxn_id]['GPR_locusNums'] = locusNums
                    
                    if len(locusNums) > 1:
                        rxn_params = params[params['Reaction Name'] == rxn_id]
                        GPRrule = rxn_params[rxn_params['Parameter Type'] == 'GPR rule']['Value'].values[0]
                        rxns_dict[rxn_id]['GPR'] = GPRrule

                        
                else:
                    print(f"Reaction {rxn_id} has no assigned proteins")

    return rxns_dict

def getChangedRxnsDict(rxns_dict, cplx_dict):
    """
    
    Get the information of metabolic reactions that changed from GPR AND/OR rule to complex
    """
    changed_rxns_dict = {}

    print("********* Checking Changed Reactions ***********")

    for rxn_id, rxn_subdict in rxns_dict.items():
        GPR_locusNums = rxn_subdict['GPR_locusNums']

        for cplx_id, cplx_subdict in cplx_dict.items():
            cplx_stoi = cplx_subdict['Stoi']
            cplx_locusTags = cplx_stoi.keys()
            cplx_locusNums = [locusTag.split('_')[1] for locusTag in cplx_locusTags]
            if cplx_id != 'ECF':
                if set(GPR_locusNums).issubset(set(cplx_locusNums)):
                    print(f"{rxn_id} Old reaction locusNums {GPR_locusNums} cplx {cplx_id} {cplx_locusNums}")
                    changed_rxns_dict[rxn_id] = rxn_subdict
                    changed_rxns_dict[rxn_id]['complex'] = cplx_id
                    
                elif set(cplx_locusNums).issubset(set(GPR_locusNums)):
                    print(f"{rxn_id} Old reaction locusNums {GPR_locusNums} cplx {cplx_id} {cplx_locusNums}")
                    changed_rxns_dict[rxn_id] = rxn_subdict
                    changed_rxns_dict[rxn_id]['complex'] = cplx_id

    print(changed_rxns_dict)

    return changed_rxns_dict


def checkComplexStoi(w, surface_doubling_times, cplx_dict, cplx, produced_prefix = 'Produced_P_'):
    """
    
    Return:
        ratio_dict: {'0009': 1, '0010': 2}

    Description: 
        Check the Stoichiometric balance among generated subunits up to SA doubling time; 
        Divide the generated subunits counts over the initial complex count
    """

    ratio_dict = {}

    Stoi = cplx_dict[cplx]['Stoi']

    init_cplx_count = cplx_dict[cplx]['init_count'] # initial count of assembled complex

    print(f"Stoi of complex {cplx}: {Stoi}")

    locusTags = list(Stoi.keys())

    produced_ptns_ids = [produced_prefix + locusTag.split('_')[1] for locusTag in locusTags ]

    produced_subunits_counts = w.get_species_traces(produced_ptns_ids)

    counts_SA = math.get_doubling_moments_value(surface_doubling_times, produced_subunits_counts)

    ratio_SA = counts_SA/init_cplx_count

    avg_ratio_SA = np.mean(ratio_SA, axis=1)

    for locusTag, ratio in zip(locusTags, avg_ratio_SA):
        ratio_dict[locusTag] = ratio

    ratio_dict['init_count'] = init_cplx_count

    return ratio_dict
 

def plotProducedPtns(w, surface_doubling_times, cplx_dict, cplx, fig_dir, fig_label, reps, produced_prefix = 'Produced_P_'):
    """

    Description: plot the traces and histogram of produced protein subunits of a complex 
    """

    Stoi = cplx_dict[cplx]['Stoi']

    init_cplx_count = cplx_dict[cplx]['init_count'] # initial count of assembled complex

    print(f"Stoi of complex {cplx}: {Stoi}")

    locusTags = list(Stoi.keys())

    produced_ptns_ids = [produced_prefix + locusTag.split('_')[1] for locusTag in locusTags ]

    produced_subunits_counts = w.get_species_traces(produced_ptns_ids)

    counts_SA = math.get_doubling_moments_value(surface_doubling_times, produced_subunits_counts)

    for i_subunit, subunit in enumerate(locusTags):

        ylabel = 'Count [\#]'
        title = produced_ptns_ids[i_subunit] + f' in {cplx} with init count {init_cplx_count} and Stoi {Stoi[subunit]}'
        w.plot_in_replicates_single(fig_dir,fig_label,
                               '.png',
                              produced_subunits_counts[i_subunit,:,:], reps, ylabel, title,
                               True, True)
        
        xlabel = produced_ptns_ids[i_subunit] + ' at SA doubling time'
        produced_subunits_SA = counts_SA[i_subunit,:]
        w.plot_hist(fig_dir, fig_label, '.png', produced_subunits_SA, xlabel, ylabel, title, bins=50)

    return None

def plotIntermediatesCplx(w, ):


    return None
# def compareCplx(cplx_dict, rxns_dict):

#     compared_cplx_dict = {}

    
    
#     return compared_cplx_dict