"""
Authors: Zane Thornburg

Define all reactions that are integrated in a system of ODEs
"""

from collections import defaultdict, OrderedDict

import odecell

# import lipid_patch_Zane_polysomes as lipid_patch #_polysomes
# import PPP_patch_Zane as ppp_patch

import pandas as pd
import numpy as np

import libsbml

### Constants
NA = 6.022e23 # Avogadro's
r_cell = 200e-9 # 200 nm radius, 400 nm diameter
V = ((4/3)*np.pi*(r_cell)**3)*(1000) # for a spherical cell

countToMiliMol = 1000/(NA*V)

#########################################################################################
def initModel(sim_properties):
    """
    Input:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """

    # Initialize the ODECell model
    model = odecell.modelbuilder.MetabolicModel()
    
    zeroOrderOnOff = '$onoff * $K'
    model.zeroOrderOnOff = odecell.modelbuilder.RateForm(zeroOrderOnOff)
    model.updateAvailableForms()

    # Set verbosity outputs to zero for now to improve performance
    model.setVerbosity(0)

    # Define Rxns and pass in the Particle Map containing enzyme concentrations
    model = defineRxns(model, sim_properties)

    return model
#########################################################################################


#########################################################################################
def defineRxns(model, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    model = addProteinMetabolites(model, sim_properties)
    
    model = defineRandomBindingRxns(model, sim_properties)
    
    model = defineNonRandomBindingRxns(model, sim_properties)
    
#     model = defineOtherRandomBindingReactions(model, sim_properties)
    
    return model
#########################################################################################


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


#########################################################################################
def partToSaC(particles, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    conc = particles/sim_properties['SA_total']

    return conc
#########################################################################################


#########################################################################################
def reptModel(model):
    """
    Report on the constructed hybrid model - but probably would only want to do after the first time step

    Arguments: 
    model (model obj.): The ODE Kinetic Model

    Returns:

    None
    """

    dictTypes = defaultdict(int)
    typeList = ["Transcription","Translation","Degradation"]

    for rxn in model.getRxnList():
        
        if rxn.getResult():
            # If an explicit result has been set, this is a dependent reaction.
            dictTypes["Dependent reactions"] += 1
            continue
        
        for rxntype in typeList:
            if rxntype in rxn.getID():
                dictTypes[rxntype] += 1

                
    #print( "There are {} ODEs in the model:".format(len(model.getRxnList())) )

    outList = list(dictTypes.items())
    outList.sort(key=lambda x: x[1], reverse=True)
    for key,val in outList:
        print("{:>20} :   {}".format(key,val) )
        return 0
    
    return None
#########################################################################################


#########################################################################################
def defineRandomBindingRxns(model, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    params_file = sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx'
    
    central_params = pd.read_excel(params_file, sheet_name='Central')
    nucleotide_params = pd.read_excel(params_file, sheet_name='Nucleotide')
    lipid_params = pd.read_excel(params_file, sheet_name='Lipid')
    cofactor_params = pd.read_excel(params_file, sheet_name='Cofactor')
    transport_params = pd.read_excel(params_file, sheet_name='Transport')
    
    metabolism_params = pd.concat([central_params, nucleotide_params, lipid_params, cofactor_params, transport_params], ignore_index=True) 
    
    reaction_list = []

    for row, item in metabolism_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])
            
    sbmlFile = sim_properties['head_directory'] + "input_data/Syn3A_updated.xml"

    docSBML = libsbml.readSBMLFromFile(sbmlFile)
    modelSBML = docSBML.getModel()

    speciesNames = [spc.name for spc in modelSBML.getListOfSpecies()]
    speciesNamesLower = [x.lower() for x in speciesNames]
    speciesIDs = [spc.id for spc in modelSBML.getListOfSpecies()]

    rxnNamesSBML = [ x.name for x in modelSBML.getListOfReactions()]
    
    for rxnID in reaction_list:

        rxn_info = getSpecIDs(rxnID, modelSBML, rxnNamesSBML)

        rxn_params = metabolism_params.loc[ metabolism_params["Reaction Name"] == rxnID ]

        substrates_list = rxn_info[0][0]
        substrates_stoich = rxn_info[0][1]
        products_list = rxn_info[1][0]
        products_stoich = rxn_info[1][1]
        
        substrate_count = int(-np.sum(substrates_stoich))
        product_count = int(np.sum(products_stoich))
        
        
        rateLaw = Enzymatic(substrate_count, product_count)
        
        
        rateName = rxnID+'_rate'
        
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)

        kcatF = rxn_params.loc[ rxn_params["Parameter Type"] == "Substrate Catalytic Rate Constant" ]["Value"].values[0]
        kcatR = rxn_params.loc[ rxn_params["Parameter Type"] == "Product Catalytic Rate Constant" ]["Value"].values[0]
        
        model.addParameter(rxnIndx, 'kcatF', kcatF)
        model.addParameter(rxnIndx, 'kcatR', kcatR)


        rxn_KMs = rxn_params.loc[ rxn_params["Parameter Type"] == "Michaelis Menten Constant" ]
    
        # Add substrates to the reaction
        sub_rxn_indx_counter = 0
        for i in range(len(substrates_list)):

            metID = substrates_list[i]
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID], sim_properties)
                
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(-substrates_stoich[i])
            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
            for j in range(stoichiometry):
                
                sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
                rateFormID = 'Sub' + str(sub_rxn_indx_counter)
            
                if metID.endswith('_e'):
                    
                    model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])

                else:

                    if metID == 'M_o2_c':
                        model.addParameter(rxnIndx, rateFormID, metID)
                    else:
                        model.addSubstrate(rxnIndx, rateFormID, metID)
                
                KM_ID = 'KmSub' + str(sub_rxn_indx_counter)
                
                model.addParameter(rxnIndx, KM_ID, met_KM)
                

        # Add products to the reaction
        prod_rxn_indx_counter = 0
        for i in range(len(products_list)):

            metID = products_list[i]
            
            if metID not in list(model.getMetDict().keys()):
                
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID], sim_properties)
                
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(products_stoich[i])

            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
            for j in range(stoichiometry):
                
                prod_rxn_indx_counter = prod_rxn_indx_counter + 1
                
                rateFormID = 'Prod' + str(prod_rxn_indx_counter)
            
                if metID.endswith('_e'):

                    model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])

                else:

                    model.addProduct(rxnIndx, rateFormID, metID)
                
                KM_ID = 'KmProd' + str(prod_rxn_indx_counter)
                
                model.addParameter(rxnIndx, KM_ID, met_KM)

            
        EnzymeConc = getEnzymeConc(rxn_params, sim_properties)
        
        model.addParameter(rxnIndx, "Enzyme", EnzymeConc)
        
        model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1)
        

    reptModel(model)

    return model
#########################################################################################


#########################################################################################
def defineOtherRandomBindingReactions(model, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    params_file = sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx'
    
    RXNS_params = pd.read_excel(params_file, sheet_name='Other-Random-Binding')
    
    reaction_list = []

    for row, item in RXNS_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])
#             print(item['Reaction Name'])
            
    for rxnID in reaction_list:
    
#         print(rxnID)

        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]
        
        substrate_count = int(rxn_params.loc[ rxn_params["Parameter Type"] == "Substrates" ]["Value"].values[0])
        product_count = int(rxn_params.loc[ rxn_params["Parameter Type"] == "Products" ]["Value"].values[0])
        
        rateLaw = Enzymatic(substrate_count, product_count)
        
        rateName = rxnID+'_rate'
        
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)

        kcatF = rxn_params.loc[ rxn_params["Parameter Type"] == "Substrate Catalytic Rate Constant" ]["Value"].values[0]
        kcatR = rxn_params.loc[ rxn_params["Parameter Type"] == "Product Catalytic Rate Constant" ]["Value"].values[0]
        
        model.addParameter(rxnIndx, 'kcatF', kcatF)
        model.addParameter(rxnIndx, 'kcatR', kcatR)
        
#         print(kcatF)
#         print(kcatR)

        rxn_KMs = rxn_params.loc[ rxn_params["Parameter Type"] == "Michaelis Menten Constant" ]
        
        
        for i in range(1, substrate_count+1):

            metID = rxn_params.loc[ rxn_params["Parameter Type"] == "Sub" + str(i) ]["Value"].values[0]
            
#             if spcID.endswith("_e"):
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID], sim_properties)
                
#                 print('Added metabolite: ', metID, spcConc)
                
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(1)
            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
#             for j in range(stoichiometry):
                
#                 sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
            rateFormID = 'Sub' + str(i)
            
            if metID.endswith('_e'):

#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

            else:
#                     print(rxnIndx, rateFormID, metID)
                model.addSubstrate(rxnIndx, rateFormID, metID)

            KM_ID = 'KmSub' + str(i)

            model.addParameter(rxnIndx, KM_ID, met_KM)
                

        for i in range(1, product_count+1):

            metID = rxn_params.loc[ rxn_params["Parameter Type"] == "Prod" + str(i) ]["Value"].values[0]
            
#             if spcID.endswith("_e"):
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID], sim_properties)
                
#                 print('Added metabolite: ', metID, spcConc)
                
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(1)
            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
#             for j in range(stoichiometry):
                
#                 sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
            rateFormID = 'Prod' + str(i)
            
            if metID.endswith('_e'):

#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

            else:
#                     print(rxnIndx, rateFormID, metID)
                model.addProduct(rxnIndx, rateFormID, metID)

            KM_ID = 'KmProd' + str(i)

            model.addParameter(rxnIndx, KM_ID, met_KM)

            
        EnzymeConc = getEnzymeConc(rxn_params, sim_properties)
            
#         EnzymeConc = partTomM(rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0], sim_properties)
        
        model.addParameter(rxnIndx, "Enzyme", EnzymeConc)
        
        model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1)
        
    return model
#########################################################################################


#########################################################################################
def defineNonRandomBindingRxns(model, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    params_file = sim_properties['head_directory'] + 'input_data/kinetic_params.xlsx'
    
    RXNS_params = pd.read_excel(params_file, sheet_name='Non-Random-Binding Reactions')
    
    reaction_list = []

    for row, item in RXNS_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])
#             print(item['Reaction Name'])
            
    for rxnID in reaction_list:
    
#         print(rxnID)

        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]
        
        rateLaw = str(rxn_params.loc[ rxn_params["Parameter Type"] == "Kinetic Law" ]["Value"].values[0])
        
        rateName = rxnID+'_rate'
        
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)
        
        for index, row in rxn_params.iterrows():
            
            param = row['Parameter Type']
            
            if (param != "Reaction Formula") and (param != "Kinetic Law"):

                if param.startswith('Sub'):
                    
                    metID = row['Value']
                    
                    if metID not in list(model.getMetDict().keys()):
            
                        if metID.endswith('_e'):

                            spcConc = sim_properties['medium'][metID]

                        else:

                            spcConc = partTomM(sim_properties['counts'][metID], sim_properties)

#                         print('Added metabolite: ', metID, spcConc)

                            model.addMetabolite(metID, metID, spcConc)

                    if metID.endswith('_e'):

                        model.addParameter(rxnIndx, param, sim_properties['medium'][metID])

                    else:

                        model.addSubstrate(rxnIndx, param, metID)

                elif param.startswith('Prod'):
                    
                    metID = row['Value']
                    
                    if metID not in list(model.getMetDict().keys()):
            
                        if metID.endswith('_e'):

                            spcConc = sim_properties['medium'][metID]

                        else:

                            spcConc = partTomM(sim_properties['counts'][metID], sim_properties)

#                         print('Added metabolite: ', metID, spcConc)

                            model.addMetabolite(metID, metID, spcConc)

                    if metID.endswith('_e'):

                        model.addParameter(rxnIndx, param, sim_properties['medium'][metID])

                    else:

                        model.addProduct(rxnIndx, param, metID)

                else:

                    model.addParameter(rxnIndx, param, row['Value'])

    return model
#########################################################################################


#########################################################################################
def addProteinMetabolites(model, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    data_file = sim_properties['head_directory'] + 'input_data/protein_metabolites.xlsx'
    
    ptnMets = pd.read_excel(data_file, sheet_name='protein metabolites')
    
    for index, row in ptnMets.iterrows():
        
        ptnID = row['Protein']
        
        metabolites = row['Metabolite IDs'].split(',')
        
#         print(ptnID)
#         print(metabolites)
        
        ptnCount = sim_properties['counts'][ptnID]
        
        formsCount = 0
        
        for metID in metabolites:
            
            if metID != metabolites[0]:
                
                formsCount = formsCount + sim_properties['counts'][metID]
                
#                 print(metID, formsCount)
                
                model.addMetabolite(metID, metID, partTomM(sim_properties['counts'][metID], sim_properties))
                
        baseFormCount = int(ptnCount - formsCount)
        
#         print(metabolites[0], baseFormCount)
        
        baseFormConc = partTomM(baseFormCount, sim_properties)

        model.addMetabolite(metabolites[0], metabolites[0], baseFormConc)
        
    return model
#########################################################################################


#########################################################################################
def Enzymatic(subs, prods):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
        
    def numerator(subs,prods):
        
        subterm = [ '( $Sub' + str(i) + ' / $KmSub' + str(i) + ' )' for i in range(1,subs+1)]
        subNumer = ' * '.join(subterm)
        
        prodterm = [ '( $Prod' + str(i) + ' / $KmProd' + str(i) + ' )' for i in range(1,prods+1)]
        prodNumer = ' * '.join(prodterm)
        
        numerator = '( ' + '$kcatF * ' + subNumer + ' - ' + '$kcatR * ' + prodNumer + ' )'
        return numerator
    
    def denominator(subs,prods):
        
        subterm = [ '( 1 + $Sub' + str(i) + ' / $KmSub' + str(i) + ' )' for i in range(1,subs+1)]
        subDenom = ' * '.join(subterm)
        
        prodterm = [ '( 1 + $Prod' + str(i) + ' / $KmProd' + str(i) + ' )' for i in range(1,prods+1)]
        prodDenom = ' * '.join(prodterm)
        
        denominator = '( ' + subDenom + ' + ' + prodDenom + ' - 1 )'
        return denominator
        
    rate = '$onoff * $Enzyme * ( ' + numerator(subs,prods) + ' / ' + denominator(subs,prods) + ' )'
    
    return rate
#########################################################################################


#########################################################################################
def getEnzymeConc(rxn_params, sim_properties):
    """
    Inputs:
    sim_properties - Dictionary of simulation variables and state trackers
    
    Returns:
    Called by:
    Description:
    """
    
    EnzymeStr = rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0]
    
    Enzymes = EnzymeStr.split('-')
#     print(Enzymes)
    
    if len(Enzymes) == 1:
        
        if Enzymes[0] == 'default':
            
#             EnzymeConc = partTomM(10, sim_properties)
            EnzymeConc = 0.001
            
            return EnzymeConc
        
        else:
            
            ptnID = Enzymes[0]
            
#             print(sim_properties['counts'][ptnID])
            
            EnzymeConc = partTomM(sim_properties['counts'][ptnID], sim_properties)
            
            return EnzymeConc
    
    else:
        
        GPRrule = rxn_params.loc[ rxn_params["Parameter Type"] == "GPR rule" ]["Value"].values[0]

        if GPRrule == 'or':
            
            ptnCount = 0
            
            for ptnID in Enzymes:
                
                ptnCount = ptnCount + sim_properties['counts'][ptnID]
                
            EnzymeConc = partTomM(ptnCount, sim_properties)
                
            return EnzymeConc
        
        elif GPRrule == 'and':
            
            ptnCounts = []
            
            for ptnID in Enzymes:
                
                ptnCounts.append(sim_properties['counts'][ptnID])
                
            ptnCount = min(ptnCounts)
            
            EnzymeConc = partTomM(ptnCount, sim_properties)
                
            return EnzymeConc
        
    print('Something went wrong getting enzyme count')
        
    return None
#########################################################################################

    
#########################################################################################
def getSpecIDs(rxnName, modelSBML, rxnNamesSBML):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    returnList = []
    
    rxnObj = modelSBML.getReaction( rxnNamesSBML.index(rxnName) )
    
    # Use model SBML to get IDs, names, and stoichiometries for reactants
    specIDs = [ x.getSpecies() for x in rxnObj.getListOfReactants() ]
    spcStoich = [ -1*float(x.getStoichiometry()) for x in rxnObj.getListOfReactants() ]
    spcNames = [ modelSBML.getSpecies( spcID ).name for spcID in specIDs]
    
    specIDs_noH = []
    spcStoich_noH = []
    
    for i in range(len(specIDs)):
        
        metID = specIDs[i]
        
        if (metID != 'M_h_c') and (metID != 'M_h_e') and (metID != 'M_h2o_c')  and (metID != 'M_h2o_e'):
            
            specIDs_noH.append(metID)
            stoich = spcStoich[i]
            spcStoich_noH.append(stoich)
    
    if np.any( np.isnan( spcStoich ) ):
        raise Exception('Invalid stoichiometry for reaction: {}'.format(rxnName)) 
    
#     returnList.append( [specIDs, spcStoich] )
#     returnList.append( [spcNames, specIDs, spcStoich] )
    returnList.append( [specIDs_noH, spcStoich_noH] )
    
    # Now do the same for products
    specIDs = [ x.getSpecies() for x in rxnObj.getListOfProducts() ]
    spcStoich = [ float(x.getStoichiometry()) for x in rxnObj.getListOfProducts() ]
    spcNames = [ modelSBML.getSpecies( spcID ).name for spcID in specIDs]
    
    specIDs_noH = []
    spcStoich_noH = []
    
    for i in range(len(specIDs)):
        
        metID = specIDs[i]
        
        if (metID != 'M_h_c') and (metID != 'M_h_e') and (metID != 'M_h2o_c')  and (metID != 'M_h2o_e'):
            
            specIDs_noH.append(metID)
            stoich = spcStoich[i]
            spcStoich_noH.append(stoich)
    
    if np.any( np.isnan( spcStoich ) ):
        raise Exception('Invalid stoichiometry for reaction: {}'.format(rxnName)) 
    
#     returnList.append( [specIDs, spcStoich] )
#     returnList.append( [spcNames, specIDs, spcStoich] )
    returnList.append( [specIDs_noH, spcStoich_noH] )
    
    return returnList
#########################################################################################


