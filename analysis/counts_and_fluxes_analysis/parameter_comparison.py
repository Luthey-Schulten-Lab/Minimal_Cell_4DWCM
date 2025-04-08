import pandas as pd
import numpy as np

new_params_file = './parameter_comparison/kinetic_params.xlsx'
old_params_file = './parameter_comparison/kinetic_params_old.xlsx'

parameter_report = './parameter_report.txt'

def read_params_file(params_file):


    central_params = pd.read_excel(params_file,sheet_name='Central')
    nucleotide_params = pd.read_excel(params_file,sheet_name='Nucleotide')
    lipid_params = pd.read_excel(params_file,sheet_name='Lipid')
    AA_params = pd.read_excel(params_file,sheet_name='Amino Acid')
    cofactor_params = pd.read_excel(params_file,sheet_name='Cofactor')
    transport_params = pd.read_excel(params_file,sheet_name='Transport')

    params = pd.concat([central_params,nucleotide_params,\
                        lipid_params,AA_params,\
                        cofactor_params,transport_params],
                       ignore_index=True)

    # params = pd.concat([central_params,nucleotide_params,\
    #                     lipid_params,AA_params,\
    #                     cofactor_params],
    #                    ignore_index=True)

    params = params[params['Parameter Type'] != 'Eff Enzyme Count']
    params = params[params['Parameter Type'] != 'GPR rule']

    return params

new_params = read_params_file(new_params_file)
old_params = read_params_file(old_params_file)

new_rxns = new_params['Reaction Name'].unique()
old_rxns = new_params['Reaction Name'].unique()

intersect_rxns = [x for x in new_rxns if x in set(old_rxns)]
diff_rxns = [x for x in new_rxns if x not in set(old_rxns)]

with open(parameter_report, 'w') as f:

    file_header = '=====[  Parameter Files  ]====='
    f.write(len(file_header)*'='+'\n')
    f.write(file_header+'\n')
    f.write(len(file_header)*'='+'\n\n')
    
    f.write('  new parameter file: {}\n'.format(new_params_file))
    f.write('  old parameter file: {}\n'.format(old_params_file))

    f.write('\n\n')

    f.write('Missing Reactions in Old:')

    for rxn in diff_rxns:
        f.write('{rxn}\n')

    f.write('\n\n')

    f.write('Overlapping Reactions in Old:\n\n')
    
    for rxn in intersect_rxns:

        f.write('=====[  {}  ]====='.format(rxn))

        new_temp = new_params[new_params['Reaction Name'] == rxn]
        old_temp = old_params[old_params['Reaction Name'] == rxn]

        for index, row in new_temp.iterrows():

            new_param_type = row['Parameter Type']

            old_entry = old_temp[old_temp['Parameter Type'] == new_param_type]

            need_species = False

            if (new_param_type not in ['Substrate Catalytic Rate Constant',
                                       'Product Catalytic Rate Constant',
                                       'Catalytic Rate Constant']):

                need_species = True
                new_species = row['Related Species']
                old_entry = old_entry[old_entry['Related Species'] == new_species]

            if old_entry.shape[0] == 0:

                missing_message = '\n  Missing Entry in Old:\t{} - ({})'.format(rxn,new_param_type)

                if need_species == True:
                    missing_message += ' - ({})'.format(new_species)

                f.write(missing_message)
                continue

            old_val = float(old_entry['Value'].iloc[0])
            new_val = float(row['Value'])

            if old_val != new_val:

                mismatch_message = '\n  Value Mismatch:\t{} - ({})'.format(rxn,new_param_type)

                if need_species == True:
                    mismatch_message += ' - ({})'.format(new_species)

                mismatch_message += '\n\t{:f} (new) vs. {:f} (old)'.format(new_val,old_val)

                f.write(mismatch_message)

        f.write('\n\n')
