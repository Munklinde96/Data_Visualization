import pandas as pd
import requests
import re
import numpy as np
import matplotlib.patches as patches


######################################
############ PREPARE DATA ############
######################################

uniprot = 'https://www.uniprot.org/uniprot/'

def count_no_of_modifications(ptm_str):
    #check if NaN value
    if pd.isnull(ptm_str):
        return None
    return 1 + ptm_str.count(';')

def get_data_and_remove_unwanted_columns():
    # df = pd.read_csv('UHT milk P036.csv')
    df = pd.read_csv('UHT milk P036.csv')
    df.drop('Protein ID', inplace=True, axis=1)
    df.drop('Unique', inplace=True, axis=1)
    df.drop('m/z', inplace=True, axis=1)
    df.drop('Scan', inplace=True, axis=1)
    df.drop('Source File', inplace=True, axis=1)
    df.drop('Found By', inplace=True, axis=1)
    #apply count_no_of_modifications to each PTM column
    df['#modifications'] = df['PTM'].apply(count_no_of_modifications)

    return df

def clean_peptide(peptide):
    mod_start_list = [m.start() for m in re.finditer('\(', peptide)]
    mod_end_list = [m.start() for m in re.finditer('\)', peptide)]
    clean_peptide = peptide
    if len(mod_start_list) is not 0:
        clean_peptide = peptide[0 : mod_start_list[0]]
        for i in range(len(mod_start_list) - 1):
            clean_peptide = clean_peptide + peptide[mod_end_list[i] + 1 : mod_start_list[i+1]]
        clean_peptide = clean_peptide + peptide[mod_end_list[-1] +1 :]
    if peptide[1] is '.':
        clean_peptide = clean_peptide[2:]
    if peptide[len(peptide) - 2] is '.':
        clean_peptide = clean_peptide[0:len(clean_peptide)-2]
    return clean_peptide

def get_protein_sequence(protein):
    url = uniprot + protein + '.fasta'
    response = requests.get(url).text
    if(response):
        first_new_line = response.find('\n')
        sequence = response[first_new_line:].replace('\n', '')
        return sequence
    else:
        return ''

def sanitize_data(df):
    protein_sequence = ''
    current_protein = ''
    for index, row in df.iterrows():
        peptide = row['Peptide']
        peptide = clean_peptide(peptide)
        if row['Length'] is not len(peptide) and row['End'] - row['Start'] is not len(peptide):
            print('length mismatch')
            df.drop(index)
        if(current_protein != row['Protein Accession']):
            current_protein = row['Protein Accession']
            protein_sequence = get_protein_sequence(current_protein)
        #The start and end are 1-indexed
        protein_sequence_substring = protein_sequence[row['Start'] - 1 : int(row['End'])]
        if protein_sequence_substring != peptide:
            print('peptide mismatch')
            df.drop(index)
        return df

def get_and_prepare_data():
    df = get_data_and_remove_unwanted_columns()
    df = sanitize_data(df)
    return df


##############################################
############ USEFUL UTILS METHODS ############
##############################################

# nomalize peptide intensity(per sample) over protetin intensity
def normalize_intensities_by_protein_intensity(df, sample_column_id='Area'):
    protein_start = [0]
    protein_end = []
    protein_id = ""

    for i, row in df.iterrows():
        if i == 0:
            protein_id = row['Protein Accession']

        if(row['Protein Accession'] != protein_id):
            protein_start.append(i)
            protein_end.append(i)
            protein_id = row['Protein Accession']

    protein_start.pop()

    dataframes = []
    for i in range(len(protein_start)):
        protein_df = df.iloc[protein_start[i] : protein_end[i]]
        protein_df = protein_df.copy()
        area_cols = [col for col in protein_df.columns if sample_column_id in col]
        for col_name in area_cols:
            intensity_sum = protein_df[col_name].sum()
            protein_df[col_name] = protein_df[col_name].divide(intensity_sum)

        dataframes.append(protein_df)
    
    return pd.concat(dataframes, axis=0, ignore_index=True)

# indicies have sign in front
def get_position_of_mass_shift(input_string):
    # get charachters after and before mass shift
    ls = input_string.split('(')
    if len(ls) == 0:
        return np.nan
    modified_string = re.sub(r"\([^()]*\)", "#", input_string) #replace everythting witthin parenthesis with "#"
    modified_string = modified_string[2:-2] #remove first and last splice sites
    indices = [m.start()-1 for m in re.finditer("#", modified_string)] #get indices of "#"
    for i in range(len(indices)):
        indices[i] -= i
    return indices

def get_color_palette_for_modifications ():
    mod_type_map = {'Oxidation (M)': '#E6194B' , 'Phosphorylation (STY)': '#F58231', 'Deamidation (NQ)': '#FFE119', 'lal': '#BFEF45',
                      'Lactosylation': '#3CB44B', 'Pyro-glu from Q': '#42D4F4', 'Glycosylation type b': '#4363D8', 'Dioxidation (M)': '#911EB4',
                    'Glycosylation type e': '#F032E6', 'Glycosylation type a': '#000000','Glycosylation type c/d': '#800000', 'Carbamidomethylation': '#FABED4', 'lan': '#808000'}
    return mod_type_map

def get_color_legend_mods(modification_types_to_color_map):
    handles = []
    for mod_type in modification_types_to_color_map:
        handles.append(patches.Patch(color=modification_types_to_color_map[mod_type], label=mod_type))
    return handles
