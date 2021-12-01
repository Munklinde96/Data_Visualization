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
    

def get_data_and_remove_unwanted_columns(wanted_columns, data_path = 'UHT milk P036.csv', sample_column_id = "Area"):
    # df = pd.read_csv('UHT milk P036.csv')
    df = pd.read_csv(data_path)
    sample_columns = [col for col in df.columns if sample_column_id in col]
    # only use wanted columns and sample colulmns 
    df = df[wanted_columns + sample_columns]
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
    if '|' in protein:
        protein = protein.split('|')[0]
    url = uniprot + protein + '.fasta'
    response = requests.get(url).text
    if(response):
        first_new_line = response.find('\n')
        sequence = response[first_new_line:].replace('\n', '')
        return sequence
    else:
        return ''

def sanitize_data(df, sample_column_id='Area'):
    selected_sample_columns = [col for col in df.columns if sample_column_id in col]
    for col in selected_sample_columns:  #make Area samples Nan values 0
        df[col] = df[col].fillna(0)
    df = df[df[selected_sample_columns].sum(axis=1) > 0]
    df["PTM"] = df["PTM"].fillna("Unmodified")

    protein_sequence = ''
    current_protein = ''
    dropped_indices = []
    for index, row in df.iterrows():
        peptide = row['Peptide']
        peptide = clean_peptide(peptide)
        if row['Length'] is not len(peptide) and row['End'] - row['Start'] is not len(peptide):
            print('length mismatch')
            dropped_indices.append(index)
            continue
        if(current_protein != row['Protein Accession']):
            current_protein = row['Protein Accession']
            protein_sequence = get_protein_sequence(current_protein)
        #The start and end are 1-indexed
        protein_sequence_substring = protein_sequence[row['Start'] - 1 : int(row['End'])]
        if protein_sequence_substring != peptide:
            print('peptide mismatch')
            dropped_indices.append(index)
            continue
        # if sum of selected_sample_columns in row is zero, drop the row
        if sum(row[selected_sample_columns]) == 0:
            dropped_indices.append(index)
            continue
    #drop rows in dropped_indices
    df = df.drop(dropped_indices)
    return df

def get_and_prepare_data(data_path = 'UHT milk P036.csv', wanted_columns = ['Protein Accession', 'Peptide', 'PTM', 'Start', 'End', 'Length'], sample_column_id = 'Area'):
    df = get_data_and_remove_unwanted_columns(wanted_columns, data_path, sample_column_id = sample_column_id)
    df = sanitize_data(df, sample_column_id = sample_column_id)
    return df

##############################################
############ USEFUL UTILS METHODS ############
##############################################

def normalize_intensities(df, sample_column_id = 'Area'):
    selected_sample_columns = [col for col in df.columns if sample_column_id in col]
    for col in selected_sample_columns:
        df[col] = df[col] / df[col].sum()
    return df

# nomalize peptide intensity(per sample) over protetin intensity
def normalize_intensities_by_protein_intensity(df, sample_column_id='Area'):
    print(df.head())
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
    print(len(dataframes))
    
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

def get_color_palette_for_modifications(modification_types = []):
    if modification_types == []:
        # use keys from mod_type map harcoded 
        modification_types = ['lal', 'Glycosylation type e', 'lan', 'Glycosylation type a', 'Glycosylation type b', 'Glycosylation type c/d', 'Phosphorylation (STY)', 'Dioxidation (M)', 'Oxidation (M)', 'Lactosylation', 'Carbamidomethylation', 'Deamidation (NQ)', 'Pyro-glu from Q']
    # COLORS = ['#3CB44B',  '#FFE119', '#F032E6', '#808000', '#000000','#BFEF45','#42D4F4', '#4363D8','#FABED4', '#800000', '#F58231', '#911EB4', '#E6194B']    
    # remove navy, blue, cyan and grey from colors from: https://sashamaps.net/docs/resources/20-colors/
    COLORS = '#e6194b', '#3cb44b', '#ffe119', '#f58231', '#911eb4',  '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000000'
    # show color palette
    
    mod_type_map = {}
    for i in range(len(modification_types)):
        mod_type_map[modification_types[i]] = COLORS[i]

    return mod_type_map

def get_color_legend_mods(modification_types_to_color_map):
    handles = []
    for mod_type in modification_types_to_color_map:
        handles.append(patches.Patch(color=modification_types_to_color_map[mod_type], label=mod_type))
    return handles

def get_selected_sample_columns(sample_columns, selected_sample_indices):
    samples = []
    for i in selected_sample_indices:
        samples.append(sample_columns[i-1])
    return samples