from numpy.lib.function_base import append
import requests
import pandas as pd
import re
from matplotlib import pyplot as plt
import numpy as np
import json

uniprot = 'https://www.uniprot.org/uniprot/'

def heatmap2d(arr: np.ndarray):
    plt.imshow(arr, cmap='viridis')
    plt.colorbar()
    plt.show()

def get_protein_sequence(protein):
    url = uniprot + protein + '.fasta'
    response = requests.get(url).text
    if(response):
        first_new_line = response.find('\n')
        sequence = response[first_new_line:].replace('\n', '')
        return sequence
    else:
        return ''

def get_data_and_remove_unwanted_columns():
    df = pd.read_csv('UHT milk P036.csv')
    df.drop('Protein ID', inplace=True, axis=1)
    df.drop('Unique', inplace=True, axis=1)
    df.drop('m/z', inplace=True, axis=1)
    df.drop('Scan', inplace=True, axis=1)
    df.drop('Source File', inplace=True, axis=1)
    df.drop('Found By', inplace=True, axis=1)
    #apply count_no_of_modifications to each PTM column
    df['#modifications'] = df['PTM'].apply(count_no_of_modifications)

    print(df.columns)
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

def get_cleavage_sites_for_peptides(df):
    start_cleavage = []
    end_cleavage = []
    protein_sequence = ''
    current_protein = ''
    for index, row in df.iterrows():
        if(current_protein != row['Protein Accession']):
            current_protein = row['Protein Accession']
            protein_sequence = get_protein_sequence(current_protein)
        if protein_sequence[row['Start']-5 : row['Start'] -1] != "":
            start_cleavage.append(f"{protein_sequence[row['Start']-5 : row['Start'] -1]}.{protein_sequence[row['Start'] -1 : row['Start'] + 3]}")
        else:
            start_cleavage.append("")
        if protein_sequence[row['End'] : row['End'] + 4] != "":
            end_cleavage.append(f"{protein_sequence[row['End']-4 : row['End']]}.{protein_sequence[row['End'] : row['End'] + 4]}")
        else: 
            end_cleavage.append("")
    df['Start Cleavage'] = start_cleavage
    df['End Cleavage'] = end_cleavage


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

import numpy as np
def count_no_of_modifications(ptm_str):
    #check if NaN value
    if pd.isnull(ptm_str):
        return 0
    return 1 + ptm_str.count(';')

def mass_div_len_column(mass, length):
    return mass/length

def split_data_in_samples(df):
    df1 = []
    df2 = []
    df3 = []
    df4 = []
    for i, row in df.iterrows():
        if row['Area Sample 1'] != 0 and row['Area Sample 1'] != None:
            df1.append(row)
        if row['Area Sample 2'] != 0 and row['Area Sample 2'] != None:
            df2.append(row)
        if row['Area Sample 3'] != 0 and row['Area Sample 3'] != None:
            df3.append(row)
        if row['Area Sample 4'] != 0 and row['Area Sample 4'] != None:
            df4.append(row)
    df1 = pd.DataFrame(df1, columns=df.columns)
    df2 = pd.DataFrame(df2, columns=df.columns)
    df3 = pd.DataFrame(df3, columns=df.columns)
    df4 = pd.DataFrame(df4, columns=df.columns)
    return df1, df2, df3, df4

def combine_and_aggregate_sample_PTM_in_dataframe(df1,df2,df3,df4):
    df1['PTM'] = df1['PTM'].str.split(';').str[0]
    df1_PTM_count = df1['PTM'].value_counts()
    df1_new = pd.DataFrame()
    df1_new['PTM'] = df1_PTM_count.index
    df1_new['#PTM'] = df1_PTM_count.values

    df2['PTM'] = df2['PTM'].str.split(';').str[0]
    df2_PTM_count = df2['PTM'].value_counts()
    df2_new = pd.DataFrame()
    df2_new['PTM'] = df2_PTM_count.index
    df2_new['#PTM'] = df2_PTM_count.values

    df3['PTM'] = df3['PTM'].str.split(';').str[0]
    df3_PTM_count = df3['PTM'].value_counts()
    df3_new = pd.DataFrame()
    df3_new['PTM'] = df3_PTM_count.index
    df3_new['#PTM'] = df3_PTM_count.values

    df4['PTM'] = df4['PTM'].str.split(';').str[0]
    df4_PTM_count = df4['PTM'].value_counts()
    df4_new = pd.DataFrame()
    df4_new['PTM'] = df4_PTM_count.index
    df4_new['#PTM'] = df4_PTM_count.values

    df1_new['Sample'] = 1
    df2_new['Sample'] = 2
    df3_new['Sample'] = 3
    df4_new['Sample'] = 4

    combined = pd.concat([df1_new[['PTM', '#PTM', 'Sample']],
                          df2_new[['PTM', '#PTM', 'Sample']],
                          df3_new[['PTM', '#PTM', 'Sample']],
                          df4_new[['PTM', '#PTM', 'Sample']]], axis=0)
    return combined

def get_protein_length_from_uniprot(protein):
    url = "https://www.ebi.ac.uk/proteins/api/proteins/"+ protein
    response = requests.get(url).text
    if(response):
        proteinJson = json.loads(response)
        proteinLength = int(proteinJson["sequence"]["length"])
        if proteinLength is 0:
            print("protein lenth was 0")
        return proteinLength
    else:
        print("no response")
        return 0

def get_protein_mass_from_uniprot(protein):
    url = "https://www.ebi.ac.uk/proteins/api/proteins/"+ protein
    response = requests.get(url).text
    if(response):
        proteinJson = json.loads(response)
        proteinMass = int(proteinJson["sequence"]["mass"])
        if proteinMass is 0:
            print("protein mass was 0")
        return proteinMass
    else:
        print("no response")
        return 0

def get_modification_count_per_protein(df, countFilter, normalize):
    df_protein_mods = df[["PTM", "Protein Accession"]]
    modificationCountByProtein = {}
    totalProteinModCount = {}
    for modString, proteinName in df_protein_mods.itertuples(index=False):
        if pd.isnull(modString):
            continue
        proteinName = proteinName.strip()
        modString = modString.strip()
        if proteinName not in modificationCountByProtein:
            modificationCountByProtein[proteinName] = {}
            totalProteinModCount[proteinName] = 0
        mods = modString.split(";")
        for mod in mods:
            mod = mod.strip()
            if mod not in modificationCountByProtein[proteinName]:
                 modificationCountByProtein[proteinName][mod] = 1
            else:
                 modificationCountByProtein[proteinName][mod] += 1
                 totalProteinModCount[proteinName] += 1
    modificationCountByProteinFiltered = {}       
    for proteinName, mods in modificationCountByProtein.items():
        if totalProteinModCount[proteinName] > countFilter:
            modificationCountByProteinFiltered[proteinName] = mods
    
    if normalize is "protein_total_mod_count":
        for protein, mods in modificationCountByProteinFiltered.items():
            print(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / totalProteinModCount[protein]
            modificationCountByProteinFiltered[protein] = updateMods
    elif normalize is "protein_amino_acid_length":
        for protein, mods in modificationCountByProteinFiltered.items():
            sequence = get_protein_sequence(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / len(sequence)
            modificationCountByProteinFiltered[protein] = updateMods
    elif normalize is "protein_length":
        for protein, mods in modificationCountByProteinFiltered.items():
            proteinLength = get_protein_length_from_uniprot(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / proteinLength
            modificationCountByProteinFiltered[protein] = updateMods
    elif normalize is "protein_mass":
        for protein, mods in modificationCountByProteinFiltered.items():
            proteinMass = get_protein_mass_from_uniprot(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / proteinMass
            modificationCountByProteinFiltered[protein] = updateMods
    return modificationCountByProteinFiltered

def get_modification_count_per_protein_reverse(df, countFilter, normalize):
    df_protein_mods = df[["PTM", "Protein Accession"]]
    modificationCountByProtein = {}
    totalProteinModCount = {}
    for modString, proteinName in df_protein_mods.itertuples(index=False):
        if pd.isnull(modString):
            continue
        proteinName = proteinName.strip()
        modString = modString.strip()
        if proteinName not in modificationCountByProtein:
            modificationCountByProtein[proteinName] = {}
            totalProteinModCount[proteinName] = 0
        mods = modString.split(";")
        for mod in mods:
            mod = mod.strip()
            if mod not in modificationCountByProtein[proteinName]:
                 modificationCountByProtein[proteinName][mod] = 1
            else:
                 modificationCountByProtein[proteinName][mod] += 1
                 totalProteinModCount[proteinName] += 1
    modificationCountByProteinFiltered = {}       
    for proteinName, mods in modificationCountByProtein.items():
        if totalProteinModCount[proteinName] < countFilter:
            modificationCountByProteinFiltered[proteinName] = mods
    # different modifications
    if normalize is "protein_total_mod_count":
        for protein, mods in modificationCountByProteinFiltered.items():
            print(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / totalProteinModCount[protein]
            modificationCountByProteinFiltered[protein] = updateMods
    elif normalize is "protein_amino_acid_length":
        for protein, mods in modificationCountByProteinFiltered.items():
            sequence = get_protein_sequence(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / len(sequence)
            modificationCountByProteinFiltered[protein] = updateMods
    elif normalize is "protein_length":
        for protein, mods in modificationCountByProteinFiltered.items():
            proteinLength = get_protein_length_from_uniprot(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / proteinLength
            modificationCountByProteinFiltered[protein] = updateMods
    elif normalize is "protein_mass":
        for protein, mods in modificationCountByProteinFiltered.items():
            proteinMass = get_protein_mass_from_uniprot(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / proteinMass
            modificationCountByProteinFiltered[protein] = updateMods
    return modificationCountByProteinFiltered