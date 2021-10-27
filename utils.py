from textwrap import fill
import requests
import pandas as pd
import re
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.patches as patches
import seaborn as sns
import json
from refactored_utils import get_protein_sequence, clean_peptide

from seaborn.palettes import color_palette


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
            start_cleavage.append(f" {protein_sequence[row['Start']-5 : row['Start'] -1]}-{protein_sequence[row['Start'] -1 : row['Start'] + 3]} ")
        else:
            start_cleavage.append("")
        if protein_sequence[row['End'] : row['End'] + 4] != "":
            end_cleavage.append(f" {protein_sequence[row['End']-4 : row['End']]}-{protein_sequence[row['End'] : row['End'] + 4]} ")
        else: 
            end_cleavage.append("")
    df['Start Cleavage'] = start_cleavage
    df['End Cleavage'] = end_cleavage
    return df, start_cleavage, end_cleavage

def mass_div_len_column(mass, length):
    return mass/length


def check_if_string_contains_substring_x_times(string, substring, no_times, exact = False):
    if exact:
        return string.count(substring) == no_times
    elif not exact and  string.count(substring) > no_times:
        return True
    else:
        return False

def get_all_rows_with_at_least_x_modifications(df, no_modifications):
    df = df[df['PTM'].notna()]
    return df[df['PTM'].apply(check_if_string_contains_substring_x_times, args =(';', no_modifications-1, False) )] #-1 because of the ";" in the PTM column

def get_all_rows_with_exactly_x_modifications(df, no_modifications):
    df = df[df['PTM'].notna()]
    return df[df['PTM'].apply(check_if_string_contains_substring_x_times, args =(';', no_modifications-1,True) )] #-1 because of the ";" in the PTM column

def get_mass_shift_per_peptide(string):
    ls = string.split('(')
    if len(ls) == 0:
        return np.nan
    res = []
    for chunch in ls:
        res.append(chunch.split(')')[0])
    return ",".join(res)

def create_mass_shift_column(df):
    #remove nan values
    df['MassShift'] = df['Peptide'].apply(get_mass_shift_per_peptide)
    return df

#Used for boxplots
def add_value_labels(ax, spacing=1):
    for rect in ax.patches:
        y_value = rect.get_height()
        x_value = rect.get_x() + rect.get_width() / 2
        space = spacing
        va='bottom'
        if y_value < 0:
            space *= -1
            va='top'
        label='{:.2f}'.format(y_value)
        #creat annotation
        ax.annotate(label,(x_value,y_value),xytext=(0,space),textcoords='offset points',ha='center',va=va)
        ax.axhline(y=0.0, color='black', linestyle='-', linewidth=2)


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


#OLD VERSION
def get_position_of_mass_shifts(input_string):
    # get charachters after and before mass shift
    ls = input_string.split('(')
    if len(ls) == 0:
        return np.nan
    firs_chunch = [ls[0]]
    ls = firs_chunch + [l.split(")")[1] for l in ls[1:]]
    before_ms = [l for l in ls[::2]]# get every second element
    after_ms= [l for l in ls[1::2]]# get every second element starting from 1
    #replace everythting witthin parenthesis with "#"
    modified_string = re.sub(r"\([^()]*\)", "#", input_string)
    modified_string = modified_string[2:-2] #remove first and last splice sites
    #get indices of "#"
    indices = [m.start() for m in re.finditer("#", modified_string)]
    return indices, before_ms, after_ms

# get colour palette from y-value distribution
def colors_from_values(values, palette_name):
    # normalize the values to range [0, 1]
    normalized = (values - min(values)) / (max(values) - min(values))
    # convert to indices
    indices = np.round(normalized * (len(values) - 1)).astype(np.int32)
    # use the indices to get the colors
    palette = sns.color_palette(palette_name, len(values))
    return np.array(palette).take(indices, axis=0)

def plot_overlap_barchart(df, selected_protein= "P02666"):
    df_new = df[df["Protein Accession"] == selected_protein]
    seq_list = list(get_protein_sequence(selected_protein))
    _len = len(seq_list)
    num_overlpas_list = [0]*_len # init as zeroes
    # add 1 into all positions where there is an overlap
    for i in range(len(df_new)):
        for j in range(df_new.iloc[i]['Start'], df_new.iloc[i]['End']):
            num_overlpas_list[j] += 1
    df_overlaps = pd.DataFrame(list(zip(range(_len), num_overlpas_list)), columns=['Position', 'Overlaps'])
    plt.figure(figsize=(30,10))
    g= sns.barplot(x="Position", y= "Overlaps", data= df_overlaps, palette=colors_from_values(np.asarray(num_overlpas_list), "YlOrRd"))
    g.set_xticklabels(seq_list)
    g.set_title(f"Number of overlaps per position - for protein: {selected_protein}")
    return g    


def add_trailing_white_spaces_to_chars(seq_list):
    # equal to number of time char is seen in the sequence
    res_list = seq_list
    char_count_dict = {}
    counter = 0
    for char in seq_list:
        if char in char_count_dict:
            res_list[counter] = char + " "*char_count_dict[char]
            char_count_dict[char] += 1
        else:
            char_count_dict[char] = 1
        counter += 1
    return res_list


#test get_position_of_mass_shift_and_sign
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
    

def get_protein_total_intensity(df, protein):
    intensity = 0
    hasSeen = False
    df_protein_intensity = df[["Protein Accession", "Area Sample 1"]]
    df_protein_intensity.sort_values(by="Protein Accession", ascending=False)
    for proteinName, area1 in df_protein_intensity.itertuples(index=False):
        if hasSeen and proteinName != protein:
            break
        if proteinName == protein:
            hasSeen = True
            if not pd.isnull(area1):
                intensity += area1
    hasSeen = False
    df_protein_intensity = df[["Protein Accession", "Area Sample 2"]]
    df_protein_intensity.sort_values(by="Protein Accession", ascending=False)
    for proteinName, area2 in df_protein_intensity.itertuples(index=False):
        if hasSeen and proteinName != protein:
            break
        if proteinName == protein:
            hasSeen = True
            if not pd.isnull(area2):
                intensity += area2
    df_protein_intensity = df[["Protein Accession", "Area Sample 3"]]
    df_protein_intensity.sort_values(by="Protein Accession", ascending=False)
    for proteinName, area3 in df_protein_intensity.itertuples(index=False):
        if hasSeen and proteinName != protein:
            break
        if proteinName == protein:
            hasSeen = True
            if not pd.isnull(area3):
                intensity += area3
    df_protein_intensity = df[["Protein Accession", "Area Sample 4"]]
    df_protein_intensity.sort_values(by="Protein Accession", ascending=False)
    for proteinName, area4 in df_protein_intensity.itertuples(index=False):
        if hasSeen and proteinName != protein:
            break
        if proteinName == protein:
            hasSeen = True
            if not pd.isnull(area4):
                intensity += area4
    
    return intensity

def get_modification_count_per_protein(df, countFilter, normalize):
    df_protein_mods = df[["PTM", "Protein Accession"]]
    print("normalization: "+normalize)
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
    if "protein_total_mod_count" in '{0}'.format(normalize):
        print("norm is: protein_total_mod_count")
        for protein, mods in modificationCountByProteinFiltered.items():
            print(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / totalProteinModCount[protein]
            modificationCountByProteinFiltered[protein] = updateMods
    elif "protein_length" in '{0}'.format(normalize):
        print("norm is: protein_length")
        for protein, mods in modificationCountByProteinFiltered.items():
            proteinLength = get_protein_length_from_uniprot(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / proteinLength
            modificationCountByProteinFiltered[protein] = updateMods
    elif "protein_mass" in '{0}'.format(normalize):
        print("norm is: protein_mass")
        for protein, mods in modificationCountByProteinFiltered.items():
            proteinMass = get_protein_mass_from_uniprot(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / proteinMass
            modificationCountByProteinFiltered[protein] = updateMods
    elif "protein_intensity" in '{0}'.format(normalize):
        print("norm is: protein_intensity")
        for protein, mods in modificationCountByProteinFiltered.items():
            proteinIntensity =  get_protein_total_intensity(df, "P02666")
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / proteinIntensity
            modificationCountByProteinFiltered[protein] = updateMods
    else:
        print("norm is: no normalization")
    return modificationCountByProteinFiltered

#test commit

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
    elif normalize is "protein_intensity":
        for protein, mods in modificationCountByProteinFiltered.items():
            proteinMass = get_protein_mass_from_uniprot(protein)
            updateMods = {}
            for mod, count in mods.items():
                updateMods[mod] = count / proteinMass
            modificationCountByProteinFiltered[protein] = updateMods
    return modificationCountByProteinFiltered


def get_overlap_overlaps_by_intensity_and_sample(df, selected_protein= "P02666"):
    df_new = df[df["Protein Accession"] == selected_protein]
    seq_list = list(get_protein_sequence(selected_protein))
    _len = len(seq_list)
    overlaps_list_1 = [0]*_len
    overlaps_list_2 = [0]*_len
    overlaps_list_3 = [0]*_len
    overlaps_list_4 = [0]*_len

    # add intensity into all positions where there is an overlap
    for i in range(len(df_new)):
        row = df_new.iloc[i]
        for j in range(row['Start'], row['End']):
            if not pd.isnull(row["Area Sample 1"]):
                overlaps_list_1[j] += row['Area Sample 1']
            if not pd.isnull(row["Area Sample 2"]):
                overlaps_list_2[j] += row['Area Sample 2']
            if not pd.isnull(row["Area Sample 3"]):
                overlaps_list_3[j] += row['Area Sample 3']            
            if not pd.isnull(row["Area Sample 4"]):
                overlaps_list_4[j] += row['Area Sample 4']

    df_overlaps1 = pd.DataFrame(list(zip(range(_len), overlaps_list_1)), columns=['Position', 'Overlaps'])
    df_overlaps2 = pd.DataFrame(list(zip(range(_len), overlaps_list_2)), columns=['Position', 'Overlaps'])
    df_overlaps3 = pd.DataFrame(list(zip(range(_len), overlaps_list_3)), columns=['Position', 'Overlaps'])
    df_overlaps4 = pd.DataFrame(list(zip(range(_len), overlaps_list_4)), columns=['Position', 'Overlaps'])

    overlap_lists = (overlaps_list_1, overlaps_list_2, overlaps_list_3, overlaps_list_4)
    overlap_dataframes = (df_overlaps1, df_overlaps2, df_overlaps3, df_overlaps4)

    return overlap_lists, overlap_dataframes

def get_overlap_pixel_plot(num_overlpas_lists, peptide_seq_list, protein_num, fig_size=(30,10), color_scale='YlOrRd'):
    fig, axs = plt.subplots(len(num_overlpas_lists), 1, figsize=fig_size)
    counter=0
    for ls in num_overlpas_lists:
        if counter == 0:
            axs[counter].set_title(f"Frequency of Overlaps for Protein {protein_num} - sample 1,2,3,4")
        im = axs[counter].imshow(np.asarray(ls).reshape(1, -1), cmap=color_scale, extent=[0, len(peptide_seq_list), 0, 10])
        axs[counter].set_xticks(np.arange(len(peptide_seq_list)))
        axs[counter].set_xticklabels(peptide_seq_list)
        axs[counter].set_yticks([])
        axs[counter].set_ylabel(f"Sample {counter+1}")
        
        counter = counter + 1
    fig.colorbar(im, ax=axs, label = "Percentage of Overlab")
    plt.show()


def get_gradient_plot(num_overlpas_lists, peptide_seq_list, protein_num, fig_size=(30,10), color_scale='YlOrRd'):
    fig, axs = plt.subplots(len(num_overlpas_lists), 1, figsize=fig_size)
    counter=0
    for ls in num_overlpas_lists:
        if counter == 0:
            axs[counter].set_title(f"Gradient plot for {protein_num} - Shows frequent clevage sites")
        im = axs[counter].imshow(abs(np.diff(np.asarray(ls).reshape(1, -1)[::-1])), cmap=color_scale, extent=[0, len(peptide_seq_list), 0, 10])
        axs[counter].set_xticks(np.arange(len(peptide_seq_list)))
        axs[counter].set_xticklabels(peptide_seq_list)
        axs[counter].set_yticks([])
        axs[counter].set_ylabel(f"Sample {counter+1}")
        counter = counter + 1
    #set label on colorbar
    fig.colorbar(im, ax=axs, label="Overlab Gradient")
    plt.show()

