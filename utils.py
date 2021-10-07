from textwrap import fill
import requests
import pandas as pd
import re
from matplotlib import pyplot as plt
import numpy as np
<<<<<<< HEAD
import matplotlib.patches as patches
import seaborn as sns
=======
import json

>>>>>>> cc6a8e746f0730275cea3e992e9de10fd75f4bf7
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

def count_no_of_modifications(ptm_str):
    #check if NaN value
    if pd.isnull(ptm_str):
        return None
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

# df = get_data_and_remove_unwanted_columns()
# df = sanitize_data(df)
# df1, df2, df3, df4 = split_data_in_samples(df)
# print(df1)
def combine_and_aggregate_sample_PTM_in_dataframe(df1,df2,df3,df4):
    maxs = [df1['Area Sample 1'].max(), df2['Area Sample 2'].max(), df3['Area Sample 3'].max(), df4['Area Sample 4'].max()]
    max = np.max(maxs)

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


def get_position_of_mass_shift(input_string):
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

# indicies have sign in front
def get_position_of_mass_shift_and_sign(input_string):
    # get charachters after and before mass shift
    ls = input_string.split('(')
    if len(ls) == 0:
        return np.nan
    #replace everythting witthin parenthesis with "#"
    modified_string = re.sub(r"\([^()]*\)", "#", input_string)
    modified_string = modified_string[2:-2] #remove first and last splice sites
    #get indices of "#"
    indices = [m.start()-1 for m in re.finditer("#", modified_string)]
    # add sign value in front of every index determined by value in ls 
    counter = 0
    for idx in indices:
        if input_string[idx+4] == '-' : #4 is 2 removed indices and + for the "(" and the sign
            indices[counter] = -indices[counter]
        counter += 1
    return indices 

def plot_peptide_segments_and_mass_shift(data, delta=0.5, _protein="P02666"):
    """data is a dictionary, {"Label":(low,hi), ... }
    return a drawing that you can manipulate, show, save etc"""
    yplaces = [.5+i for i in range(len(data))]
    #sort map by values, lowest first
    protein_labels = sorted(data.keys(), key=lambda x: data[x][0], reverse=True)
    fig = plt.figure(figsize=(30,25))
    ax = fig.add_subplot(111)
    ax.set_yticks(yplaces)
    ax.set_yticklabels(protein_labels)
    ax.set_ylim((0,len(data)))

    low, hi, _, ms_pos=  data[protein_labels[0]][0]
    for pos, label in zip(yplaces,protein_labels):
        start, end, peptide, ms_pos = data[label][0]
        start, end = start-1, end-1
        # add recatangle
        ax.add_patch(patches.Rectangle((start,pos-delta/2.0),end-start+1, delta, facecolor = '#add8e6', alpha=0.5))  #light blue      
        ax.text(start, pos, label, ha='right', va='center') #add protetin label to the left og the rectatngle
        for _ms in ms_pos: #add mass shift color on rectangles if present
            if _ms > 0:
                ax.add_patch(patches.Rectangle((start+_ms,pos-delta/2.0),1,delta, color="green"))
            else:
                _ms = (-1)*_ms
                ax.add_patch(patches.Rectangle((start+_ms,pos-delta/2.0),1,delta, color="red"))
        if start<low : low=start
        if end>hi : hi=end

    ax.plot((low,hi),(0,0))
    ax.grid(axis='x')
    seqq = get_protein_sequence(_protein)
    ax.set_xlabel("Full Protein Sequence")
    ax.set_xticks(range(0,len(seqq)))
    # create list of chars from string 
    protein_seq_list = list(seqq)
    ax.set_xticklabels(protein_seq_list)

    plt.show()
    return ax

def preprocess_data_for_peptide_segment_plot(df, _protein="P02666", size=50):
    # get position of mass shift in "peptide" for each row
    df["Position of Mass Shift"] = df["Peptide"].apply(get_position_of_mass_shift_and_sign)
    # this is the main script, note that we have imported pyplot as plt
    start_end_df = df[["Start", "End", "Protein Accession", "Peptide", 'Position of Mass Shift']]
    #only look at values for protein : P02666
    start_end_df = start_end_df[start_end_df["Protein Accession"] == _protein]
    start_end_df.sort_values('Start', inplace=True)
    start_end_df['index1'] = start_end_df.index
    #concat index1 and protein accession
    start_end_df['Protein_Accession_idx'] = start_end_df['Protein Accession'] +"_" + start_end_df['index1'].astype(str) 
    start_end_df["(start,end,peptide,pos_ms)"] = start_end_df[["Start", "End", "Peptide",'Position of Mass Shift']].apply(tuple, axis=1)
    start_end_df.drop(["Start", "End", "index1"], axis=1, inplace=True)
    start_end_df.sort_values('Protein_Accession_idx', inplace=True)
    new = start_end_df.head(size)
    # make dictionary with index as keys and (Start,End) as values
    data = new.groupby("Protein_Accession_idx").apply(lambda x: x["(start,end,peptide,pos_ms)"].tolist())
    return data 

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
print(get_position_of_mass_shift_and_sign("K.jkfnekj(-12)8787"))
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

def normalize_intensities_by_protein_intensity(df):
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
        intensity_sum1 = protein_df['Area Sample 1'].sum()
        intensity_sum2 = protein_df['Area Sample 2'].sum()
        intensity_sum3 = protein_df['Area Sample 3'].sum()
        intensity_sum4 = protein_df['Area Sample 4'].sum()
        protein_df['Area Sample 1'] = protein_df['Area Sample 1'].divide(intensity_sum1)
        protein_df['Area Sample 2'] = protein_df['Area Sample 2'].divide(intensity_sum2)
        protein_df['Area Sample 3'] = protein_df['Area Sample 3'].divide(intensity_sum3)
        protein_df['Area Sample 4'] = protein_df['Area Sample 4'].divide(intensity_sum4)

        dataframes.append(protein_df)
    
    return dataframes

def combine_and_aggregate_intensity(df1, df2, df3, df4):
    df1 = df1[df1["PTM"].notnull()]
    df1 = df1.copy()
    df1['Intensity'] = df1['Area Sample 1'].divide(df1['#modifications'])
    df1['PTM'] = df1['PTM'].str.split(';').str[0]
    df1_new = df1[['PTM', 'Intensity']]
    df1_new = df1_new.groupby(['PTM']).sum()
    df1_new = df1_new.reset_index()

    df2 = df2[df2["PTM"].notnull()]
    df2 = df2.copy()
    df2['Intensity'] = df2['Area Sample 2'].divide(df2['#modifications'])
    df2['PTM'] = df2['PTM'].str.split(';').str[0]
    df2_new = df2[['PTM', 'Intensity']]
    df2_new = df2_new.groupby(['PTM']).sum()
    df2_new = df2_new.reset_index()

    df3 = df3[df3["PTM"].notnull()]
    df3 = df3.copy()
    df3['Intensity'] = df3['Area Sample 3'].divide(df3['#modifications'])
    df3['PTM'] = df3['PTM'].str.split(';').str[0]
    df3_new = df3[['PTM', 'Intensity']]
    df3_new = df3_new.groupby(['PTM']).sum()
    df3_new = df3_new.reset_index()

    df4 = df4[df4["PTM"].notnull()]
    df4 = df4.copy()
    df4['Intensity'] = df4['Area Sample 4'].divide(df4['#modifications'])
    df4['PTM'] = df4['PTM'].str.split(';').str[0]
    df4_new = df4[['PTM', 'Intensity']]
    df4_new = df4_new.groupby(['PTM']).sum()
    df4_new = df4_new.reset_index()

    df1_new['Sample'] = 1
    df2_new['Sample'] = 2
    df3_new['Sample'] = 3
    df4_new['Sample'] = 4

    combined = pd.concat([df1_new[['PTM', 'Intensity', 'Sample']],
                          df2_new[['PTM', 'Intensity', 'Sample']],
                          df3_new[['PTM', 'Intensity', 'Sample']],
                          df4_new[['PTM', 'Intensity', 'Sample']]], axis=0)
    return combined
    

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
