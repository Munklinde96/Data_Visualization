from textwrap import fill
import requests
import pandas as pd
import re
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.patches as patches
import seaborn as sns
import json

from seaborn.palettes import color_palette

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


# create modification_types to color mapping
def get_color_palette_for_modifications ():
    mod_type_map = {'Oxidation (M)': '#E6194B' , 'Phosphorylation (STY)': '#F58231', 'Deamidation (NQ)': '#FFE119', 'lal': '#BFEF45',
                      'Lactosylation': '#3CB44B', 'Pyro-glu from Q': '#42D4F4', 'Glycosylation type b': '#4363D8', 'Dioxidation (M)': '#911EB4',
                    'Glycosylation type e': '#F032E6', 'Glycosylation type a': '#000000','Glycosylation type c/d': '#800000', 'Carbamidomethylation': '#FABED4', 'lan': '#808000'}
    return mod_type_map

def get_rectangles_for_peptides_and_mods(data):
    """data is a list of tuples on the form (low,hi, [modifications], [modtypes], [agg intensity], [quintile])"""
    data = sorted(data, key=lambda x: (x[0], -x[1]), reverse=False)
    modification_types_to_color_map = get_color_palette_for_modifications()
    rectangles = []
    for i in range (len(data)):
        modifications = []
        low, hi, mod_positions, mod_types, intensity, quintile = data[i]
        width = hi-low+1
        rec = (low, width)
        if len(mod_positions) > 0:
            for _ms_pos, mod_type in zip(mod_positions, mod_types): #add mass shift color on rectangles if present
                ms_color = modification_types_to_color_map[mod_type]
                modifications.append((_ms_pos, ms_color))
        rectangles.append((rec, modifications))
    return rectangles


#Adjsut height på baggrund af intensitet - Lav hver index i y_spaces til min højden,
# hvis et peptid fylder mere, indsæt x + width + spacing i 1, 3, 5, 7, 9 indices af y_spaces
def get_stacking_patches(rectangles):
    last_x = 0
    last_y = 0
    last_rec_width = 0
    y_spaces = [0]*5*len(rectangles)
    peptide_patches = []
    mod_patches = []
    for rect, mods in rectangles:
        if(rect[0] != last_x and rect[1] != last_rec_width):
            x, width = rect
            patch = None
            for i in range(len(y_spaces)):
                if(y_spaces[i] <= x):
                    patch = patches.Rectangle((x, i), width, 0.5, facecolor = '#add8e6', ec='black', lw=3, alpha=0.5)
                    y_spaces[i] = x + width
                    break

            peptide_patches.append(patch)
            last_x, last_y = patch.get_xy()
            last_rec_width = width

        for m_pos, m_color in mods:
            mod_rec = patches.Rectangle((last_x+m_pos,last_y),1,0.5, facecolor= m_color, edgecolor="black")
            mod_patches.append(mod_rec)
    height = y_spaces.index(0)
    return peptide_patches, mod_patches, height

def plot_peptide_segments(peptide_patches, mod_patches, height, _protein="P02666"):
    fig = plt.figure(figsize=(30,25))
    ax = fig.add_subplot(111)
    ax.set_ylim((0,height))

    for patch in peptide_patches:
        ax.add_patch(patch)

    for patch in mod_patches:
        ax.add_patch(patch)

    ax.grid(axis='x')
    seqq = get_protein_sequence(_protein)
    ax.set_xlabel("Full Protein Sequence")
    ax.set_xticks(range(0,len(seqq)))
    # create list of chars from string 
    protein_seq_list = list(seqq)
    ax.set_xticklabels(protein_seq_list)
    ax.get_yaxis().set_visible(False)
    modification_types_to_color_map = get_color_palette_for_modifications()
    # create legend based on modification_types_to_color_map
    handles = get_color_legend_mods(modification_types_to_color_map)
    ax.legend(handles=handles)
    # plt.colorbar(modification_types_to_color_map, ticks=list(modification_types_to_color_map.keys()))
    #plt.gca().invert_yaxis()
    plt.show()
    return ax

def get_color_legend_mods(modification_types_to_color_map):
    handles = []
    for mod_type in modification_types_to_color_map:
        handles.append(patches.Patch(color=modification_types_to_color_map[mod_type], label=mod_type))
    return handles
    
def preprocess_data_for_peptide_segment_plot(df, _protein="P02666", size=50):
    # get position of mass shift in "peptide" for each row
    df["Position of Mass Shift"] = df["Peptide"].apply(get_position_of_mass_shift)
    # get list of modification for each PTM
    df["Modification_types"] = df["PTM"].apply(lambda x: x if pd.isnull(x) else [s.strip() for s in x.split(";")])

    # this is the main script, note that we have imported pyplot as plt
    start_end_df = df[["Start", "End", "Protein Accession", "Peptide", 'Position of Mass Shift', 'PTM', 'Modification_types', 'Area Sample 1', 'Area Sample 2', 'Area Sample 3', 'Area Sample 4']]
    #only look at values for protein : P02666
    start_end_df = start_end_df[start_end_df["Protein Accession"] == _protein]
    # start_end_df.sort_values(['Start', 'End'], ascending=[True, False], inplace=True)
    start_end_df['index1'] = start_end_df.index

    #Aggregate sample intensity, and normalize it
    start_end_df['Agg Intensity'] = start_end_df['Area Sample 1'] + start_end_df['Area Sample 2'] + start_end_df['Area Sample 3'] + start_end_df['Area Sample 4']
    start_end_df['Agg Intensity'] = start_end_df['Agg Intensity'] / start_end_df['Agg Intensity'].sum()

    start_end_df['quintile'] = pd.qcut(start_end_df['Agg Intensity'], q=5)
    intervals = pd.unique(start_end_df['quintile'])
    print(intervals)

    #concat index1 and protein accession
    start_end_df['Protein_Accession_idx'] = start_end_df['Protein Accession'] +"_" + start_end_df['index1'].astype(str) 
    start_end_df["(start, end, pos_ms, mod_types, agg_intensity, quintile)"] = start_end_df[["Start", "End", 'Position of Mass Shift', 'Modification_types', 'Agg Intensity', 'quintile']].apply(tuple, axis=1)
    start_end_df.drop(["Start", "End", "index1", 'PTM','Modification_types', 'Area Sample 1', 'Area Sample 2', 'Area Sample 3', 'Area Sample 4'], axis=1, inplace=True)
    start_end_df.sort_values('Protein_Accession_idx', inplace=True)
    new = start_end_df.head(size)

    start_end_ms_modtype_list = new['(start, end, pos_ms, mod_types, agg_intensity, quintile)'].tolist()
    return start_end_ms_modtype_list

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

# nomalize peptide intensity(per sample) over protetin intensity
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

def get_overlap_heapmap(num_overlpas_lists, peptide_seq_list, protein_num, fig_size=(30,10), color_scale='YlOrRd'):
    plt.figure(figsize=fig_size)
    plt.title(f"Frequency of Overlaps for Protein {protein_num} - sample 1,2,3,4")
    ax = sns.heatmap(num_overlpas_lists, cmap=color_scale)
    plt.xticks(np.arange(len(peptide_seq_list)), peptide_seq_list, rotation = 0)
    ylabels = ["Sample 1", "Sample 2", "Sample 3", "Sample 4"]
    ax.set_yticklabels(ylabels)
    plt.show()

def get_overlap_gradient_heapmap(num_overlpas_lists, peptide_seq_list, protein_num, fig_size=(30,10), color_scale='YlOrRd'):
    gradient_list = []
    for i in range(len(num_overlpas_lists)):
        gradient_list.append( abs(np.diff(np.asarray(num_overlpas_lists[i]).reshape(1, -1)[::-1])))
    gradient_list = np.asarray(gradient_list).reshape(4, -1)
    plt.figure(figsize=fig_size)
    plt.title(f"Gradient plot for {protein_num} - Shows frequent clevage sites")
    ax = sns.heatmap(gradient_list , cmap=color_scale)
    plt.xticks(np.arange(len(peptide_seq_list)), peptide_seq_list, rotation = 0)
    ylabels = ["Sample 1", "Sample 2", "Sample 3", "Sample 4"]
    ax.set_yticklabels(ylabels)
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

