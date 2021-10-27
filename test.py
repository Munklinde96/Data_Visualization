from utils import get_data_and_remove_unwanted_columns, sanitize_data, get_protein_sequence, get_position_of_mass_shift
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.patches as patches



# create modification_types to color mapping
def get_color_palette_for_modifications ():
    modification_types = ['Oxidation (M)', 'Phosphorylation (STY)', 'Deamidation (NQ)', 'lal',
                      'Lactosylation', 'Pyro-glu from Q', 'Glycosylation type b', 'Dioxidation (M)',
                    'Glycosylation type e', 'Glycosylation type a','Glycosylation type c/d', 'Carbamidomethylation', 'lan']
    modification_types_to_color = {}
    color_palette = sns.color_palette("tab10", n_colors=len(modification_types)) 
    for i in range(len(modification_types)):
        modification_types_to_color[modification_types[i]] = color_palette[i]
    return modification_types_to_color

def get_peptide_segments_and_modifications(data, delta=0.5, _protein="P02666"):
    """data is a list of tuples on the form (low,hi, [modifications], [modtypes])"""
    yplaces = [.5+i for i in range(len(data))]
    data = sorted(data, key=lambda x: x[0], reverse=False)
    modification_types_to_color_map = get_color_palette_for_modifications()
    rectatngles = []
    modifications = []
    highest_values_for_index = [0]*len(data)
    position_idx = 0 #value to keep track of where we should place next rectangle
    for i in range (len(data)):
        low, hi, mod_positions, mod_types = data[i]
        low, hi = low-1, hi-1 #convert to zero based index
        position_idx = i 
        if i == 12:
            print("")
        for j in range (len(data)):
            if j >= i:
                break
            hi_j = highest_values_for_index[j]
            if hi_j < low:
                position_idx = j
                break
        pos = yplaces[position_idx]
        highest_values_for_index[position_idx] = hi
        rectatngles.append(patches.Rectangle((low,pos-delta/2.0),hi-low+1, delta, facecolor = '#add8e6', ec='black', lw=3, alpha=0.5))
        #add modification type-color to mass shift position
        if len(mod_positions) > 0:
            for _ms_pos, mod_type in zip(mod_positions, mod_types): #add mass shift color on rectangles if present
                ms_color = modification_types_to_color_map[mod_type]
                modifications.append(patches.Rectangle((low+_ms_pos,pos-delta/2.0),1,delta, facecolor= ms_color, edgecolor="black"))
    #remove 0's from highest_values_for_index
    highest_values_for_index = [x for x in highest_values_for_index if x != 0]
    height = len(highest_values_for_index)
    return rectatngles, modifications, height

def plot_peptide_segments(segments_patches, modifications_patches, height, _protein="P02666"):
    fig = plt.figure(figsize=(40,25))
    ax = fig.add_subplot(111)
    ax.set_ylim((0,height))

    for rect in segments_patches:
        ax.add_patch(rect)
    for mod in modifications_patches:
        ax.add_patch(mod)
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
    plt.show()
    return ax

def get_color_legend_mods(modification_types_to_color_map):
    handles = []
    for mod_type in modification_types_to_color_map:
        handles.append(patches.Patch(color=modification_types_to_color_map[mod_type], label=mod_type))
    return handles
    
def preprocess_data_for_peptide_segment_plot_new(df, _protein="P02666", size=50):
    # get position of mass shift in "peptide" for each row
    df["Position of Mass Shift"] = df["Peptide"].apply(get_position_of_mass_shift)
    # get list of modification for each PTM
    df["Modification_types"] = df["PTM"].apply(lambda x: x if pd.isnull(x) else [s.strip() for s in x.split(";")])

    # this is the main script, note that we have imported pyplot as plt
    start_end_df = df[["Start", "End", "Protein Accession", "Peptide", 'Position of Mass Shift', 'PTM', 'Modification_types']]
    #only look at values for protein : P02666
    start_end_df = start_end_df[start_end_df["Protein Accession"] == _protein]
    start_end_df.sort_values('Start', inplace=True)
    start_end_df['index1'] = start_end_df.index
    #concat index1 and protein accession
    start_end_df['Protein_Accession_idx'] = start_end_df['Protein Accession'] +"_" + start_end_df['index1'].astype(str) 
    start_end_df["(start,end,pos_ms,mod_types)"] = start_end_df[["Start", "End", 'Position of Mass Shift', 'Modification_types']].apply(tuple, axis=1)
    start_end_df.drop(["Start", "End", "index1", 'PTM','Modification_types'], axis=1, inplace=True)
    start_end_df.sort_values('Protein_Accession_idx', inplace=True)
    new = start_end_df.head(size)
    # make dictionary with index as keys and (Start,End) as values
    #data = new.groupby("Protein_Accession_idx").apply(lambda x: x["(start,end,peptide,pos_ms)"].tolist())
    start_end_ms_modtype_list = new['(start,end,pos_ms,mod_types)'].tolist()
    return start_end_ms_modtype_list


from utils import get_data_and_remove_unwanted_columns, sanitize_data, preprocess_data_for_peptide_segment_plot, get_rectangles_for_peptides_and_mods, plot_peptide_segments, normalize_intensities_by_protein_intensity, get_protein_sequence, get_overlap_overlaps_by_intensity_and_sample, get_overlap_pixel_plot, get_gradient_plot, get_stacking_patches, stack_recs
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

df = get_data_and_remove_unwanted_columns()
df = sanitize_data(df)
dfs = normalize_intensities_by_protein_intensity(df)

data = preprocess_data_for_peptide_segment_plot(df, size=df.shape[0])
rectangles_and_mods = get_rectangles_for_peptides_and_mods(data)
rectangles_and_mods = stack_recs(rectangles_and_mods)
peptide_patches, mod_patches, height = get_stacking_patches(rectangles_and_mods)
plot_peptide_segments(peptide_patches, mod_patches, height)

