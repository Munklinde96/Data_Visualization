import pandas as pd
import matplotlib.patches as patches
from matplotlib import pyplot as plt
from refactored_utils import get_position_of_mass_shift, get_color_palette_for_modifications, get_protein_sequence, get_color_legend_mods

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

    #make Area samples Nan values 0
    start_end_df['Area Sample 1'] = start_end_df['Area Sample 1'].fillna(0) #TODO refactor this
    start_end_df['Area Sample 2'] = start_end_df['Area Sample 2'].fillna(0)
    start_end_df['Area Sample 3'] = start_end_df['Area Sample 3'].fillna(0)
    start_end_df['Area Sample 4'] = start_end_df['Area Sample 4'].fillna(0)
    # drop rows where all area samples are 0
    start_end_df = start_end_df[(start_end_df['Area Sample 1'] != 0) | (start_end_df['Area Sample 2'] != 0) | (start_end_df['Area Sample 3'] != 0) | (start_end_df['Area Sample 4'] != 0)]

    #Aggregate sample intensity, and normalize it
    start_end_df['Agg Intensity'] = start_end_df['Area Sample 1'] + start_end_df['Area Sample 2'] + start_end_df['Area Sample 3'] + start_end_df['Area Sample 4']
    start_end_df['Agg Intensity'] = start_end_df['Agg Intensity'] / start_end_df['Agg Intensity'].sum()

    #concat index1 and protein accession
    start_end_df['Protein_Accession_idx'] = start_end_df['Protein Accession'] +"_" + start_end_df['index1'].astype(str) 
    start_end_df["(start, end, pos_ms, mod_types, agg_intensity)"] = start_end_df[["Start", "End", 'Position of Mass Shift', 'Modification_types', 'Agg Intensity']].apply(tuple, axis=1)
    start_end_df.drop(["Start", "End", "index1", 'PTM','Modification_types', 'Area Sample 1', 'Area Sample 2', 'Area Sample 3', 'Area Sample 4'], axis=1, inplace=True)
    start_end_df.sort_values('Protein_Accession_idx', inplace=True)
    new = start_end_df.head(size)

    start_end_ms_modtype_list = new['(start, end, pos_ms, mod_types, agg_intensity)'].tolist()
    return start_end_ms_modtype_list

def get_rectangles_for_peptides_and_mods(data):
    """data is a list of tuples on the form (low,hi, [modifications],[modtypes], [agg_intensity])"""
    data = sorted(data, key=lambda x: (x[0], -x[1]), reverse=False)
    modification_types_to_color_map = get_color_palette_for_modifications()
    rectangles_and_mods = []
    for i in range (len(data)):
        modifications = []
        low, hi, mod_positions, mod_types, agg_intensity = data[i]
        width = hi-low+1
        rec = ((low, width), agg_intensity)
        if len(mod_positions) > 0:
            for _ms_pos, mod_type in zip(mod_positions, mod_types): #add mass shift color on rectangles if present
                ms_color = modification_types_to_color_map[mod_type]
                modifications.append((_ms_pos, ms_color))
        rectangles_and_mods.append((rec, modifications))
        
    return rectangles_and_mods

def get_intervals(agg_intensities: list, no_intervals = 5):
    intervals = pd.qcut(agg_intensities, q=no_intervals)
    intervals = pd.unique(intervals)
    intervals = intervals[::-1] #reverse intervals 
    interval_len = len(intervals)
    intervals = {intervals[i]: (i+1)/interval_len for i in range(interval_len)}
    return intervals

def map_to_intervals(res_intensities, res_rectangles_and_mods, predefined_intervals: dict = None):
    if predefined_intervals is None:
        intervals = get_intervals(res_intensities)
    
    rects = [r[0] for r in res_rectangles_and_mods]
    res_intensities = pd.Series(res_intensities)
    #map res_intensities to quantiles
    res_intensities = res_intensities.map(intervals)
    rects = [(r[0],r[1]) for r in rects]
    rects_and_quantiles = list(zip(rects, res_intensities))
    mods = [r[1] for r in res_rectangles_and_mods]
    rects_and_quantiles_mods = list(zip(rects_and_quantiles, mods))
    return rects_and_quantiles_mods

def stack_recs(rectangles_and_mods: list):
    """
    rectangles_and_mods is a tuple list of (rectangle, modifications)
    a rectangle is tuples on the form (low, width, agg_intensity)
    rectabgles are sorted after both start_pos (low) and end (width)
    """
    last_x = 0
    last_rec_width = 0
    res_intensities = []
    res_rectangles_and_mods = []

    for rect, mods in rectangles_and_mods:
        rect_info, intensity = rect
        low, width = rect_info
        if low != last_x and width != last_rec_width:
            res_rectangles_and_mods.append(((low, width, intensity), mods))
            res_intensities.append(intensity)
            last_x = low
            last_rec_width = width
        else:
            last_rect, last_mods = res_rectangles_and_mods[-1]
            last_rect_low, last_rect_width, last_intensity = last_rect
            new_intensity = last_intensity + intensity
            new_mods = last_mods + mods
            res_rectangles_and_mods[-1] = ((last_rect_low, last_rect_width, new_intensity), new_mods)
            res_intensities[-1] = new_intensity

    rects_and_quantiles_mods = map_to_intervals( res_intensities, res_rectangles_and_mods)
    return rects_and_quantiles_mods

def get_stacking_patches(rectangles, spacing=0.2):
    last_x = 0
    last_y = 0
    y_spaces = [0]*10*len(rectangles) #10 indicates a division into heights of 0.1
    peptide_patches = []
    mod_patches = []
    for rect, mods in rectangles:
        rect_info, quantile = rect
        x, width = rect_info
        for i in range(len(y_spaces)):
            is_available = True
            for k in range(int(quantile*10)): 
                if y_spaces[i+k] > x:
                    is_available = False
                    break
            if is_available: 
                patch = patches.Rectangle((x, i/10), width, quantile, facecolor = '#add8e6', ec='black', lw=1, alpha=0.5)
                for j in range(i, i + int(quantile * 10) + int(spacing*10)):
                    y_spaces[j] = x + width
                break

        peptide_patches.append(patch)
        last_x, last_y = patch.get_xy()

        for m_pos, m_color in mods:
            mod_rec = patches.Rectangle((last_x+m_pos,last_y), 1, quantile, facecolor= m_color, edgecolor="black")
            mod_patches.append(mod_rec)
    height = y_spaces.index(0)
    return peptide_patches, mod_patches, height/10


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


def create_and_plot_segment_plot(df, _protein="P02666", size=50, spacing=0):
    data = preprocess_data_for_peptide_segment_plot(df, _protein, size)
    rectangles_and_mods = get_rectangles_for_peptides_and_mods(data)
    rectangles_and_mods = stack_recs(rectangles_and_mods)
    peptide_patches, mod_patches, height = get_stacking_patches(rectangles_and_mods, spacing)
    plot_peptide_segments(peptide_patches, mod_patches, height)

