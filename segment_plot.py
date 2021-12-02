import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.patches as patches
from matplotlib import pyplot as plt
from matplotlib import colors as cls
from refactored_utils import get_position_of_mass_shift, get_color_palette_for_modifications, get_protein_sequence, get_color_legend_mods, get_selected_sample_columns
from utils import colors_from_values
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sklearn.preprocessing import normalize as norm
import time

######################################
############ SEGMENT PLOT ############
######################################
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = cls.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def preprocess_data_for_peptide_segment_plot(df, _protein="P02666", sample_column_id = 'Area', selected_sample_indices = [], start_end_indicies =None):
    df = df.copy()
    df = df[df["Protein Accession"] == _protein]
    if start_end_indicies is not None:
        df = df[(df.Start == start_end_indicies[0]) & (df.End == start_end_indicies[1])]

    df["Position of Mass Shift"] = df["Peptide"].apply(get_position_of_mass_shift)
    # get list of modification for each PTM
    df["Modification_types"] = df["PTM"].apply(lambda x: x if pd.isnull(x) or x == 'Unmodified' else [s.strip() for s in x.split(";")])    
    selected_sample_columns = [col for col in df.columns if sample_column_id in col]

    if len(selected_sample_indices) > 0:
        selected_sample_columns = get_selected_sample_columns(selected_sample_columns, selected_sample_indices)
    
    ##SHOULD be deleted
    selected_columns = ["Start", "End", "Protein Accession", 'Position of Mass Shift', 'Modification_types'] + selected_sample_columns
    df = df[selected_columns]
    
    unique_mod_types = get_unique_modification_types(df)

    #Aggregate sample intensity, and normalize it
    df["Agg Intensity"] = df[selected_sample_columns].sum(axis=1)
    df['Agg Intensity'] = df['Agg Intensity'] / df['Agg Intensity'].sum()

    #remove all rows where Agg Intensity is 0
    df = df[df['Agg Intensity'] > 0]

    df["tuple"] = df[["Start", "End", 'Position of Mass Shift', 'Modification_types', 'Agg Intensity']].apply(tuple, axis=1)

    start_end_ms_modtype_list = df['tuple'].tolist()
    return start_end_ms_modtype_list, unique_mod_types

def get_unique_modification_types(df):
    ls = df.Modification_types.tolist()
    mod_types = []
    for l in ls :
        if type(l) is list:
            mod_types.extend(l)
    return list(set(mod_types))

def get_rectangles_for_peptides_and_mods(data, modtypes_color_map):
    """data is a list of tuples on the form (low,hi, [modifications],[modtypes], [agg_intensity])"""
    data = sorted(data, key=lambda x: (x[0], -x[1]), reverse=False)
    
    rectangles_and_mods = []
    res_intensities = []
    for i in range (len(data)):
        modifications = []
        low, hi, mod_positions, mod_types, agg_intensity = data[i]
        width = hi-low + 1
        low = low - 1
        rec = (low, width, agg_intensity)
        res_intensities.append(agg_intensity)
        if len(mod_positions) > 0:
            for _ms_pos, mod_type in zip(mod_positions, mod_types): #add mass shift color on rectangles if present
                ms_color = modtypes_color_map[mod_type]
                modifications.append(((_ms_pos), ms_color, agg_intensity, mod_type))
        rectangles_and_mods.append((rec, modifications))
        
    return res_intensities, rectangles_and_mods

def get_intervals(agg_intensities: list, no_intervals = 5):
    intervals = pd.qcut(agg_intensities, q=no_intervals)
    intervals = pd.unique(intervals)
    intervals = intervals[::-1] #reverse intervals 
    interval_len = len(intervals)
    intervals = {intervals[i]: (i+1)/interval_len for i in range(interval_len)}
    return intervals

def make_continous_color_scale(intervals: dict):
    """
    intervals is a dict with keys as quantiles and values as colors
    """
    #make continuous grey color scale

def normalize(res_intensities, is_log_scaled = False,  min_val=0, intensity_value= None):
    values = np.array(res_intensities)
    print("values",values)
    if is_log_scaled:
        values = np.log(values)
    normalized = (values - min(values)) / (max(values) - min(values)) # normalize
    if intensity_value is not None:
        normalized = normalized * intensity_value       
    normalized = normalized * (1 - min_val) + min_val
    print("normalized",normalized)
    return normalized

def colors_(values: list , color_scale = 'Blues', is_log_scaled = True, is_normalized = True, intensity_value = None):
    normalized = np.asarray(values)
    cmap = plt.cm.get_cmap(color_scale)
    if len(values) == 1:
        if intensity_value is not None:
            return [cls.rgb2hex(cmap(intensity_value))]
        return [cls.rgb2hex(cmap(1.0))]
    elif is_normalized:
        normalized = normalize(values, is_log_scaled,  min_val=0.2, intensity_value=intensity_value)
    #get max value from cmap
    color_list = [cls.rgb2hex(cmap(i)) for i in normalized]
    return color_list


def map_to_colors(res_intensities, res_rectangles_and_mods, color_scale = 'Blues', is_log_scaled = True, intensity_value = None):
    rects = [r[0] for r in res_rectangles_and_mods]
    color_values = colors_(res_intensities, color_scale=color_scale, is_log_scaled=is_log_scaled, intensity_value=intensity_value)
    rects_and_colors = list(zip(rects, color_values))
    mods = [r[1] for r in res_rectangles_and_mods]
    rects_and_colors_mods = list(zip(rects_and_colors, mods))
    return rects_and_colors_mods

#documentation:
"""
res_intensities: list of intensities
res_rectangles_and_mods: list of tuples of the form (((low, width), agg_intensity), [(mod_pos, mod_color, agg_intensity), ...])
"""
def map_to_norm_intensities(res_intensities, res_rectangles_and_mods, is_log_scaled = True):
    res_intensities = normalize(res_intensities, is_log_scaled)
    rects = [r[0] for r in res_rectangles_and_mods]
    rects_and_intensities = list(zip(rects, res_intensities))
    mods = [r[1] for r in res_rectangles_and_mods]
    rects_and_intensities = list(zip(rects_and_intensities, mods))
    return rects_and_intensities

def map_to_attribute(colors, color_scale, is_log_scaled, res_intensities, rectangles_and_mods, intensity_value = None):
    if colors:
        rects_and_attribute = map_to_colors(res_intensities, rectangles_and_mods, color_scale = color_scale, is_log_scaled=is_log_scaled, intensity_value=intensity_value)
    else:
        rects_and_attribute = map_to_norm_intensities(res_intensities, rectangles_and_mods)
    return rects_and_attribute

def stack_recs(rectangles_and_mods: list, colors = True, color_scale = 'Blues', is_log_scaled = True):
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
        low, width, intensity = rect
        if low == last_x and width == last_rec_width:
            last_rect, last_mods = res_rectangles_and_mods[-1]
            last_rect_low, last_rect_width, last_intensity = last_rect
            new_intensity = last_intensity + intensity
            new_mods = last_mods + mods
            res_rectangles_and_mods[-1] = ((last_rect_low, last_rect_width, new_intensity), new_mods)
            res_intensities[-1] = new_intensity
        else: 
            res_rectangles_and_mods.append(((low, width, intensity), mods))
            res_intensities.append(intensity)
            last_x = low
            last_rec_width = width
    
    return res_intensities, res_rectangles_and_mods


def get_mods_height_distribution(mods):
    #sort by first pos then color
    mods = sorted(mods, key=lambda x: (x[0], x[1]))
    mods_list_list = []
    mods_list = [mods[0]]
    mods_list_list.append(mods_list)
    pos = mods[0][0]
    color = mods[0][1]
    for i in range(1, len(mods)):
        mod = mods[i]
        if mod[0] == pos and mod[1] == color:
            mods_list[-1] = (mods_list[-1][0], mods_list[-1][1], mods_list[-1][2] + mod[2])
        elif mod[0] == pos and mod[1] != color:
            mods_list.append(mod)
            color = mod[1]
        else:
            mods_list = []
            mods_list.append(mod)
            mods_list_list.append(mods_list)
            pos = mod[0]
            color = mod[1]
    return mods_list_list

# standard height only used if colors is True
def get_patch_attributes(rectangles, spacing=0.2, colors = True, standard_height= 0.5):
    y_spaces = [0]*int(10*len(rectangles)*standard_height + len(rectangles)*10*spacing + 1) #10 indicates a division into heights of 0.1
    peptide_patches = []
    mod_patches = []
    for rect, mods in rectangles:
        rect_info, attribute = rect #attribute is either color or height depepnding on the colors parameter
        x, width, intensity = rect_info
        for i in range(len(y_spaces)):
            if colors: # Encode intensities as COLOR
                is_available = True
                for k in range(int(standard_height*10)): 
                    if y_spaces[i+k] > x:
                        is_available = False
                        break
                if is_available:
                    patch_attributes = (x, i/10, width, standard_height, attribute, intensity) # add in d3.js ec='black', lw=1, alpha=1
                    for j in range(i, i + int(standard_height * 10) + int(spacing*10)):
                        y_spaces[j] = x + width
                    break
            else: # Encode intensities as HEIGHT
                is_available = True
                for k in range(int((attribute + spacing)*10)): 
                    if y_spaces[i+k] > x:
                        is_available = False
                        break
                if is_available: 
                    patch_attributes = (x, i/10, width, attribute, '#add8e6', intensity) # add in d3.js ec='black', lw=1, alpha=1
                    for j in range(i, i + int(attribute * 10) + int(spacing*10)):
                        y_spaces[j] = x + width
                    break
        
        peptide_patches.append(patch_attributes)
        last_x = patch_attributes[0]
        last_y = patch_attributes[1]
        if mods != []:
            mods_list_list = get_mods_height_distribution(mods)
        else: 
            mods_list_list = []
        for mods_list in mods_list_list:
            intensities = [m[2] for m in mods_list]
            normalized_intensities_mods = np.asarray(intensities) / intensity
            y_pos = last_y
            for i in range(len(mods_list)):
                m_pos = None
                m_color = None 
                mod_type = ""
                if len(mods_list[i]) == 4:
                    m_pos, m_color, _, mod_type = mods_list[i]
                else:
                    m_pos, m_color, _ = mods_list[i]
                if colors:
                    mod_height = standard_height * normalized_intensities_mods[i]
                    mod_attributes = (last_x+m_pos,y_pos, 1, mod_height, m_color, intensities[i], mod_type) # add in d3.js ec='black', lw=1, alpha=1
                else:
                    mod_height =  attribute * normalized_intensities_mods[i]
                    mod_attributes = (last_x+m_pos, y_pos, 1, mod_height, m_color, intensities[i], mod_type) # add in d3.js ec='black', lw=1, alpha=1
                y_pos = y_pos + mod_height
                mod_patches.append(mod_attributes)
    height = y_spaces.index(0)
    return peptide_patches, mod_patches, height/10

def get_patches_from_patch_attributes(peptide_patches, mod_patches):
    peptide_patches_list = []
    mod_patches_list = []
    for i in range(len(peptide_patches)):
        x, y, width, height, color, intensity = peptide_patches[i]
        patch = patches.Rectangle((x, y), width, height, facecolor = color, ec='black', lw=1, alpha=1)
        peptide_patches_list.append(patch)
    for i in range(len(mod_patches)):
        x, y, width, height, color, intensity = mod_patches[i]
        patch = patches.Rectangle((x, y), width, height, facecolor = color, ec='black', lw=1, alpha=1)
        mod_patches_list.append(patch)
    return peptide_patches_list, mod_patches_list

def plot_peptide_segments(peptide_patches, mod_patches, height, modification_color_map, _protein="P02666", colors = True, color_scale = 'Blues', is_log_scaled = True):
    fig = plt.figure(figsize=(30,25))
    ax = fig.add_subplot(111)
    ax.set_ylim((0,height))

    for patch in peptide_patches:
        ax.add_patch(patch)

    for patch in mod_patches:
        ax.add_patch(patch)

    seqq = get_protein_sequence(_protein)
    ax.set_xlabel("Full Protein Sequence")
    ax.set_xticks(range(0,len(seqq)))
    # create list of chars from string 
    protein_seq_list = list(seqq)
    ax.set_xticklabels(protein_seq_list)
    ax.get_yaxis().set_visible(False)
    modification_types_to_color_map = modification_color_map
    
    handles = get_color_legend_mods(modification_types_to_color_map)
    cmap = plt.cm.get_cmap(color_scale)
    new_cmap = truncate_colormap(cmap, 0.2, 1)
    if colors:
        if is_log_scaled:
            sm = plt.cm.ScalarMappable(cmap=new_cmap, norm=cls.LogNorm(vmin=0.00001, vmax=1))
            sm._A = []
            label_ = 'Normalized and logscaled intensity'
        else: 
            sm = plt.cm.ScalarMappable(cmap=new_cmap, norm=cls.Normalize(vmin=0, vmax=1))
            sm._A = []
            label_ = 'Normalized intensity'
        cbaxes = inset_axes(ax, width="2%", height="15%", loc=2)
        cbar = plt.colorbar(sm, cax=cbaxes, orientation='vertical')
        cbar.set_label(label_)
        leg1 = ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.1, 1), ncol=1)
    else: 
        leg1 = ax.legend(handles=handles, loc='upper left')
        
    plt.show()
    return ax


def create_histogram_over_mod_positions(mod_patches, modtypes_color_map, peptide_seq):
    mod_positions = []
    # make map from color to modtype
    color_to_modtype = {}
    mod_types = list(modtypes_color_map.keys())
    colors = list(modtypes_color_map.values())
    for i in range(len(colors)): # make map from colors to modtypes
        color_to_modtype[colors[i]] = mod_types[i]

    mod_positions_df = pd.DataFrame(mod_patches)
    columns = ['x', 'y', 'width', 'height', 'color', 'intensity', 'mod_type']
    mod_positions_df.columns = columns
    # apply color_to_modtype to color column
    mod_positions_df['mod_type'] = mod_positions_df['color'].apply(lambda x: color_to_modtype[x])
    # drop all other than x, color and modtype
    mod_positions_df = mod_positions_df[['x', 'mod_type', 'color']]
    # set x as index
    # plt.figure(figsize=(20,10))
    # sns.histplot(data=mod_positions_df, x = 'x',  hue = 'mod_type', multiple='stack', bins=50)
    # plt.show()
    return mod_positions_df

"""
   Poetry. slam poetry. Slam Poetry is a collection of poems written by the poet, slam poet.
    The poems are collected in the slam poetry collection.
"""
def create_and_plot_segment_plot(df, _protein="P02666", spacing=0.2, colors = True, color_scale='Blues', is_log_scaled = True, standard_height = 2, start_end_indicies = None, is_stacked = True, sample_column_id = 'Area' , selected_sample_indices = []):
    peptide_patches, mod_patches, height, seqq, modification_color_map, min_ind, max_ind = create_data_for_segment_plot(df, _protein, spacing=spacing, colors=colors, color_scale=color_scale, is_log_scaled=is_log_scaled, standard_height=standard_height, start_end_indicies=start_end_indicies, is_stacked=is_stacked, sample_column_id=sample_column_id, selected_sample_indices=selected_sample_indices)
    peptide_patches, mod_patches = get_patches_from_patch_attributes(peptide_patches, mod_patches)
    plot_peptide_segments(peptide_patches, mod_patches, height, modification_color_map, colors = colors, color_scale=color_scale, is_log_scaled=is_log_scaled, _protein = _protein)


def create_data_for_segment_plot(df, _protein="P02666", spacing=0.2, colors = True, color_scale='Blues', is_log_scaled = True, standard_height = 2, start_end_indicies = None, is_stacked = True, sample_column_id = 'Area' , selected_sample_indices = [], intensity_value = None):
    # alterntive protein sequence P80457|XDH_BOVIN
    data, unique_mod_types= preprocess_data_for_peptide_segment_plot(df, _protein, sample_column_id=sample_column_id, selected_sample_indices=selected_sample_indices, start_end_indicies=start_end_indicies)
    modtypes_color_map = get_color_palette_for_modifications(unique_mod_types)
    res_intensities, rectangles_and_mods = get_rectangles_for_peptides_and_mods(data, modtypes_color_map)
    
    if(is_stacked):
        res_intensities, rectangles_and_mods = stack_recs(rectangles_and_mods)
    rects_and_attribute = map_to_attribute(colors, color_scale, is_log_scaled, res_intensities, rectangles_and_mods, intensity_value=intensity_value)
    print("here")
    print(rectangles_and_mods)
    print("here")
    peptide_patches, mod_patches, height = get_patch_attributes(rects_and_attribute, spacing = spacing, standard_height=standard_height)
    seqq = get_protein_sequence(_protein)

    histogram_df = create_histogram_over_mod_positions(mod_patches, modtypes_color_map, seqq)

    min_ind = np.argmin(res_intensities)
    max_ind = np.argmax(res_intensities)
    min_intensity = res_intensities[min_ind]
    max_intensity = res_intensities[max_ind]
    min_color = peptide_patches[min_ind][4]
    max_color = peptide_patches[max_ind][4]

    return peptide_patches, mod_patches, height, seqq, modtypes_color_map, (min_intensity, min_color), (max_intensity, max_color), histogram_df

