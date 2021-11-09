import pandas as pd
import numpy as np
import matplotlib.patches as patches
from matplotlib import pyplot as plt
from matplotlib import colors as cls
from refactored_utils import get_position_of_mass_shift, get_color_palette_for_modifications, get_protein_sequence, get_color_legend_mods 
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

def preprocess_data_for_peptide_segment_plot(df, _protein="P02666", size=50, sample_column_id = 'Area', selected_samples = ['Area Sample 1', 'Area Sample 2']):
    df = df.copy()
    df["Position of Mass Shift"] = df["Peptide"].apply(get_position_of_mass_shift)
    # get list of modification for each PTM
    df["Modification_types"] = df["PTM"].apply(lambda x: x if pd.isnull(x) else [s.strip() for s in x.split(";")])
    if len(selected_samples) == 0:
        selected_samples = [col for col in df.columns if sample_column_id in col]
    
    selected_columns = ["Start", "End", "Protein Accession", 'Position of Mass Shift', 'Modification_types'] + selected_samples
    df = df[selected_columns]
    df = df[df["Protein Accession"] == _protein] #only look at values for protein : P02666
   
    for col in selected_samples:  #make Area samples Nan values 0
        df[col] = df[col].fillna(0)
    df = df[df[selected_samples].sum(axis=1) > 0]   # drop rows where all sample_columns are 0

    #Aggregate sample intensity, and normalize it
    df["Agg Intensity"] = df[selected_samples].sum(axis=1)
    df['Agg Intensity'] = df['Agg Intensity'] / df['Agg Intensity'].sum()

    df["tuple"] = df[["Start", "End", 'Position of Mass Shift', 'Modification_types', 'Agg Intensity']].apply(tuple, axis=1)
    new = df.head(size)

    start_end_ms_modtype_list = new['tuple'].tolist()
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
                modifications.append((_ms_pos, ms_color, agg_intensity))
        rectangles_and_mods.append((rec, modifications))
        
    return rectangles_and_mods

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

def normalize(res_intensities, is_log_scaled = False,  min_val=0):
    values = np.array(res_intensities)
    if is_log_scaled:
        values = np.log(values)
    values = (values - min(values)) / (max(values) - min(values)) # normalize
    normalized = values * (1 - min_val) + min_val
    return normalized

def colors_(values: list , color_scale = 'Greys', is_log_scaled = True):
    normalized = normalize(values, is_log_scaled,  min_val=0.2)
    cmap = plt.cm.get_cmap(color_scale)
    color_list = [cls.rgb2hex(cmap(i)) for i in normalized]
    return color_list
    
def map_to_colors(res_intensities, res_rectangles_and_mods, color_scale = 'Greys', is_log_scaled = True):
    rects = [r[0] for r in res_rectangles_and_mods]
    color_values = colors_(res_intensities, color_scale=color_scale, is_log_scaled=is_log_scaled)
    rects_and_colors = list(zip(rects, color_values))
    mods = [r[1] for r in res_rectangles_and_mods]
    rects_and_colors_mods = list(zip(rects_and_colors, mods))
    return rects_and_colors_mods

def map_to_norm_intensities(res_intensities, res_rectangles_and_mods, is_log_scaled = True):
    res_intensities = normalize(res_intensities, is_log_scaled)
    rects = [r[0] for r in res_rectangles_and_mods]
    rects_and_intensities = list(zip(rects, res_intensities))
    mods = [r[1] for r in res_rectangles_and_mods]
    rects_and_intensities = list(zip(rects_and_intensities, mods))
    return rects_and_intensities

def stack_recs(rectangles_and_mods: list, colors = True, color_scale = 'Greys', is_log_scaled = True):
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

    if colors:
        rects_and_attribute = map_to_colors(res_intensities, res_rectangles_and_mods, color_scale = color_scale, is_log_scaled=is_log_scaled)
    else:
        rects_and_attribute = map_to_norm_intensities(res_intensities, res_rectangles_and_mods)

    return rects_and_attribute


# standard height only used if colors is True
def get_stacking_patches(rectangles, spacing=0.2, colors = True, standard_height= 0.5):
    y_spaces = [0]*10*len(rectangles) #10 indicates a division into heights of 0.1
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
                    patch = patches.Rectangle((x, i/10), width, standard_height, facecolor = attribute, ec='black', lw=1, alpha=1)
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
                    patch = patches.Rectangle((x, i/10), width, attribute, facecolor = '#add8e6', ec='black', lw=1, alpha=0.5)
                    for j in range(i, i + int(attribute * 10) + int(spacing*10)):
                        y_spaces[j] = x + width
                    break
        
        peptide_patches.append(patch)
        last_x, last_y = patch.get_xy()
        if mods != []:
            mods_list_list = get_mods_height_distribution(mods)
        else: 
            mods_list_list = []
        for mods_list in mods_list_list:
            intensities = [m[2] for m in mods_list]
            #normalize intensities with numpy
            normalized_intensities_mods = np.asarray(intensities) / intensity
            y_pos = last_y
            for i in range(len(mods_list)):
                m_pos, m_color, _ = mods_list[i]
                if colors:
                    mod_height = standard_height * normalized_intensities_mods[i]
                    mod_rec = patches.Rectangle((last_x+m_pos,y_pos), 1, mod_height, facecolor= m_color, edgecolor="black")
                else:
                    mod_height =  attribute * normalized_intensities_mods[i]
                    mod_rec = patches.Rectangle((last_x+m_pos,y_pos), 1, mod_height, facecolor= m_color, edgecolor="black")
                y_pos = y_pos + mod_height
                mod_patches.append(mod_rec)
    height = y_spaces.index(0)
    return peptide_patches, mod_patches, height/10

# standard height only used if colors is True
def get_stacking_patch_attributes(rectangles, spacing=0.2, colors = True, standard_height= 0.5):
    y_spaces = [0]*10*len(rectangles) #10 indicates a division into heights of 0.1
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
                    patch_attributes = (x, i/10, width, standard_height, attribute) # add in d3.js ec='black', lw=1, alpha=1
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
                    patch_attributes = (x, i/10, width, attribute, '#add8e6') # add in d3.js ec='black', lw=1, alpha=1
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
            #normalize intensities with numpy
            normalized_intensities_mods = np.asarray(intensities) / intensity
            y_pos = last_y
            for i in range(len(mods_list)):
                m_pos, m_color, _ = mods_list[i]
                if colors:
                    mod_height = standard_height * normalized_intensities_mods[i]
                    mod_attributes = (last_x+m_pos,y_pos, 1, mod_height, m_color) # add in d3.js ec='black', lw=1, alpha=1
                else:
                    mod_height =  attribute * normalized_intensities_mods[i]
                    mod_attributes = (last_x+m_pos, y_pos, 1, mod_height, m_color) # add in d3.js ec='black', lw=1, alpha=1
                y_pos = y_pos + mod_height
                mod_patches.append(mod_attributes)
    height = y_spaces.index(0)
    return peptide_patches, mod_patches, height/10

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



def get_mods_heights(mods: list, total_height: float, use_intensity = False):
    """
    mods is a list of tuples of the form (mod_pos, mod_color, intensity)
    total_height is the total height of the peptide
    use_intensity is a boolean that indicates whether the intensity should be used to determine the height of the mod
    """
    height_list = []
    if use_intensity:
        pass
    else:
    
        mods_set = set(mods)
        # get size of mods_set
        num_unique = len(mods_set)
        # get count for each color in mods
        color_counts = [mods.count(tup[1]) for tup in mods]


def plot_peptide_segments(peptide_patches, mod_patches, height, _protein="P02666", colors = True, color_scale = 'Greys', is_log_scaled = True):
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
    # create 2 legends and place them next to each other
    
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


def create_and_plot_segment_plot(df, _protein="P02666", size=50, spacing=0, colors = True, color_scale='Greys', is_log_scaled = True):
    start = time.time()
    data = preprocess_data_for_peptide_segment_plot(df, _protein, size)
    end = time.time()
    print("preprocessing time: ", end - start)
    start = time.time()
    rectangles_and_mods = get_rectangles_for_peptides_and_mods(data)
    end = time.time()
    print("get_rectangles_for_peptides_and_mods time: ", end - start)
    start = time.time()
    rectangles_and_mods = stack_recs(rectangles_and_mods, colors=colors, color_scale=color_scale, is_log_scaled=is_log_scaled)
    end = time.time()
    print("stack_recs time: ", end - start)
    start = time.time()
    peptide_patches, mod_patches, height = get_stacking_patches(rectangles_and_mods, spacing = spacing, colors = colors)
    end = time.time()
    print("get_stacking_patches time: ", end - start)
    start = time.time()
    plot_peptide_segments(peptide_patches, mod_patches, height, colors = colors, color_scale=color_scale, is_log_scaled=is_log_scaled)
    end = time.time()
    print("plot_peptide_segments time: ", end - start)
