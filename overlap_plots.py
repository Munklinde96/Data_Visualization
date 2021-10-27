import pandas as pd
from refactored_utils import get_protein_sequence, normalize_intensities_by_protein_intensity
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

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

def create_and_plot_overlap_plot(df, _protein="P02666"):
    df = normalize_intensities_by_protein_intensity(df)
    lists, dataframes = get_overlap_overlaps_by_intensity_and_sample(df, _protein)

    seq_list = list(get_protein_sequence(_protein))
    fig, ax = plt.subplots(figsize = (30, 10))
    ax.set_xticks(range(len(lists[0])))
    ax.set_xticklabels(seq_list)

    ax.bar(dataframes[0]['Position'], dataframes[0]['Overlaps'], label='Sample 1', color = '#A6CEE3')
    ax.bar(dataframes[1]['Position'], dataframes[1]['Overlaps'], bottom=dataframes[0]['Overlaps'], label='Sample 2', color='#1F78B4')
    ax.bar(dataframes[2]['Position'], dataframes[2]['Overlaps'], bottom=dataframes[1]['Overlaps'] + dataframes[0]['Overlaps'], label="Sample 3", color='#B2DF8A')
    ax.bar(dataframes[3]['Position'], dataframes[3]['Overlaps'], bottom=dataframes[2]['Overlaps'] + dataframes[1]['Overlaps'] + dataframes[0]['Overlaps'], label="Sample 4", color='#33A02C')

    ax.set_title('Overlaps by Position and Sample')
    ax.legend()

def create_and_plot_overlap_heat_map_and_gradient_heatmap(df, _protein="P02666"):
    seq_list = list(get_protein_sequence(_protein))
    df = normalize_intensities_by_protein_intensity(df)
    lists, dataframes = get_overlap_overlaps_by_intensity_and_sample(df, _protein)
    get_overlap_heapmap(lists, seq_list, _protein, fig_size=(30,4))
    get_overlap_gradient_heapmap(lists, seq_list, _protein, fig_size=(30,4))