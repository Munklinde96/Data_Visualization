import pandas as pd
from refactored_utils import get_protein_sequence, normalize_intensities_by_protein_intensity
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

def get_overlap_overlaps_by_intensity_and_sample(df, selected_protein= "P02666", sample_column_id='Area'):
    df_new = df[df["Protein Accession"] == selected_protein]
    seq_list = list(get_protein_sequence(selected_protein))
    _len = len(seq_list)
    overlap_lists = []
    overlap_dataframes = []
    sample_columns = [col for col in df.columns if sample_column_id in col]

    for _ in sample_columns:
        overlap_lists.append([0]*_len)

    for i in range (len(df_new)):
        row = df_new.iloc[i]
        for j in range(row['Start'], row['End']):
            for k in range(len(sample_columns)):
                if not pd.isnull(row[sample_columns[k]]):
                    overlap_lists[k][j] += row[sample_columns[k]]

    for overlap_list in overlap_lists:
        df_overlaps = pd.DataFrame(list(zip(range(_len), overlap_list)), columns=['Position', 'Overlaps'])
        overlap_dataframes.append(df_overlaps)

    return overlap_lists, overlap_dataframes

def get_overlap_heapmap(num_overlpas_lists, peptide_seq_list, protein_num, fig_size=(30,10), color_scale='YlOrRd'):
    # make string, with 1,2,3,... 
    title_samples_ennumeration = [str(i) for i in range(1, len(num_overlpas_lists)+1)]
    title_samples_ennumeration = ','.join(title_samples_ennumeration)

    plt.figure(figsize=fig_size)
    plt.title(f"Frequency of Overlaps for Protein {protein_num} - sample "+title_samples_ennumeration)
    ax = sns.heatmap(num_overlpas_lists, cmap=color_scale)
    plt.xticks(np.arange(len(peptide_seq_list)), peptide_seq_list, rotation = 0)
    ylabels = [ 'Sample '+ str(i) for i in range(1, len(num_overlpas_lists)+1)]
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
    ylabels = [ 'Sample '+ str(i) for i in range(1, len(num_overlpas_lists)+1)]    
    ax.set_yticklabels(ylabels)
    plt.show()

def create_and_plot_overlap_plot(df, _protein="P02666", sample_column_id='Area'):
    df = normalize_intensities_by_protein_intensity(df, sample_column_id=sample_column_id)
    lists, dataframes = get_overlap_overlaps_by_intensity_and_sample(df, _protein, sample_column_id=sample_column_id)

    seq_list = list(get_protein_sequence(_protein))
    fig, ax = plt.subplots(figsize = (30, 10))
    ax.set_xticks(range(len(lists[0])))
    ax.set_xticklabels(seq_list)
    
    labels = [ 'Sample '+ str(i) for i in range(1, len(dataframes)+1)]
    color_list = sns.color_palette("hls", len(labels))

    bottom = dataframes[0]['Overlaps']*0
    for i in range(len(dataframes)):
        ax.bar(dataframes[i]['Position'], dataframes[i]['Overlaps'], bottom = bottom, label=labels[i], color = color_list[i])
        bottom = bottom + dataframes[i]['Overlaps']        

    ax.set_title('Overlaps by Position and Sample')
    ax.legend()

def create_and_plot_overlap_heat_map_and_gradient_heatmap(df, _protein="P02666", sample_column_id='Area'):
    seq_list = list(get_protein_sequence(_protein))
    df = normalize_intensities_by_protein_intensity(df, sample_column_id=sample_column_id)
    lists, dataframes = get_overlap_overlaps_by_intensity_and_sample(df, _protein, sample_column_id = sample_column_id)
    get_overlap_heapmap(lists, seq_list, _protein, fig_size=(30,4))
    get_overlap_gradient_heapmap(lists, seq_list, _protein, fig_size=(30,4))