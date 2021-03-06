import pandas as pd
import seaborn as sns
from refactored_utils import normalize_intensities_by_protein_intensity, normalize_intensities

def split_data_in_samples(df, sample_column_id='Area'):
    sample_columns = [col for col in df.columns if sample_column_id in col]
    data_lists = []
    #Make list containing a list for each sample
    for _ in range(len(sample_columns)):
        data_lists.append([])

    for i, row in df.iterrows():
        for j, col in enumerate(sample_columns):
            if row[col] is not None and row[col] > 0:
                data_lists[j].append(row)

    dataframes = [pd.DataFrame(data_list) for data_list in data_lists]
    return dataframes

def combine_and_aggregate_intensity(dataframes: list, sample_column_id='Area'):
    sample_columns = [col for col in dataframes[0].columns if sample_column_id in col]

    df_collection = []
    for i, df in enumerate(dataframes):
        if(df.empty):
            df = pd.DataFrame(columns=[['PTM', 'Intensity', 'Sample']])
            df['Sample'] = i+1
            df['PTM'] = 'Unmodified'
            df['Intensity'] = 0.0000001
            df_collection.append(df)
            continue
        df = df.copy()
        # df['#modifications'] = df['#modifications'].fillna(0)
        # df['#modifications'] = df['#modifications'].replace(0,1)
        df['Intensity'] = df[sample_columns[i]]
        df['Intensity'] = df['Intensity'] / df['Intensity'].sum()
        df['PTM'] = df['PTM'].str.split(';').str[0]
        df_new = df[['PTM', 'Intensity']]
        df_new = df_new.groupby(['PTM']).sum()
        df_new = df_new.reset_index()
        df_new['Sample'] = i+1
        df_collection.append(df_new[['PTM', 'Intensity', 'Sample']])

    combined = pd.concat(df_collection, axis=0) # concat each df in df_collection on axis 0
    return combined


def create_and_plot_sample_heatmap(df, _protein='P02666'):
    data = get_sample_heatmap_data(df, _protein)
    sns.heatmap(data, cmap='viridis')

def get_sample_heatmap_data(df, _protein='P02666', sample_column_id='Area'):
    df = df[df['Protein Accession'] == _protein]
    # df = normalize_intensities(df, sample_column_id=sample_column_id)
    df_list = split_data_in_samples(df, sample_column_id=sample_column_id)
    combined = combine_and_aggregate_intensity(df_list, sample_column_id=sample_column_id).sort_values(by=['PTM'])

    data = pd.pivot_table(data = combined, index = 'PTM', values = 'Intensity', columns='Sample')
    return data
