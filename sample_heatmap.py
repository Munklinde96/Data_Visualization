import pandas as pd
import seaborn as sns
from refactored_utils import normalize_intensities_by_protein_intensity

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


def create_and_plot_sample_heatmap(df, _protein='P02666'):
    df = normalize_intensities_by_protein_intensity(df)
    df = df[df['Protein Accession'] == _protein]
    df1, df2, df3, df4 = split_data_in_samples(df)
    combined = combine_and_aggregate_intensity(df1, df2, df3, df4).sort_values(by=['PTM'])

    data = pd.pivot_table(data = combined, index = 'PTM', values = 'Intensity', columns='Sample')
    sns.heatmap(data, cmap='viridis')



        
