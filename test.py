
from refactored_utils import get_and_prepare_data
from segment_plot import create_and_plot_segment_plot
from sample_heatmap import create_and_plot_sample_heatmap
from overlap_plots import create_and_plot_overlap_plot, create_and_plot_overlap_heat_map_and_gradient_heatmap

data_path = r'Data_Visualization/protein-peptides.csv' #FOR DEBUG
data_path = r'protein-peptides.csv'
df = get_and_prepare_data(data_path)
# get number of different proteins
protein_ids = df['Protein Accession'].unique()
print(f'Number of different proteins: {len(protein_ids)}')

#create_and_plot_segment_plot(df, _protein = 'Q2UVX4|CO3_BOVIN' ,spacing = 0.2, color_scale ='Blues', is_log_scaled=True)
create_and_plot_segment_plot(df, _protein = 'P02663|CASA2_BOVIN' ,spacing = 0.2, color_scale ='Blues', is_log_scaled=True)