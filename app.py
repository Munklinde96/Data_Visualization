import dash
import dash_html_components as html
import plotly.graph_objects as go
import dash_core_components as dcc
import plotly.express as px
from dash.dependencies import Input, Output
from utils import get_modification_count_per_protein
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

path = r"UHT milk P036.csv"
 
df = pd.read_csv(path)

def count_no_of_modifications(ptm_str):
    #check if NaN value
    if pd.isnull(ptm_str):
        return 0
    return 1 + ptm_str.count(';')

#apply count_no_of_modifications to each PTM column
df['#modifications'] = df['PTM'].apply(count_no_of_modifications)
#print non zero values in modifcations column
df[df['#modifications'] > 0]

app = dash.Dash()

def df_to_plotly(df):
    return {'z': df.values.tolist(),
            'x': df.columns.tolist(),
            'y': df.index.tolist()}

app.layout = html.Div(id = 'parent', children = [
    html.H1(id = 'H1', children = 'Styling using html components', style = {'textAlign':'center',\
                                            'marginTop':40,'marginBottom':40}),
        dcc.Dropdown( id = 'norm_dropdown',
        options = [
            {'label':'No normalization', 'value':'' },
            {'label': 'Total protein modification count', 'value':'protein_total_mod_count'},
            {'label': 'Protein length', 'value':'protein_length'},
            {'label': 'Protein mass', 'value':'protein_mass'},
            ],
        value = 'protein_total_mod_count'),
        dcc.Graph(id = 'mod_heatmap_plot')
    ])
    
    
@app.callback(Output(component_id='mod_heatmap_plot', component_property= 'figure'),
              [Input(component_id='norm_dropdown', component_property= 'value')])
def graph_update(dropdown_value):
    print(dropdown_value)
    modPd = pd.DataFrame(get_modification_count_per_protein(df, 200, '{}'.format(dropdown_value)))
    fig = go.Figure(data=go.Heatmap(df_to_plotly(modPd)))
    fig.update_layout(title = 'Distribution of modifications over proteins',
                      xaxis_title = 'Name of protein',
                      yaxis_title = 'Name of modification'
                      )
    return fig  



if __name__ == '__main__': 
    app.run_server()