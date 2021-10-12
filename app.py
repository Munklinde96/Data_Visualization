import dash
import dash_html_components as html
import plotly.graph_objects as go
import dash_core_components as dcc
from dash.dependencies import Input, Output
from utils import get_modification_count_per_protein, split_data_in_samples
import pandas as pd
import plotly.express as px

path = r"UHT milk P036.csv"


def getdata(sample=0):
    df = pd.read_csv(path)
    df1, df2, df3, df4 = split_data_in_samples(df)
    if sample == 1:
        df = df1
    elif sample == 2:
        df = df2
    elif sample == 3:
        df = df3
    elif sample == 4:
        df = df4
    #apply count_no_of_modifications to each PTM column
    df['#modifications'] = df['PTM'].apply(count_no_of_modifications)
    #print non zero values in modifcations column
    df[df['#modifications'] > 0]
    return df


def count_no_of_modifications(ptm_str):
    #check if NaN value
    if pd.isnull(ptm_str):
        return 0
    return 1 + ptm_str.count(';')

app = dash.Dash()

colorscales = px.colors.named_colorscales()

def df_to_plotly(df):
    return {'z': df.values.tolist(),
            'x': df.columns.tolist(),
            'y': df.index.tolist()}


app.layout = html.Div(
    children=[
        html.Div(className='row',
                 children=[
                    html.Div(className='three columns div-user-controls',
                             children=[
                              html.H1(
        id = 'H1', children = 'Distribution of modifications over different proteins', style = {'textAlign':'center',\
                                            'marginTop':40,'marginBottom':40}),
        dcc.Dropdown( id = 'norm_dropdown',
        style= {'padding': 5,},
        options = [
            {'label':'Normalization: nothing', "value":"" },
            {'label': 'Normalization: total protein modification count', 'value':"protein_total_mod_count"},
            {'label': 'Normalizaton: protein length', 'value':"protein_length"},
            {'label': 'Normalization: protein mass', 'value':"protein_mass"},
            {'label': 'Normalization: protein intensity', 'value':"protein_intensity"},
            ],
        value = ''),
        dcc.Dropdown( id = 'count_dropdown',
        style= {'padding': 5,},
        options = [
            {'label':'Proteins with modification count OVER 10', "value":10 },
            {'label': 'Proteins with modification count OVER 20', 'value':20},
            {'label': 'Proteins with modification count OVER 40', 'value':40},
            {'label': 'Proteins with modification count OVER 80', 'value':80},
            {'label': 'Proteins with modification count OVER 160', 'value':160},
            {'label': 'Proteins with modification count OVER 320', 'value':320},
            {'label': 'Proteins with modification count OVER 640', 'value':640},
            ],
        value = 80),
        dcc.Dropdown(
        id='colorscale', 
        style= {'padding': 5,},
        options=[{"value": x, "label": x} 
                 for x in colorscales],
        value='viridis'
        ),
        dcc.Dropdown( id = 'data_dropdown',
        style= {'padding': 5,},
        options = [
            {'label': 'All samples', "value":0 },
            {'label': 'Sample 1', 'value':1},
            {'label': 'Sample 2', 'value':2},
            {'label': 'Sample 3', 'value':3},
            {'label': 'Sample 4', 'value':4},
            ],
        value = 0),
                                ]
                             ),
                    html.Div(className='nine columns div-for-charts bg-grey',
                             children=[
                                 dcc.Graph(style= {'padding': 5,}, id='mod_heatmap_plot')
                             ])
                              ])
        ]

)    
    
@app.callback(Output(component_id='mod_heatmap_plot', component_property= 'figure'),
              [
                  Input(component_id='norm_dropdown', component_property= 'value'), 
                  Input(component_id='count_dropdown', component_property= 'value'), 
                  Input("colorscale", "value"),
                  Input("data_dropdown", component_property='value')
                  ])
def graph_update(dropdown_value_norm, dropdown_value_count, scale, data):
    print(dropdown_value_norm)
    print(dropdown_value_count)
    modPd = pd.DataFrame(get_modification_count_per_protein(getdata(data), dropdown_value_count, '{0}'.format(dropdown_value_norm)))
    fig = go.Figure(data=go.Heatmap(df_to_plotly(modPd), colorscale = scale,))
    fig.update_layout(xaxis_title = 'Name of protein',
                      yaxis_title = 'Name of modification',
                      )
    return fig  


if __name__ == '__main__':
    app.run_server(host='0.0.0.0')