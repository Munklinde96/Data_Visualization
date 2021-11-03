from flask import Flask, flash, redirect, render_template, request, session, abort,send_from_directory,send_file,jsonify
import pandas as pd
from utils import combine_and_aggregate_intensity, get_modification_count_per_protein, sanitize_data, preprocess_data_for_peptide_segment_plot, plot_peptide_segments, normalize_intensities_by_protein_intensity, split_data_in_samples, get_protein_sequence, get_overlap_overlaps_by_intensity_and_sample, get_overlap_pixel_plot, get_gradient_plot
from flask_cors import CORS, cross_origin
import json
import mpld3
from mpld3 import plugins

#0. Helpers
def df_to_map(df):
    return {'z': df.values.tolist(),
            'x': df.columns.tolist(),
            'y': df.index.tolist()}

def count_no_of_modifications(ptm_str):
    #check if NaN value
    if pd.isnull(ptm_str):
        return 0
    return 1 + ptm_str.count(';')


#1. Declare application
application= Flask(__name__)
cors = CORS(application)
application.config['CORS_HEADERS'] = 'Content-Type'

#2. Declare data stores
class DataStore():
    Normalization=None
    MinCount= None
    ProteinModData=None
    SampleModData=None
data=DataStore()

def buildProteinModData(df, count, normalization):
    # Count modifications
    df['#modifications'] = df['PTM'].apply(count_no_of_modifications)
    df[df['#modifications'] > 0]

    # Get data for first viz
    modData = pd.DataFrame(get_modification_count_per_protein(df, count,  normalization))
    
    # Convert to json
    return modData.to_json()

def buildSampleModData(df):
    # split data into samples
    df1, df2, df3, df4 = split_data_in_samples(df)
    combined = combine_and_aggregate_intensity(df1, df2, df3, df4).sort_values(by=['PTM'])
    data = pd.pivot_table(data = combined, index = 'PTM', values = 'Intensity', columns='Sample')
    return data.to_json()

@application.route("/main",methods=["GET","POST"])
#3. Define main code
@application.route("/",methods=["GET","POST"])
def homepage():
    # Get data and build dataframe
    path = r"UHT milk P036.csv"
    df = pd.read_csv(path)
    #df1, df2, df3, df4 = split_data_in_samples(df)
    
    # Create initial configuration
    data.Normalization = request.form.get('Normalization_field','protein_intensity')
    Normalization=data.Normalization
    data.MinCount = request.form.get('Min_count_field', 100)
    MinCount=data.MinCount

    # Get JSON data for protein mod viz
    proteinModJson = buildProteinModData(df, MinCount, Normalization)

    # Save to datastore
    data.ProteinModData = json.loads(proteinModJson)

    # Save to temporary variable
    ProteinModData=data.ProteinModData

    # Get JSON data for protein mod viz
    sampleModJson = buildSampleModData(df)

    # Save to datastore
    data.SampleModData = json.loads(sampleModJson)

    # Save to temporary variable
    SampleModData=data.SampleModData

    # Return data to frontend
    return render_template(
        "index.html", 
        Normalization=Normalization, 
        MinCount=MinCount, 
        ProteinModData=ProteinModData, 
        SampleModData=SampleModData,
        )


@application.route("/get-protein-mod-data",methods=["GET","POST"])
@cross_origin()
def returnProteinModData():
    f=data.ProteinModData
    return f

@application.route("/get-sample-mod-data",methods=["GET","POST"])
@cross_origin()
def reutrnSampleModData():
    f=data.SampleModData
    return f

if __name__ == "__main__":
    application.run(debug=True)



