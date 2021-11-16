from flask import Flask, flash, redirect, render_template, request, session, abort,send_from_directory,send_file,jsonify
import pandas as pd
import sys
sys.path.append("./..")
from utils import get_modification_count_per_protein
from sample_heatmap import get_sample_heatmap_data
from segment_plot import create_data_for_segment_plot

from flask_cors import CORS, cross_origin
import json

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
    SegmentPlotData=None
    DataFrame=None
data=DataStore()

def buildProteinModData(df, count, normalization):
    # Count modifications
    df['#modifications'] = df['PTM'].apply(count_no_of_modifications)
    df[df['#modifications'] > 0]

    # Get data for first viz
    modData = pd.DataFrame(get_modification_count_per_protein(df, count,  normalization))
    
    # Convert to json
    return modData.to_json()

def buildSampleModData(df, count):
    df['#modifications'] = df['PTM'].apply(count_no_of_modifications)
    df[df['#modifications'] > count]
    data = get_sample_heatmap_data(df)
    return data.to_json()

@application.route("/main",methods=["GET","POST"])
#3. Define main code
@application.route("/",methods=["GET","POST"])
def homepage():
    # Get data and build dataframe
    path = r"UHT milk P036.csv"
    df = pd.read_csv(path)
    data.DataFrame = df
    #df1, df2, df3, df4 = split_data_in_samples(df)
    
    # Create initial configuration
    data.Normalization = "protein_intensity"
    data.MinCount = 0
    Normalization=data.Normalization
    MinCount=data.MinCount

    # Get JSON data for protein mod viz
    proteinModJson = buildProteinModData(df, MinCount, Normalization)

    # Save to datastore
    data.ProteinModData = json.loads(proteinModJson)

    # Get JSON data for protein mod viz
    sampleModJson = buildSampleModData(df, data.MinCount)

    # Save to datastore
    data.SampleModData = json.loads(sampleModJson)

    # Create SegmentPlotData
    peptide_patches, mod_patches, height, seqq, modification_color_map, min_peptide, max_peptide = create_data_for_segment_plot(df, start_end_indices=(16, 40))
    segmentObject = {
        'peptide_patches': peptide_patches,
        'mod_patches': mod_patches,
        'height': height,
        'seqq': seqq,
        'modification_color_map': modification_color_map,
        'min_peptide': min_peptide,
        'max_peptide': max_peptide
    }
    segmentPlotJson = json.dumps(segmentObject)
    data.SegmentPlotData = segmentPlotJson


    # Return data to frontend
    return ("", 200)


@application.route("/get-protein-mod-data",methods=["GET","POST"])
@cross_origin()
def returnProteinModData():
    print("received get-protein-mode-data request")
    minModCount = request.args.get('min_mod_count', default = 0, type = int)
    normalization = request.args.get('normalization_type', default = "", type = str)
    proteinModJson = buildProteinModData(data.DataFrame, minModCount, normalization)
    data.ProteinModData = json.loads(proteinModJson)
    f=data.ProteinModData
    return f

@application.route("/get-sample-mod-data",methods=["GET","POST"])
@cross_origin()
def reutrnSampleModData():
    print("received get-sample-mod-data request")
    minModCount = request.args.get('min_mod_count', default = 0, type = int)
    # Get JSON data for protein mod viz
    sampleModJson = buildSampleModData(data.DataFrame, minModCount)
    # Save to datastore
    data.SampleModData = json.loads(sampleModJson)
    f=data.SampleModData
    return f

@application.route("/get-segment-data",methods=["GET","POST"])
@cross_origin()
def returnSegmentData():
    print("received get-segment-data request")
    f=data.SegmentPlotData
    return f


if __name__ == "__main__":
    application.run(debug=True)