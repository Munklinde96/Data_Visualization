from flask import Flask, flash, redirect, render_template, request, session, abort,send_from_directory,send_file,jsonify
import pandas as pd
import sys
from refactored_utils import get_and_prepare_data
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
    HeatMapSortingMods=None
    HeatMapSortingProteins=None
data=DataStore()

def buildProteinModData(df, count, normalization, _samples=[]):

    # Get data for first viz
    modData = pd.DataFrame(get_modification_count_per_protein(df, count,  normalization))
    # data.HeatMapSortingMods = 
    # Convert to json
    return modData.to_json()

# def buildSampleModData(df, protein):
#     data = get_sample_heatmap_data(df, _protein=protein)
#     return data.to_json()

def buildSampleModData(df, protein):
    data = None
    if protein == "":
        data = get_sample_heatmap_data(df)
    else:
        data = get_sample_heatmap_data(df, _protein=protein)
    return data.to_json()

def buildSegmentData(df, protein, _samples=[], start_pos=0, end_pos=0, _stacked=True, intensity_value=None):
    if len(_samples) == 1:
        if _samples[0] == '':
            _samples = []
    intSamples = []
    for sample in _samples:
        intSamples.append(int(sample))
    # Create start and end indicies
    start_end_indicies = None
    if start_pos != 0 and end_pos != 0:
        start_end_indicies = (start_pos, end_pos)
    # Create SegmentPlotData
    if protein == "" or protein == "null":
        peptide_patches, mod_patches, height, seqq, modification_color_map, min_peptide, max_peptide, histogram_df = create_data_for_segment_plot(df, spacing=0.0, start_end_indicies =start_end_indicies, is_stacked=_stacked, intensity_value=intensity_value)
        histogram_data_json = json.loads(histogram_df.to_json())
        segmentObject = {
        'peptide_patches': peptide_patches,
        'mod_patches': mod_patches,
        'height': height,
        'seqq': seqq,
        'modification_color_map': modification_color_map,
        'min_peptide': min_peptide,
        'max_peptide': max_peptide,
        'histogram_data': histogram_data_json
        }
        segmentPlotJson = json.dumps(segmentObject)
        return segmentPlotJson
    else:
        peptide_patches, mod_patches, height, seqq, modification_color_map, min_peptide, max_peptide, histogram_df = create_data_for_segment_plot(df, _protein=protein, selected_sample_indices=intSamples, spacing = 0.0, start_end_indicies =start_end_indicies, is_stacked=_stacked, intensity_value=intensity_value)
        histogram_data_json = json.loads(histogram_df.to_json())
        segmentObject = {
        'peptide_patches': peptide_patches,
        'mod_patches': mod_patches,
        'height': height,
        'seqq': seqq,
        'modification_color_map': modification_color_map,
        'min_peptide': min_peptide,
        'max_peptide': max_peptide,
        'histogram_data': histogram_data_json 
        }
        segmentPlotJson = json.dumps(segmentObject)
        return segmentPlotJson
    
@application.route("/main",methods=["GET","POST"])
#3. Define main code
@application.route("/",methods=["GET","POST"])
def homepage():
    # Get data and build dataframe
    # path = r"protein-peptides.csv"
    path = r"UHT milk P036.csv"
    #path = r"protein-peptides.csv"
    df = get_and_prepare_data(path)
    data.DataFrame = df
    
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
    # sampleModJson = buildSampleModData(df, "")

    # Save to datastore
    # data.SampleModData = json.loads(sampleModJson)

    # segmentPlotJson = buildSegmentData(df, "")
    # data.SegmentPlotData = segmentPlotJson


    # Return data to frontend
    return ("", 200)


@application.route("/get-protein-mod-data",methods=["GET","POST"])
@cross_origin()
def returnProteinModData():
    print("received get-protein-mode-data request")
    minModCount = request.args.get('min_mod_count', default = 0, type = int)
    normalization = request.args.get('normalization_type', default = "", type = str)
    samples = request.args.get('samples', default = "", type = str).split(",")
    proteinModJson = buildProteinModData(data.DataFrame, minModCount, normalization)
    data.ProteinModData = json.loads(proteinModJson)
    f=data.ProteinModData
    return f

@application.route("/get-sample-mod-data",methods=["GET","POST"])
@cross_origin()
def reutrnSampleModData():
    print("received get-sample-mod-data request")
    protein = request.args.get('protein', default = "", type = str)
    # Get JSON data for protein mod viz
    sampleModJson = buildSampleModData(data.DataFrame, protein)
    # Save to datastore
    data.SampleModData = json.loads(sampleModJson)
    f=data.SampleModData
    return f

@application.route("/get-segment-data",methods=["GET","POST"])
@cross_origin()
def returnSegmentData():
    print("received get-segment-data request")
    protein = request.args.get('protein', default = "", type = str)
    samples = request.args.get('samples', default = "", type = str).split(",")
    segmentPlotJson = buildSegmentData(data.DataFrame, protein, _samples=samples)
    data.SegmentPlotData = segmentPlotJson
    f=data.SegmentPlotData
    return f

@application.route("/get-segment-protein-data",methods=["GET","POST"])
@cross_origin()
def returnSegmentProteinData():
    print("received get-segment-protein-data request")
    startPos = request.args.get('start_pos', default = 0, type = int)
    endPos = request.args.get('end_pos', default = 0, type = int)
    protein = request.args.get('protein', default = "", type = str)
    intensity = request.args.get('intensity', default = None, type = float)
    samples = request.args.get('samples', default = "", type = str).split(",")
    segmentPlotJson = buildSegmentData(data.DataFrame, protein, _samples=samples, start_pos=startPos, end_pos=endPos, _stacked=False, intensity_value=intensity)
    data.SegmentPlotData = segmentPlotJson
    f=data.SegmentPlotData
    return f


if __name__ == "__main__":
    application.run(debug=True)