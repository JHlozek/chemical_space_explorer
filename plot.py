import numpy as np
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

import bokeh
from bokeh.models import PanTool, ResetTool, HoverTool, WheelZoomTool, ColumnDataSource
from bokeh import plotting as bplot
from bokeh.plotting import figure, show

import umap

from processing import clean_df, featurize_morgan, b64_image_files

def plot(df):
    df = process(df)
    
    TOOLS = [PanTool, WheelZoomTool, ResetTool]
    tooltip = """
        <div>
            <div>
                <img
                src="@image_files" height="200" alt="image" width="200"
                style="float: left; margin: 0px 15px 15px 0px;"
                border="2"
                ></img>
            </div>
            <div>
                <span style="font-size: 17px;">@source_filenames</span>
            </div>
        </div>
              """

    hover = HoverTool(tooltips=tooltip)
    tools = [t() for t in TOOLS] + [hover]

    p = figure(width=600, height=600, tools=tools)
    p.scatter(x='x', y='y', size=8, color="navy", source=ColumnDataSource(df.drop(["fp"], axis=1)))
    return p
    
def process(df):
    df["fp"] = featurize_morgan(df["SMILES"].tolist())
    df = clean_df(df)
    
    transformed_data = create_umap(df["fp"].tolist())
    df['x'] = transformed_data.T[0]
    df['y'] = transformed_data.T[1]

    df['image_files'] = b64_image_files(df['SMILES'].tolist())
    df["source_filenames"] = df.index #df["ID"]
    return df

def create_umap(fp_list):
    umap_transformer = umap.UMAP(a=0.001, b=1.5, min_dist=0.1)
    return umap_transformer.fit_transform(fp_list)
