import numpy as np
import pandas as pd
import os
from io import BytesIO
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from skimage import img_as_ubyte
import base64

def clean_df(df):
    """
    df: a Pandas dataframe with a column "SMILES"
    """
    
    #df["CAN_SMILES"] = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) if s is not np.nan else None for s in df["SMILES"]]

    df.drop_duplicates(subset ="CAN_SMILES", keep = 'first', inplace = True)
    df = df.dropna()
    df.reset_index(level=None, drop=True, inplace=True, col_level=0, col_fill='')

    return df

def featurize_morgan(smiles_list):
    X = []
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
            if fp is not None:
                X.append(fp)
            else:
                X.append(None)
        else:
            X.append(None)
    return X

def to_png(arr):
    out = BytesIO()
    im = Image.fromarray(arr)
    im.save(out, format='png')
    return out.getvalue()

def b64_image_files(smiles_list):
    urls = []
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        image = Draw.MolToImage(mol, size=(150,150))
        png = to_png(img_as_ubyte(image))
        url = 'data:image/png;base64,' + base64.b64encode(png).decode('utf-8')
        urls.append(url)
    return urls