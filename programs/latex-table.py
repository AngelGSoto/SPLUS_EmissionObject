'''
Create latex table
'''
from __future__ import print_function
import numpy as np
from astropy.io import fits
import os
import glob
import json
import matplotlib.pyplot as plt
import pandas as pd
#import StringIO
from astropy.table import Table, vstack
import seaborn as sns
import sys
from scipy.optimize import fsolve
import colours
from pathlib import Path
ROOT_PATH = Path("../paper")

# Read the files
df = pd.read_csv("simbad.csv")

df_blue = pd.read_csv("Simbad-Blue0-Good-LD-Halpha-DR3_noFlag_merge-takeoutbad-Final.csv")

df_red = pd.read_csv("Simbad-Red1-Good-LD-Halpha-DR3_noFlag_merge-takeoutbad-Final.csv")

# Replace values on Label
df_blue["Label_hier"] = df_blue["Label_hier"].replace(0, "Blue")
df_red["Label_hier"] = df_red["Label_hier"].replace(1, "Red")

# Join the other two tables
df_cont = pd.concat([df_blue, df_red])

# Converting in astropy table
tab = Table.from_pandas(df)
tab_cont = Table.from_pandas(df_cont)

# Column in main table
tab['Label_hier'] = '-'

#a = np.array(tab_cont['Label_hier'], dtype=str)
#tab_cont['Label_hier'] = a

# Masking and applying
id1 = tab["ID"]
id2 = tab_cont["ID"]
mask = np.array([not source in id2 for source in id1])

tab_drop = tab[mask]

tab_final = vstack([tab_cont, tab_drop])

Col = ["main_type", "RA","DEC", "Label_hier"]

tab_final.sort('RA')
tab_final[Col].write(ROOT_PATH / 'table-simbad-sort.tex', format = "ascii.latex") 


