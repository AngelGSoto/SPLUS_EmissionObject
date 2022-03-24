'''
02/23/22
Create final table with all information
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
from astropy.table import Table, vstack, Column, MaskedColumn
import seaborn as sns
import sys
from scipy.optimize import fsolve
import colours
#from PyAstronomy import pyasl
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from pathlib import Path
ROOT_PATH = Path("../paper")

# Read the files
df = pd.read_csv("Final-list-emitters-allparam-unique.csv")

tab_blue = Table.read("Blue0-Good-LD-Halpha-DR3_noFlag_merge-takeoutbad-Final.ecsv", format="ascii.ecsv")
df_blue = tab_blue.to_pandas()
tab_red = Table.read("Red1-Good-LD-Halpha-DR3_noFlag_merge-takeoutbad-Final.ecsv", format="ascii.ecsv")
df_red = tab_red.to_pandas()

df_hdbscan = pd.read_csv("Final-list-emitters-allparam-unique-hdbscan.csv")

# Join the other two tables
tab_hier = vstack([tab_blue, tab_red])

#coverting pandas to astropy tables
tab = Table.from_pandas(df)
tab_hdbscan = Table.from_pandas(df_hdbscan)

# Sorting
tab_hier.sort('ID')
tab_hdbscan.sort('ID')

# Additing columns with the probabilities on table of hier
tab_hier['Label_hdbscan'] = tab_hdbscan["Label"]
tab_hier['P(Blue)'] = tab_hdbscan["P(Blue)"]
tab_hier['P(Red)'] = tab_hdbscan["P(Red)"]

# Column in main table
tab['Label_hier'] = np.NaN
tab['Label_hdbscan'] = np.NaN
tab['P(Blue)'] = np.NaN
tab['P(Red)'] = np.NaN

# Masking and applying
id1 = tab["ID"]
id2 = tab_hier["ID"]
mask = np.array([not source in id2 for source in id1])

tab_drop = tab[mask]

tab_final = vstack([tab_hier, tab_drop])

#Removing columns
tab_final.remove_column('Label_k')
tab_final.remove_column('Label_sp')

tab_final.sort('ID')

#Saving the results
tab_final.write("Halpha-emiters-splusdr3.dat", format="ascii")
#Pandas
df_final = tab_final.to_pandas()
df_final.to_csv("Halpha-emiters-splusdr3.csv", index=False)
