'''
02/07/21
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
#from PyAstronomy import pyasl
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from pathlib import Path
ROOT_PATH = Path("../paper")

def format_RA(ra):
    return coord.Angle(ra, unit=u.deg).to_string(u.hour, sep=':', precision=2, pad=True)

def format_DEC(dec):
    s = coord.Angle(dec, unit=u.deg).to_string(sep=':', precision=1, pad=True)
    if s.startswith('-'):
        return r'$-$' + s[1:]
    else:
        return s

def replace_1(colum):
    if '[' in colum:
        return colum.replace('[', '$[$')
    else:
        return colum

def replace_2(colum1):
    if ']' in colum1:
        return colum1.replace(']', '$]$')
    else:
        return colum1

def replace_3(colum):
    if '_' in colum:
        return colum.replace('_', ' ')
    else:
        return colum

def ra_fmt(x):
    """Write RA  to accuracy of 0.2 deg"""
    return "{:.2f}".format(x)

def dec_fmt(x):
    """Write DEC to accuracy of 0.1 deg"""
    return "{:.1f}".format(x)

def prob_fmt(p):
    """Write proybablity to accuracy of 0.01"""
    return "{:.2f}".format(p)

# Read the files
df = pd.read_csv("simbad.csv")

df_blue = pd.read_csv("Simbad-Blue0-Good-LD-Halpha-DR3_noFlag_merge-takeoutbad-Final.csv")

df_red = pd.read_csv("Simbad-Red1-Good-LD-Halpha-DR3_noFlag_merge-takeoutbad-Final.csv")

df_hdbscan = pd.read_csv("Simbad-hdbscan-Good-LD-Halpha-DR3_noFlag_merge-takeoutbad-Final.csv")


# Replace values on Label
df_blue["Label_hier"] = df_blue["Label_hier"].replace(0, "Blue")
df_red["Label_hier"] = df_red["Label_hier"].replace(1, "Red")

# Join the other two tables
df_cont = pd.concat([df_blue, df_red])

# Converting in astropy table
tab = Table.from_pandas(df)
tab_cont = Table.from_pandas(df_cont)
tab_hdbscan = Table.from_pandas(df_hdbscan)

tab_cont.sort('RA')
tab_hdbscan.sort('RA')


blue = [prob_fmt(a) for a in tab_hdbscan['P(Blue)']]
red = [prob_fmt(a) for a in tab_hdbscan['P(red)']]

# tab_hdbscan['P(Blue)'].format = "%.2f"
# tab_hdbscan['P(red)'].format = "%.2f"

# Additing columns with the probabilities on table of HBDSCAN
tab_cont['P(Blue)'] = blue
tab_cont['P(red)'] = red

# print(tab_cont['P(red)'])
# sys.exit()
# Column in main table
tab['Label_hier'] = '--'
tab['P(Blue)'] = '--'
tab['P(red)'] = '--'

#a = np.array(tab_cont['Label_hier'], dtype=str)
#tab_cont['Label_hier'] = a

# Masking and applying
id1 = tab["ID"]
id2 = tab_cont["ID"]
mask = np.array([not source in id2 for source in id1])

tab_drop = tab[mask]

tab_final = vstack([tab_cont, tab_drop])
#c = SkyCoord(ra=tab_final["RA"]*u.degree, dec=tab_final["DEC"]*u.degree)

#ra = c.ra.to_string(u.hour, sep=':',  precision=2, pad=True)
# dec = c.dec.to_string(u.degree, alwayssign=True, sep=':',  precision=1, pad=True)

# Selected columns
Col = ["main_id", "RA", "DEC", "main_type", "Label_hier", 'P(Blue)', 'P(red)']

latex_columns = ['Id Simbad', 'RA', 'DEC',
                 'Type' ,
                 r'Group -- {\sc hac}',
                 r'P(Blue) -- {\sc hdbscan}',
                 r'P(Red) -- {\sc hdbscan}'
                 ]

column_formats = {}

column_formats['RA'] = format_RA
column_formats['DEC'] = format_DEC
column_formats['Id Simbad'] = replace_1
column_formats['Id Simbad'] = replace_2
column_formats['Type'] = replace_3

tab_final.sort('RA')
tab_final[Col].write(ROOT_PATH / 'table-simbad-sort.tex', format = "ascii.latex",
                     col_align='r'*len(latex_columns),
                     names=latex_columns,
                     formats=column_formats,
                     fill_values=[('nan', r'0.00')],
                     overwrite=True) 


