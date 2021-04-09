'''
Create file.tex with several figures (table)
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
from astropy.table import Table
import seaborn as sns
import sys
from scipy.optimize import fsolve
import colours


#Read de files
pattern = "spec*.pdf"
file_list = glob.glob(pattern)

# Number of objects
n = len(file_list)
n_col = int(len(file_list)/2.)
count = 0

# Lists
list1, list2 = [], []

for i, a in zip(range(n), file_list):
    if count < n_col:
        count +=1
        list1.append("\includegraphics[width=0.5\linewidth, clip]{"+a+"}")
    else:
        list2.append("\includegraphics[width=0.5\linewidth, clip]{"+a+"}")

list1.sort()
list2.sort()
table_fig = Table([list1, list2],  names=('Spectra 1', 'Spectra 2'), meta={'name': 'first table'})
table_fig.write('table-apectra-lamost.tex', format = "ascii.latex", latexdict=dict(tabletype='longtable*'), caption='Emission line objects', fill_values=[('nan', r'\nodata'), ('0.00', r'\nodata')], overwrite=True) 
