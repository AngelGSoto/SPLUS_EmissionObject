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
tab = Table.read("hydra_Halpha_simbad.dat", format="ascii")

# Number of objects
n = len(tab)
n_col = int(len(tab)/4.)
count = 0

# Lists
list1, list2, list3, list4 = [], [], [], []

for i in range(n):
    if count < n_col:
        count +=1
        spec = "Fig_spectra/Simbad/photopectrum_splus_"+str(tab["ID"][i].split("R3.")[-1]).replace(".", "-")+"_hydra_Halpha_simbad_PStotal.pdf"
        list1.append("\includegraphics[width=0.2\linewidth, clip]{"+spec+"}")
    elif count < 2*n_col:
        count +=1
        spec = "Fig_spectra/Simbad/photopectrum_splus_"+str(tab["ID"][i].split("R3.")[-1]).replace(".", "-")+"_hydra_Halpha_simbad_PStotal.pdf"
        list2.append("\includegraphics[width=0.2\linewidth, clip]{"+spec+"}")
    elif count < 3*n_col:
        count +=1
        spec = "Fig_spectra/Simbad/photopectrum_splus_"+str(tab["ID"][i].split("R3.")[-1]).replace(".", "-")+"_hydra_Halpha_simbad_PStotal.pdf"
        list3.append("\includegraphics[width=0.2\linewidth, clip]{"+spec+"}")
    else:
        spec = "Fig_spectra/Simbad/photopectrum_splus_"+str(tab["ID"][i].split("R3.")[-1]).replace(".", "-")+"_hydra_Halpha_simbad_PStotal.pdf"
        list4.append("\includegraphics[width=0.2\linewidth, clip]{"+spec+"}")
        
tab.sort('RA')
table_fig = Table([list1, list2, list3, list4],  names=('S-spectra 1', 'S-spectra 2', 'S-spectra 3', 'S-spectra 4'), meta={'name': 'first table'})
    #table_fig.sort('Auto')
table_fig.write('table-hydra_Halpha_simbad.tex', format = "ascii.latex", latexdict=dict(tabletype='longtable*'),  caption='Emission line objects', fill_values=[('nan', r'\nodata'), ('0.00', r'\nodata')], overwrite=True) 


#latexdict=dict(tabletype='deluxetable*')
#longtable

    
