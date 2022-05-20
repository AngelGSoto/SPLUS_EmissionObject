'''
This scriti allows us to put in a latex file spectra and colored images of the SPLUS sources
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
from astropy.table import Table, hstack
import seaborn as sns
import argparse
import sys

parser = argparse.ArgumentParser(
    description="""Make a table from the S-PLUS catalogs """)

parser.add_argument("source", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")


cmd_args = parser.parse_args()
file_ = cmd_args.source + ".ecsv"

tab = Table.read(file_, format="ascii.ecsv")

#Read de files
pattern1 = "*-emitters/images/*.pdf"
file_list1 = glob.glob(pattern1)

pattern2 = "*-emitters/S-spectra/*.pdf"
file_list2 = glob.glob(pattern2)


pattern3 = "remaindOld-images/*.pdf"
file_list3 = glob.glob(pattern3)

pattern4 = "remaindOld-S-spectra/*.pdf"
file_list4 = glob.glob(pattern4)

file_list_img = file_list1 + file_list3
file_list_spec = file_list2 + file_list4

fig_template = r'\includegraphics[width=0.4\linewidth, clip]{{{:s}}}'
fig_template_ = r'\includegraphics[width=0.3\linewidth, clip]{{{:s}}}'

id_ = []
img = []
print(len(tab))
for i in tab:
    for j in file_list_img:
        if i["ID"].split("R3.")[-1].replace(".", "-") == j.split("es/")[-1].split("_")[0]:
            idd = i["ID"]
            id_.append(idd)
            img.append(fig_template.format(j))   

newtab1 = Table([id_,  img], names=('ID', 'Image'), meta={'name': 'first table'})

id_1 = []
spec = []
for i in tab:
    for k in file_list_spec:
        if i["ID"].split("R3.")[-1].replace(".", "-") == k.split("splus_")[-1].split("_", 1)[0]:
            idd = i["ID"]
            id_1.append(idd)
            spec.append(fig_template_.format(k))
            
                
newtab2 = Table([id_1, spec], names=('ID_1', 'S-spectra'), meta={'name': 'first table'})

newtab = hstack([newtab1, newtab2])

Col = ['ID', 'S-spectra', 'Image']

tabnew_final = newtab[Col]

latexname = "table-plots-Halpha-emitters-final.tex"
tabnew_final.write(latexname, format = "ascii.latex",
                     fill_values=['--'],
                     overwrite=True)
