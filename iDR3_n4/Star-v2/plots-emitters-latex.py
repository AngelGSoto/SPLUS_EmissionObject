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
from astropy.table import Table
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
pattern1 = "images/*.pdf"
file_list1 = glob.glob(pattern1)


pattern2 = "S-spectra/*.pdf"
file_list2 = glob.glob(pattern2)

#Table
#tab = Table.read("Halpha-DR3-SPLUS-PStotal-STAR-r16-v2.ecsv", format="ascii.ecsv")

fig_template = r'\includegraphics[width=0.4\linewidth, clip]{{{:s}}}'
fig_template_ = r'\includegraphics[width=0.3\linewidth, clip]{{{:s}}}'

# file_list1.sort()
# file_list2.sort()
# tab["ID"].sort()

id_ = []
img = []
spec = []
for i, ii in zip(file_list1, file_list2):
    idd1 = "iDR3." + i.split("s/")[-1].rsplit("-", 1)[0]
    idd2 = i.rsplit("-", 1)[-1].split("_")[0]
    idd = idd1 + "." + idd2
    id_.append(idd)
    img.append(fig_template.format(i))
    spec.append(fig_template_.format(ii))

id_.sort()
img.sort()
spec.sort()    
newtab = Table([id_, spec, img], names=('ID', 'S-spectra', 'Image'),
         meta={'name': 'first table'})

latexname = "table-plots-{}.tex".format(file_.split("STAR-")[-1].split("-v2")[0])
newtab.write(latexname, format = "ascii.latex",
                     fill_values=['--'],
                     overwrite=True)

