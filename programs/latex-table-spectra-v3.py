'''
Create file.tex with several figures (table).
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
pattern = "Figs/Figs-lamost/spec*.pdf"
file_list = glob.glob(pattern)

# Number of objects
nrows = len(file_list)

fig_template = r'\includegraphics[width=0.3\linewidth, clip]{{{:s}}}'

NCOLS = 6
thefiles = []
thistable = []
for filename in file_list:
    thefiles.append(fig_template.format(filename))
    if len(thefiles) == NCOLS:
        thistable.append(thefiles)
        thefiles = []

if thefiles:
    # catch any partial row at the end
    thistable.append(thefiles)
             

def output_row(row):
    return " & ".join(row) + r" \\"


def output_figtable(table):
    result = r"\begin{center}" + "\n"
    result += r"  \begin{longtable}{" + "l "*NCOLS + "}" + "\n"
    result += r"  \caption{Espectra from LAMOST DR6" + " " + "\label{tab:spec-lamost}}" + "\\" + "\n"
    result += r"  \endfirsthead" + "\n"
    result += r"  \caption[]{--continued}\\" + "\n"
    result += r"  \endhead" + "\n" 
    result += r"  \hline \endfoot" + "\n"
    for row in table:
        result += "    " + output_row(row) + "\n"
    result += r"  \end{longtable}" + "\n"
    result += r"\end{center}" + "\n"
    return result

with open("S-spectra-lamost.tex", 'w') as f:
    f.write(output_figtable(thistable))
