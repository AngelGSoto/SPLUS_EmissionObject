import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import pandas as pd
import numpy as np
from astropy.table import Table
import seaborn as sns
import argparse
import sys
import os
import glob
import json
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from pathlib import Path
from density_scatter import density_scatter
sns.set_color_codes()
ROOT_PATH = Path("..")

# Reading the json files with synthectic photometry of the star library Pickles, A. J. (1998)
def filter_mag(e, s, f1, f2, f3):
    '''
    Calculate the colors using any of set of filters
    '''
    col, col0 = [], []
    if data['id'].endswith(e):
        if data['id'].startswith(str(s)):
            filter1 = data[f1]
            filter2 = data[f2]
            filter3 = data[f3]
            diff = filter1 - filter2
            diff0 = filter1 - filter3
            col.append(diff)
            col0.append(diff0)
    
    return col, col0

def plot_mag(f1, f2, f3):
    x, y = filter_mag("Star", "", f1, f2, f3)
    for a, b in zip(x, y):
        A1.append(a)
        B1.append(b)

# Read the file
parser = argparse.ArgumentParser(
    description="""Make a table from the S-PLUS catalogs """)

parser.add_argument("fileName", type=str,
                    default="teste-program",
                    help="Name of table, taken the prefix ")

cmd_args = parser.parse_args()
file_ = cmd_args.fileName + ".ecsv"

table = Table.read(file_, format="ascii.ecsv")

ra = table["RA"]
dec = table["DEC"]

icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
gal = icrs.galactic  

l_rad = gal.l.radian
l_rad[l_rad > np.pi] -= 2. * np.pi
b_rad = gal.b.radian

# Sintectin MS track
A1, B1 = [], []

pattern = "../../MS_stars/*.json"
file_list = glob.glob(pattern)

for file_name in file_list:
    with open(file_name) as f:
        data = json.load(f)
        plot_mag("F0626_rSDSS", "F0660", "F0769_iSDSS")

# Plots
color_map = plt.cm.Spectral_r
color_palette = sns.color_palette('Paired', 55)
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize=(15, 11))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    plt.xlabel(r"$r - i$", fontsize=35)
    plt.ylabel(r"$r - J0660$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    scat = ax.scatter(table['r - i'], table['r - J0660'], s=15*table["FWHM"], edgecolor='black',
                             c=table["R_PStotal"], alpha=0.7, zorder = 2, cmap='RdBu_r')
    #pal = sns.dark_palette("magma", as_cmap=True)
    #pal = sns.cubehelix_palette(as_cmap=True)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    ax = sns.kdeplot(B1, A1, zorder = 3, cmap=pal);
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    ax.set(
      xlim=[-3.5, 5.],
      ylim=[-2.0, 6.])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    cb = fig.colorbar(scat,extend='both', ax=ax).set_label("$r-band$", fontsize=35)
    #Symbol size indicates outer shell radius
    plt.text(0.01, 0.95, 'Symbol size indicates FWHM',
             transform=ax.transAxes, fontsize=20)

    # main sequence and giant stars loci
    x1, y1 = 0.3, 0.3
    el = mpatches.Ellipse((x1, y1), 0.3, 0.4, angle=30, alpha=0.3)
    ax.annotate("Contour indicates main-sequence and giant stars loci", xy=(0.1, -0.1), xytext=(-1.5, -1.), color='black', size=20, zorder= 111, arrowprops=dict(arrowstyle="fancy",
                            color="0.5",
                            patchB=el,
                            shrinkB=5,
                            connectionstyle="arc3,rad=0.3",
                                                                                                                                                                 ))
    
    plt.savefig("../paper/Figs/final-emitters.pdf")

#Distribution of Halpha emitters
with sns.axes_style("ticks"):
    fig = plt.figure(figsize=(14,7))
    ax = fig.add_subplot(1,1,1, projection='aitoff')
    plt.xlabel(r'$l (Gal)$')
    plt.ylabel(r'$b (Gal)$')
    ax.xaxis.label.set_fontsize(23)
    ax.yaxis.label.set_fontsize(23)
    plt.tick_params(axis='x', labelsize=23) 
    plt.tick_params(axis='y', labelsize=23)
    #ax.scatter(l_rad, b_rad, s=1, color='black', alpha=0.2)
    density_scatter(l_rad, b_rad, ax=ax)
    ax.grid(True, linestyle='-.', linewidth=0.7)
    #plt.colorbar(image, spacing='uniform', extend='max')
    plt.savefig("../paper/Figs/halpha-emitters-galactic-aitoff.pdf")

    # Bar diagram
    fig1, ax1 = plt.subplots(1, 1, figsize=(10, 5), sharex=True)
    sns.distplot(table["r - J0660"], 
                 norm_hist=True, kde=True, ax=ax1,
                 bins=20, hist_kws=dict(range=[-3.0, 3.0], color='r')
                )
    #ax1.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    ymax = ax.get_ybound()[1]
    sns.despine()
    plt.savefig("../paper/Figs/distribution-Halpha.pdf")

    # Distribution r - i color
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 5), sharex=True)
    sns.distplot(table["r - i"], 
                 norm_hist=True, kde=True, ax=ax2,
                 bins=20, hist_kws=dict(range=[-3.0, 3.0], color='r')
                )
    #ax2.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    ymax = ax.get_ybound()[1]
    sns.despine()
    plt.savefig("../paper/Figs/distribution-ri.pdf")

    # Distribution  r-mag
    fig3, ax3 = plt.subplots(1, 1, figsize=(10, 5), sharex=True)
    sns.distplot(table["R_PStotal"], 
                 norm_hist=False, kde=True, ax=ax3,
                 bins=20)#, hist_kws=dict(range=[-3.0, 3.0])
                #)
    #ax3.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.savefig("../paper/Figs/distribution_r.pdf")

    # Distribution b coordinate
    fig4, ax4 = plt.subplots(1, 1, figsize=(10, 5), sharex=True)
    sns.distplot(b_rad, 
                 norm_hist=True, kde=True, ax=ax4,
                 bins=50, hist_kws=dict(range=[-3.0, 3.0])
                )
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.savefig("../paper/Figs/distribution-bgalactic.pdf")
