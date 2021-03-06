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

# Reading the json files with synthectic photometry of the star and giant library Pickles, A. J. (1998)
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
    description="""Make a table from the S-PLUS catalogs""")

parser.add_argument("fileName", type=str,
                    default="teste-program",
                    help="Name of table, taken the prefix")

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
###############################################
# New colors    ###############################
###############################################
m = (table["e_G_PStotal"] <= 0.2) & (table["e_R_PStotal"] <= 0.2) & (table["e_Z_PStotal"] <= 0.2)
m1 = (table["e_U_PStotal"] <= 0.2) & (table["e_G_PStotal"] <= 0.2) & (table["e_R_PStotal"] <= 0.2) 
zg = table['Z_PStotal'] - table['G_PStotal']
gr = table['G_PStotal'] - table['R_PStotal']
ug = table['U_PStotal'] - table['G_PStotal']

##############################################
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
    maskfw = table["FWHM"] < 240
    scat = ax.scatter(table['r - i'][maskfw], table['r - J0660'][maskfw], s=5*table["FWHM"][maskfw], edgecolor='black',
                             c=table["R_PStotal"][maskfw], alpha=0.7, zorder = 2, cmap='RdBu_r')
    #pal = sns.dark_palette("magma", as_cmap=True)
    #pal = sns.cubehelix_palette(as_cmap=True)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    axx = sns.kdeplot(B1, A1, zorder = 3, cmap=pal);
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    ax.set(
      xlim=[-3.5, 5.],
      ylim=[-2.0, 6.])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    cb = fig.colorbar(scat, extend='both', ax=ax)#
    cb.set_label("$r-mag$", fontsize=35)
    cb.ax.tick_params(labelsize=30)
    #Symbol size indicates outer shell radius
    plt.text(0.02, 0.95, 'Symbol size indicates FWHM',
             transform=ax.transAxes, fontsize=21)

    # main sequence and giant stars loci
    x1, y1 = 0.3, 0.3
    el = mpatches.Ellipse((x1, y1), 0.3, 0.4, angle=30, alpha=0.3)
    ax.annotate("Contour indicates main-sequence and giant stars loci",
                xy=(0.1, -0.1), xytext=(-1.8, -1.), color='black', size=21,
                zorder= 111, arrowprops=dict(arrowstyle="fancy",
                            color="0.5",
                            patchB=el,
                            shrinkB=5,
                            connectionstyle="arc3,rad=0.3",
                                                        ))
    
    plt.savefig("../paper/Figs/final-emitters.pdf")
    ##########################################################
    # (g - r) vs (z - g)
    fig, ax1 = plt.subplots(figsize=(15, 11))
    ax1.spines["top"].set_visible(False)  
    ax1.spines["right"].set_visible(False)
    plt.xlabel(r"$z - g$", fontsize=35)
    plt.ylabel(r"$g - r$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    ax1.set(
        xlim=[-6.8, 2.5], ylim=[-3., 5.]
      )
    #scat = ax.scatter(zg, gr, s=15*table["FWHM"], edgecolor='black',
                             #c=table["R_PStotal"], alpha=0.7, zorder = 2, cmap='RdBu_r')
    # Limiting the blue and red region
    x_new = np.linspace(-15.0, 1000, 200)
    y = 0.45*x_new + 1.48

    ax1.plot(x_new, y, color='k', zorder=100, linestyle='-.')
    density_scatter(zg[m], gr[m], ax=ax1)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    #ax1.set(
      #xlim=[-15, 15],
      #ylim=[-15, 15])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    #cb = fig.colorbar(scat,extend='both', ax=ax).set_label("$r-band$", fontsize=35)   
    plt.savefig("../paper/Figs/red-blue-colorObjects-gr.pdf")

    ##########################################################
    # (g - r) vs (u - g)
    fig, ax11 = plt.subplots(figsize=(15, 11))
    ax11.spines["top"].set_visible(False)  
    ax11.spines["right"].set_visible(False)
    plt.xlabel(r"$u - g$", fontsize=35)
    plt.ylabel(r"$g - r$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    #scat = ax.scatter(zg, gr, s=15*table["FWHM"], edgecolor='black',
                             #c=table["R_PStotal"], alpha=0.7, zorder = 2, cmap='RdBu_r')
    density_scatter(ug[m1], gr[m1], ax=ax11)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    #ax1.set(
      #xlim=[-15, 15],
      #ylim=[-15, 15])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    #cb = fig.colorbar(scat,extend='both', ax=ax).set_label("$r-band$", fontsize=35)   
    plt.savefig("../paper/Figs/red-blue-colorObjects-ug.pdf")
###################################################################################################################
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
    ##########################
    # Bar diagram
    fig1, ax1 = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
    plt.xlabel(r"$r - J0660$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    r_j0660 = [x for x in table["r - J0660"]]
    g = sns.distplot(r_j0660, 
                 norm_hist=True, kde=True, ax=ax1,
                 bins=20, hist_kws=dict(range=[-3.0, 3.0], color='r')
                )
    #ax1.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    ymax = ax.get_ybound()[1]
    sns.despine()
    plt.tight_layout()
    plt.savefig("../paper/Figs/distribution-Halpha.pdf")
    ##########################
    # Distribution r - i color
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
    plt.xlabel(r"$r - i$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    r_i = [x for x in table["r - i"]]
    sns.distplot(r_i, 
                 norm_hist=True, kde=True, ax=ax2,
                 bins=20, hist_kws=dict(range=[-3.0, 3.0], color='r')
                )
    #ax2.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    ymax = ax.get_ybound()[1]
    sns.despine()
    plt.tight_layout()
    plt.savefig("../paper/Figs/distribution-ri.pdf")
    #########################
    # Distribution  r-mag
    fig3, ax3 = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
    plt.xlabel(r"$r$", fontsize=33)
    #plt.ylabel(r"Density", fontsize=28)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    r = [x for x in table["R_PStotal"]]
    sns.distplot(r, 
                 norm_hist=False, kde=True, ax=ax3,
                 bins=20)#, hist_kws=dict(range=[-3.0, 3.0])
                #)
    #ax3.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../paper/Figs/distribution_r.pdf")
    #########################
    # Distribution b coordinate
    fig4, ax4 = plt.subplots(1, 1, figsize=(10, 5), sharex=True)
    plt.xlabel(r"$b(Gal)$", fontsize=33)
    plt.ylabel(r"# of sources", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    sns.distplot(b_rad, 
                 norm_hist=False, kde=False, ax=ax4,
                 bins=30, hist_kws=dict(range=[-3.0, 3.0],  color='y')
                )
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../paper/Figs/distribution-bgalactic.pdf")
    #########################
    # r vs b (Gal)
    fig4, ax4 = plt.subplots(1, 1, figsize=(9, 12), sharex=True)
    plt.xlabel(r"$b(Gal)$", fontsize=33)
    plt.ylabel(r"r", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    density_scatter(b_rad, table["R_PStotal"], ax=ax4)
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../paper/Figs/bvsr.pdf")
    ###############################################################
    # Distribution z - g coordinate
    fig5, ax5 = plt.subplots(1, 1, figsize=(11, 5), sharex=True)
    plt.xlabel(r"$z - g$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    ax5.set(
      xlim=[-5.0, 2.5]
      )
    zg = [x for x in zg[m]]
    sns.distplot(zg, 
                 norm_hist=True, kde=True, ax=ax5,
                 bins=50, hist_kws=dict(range=[-6.0, 6.0],  color='r')
                )
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.tight_layout()
    plt.savefig("../paper/Figs/distribution-zg.pdf")
    ###############################################################
    # Distribution z - g coordinate
    fig6, ax6 = plt.subplots(1, 1, figsize=(11, 5), sharex=True)
    plt.xlabel(r"$g - r$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    ax6.set(
      xlim=[-1.5, 2.5]
      )
    gr = [x for x in gr[m]]
    sns.distplot(gr, 
                 norm_hist=True, kde=True, ax=ax6,
                 bins=80, hist_kws=dict(range=[-6.0, 6.0],  color='r')
                )
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../paper/Figs/distribution-gr.pdf")
    
