'''
This script make the Lamost spectra overlapped the SPLUS photometry
'''
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sn
sn.set_context("poster")
import glob
import argparse
import sys
import os
from astropy.visualization import hist
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.crossmatch import crossmatch_angular

# Read the file
parser = argparse.ArgumentParser(
    description="""Make a spectras""")

parser.add_argument("fileLamost", type=str,
                    default="teste-program",
                    help="Name of file, taken the prefix")

parser.add_argument("TableSplus", type=str,
                    default="teste-program",
                    help="Name of table, taken the prefix")

cmd_args = parser.parse_args()
file_spec = cmd_args.fileLamost + ".fits"
file_table = cmd_args.TableSplus + ".dat"

hdulist = fits.open(file_spec)
datadir = "../"

try:
    table = Table.read(os.path.join(datadir, file_table), format="ascii")
except FileNotFoundError:
    table = Table.read(file_table, format="ascii")
    

# Data from the lamost spectra
hdu = hdulist[0]
nx, wav0, i0, dwav = [hdu.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CDELT1")]
wl = wav0 + (np.arange(nx) - (i0 - 1))*dwav

Flux = hdulist[0].data
# Data of the SPLUs list
mag, mag_err = [], []
wl_sp = [3485, 3785, 3950, 4100, 4300, 4803, 5150, 6250, 6600, 7660, 8610, 9110]
color = ["#CC00FF", "#9900FF", "#6600FF", "#0000FF", "#009999", "#006600", "#DD8000", "#FF0000", "#CC0066", "#990033", "#660033", "#330034"]
marker = ["s", "o", "o", "o", "o", "s", "o", "s", "o", "s", "o", "s"] ### tienen todos los filtros

mag.append(table["U_PStotal"]) 
mag.append(table["F378_PStotal"])
mag.append(table["F395_PStotal"])
mag.append(table["F410_PStotal"])
mag.append(table["F430_PStotal"])
mag.append(table["G_PStotal"])
mag.append(table["F515_PStotal"]) 
mag.append(table["R_PStotal"]) 
mag.append(table["F660_PStotal"])
mag.append(table["I_PStotal"]) 
mag.append(table["F861_PStotal"]) 
mag.append(table["Z_PStotal"])

#ERRO PStotal
mag_err.append(table["e_U_PStotal"])
mag_err.append(table["e_F378_PStotal"])
mag_err.append(table["e_F395_PStotal"])
mag_err.append(table["e_F410_PStotal"])
mag_err.append(table["e_F430_PStotal"])
mag_err.append(table["e_G_PStotal"])
mag_err.append(table["e_F515_PStotal"]) 
mag_err.append(table["e_R_PStotal"]) 
mag_err.append(table["e_F660_PStotal"]) 
mag_err.append(table["e_I_PStotal"])
mag_err.append(table["e_F861_PStotal"])
mag_err.append(table["e_Z_PStotal"])
# ff = (10**(-(table["R_PStotal"][ind] + 2.41) / 2.5)) / 6250.0**2
# print(ff)
# for i, ii in zip(wl, Flux):
#     if i> 6000 and i< 6300:
#          print(i, ii)

# Find scale factor
m = wl == 6250.
wl_part = wl[m]

flux_part = Flux[m]
Fsp = (10**(-(table["R_PStotal"] + 2.41) / 2.5)) / 6250.0**2
factor = flux_part / Fsp

# Propagation of error
err_ = []
for wll, magg, magerr in zip(wl_sp, mag, mag_err):
    c = (10**(-2.41/2.5)) / wll**2
    b = -(1 / 2.5)
    err = np.sqrt(((c*10**(b*magg))**2)*(np.log(10)*b*magerr)**2)
    err_.append(err)

# PLOTS
fig, ax = plt.subplots(figsize=(12, 9))
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)
ax.set(xlim=[3350,9300])
ax.set(ylim=[-5.,75])
#plt.ylim(ymin=-50.0,ymax=200)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Normalized flux')
ax.plot(wl, Flux, c = "gray", linewidth=0.7, alpha=0.6, zorder=5)
for wl1, mag, magErr, colors, marker_ in zip(wl_sp, mag, err_, color, marker): #
    F = (10**(-(mag + 2.41) / 2.5)) / wl1**2
    F *= 2.2*factor
    ax.scatter(wl1, F, c = colors, marker=marker_, s=80, zorder=4)
    ax.errorbar(wl1, F, yerr=magErr, marker='.', fmt='.', color=colors, ecolor=colors, elinewidth=3.9, markeredgewidth=3.2, capsize=10)
#ax.axvline(4686, color='r', linewidth=0.3, linestyle='-', zorder = 6, label="He II")
rmag = [r for r in table["R_PStotal"]]

ax.annotate("PN G006.0-41.9", xy=(8500, 29),  xycoords='data', size=13,
            xytext=(-120, -60), textcoords='offset points', 
            bbox=dict(boxstyle="round4,pad=.5", fc="0.9"),)
ax.annotate(str(table["ID"]).split("R3.")[-1].replace(".", "-"), xy=(8500, 24.5),  xycoords='data', size=13,
            xytext=(-120, -60), textcoords='offset points', 
            bbox=dict(boxstyle="round4,pad=.5", fc="0.9"),)
ax.annotate("r =" + format(float(table["R_PStotal"]), '.2f'), xy=(8500, 20),  xycoords='data', size=13,
            xytext=(-120, -60), textcoords='offset points', 
            bbox=dict(boxstyle="round4,pad=.5", fc="0.9"),)
ax.legend()
plt.tight_layout()
asciifile = file_spec.replace(".fits", 
                  "-"+(str(table["ID"]).split("R3.")[-1]).replace(".", "-")+".pdf")
plt.savefig(asciifile)
