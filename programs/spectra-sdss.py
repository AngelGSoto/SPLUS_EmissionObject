'''
This script makes SDSS spectra only
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
from pathlib import Path
import sys
import os

ROOT_PATH = Path("../../iDR3_n4/SDSS-spectra-all-again/")

# Read the file
parser = argparse.ArgumentParser(
    description="""Make a spectras""")

parser.add_argument("fileSdss", type=str,
                    default="teste-program",
                    help="Name of file, taken the prefix")

parser.add_argument("--ymin", required=False, type=float, default=None,
                    help="""Value y-axis min""")

parser.add_argument("--ymax", required=False, type=float, default=None,
                    help="""Value y-axis max""")


cmd_args = parser.parse_args()
file_spec = cmd_args.fileSdss + ".fits"


hdu = fits.open(ROOT_PATH / file_spec)

# Data from the SDSS spectra
hdudata = hdu[1].data
wl = 10**hdudata.field("loglam")
Flux = 1E-17*hdudata.field("flux")


# PLOTS
fig, ax = plt.subplots(figsize=(12, 9))
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)
ax.set(xlim=[3350,9300])
# set Y-axis range (if applicable)
if cmd_args.ymin is not None and cmd_args.ymax is not None:
    plt.ylim(cmd_args.ymin,cmd_args.ymax)
elif cmd_args.ymin is not None:
    plt.ylim(ymin=cmd_args.ymin)
elif cmd_args.ymax is not None:
    plt.ylim(ymax=cmd_args.ymax)
#plt.ylim(ymin=-50.0,ymax=200)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel=r'F$(\mathrm{10^{-15} erg\ s^{-1} cm^{-2} \AA^{-1}})$')
Flux /=1e-15
ax.plot(wl, Flux, c = "gray", linewidth=3.5, alpha=0.8, zorder=5)
#ax.axvline(4686, color='r', linewidth=0.3, linestyle='-', zorder = 6, label="He II")
plt.annotate(r"H$\alpha$", xy=(7445, 14.2),  xycoords='data', size=23,
            xytext=(-120, -60), textcoords='offset points', 
            bbox=dict(boxstyle="round4, pad=.5", fc="#CC0066", alpha=0.7),)
ax.legend()
plt.tight_layout()
save_file = file_spec.replace(".fits", "-only-spectra.pdf")
plt.savefig(save_file)
