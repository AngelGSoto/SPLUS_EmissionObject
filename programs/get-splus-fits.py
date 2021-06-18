'''
Scrit to download fits file from database
'''
# Import the necessary packages 
import splusdata 
import pandas as pd
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
from astropy.wcs import WCS
import os
import argparse
import sys

parser = argparse.ArgumentParser(
    description="""Get colored image and cut image in the r-band""")

parser.add_argument("ra", type=float,
                    default="316.473196",
                    help="RA of the object")

parser.add_argument("dec", type=float,
                    default="-37.144562",
                    help="Dec of the object")

parser.add_argument("--radi", type=float, default=None,
                    help="""Size of the images in pixel""")

parser.add_argument("--band", type=str, default=None,
                    help="""Filter to be download""")

parser.add_argument("--name", type=str, default=None,
                    help="""Name of the object""")

cmd_args = parser.parse_args()
ra = cmd_args.ra
dec = cmd_args.dec

# Radius
rad = int(cmd_args.radi)

# Filer
band = cmd_args.band

# Nome of the object if has
Name = cmd_args.name

# Connect
conn = splusdata.connect('Luis', 'plutarco*80')

# Getting the Fits image in the r-band
hdu = conn.get_cut(ra, dec, rad, band)

# Save the image, note that the output image in compress
hdu.writeto('{}_{}-{}_{}_{}.fz'.format(Name, int(ra), int(dec), rad, band), overwrite=True) # write to fits

############################################################
# Definition to decompress the images ######################
############################################################
def fz2fits(image):
    """
    It converts SPLUS images
    from .fz to .fits
    """
    datos = fits.open(image)[1].data
    heada = fits.open(image)[1].header
    imageout = image[:-2] + 'fits'
    print ('Creating file: ')
    print (imageout)
    fits.writeto(imageout, datos, heada, overwrite=True)
############################################################
# Decompress
fz2fits('{}_{}-{}_{}_{}.fz'.format(Name, int(ra), int(dec), rad, band))

# Read the FITS file
hdul = fits.open('{}_{}-{}_{}_{}.fits'.format(Name, int(ra), int(dec), rad, band))[0]
wcs = WCS(hdul.header)

print("Testing the download:")                 
print(wcs)
# Delete fz file
os.remove('{}_{}-{}_{}_{}.fz'.format(Name, int(ra), int(dec), rad, band))
print("FZ File Removed!")

