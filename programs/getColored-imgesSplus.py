'''
Scrit to download colored images from database
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

parser.add_argument("--name", type=str, default=None,
                    help="""Name of the object""")

cmd_args = parser.parse_args()
ra = cmd_args.ra
dec = cmd_args.dec

# Radius
rad = int(cmd_args.radi)

# Nome of the object if has
Name = cmd_args.name

# Connect
conn = splusdata.connect('Luis', 'plutarco*80')

# Getting the colored imge
img = conn.twelve_band_img(ra, dec, radius=rad, noise=0.15, saturation=0.15)

# Getting the Fits image in the r-band
hdu = conn.get_cut(ra, dec, rad, 'R')

# Save the image, note that the output image in compress
hdu.writeto('{}_{}-{}_{}_r.fz'.format(Name, ra, dec, rad), overwrite=True) # write to fits

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
hdufits = fz2fits('{}_{}-{}_{}_r.fz'.format(Name, ra, dec, rad))

# Read the FITS file
hdul = fits.open('{}_{}-{}_{}_r.fits'.format(Name, ra, dec, rad))[0]
wcs = WCS(hdul.header)

print(wcs)                 

f = plt.figure(figsize=(18,9))

ax1 = aplpy.FITSFigure(hdul, figure=f, subplot=(1, 1, 1))#, north=True)
plt.imshow(img, origin='lower', cmap='cividis', aspect='equal')
                 
ax1.add_scalebar(20.0/3600)
ax1.scalebar.set_label('20 arcsec')
ax1.scalebar.set(color='yellow', linewidth=4, alpha=0.9)
ax1.scalebar.set_font(size=35, weight='bold',
                      stretch='normal', family='sans-serif',
                      style='normal', variant='normal')

ax1.axis_labels.set_font(size=22, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
ax1.axis_labels.hide()
ax1.tick_labels.hide()
#ax1.axis_labels.hide_y()

ax1.tick_labels.set_font(size=22, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
#ax1.list_layers()
#ax1.show_markers(ra, dec, layer='marker', edgecolor='green', facecolor='none', marker='o', s=10, alpha=0.9, linewidths=60)#, layer='marker_set_1', edgecolor='black', facecolor='none', s=30, alpha=0.5, linewidths=20)


# ax1.axis_labels.hide_y()
# ax1.tick_labels.hide_y()

#ax2.colorbar.set_box([0.95, 0.1, 0.015, 0.8])
ax1.set_theme('publication')
#f.tight_layout()
#f.savefig("-".join([image_name, "images.pdf"]))

plt.savefig('{}_{}-{}_{}_r.fits'.format(Name, ra, dec, rad).replace(".fits", ".pdf"))
