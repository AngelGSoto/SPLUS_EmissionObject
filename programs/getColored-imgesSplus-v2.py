'''
Scrit to download colored images from database
'''
# Import the necessary packages 
import splusdata 
import pandas as pd
from astropy.table import Table
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
from astropy.wcs import WCS
import os
import argparse
import sys
from pathlib import Path
ROOT_PATH = Path("..")

parser = argparse.ArgumentParser(
    description="""Get colored image and cut image in the r-band""")

parser.add_argument("source", type=str,
                    default="known-PN-jplus-idr",
                    help="Name of source, taken the prefix ")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

args = parser.parse_args()
file_ = args.source + ".ecsv"


try:
    data = Table.read(ROOT_PATH / file_, format="ascii.ecsv")
except FileNotFoundError:
    file_ = args.source + ".dat"
    data = Table.read(ROOT_PATH / file_, format="ascii")

# Connect
conn = splusdata.connect('Luis', 'plutarco*80')

for tab in data:
    if tab["FWHM"] != 0 and tab["FWHM"] < 10:
        rad = int(100*tab["FWHM"])
    else:
        rad = 100
    ra = tab["RA"]
    dec = tab["DEC"]
    Name = tab["ID"].split("R3.")[-1].replace(".", "-")
     
    # Getting the colored imge
    try:
        img = conn.twelve_band_img(ra, dec, radius=rad, noise=0.15, saturation=0.15)
        # Getting the Fits image in the r-band
        hdu = conn.get_cut(ra, dec, rad, 'R')
    except (KeyError, ValueError):
        img = conn.twelve_band_img(195.3910993863848,-13.7146155838678, radius=100, noise=0.15, saturation=0.15)
        hdu = conn.get_cut(195.3910993863848, -13.7146155838678, 100, 'R')
        
    # Getting the Fits image in the r-band
    #hdu = conn.get_cut(ra, dec, rad, 'R')

    # Save the image, note that the output image in compress
    hdu.writeto('{}_{}_r.fz'.format(Name, rad), overwrite=True) # write to fits

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
    hdufits = fz2fits('{}_{}_r.fz'.format(Name, rad))

    # Read the FITS file
    hdul = fits.open('{}_{}_r.fits'.format(Name, rad))[0]
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

    plt.savefig('{}_{}_r.fits'.format(Name,  rad).replace("_r.fits", ".pdf"))
    f.clear()
    plt.close(f)

    # Deleting files
    os.remove('{}_{}_r.fz'.format(Name,  rad))
    print("FZ File Removed!")

    os.remove('{}_{}_r.fits'.format(Name,  rad))
    print("FITS File Removed!")
