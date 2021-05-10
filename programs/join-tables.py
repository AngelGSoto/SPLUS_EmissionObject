from astropy.table import Table, vstack
import numpy as np
import argparse
import os
import gzip
from astropy.io import fits

parser = argparse.ArgumentParser(
    description="""Firts table from the S-PLUS catalogs """)

parser.add_argument("table1", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")
parser.add_argument("table2", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")

cmd_args = parser.parse_args()
file1 = cmd_args.table1 + ".ecsv"

cmd_args = parser.parse_args()
file2 = cmd_args.table2 + ".ecsv"

tab1 = Table.read(file1, format="ascii.ecsv")

# Table with the IDs
tab2 = Table.read(file2, format="ascii.ecsv")


# Merge the tables
table_merge = vstack([tab1, tab2])

# Save the final file
asciifile = file1.replace(".ecsv", "-Final.ecsv")
fitsfile = file1.replace(".ecsv", "-Final.fits")
gzfitsfile = file1.replace(".ecsv", "-Final.fits.gz")
table_merge.write(asciifile, format="ascii.ecsv", overwrite=True)
table_merge.write(fitsfile, format='fits',  overwrite=True)

f = gzip.open(gzfitsfile, 'wb')
h = fits.PrimaryHDU()
table_merge.writeto(f)
