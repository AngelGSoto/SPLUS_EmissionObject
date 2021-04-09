'''
This is a simply script to make table with the format of Lamost for cross-match.
'''
from astropy.table import Table
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(
    description="""Make a table from the S-PLUS catalogs """)

parser.add_argument("source", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")


cmd_args = parser.parse_args()
file_ = cmd_args.source + ".ecsv"

datadir = "../"
tab = Table.read(os.path.join(datadir, file_), format="ascii.ecsv")

n = len(tab["RA"])
sep = np.linspace(2.0, 2.0, num=n)
ra = tab["RA"]
dec = tab["DEC"]
table = Table([ra, dec, sep], names=('ra', 'dec', 'radius'), meta={'name': 'first table'})

# Save the file
asciifile = file_.replace(".ecsv", 
                  "-coorLamost.dat")
table.write(asciifile, format="ascii.commented_header", delimiter=',')

