from astropy.table import Table
import numpy as np
import argparse
import os

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
file2 = cmd_args.table2 + ".dat"

tab = Table.read(file1, format="ascii.ecsv")

# Table with the IDs
tab_id = Table.read(file2, format="ascii")

# Making mask and applying
id1 = tab["ID"]
id2 = tab_id["ID"]
mask = np.array([not source in id2 for source in id1])

# Save the final file
asciifile = file1.replace(".ecsv", "-visualCleaning.ecsv")
tab[mask].write(asciifile, format="ascii.ecsv", overwrite=True)
