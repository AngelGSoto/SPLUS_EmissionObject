from astropy.table import Table
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(
    description="""Firts table from the S-PLUS catalogs """)

parser.add_argument("source", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")

cmd_args = parser.parse_args()
file_ = cmd_args.source + ".dat"

# All table with emission line objects 
#datadir = "../"
#tab = Table.read(os.path.join(datadir, file_), format="ascii.ecsv")
tab = Table.read("Halpha-DR3_noFlag_merge.ecsv", format="ascii.ecsv")

# Table with the IDs
tab_id = Table.read(file_, format="ascii")

# Making mask and applying
id1 = tab["ID"]
id2 = tab_id["ID"]
mask = np.array([source in id2 for source in id1])

# Save the final file
asciifile = file_.replace(".dat", "-trainig.ecsv")
tab[mask].write(asciifile, format="ascii.ecsv", overwrite=True)
