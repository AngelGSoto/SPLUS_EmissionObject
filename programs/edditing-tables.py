from astropy.table import Table, vstack
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

# Table
datadir = "../3filter_noflat/Claudia-objects/"
tab = Table.read(file_, format="ascii")

tab.remove_columns(['RA_1', 'DEC_1', 'P(GoodPho)',  'P(BadPho)', 'GroupID', 'GroupSize', 'Separation'])
tab.rename_column('RA_2', 'RA')
tab.rename_column('DEC_2', 'DEC')

df = tab.to_pandas()
dffile = file_.replace(".dat", "-hydra.csv")
df.to_csv(dffile)

#ASCII
asciifile = file_.replace(".dat", "-hydra.ecsv")
tab.write(asciifile, format="ascii.ecsv")
