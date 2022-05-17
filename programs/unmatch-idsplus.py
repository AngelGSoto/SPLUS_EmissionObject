'''
Scrit to do unmatch using ID 
'''
from astropy.table import Table
import numpy as np
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(
    description="""Tables from the S-PLUS catalogs """)

parser.add_argument("table1", type=str,
                    default=" teste-program",
                    help="Table 1, taken the prefix ")
parser.add_argument("table2", type=str,
                    default=" teste-program",
                    help="Table 2, taken the prefix ")

cmd_args = parser.parse_args()
file1 = cmd_args.table1 + ".ecsv"

cmd_args = parser.parse_args()
file2 = cmd_args.table2 + ".dat"

# Table 1
try:
    tab1 = Table.read(file1, format="ascii.ecsv")
except FileNotFoundError:
    file1 = cmd_args.table1 + ".csv"
    tab1 = pd.read_csv(file1)

print("#############################################")
print("Number of objects of table 1:", len(tab1))

# Table 2

tab2 = Table.read(file2, format="ascii")
print("#############################################")    
print("Number of objects of table 2:", len(tab2))

# Making mask and applying
id1 = tab1["ID"]
id2 = tab2["ID"]
mask = np.array([not source in id2 for source in id1])

print("Remain objects:", len(tab1[mask]))
print("#############################################")

# Save the final file
asciifile = file1.replace(".ecsv", "-revised.ecsv")
tab1[mask].write(asciifile, format="ascii.ecsv", overwrite=True)

# Save dataframe
asciifile_df = file1.replace(".ecsv", "-revised.csv")
df = (tab1[mask].to_pandas())
df.to_csv(asciifile_df, index=False)
