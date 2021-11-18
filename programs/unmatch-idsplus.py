'''
Scrit to do unmatch using ID 
'''
from astropy.table import Table
import numpy as np
import argparse
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
file2 = cmd_args.table2 + ".ecsv"

# Table 1
tab1 = Table.read(file1, format="ascii.ecsv")
print("#############################################")
print("Number of objects of table 1:", len(tab1))

# Table 2
tab2 = Table.read(file2, format="ascii.ecsv")
print("Number of objects of table 2:", len(tab2))

# Making mask and applying
id1 = tab1["ID"]
id2 = tab2["ID"]
mask = np.array([not source in id2 for source in id1])

print("Remain objects:", len(tab1[mask]))
print("#############################################")

# Save the final file
asciifile = file2.replace(".ecsv", "-remain.ecsv")
tab1[mask].write(asciifile, format="ascii.ecsv", overwrite=True)
