'''
This is a simply script to make table with the SDSS format for cross-match.
'''
from astropy.table import Table, vstack
import numpy as np
import argparse
import sys
import os

parser = argparse.ArgumentParser(
    description="""Make a table from the S-PLUS catalogs """)

parser.add_argument("source", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")

cmd_args = parser.parse_args()
file_ = cmd_args.source + ".ecsv"

tab = Table.read(file_, format="ascii.ecsv")

n = len(tab)

sep = np.linspace(2.0/60., 2.0/60., num=n)
ra = tab["RA"]
dec = tab["DEC"]
#table = Table([ra, dec, sep], names=('ra', 'dec', 'sep'), meta={'name': 'first table'})

#n_new = len(table)
#n_ = n_new/900.
count = 0
chunksize = 900

number = np.linspace(1, 16, num=16, dtype = int)

ra_ = [[] for _ in range(len(number))]
dec_ = [[] for _ in range(len(number))]
sep_ = [[] for _ in range(len(number))]

for i in range(n):
    if count < chunksize:
        count += 1
        ra_[0].append(tab["RA"][i])
        dec_[0].append(tab["DEC"][i])
        sep_[0].append(sep[i])
    elif count < 2*chunksize:
        count += 1
        ra_[1].append(tab["RA"][i])
        dec_[1].append(tab["DEC"][i])
        sep_[1].append(sep[i])
    elif count < 3*chunksize:
        count += 1
        ra_[2].append(tab["RA"][i])
        dec_[2].append(tab["DEC"][i])
        sep_[2].append(sep[i])
    elif count < 4*chunksize:
        count += 1
        ra_[3].append(tab["RA"][i])
        dec_[3].append(tab["DEC"][i])
        sep_[3].append(sep[i])
    elif count < 5*chunksize:
        count += 1
        ra_[4].append(tab["RA"][i])
        dec_[4].append(tab["DEC"][i])
        sep_[4].append(sep[i])
    elif count < 6*chunksize:
        count += 1
        ra_[5].append(tab["RA"][i])
        dec_[5].append(tab["DEC"][i])
        sep_[5].append(sep[i])
    elif count < 7*chunksize:
        count += 1
        ra_[6].append(tab["RA"][i])
        dec_[6].append(tab["DEC"][i])
        sep_[6].append(sep[i])
    elif count < 8*chunksize:
        count += 1
        ra_[7].append(tab["RA"][i])
        dec_[7].append(tab["DEC"][i])
        sep_[7].append(sep[i])
    elif count < 9*chunksize:
        count += 1
        ra_[8].append(tab["RA"][i])
        dec_[8].append(tab["DEC"][i])
        sep_[8].append(sep[i])
    elif count < 10*chunksize:
        count += 1
        ra_[9].append(tab["RA"][i])
        dec_[9].append(tab["DEC"][i])
        sep_[9].append(sep[i])
    elif count < 11*chunksize:
        count += 1
        ra_[10].append(tab["RA"][i])
        dec_[10].append(tab["DEC"][i])
        sep_[10].append(sep[i])
    elif count < 12*chunksize:
        count += 1
        ra_[11].append(tab["RA"][i])
        dec_[11].append(tab["DEC"][i])
        sep_[11].append(sep[i])
    elif count < 13*chunksize:
        count += 1
        ra_[12].append(tab["RA"][i])
        dec_[12].append(tab["DEC"][i])
        sep_[12].append(sep[i])
    elif count < 14*chunksize:
        count += 1
        ra_[13].append(tab["RA"][i])
        dec_[13].append(tab["DEC"][i])
        sep_[13].append(sep[i])
    elif count < 15*chunksize:
        count += 1
        ra_[14].append(tab["RA"][i])
        dec_[14].append(tab["DEC"][i])
        sep_[14].append(sep[i])
    else:
        ra_[15].append(tab["RA"][i])
        dec_[15].append(tab["DEC"][i])
        sep_[15].append(sep[i])

t = [Table() for _ in range(len(number))]
for j in range(len(number)):
    t[j]["ra"] = ra_[j]
    t[j]["dec"] = dec_[j]
    t[j]["sep"] = sep_[j]

#Save the file
datadir = "Spectra-Sdss/"
for m in range(len(number)):
    asciifile = file_.replace(".ecsv", "-"+str(number[m])+".dat")
    t[m].write(os.path.join(datadir, asciifile), format="ascii", delimiter=',', overwrite=True)
    
