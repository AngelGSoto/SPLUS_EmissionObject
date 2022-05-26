import splusdata
import getpass
import pandas as pd
from astropy.table import Table, vstack
import numpy as np
import argparse
from pathlib import Path
  
ROOT = Path("Final-list/")
  
# Connecting with SPLUS database
username = str(input("Login: "))
password = getpass.getpass("Password: ")

conn = splusdata.connect(username, password)
  
parser = argparse.ArgumentParser(
    description="""Get data form splus.cloudy""")

parser.add_argument("rmagi", type=float,
                    default="12",
                    help="r-mag; initial")

parser.add_argument("rmagf", type=float,
                    default="16",
                    help="r-mag; final")

cmd_args = parser.parse_args()
ri = cmd_args.rmagi
rf = cmd_args.rmagf


rii = np.linspace(ri, rf, 10)

riii = []
rfff = []
for i, j in enumerate(rii, start = 1):
    if j < rf:
         riii.append(j)
    if j > ri:
        rfff.append(j)

merged_table_list = []
for k, t in zip(riii, rfff):
    Query = f"""SELECT Field, ID, RA, DEC, FWHM, ISOarea, KRON_RADIUS, MU_MAX, nDet_PStotal, PhotoFlagDet, 
                  CLASS_STAR, u_PStotal, J0378_PStotal, 
                  J0395_PStotal,  J0410_PStotal, J0430_PStotal, 
                  g_PStotal, J0515_PStotal, r_PStotal, 
                  J0660_PStotal, i_PStotal, 
        	  J0861_PStotal, z_PStotal, e_u_PStotal, 
                  e_J0378_PStotal, e_J0395_PStotal, 
                  e_J0410_PStotal, e_J0430_PStotal, 
        	  e_g_PStotal, e_J0515_PStotal, e_R_PStotal, 
                  e_J0660_PStotal, e_i_PStotal, 
                  e_J0861_PStotal, e_z_PStotal 
        	  FROM "idr3"."all_idr3"
                  WHERE CLASS_STAR >= 0.5 AND r_PStotal >= """ + str(k) + """ AND r_PStotal < """  + str(t)  +  """ AND e_r_PStotal <= 0.2 AND e_J0660_PStotal <= 0.2 
                  AND e_i_PStotal <= 0.2"""
    # Executing the query
    results = conn.query(Query)
    print("Number of soures in the interval:", str(k), str(t), "=>", len(results))
    merged_table_list.append(results)
    

# Merging all result astropy tables 
merged_table = vstack(merged_table_list)

# coverting in pandas table
df = (merged_table.to_pandas())

df.to_csv("DR3-SPLUS-PStotal-STAR-{}r{}.csv".format(int(ri), int(rf)), index = False)
