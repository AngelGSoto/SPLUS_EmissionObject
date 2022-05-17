import splusdata
import getpass
import pandas as pd
from astropy.table import Table, vstack
import numpy as np
from pathlib import Path
  
ROOT = Path("Final-list/") 
  
# Connecting with SPLUS database
username = str(input("Login: "))
password = getpass.getpass("Password: ")

conn = splusdata.connect(username, password)
  
tab = Table.read(ROOT/ "Final-list-emitters-allparam-unique.ecsv", format="ascii.ecsv")

print("Number of objects:", len(tab))

Query = f"""SELECT detection.Field, detection.ID, detection.RA, detection.DEC,
		  detection.FWHM, detection.ISOarea, detection.KRON_RADIUS, 
		  detection.MU_MAX, detection.nDet_PStotal, detection.PhotoFlagDet, 
                  detection.CLASS_STAR, detection.u_PStotal, detection.J0378_PStotal, 
                  detection.J0395_PStotal, detection.J0410_PStotal, detection.J0430_PStotal, 
                  detection.g_PStotal, detection.J0515_PStotal, detection.r_PStotal, 
                  detection.J0660_PStotal, detection.i_PStotal, 
		  detection.J0861_PStotal, detection.z_PStotal, detection.e_u_PStotal, 
                  detection.e_J0378_PStotal, detection.e_J0395_PStotal, 
                  detection.e_J0410_PStotal, detection.e_J0430_PStotal, 
		  detection.e_g_PStotal, detection.e_J0515_PStotal, detection.e_R_PStotal, 
                  detection.e_J0660_PStotal, detection.e_i_PStotal, 
                  detection.e_J0861_PStotal, detection.e_z_PStotal 
		  FROM TAP_UPLOAD.upload as tap 
                  LEFT OUTER JOIN "idr3"."all_idr3" as detection ON tap.ID = detection.ID""" 
  
# coverting in pandas table
df = (tab.to_pandas())
n = int(len(df) / 1000.) + 1
  
df_ = [] # list
j = 0 # counter
d = {} # empty
for i in range(n):
    j += 1  
    df_.append(df.iloc[1000*i:1000*j])
  
# Applying query
merged_table_list = []
for a in range(n):
    results = conn.query(Query, df_[a])
    merged_table_list.append(results)
  
# Merging all result astropy tables 
merged_table = vstack(merged_table_list)

print("Number objects with match:", len(merged_table))
  
asciifile = "Final-list-emitters-allparam-unique-class.ecsv"
merged_table.write(ROOT / asciifile, format="ascii.ecsv", overwrite=True)

