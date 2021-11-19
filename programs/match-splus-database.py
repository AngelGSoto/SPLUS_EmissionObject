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
  
tab = Table.read(ROOT/ "Final-list-emitters-radec.ecsv", format="ascii.ecsv")

print("Number of objects:", len(tab))
  
Query = f"""SELECT detection.Field, detection.ID, detection.RA, detection.DEC,
		  detection.FWHM, detection.ISOarea, detection.KRON_RADIUS, 
		  detection.MU_MAX, detection.nDet_PStotal, detection.PhotoFlagDet,
		  u.U_PStotal, J0378.J0378_PStotal, J0395.J0395_PStotal,
		  J0410.J0410_PStotal, J0430.J0430_PStotal, g.G_PStotal,
		  J0515.J0515_PStotal, r.R_PStotal, J0660.J0660_PStotal, i.I_PStotal, 
		  J0861.J0861_PStotal, z.Z_PStotal, u.e_U_PStotal, J0378.e_J0378_PStotal,
		  J0395.e_J0395_PStotal, J0410.e_J0410_PStotal, J0430.e_J0430_PStotal, 
		  g.e_G_PStotal, J0515.e_J0515_PStotal, r.e_R_PStotal, J0660.e_J0660_PStotal,
		  i.e_I_PStotal, J0861.e_J0861_PStotal, z.e_Z_PStotal 
		  FROM TAP_UPLOAD.upload as tap 
                  LEFT OUTER JOIN idr3.detection_image as detection ON tap.ID= detection.ID
		  LEFT OUTER JOIN idr3.u_band as u ON tap.ID=u.ID 
		  LEFT OUTER JOIN idr3.J0378_band as J0378 ON tap.ID=J0378.ID
		  LEFT OUTER JOIN idr3.J0395_band as J0395 ON tap.ID=J0395.ID
		  LEFT OUTER JOIN idr3.J0410_band as J0410 ON tap.ID=J0410.ID
		  LEFT OUTER JOIN idr3.J0430_band as J0430 ON tap.ID=J0430.ID
		  LEFT OUTER JOIN idr3.g_band as g ON tap.ID=g.ID 
		  LEFT OUTER JOIN idr3.J0515_band as J0515 ON tap.ID=J0515.ID
		  LEFT OUTER JOIN idr3.r_band as r ON tap.ID=r.ID
		  LEFT OUTER JOIN idr3.J0660_band as J0660 ON tap.ID=J0660.ID 
		  LEFT OUTER JOIN idr3.i_band as i ON tap.ID=i.ID
		  LEFT OUTER JOIN idr3.J0861_band as J0861 ON tap.ID=J0861.ID
		  LEFT OUTER JOIN idr3.z_band as z ON tap.ID=z.ID""" 
  
# coverting in pandas table
df = (tab.to_pandas())
n = int(len(df) / 5000.) + 1
  
df_ = [] # list
j = 0 # counter
d = {} # empty
for i in range(n):
    j += 1  
    df_.append(df.iloc[5000*i:5000*j])
  
# Applying query
merged_table_list = []
for a in range(n):
    results = conn.query(Query, df_[a])
    merged_table_list.append(results)
  
# Merging all result astropy tables 
merged_table = vstack(merged_table_list)

print("Number objects with match:", len(merged_table))
  
asciifile = "Final-list-emitters-allparam.ecsv"
merged_table.write(ROOT / asciifile, format="ascii.ecsv", overwrite=True)
