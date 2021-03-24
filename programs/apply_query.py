'''
Script to download data from my database using query
''' 
import sqlite_splus
import pandas as pd
from astropy.table import Table

# Names of the catalogs in my database
catalogs = ["HYDRA", "STRIPE82", "MAIN3_1", "MAIN3_2", "MAIN3_3", "MAIN3_4", "MAIN3_5", "MAIN3_6", "MAIN3_7", "MAIN3_8" ]

# Number of catalogs in the database
n_cat = len(catalogs)

#df_final = pd.DataFrame([])
def run_query(catalog):
    df = sqlite_splus.dr3qury("Field, ID, RA, DEC, FWHM, ISOarea, KRON_RADIUS, nDet_magPStotal, PhotoFlagDet, U_PStotal, F378_PStotal, F395_PStotal, F410_PStotal, F430_PStotal, G_PStotal, F515_PStotal, R_PStotal, F660_PStotal, I_PStotal, F861_PStotal, Z_PStotal, e_U_PStotal, e_F378_PStotal, e_F395_PStotal, e_F410_PStotal, e_F430_PStotal, e_G_PStotal, e_F515_PStotal, e_R_PStotal, e_F660_PStotal, e_I_PStotal, e_F861_PStotal, e_Z_PStotal", catalog, "PhotoFlag_U <= 3.0 AND PhotoFlag_F378 <= 3.0 AND PhotoFlag_F395 <= 3.0 AND PhotoFlag_F410 <= 3.0 AND PhotoFlag_F430 <= 3.0 AND PhotoFlag_G <= 3.0 AND PhotoFlag_F515 <= 3.0 AND PhotoFlag_R <= 3.0 AND PhotoFlag_F660 <= 3.0 AND PhotoFlag_I <= 3.0 AND PhotoFlag_F861 <= 3.0 AND  PhotoFlag_Z <= 3.0 AND R_PStotal >= 20.0 AND R_PStotal <= 21.0 AND e_U_PStotal <= 0.2 AND e_F378_PStotal <= 0.2 AND e_F395_PStotal <= 0.2 AND e_F410_PStotal <= 0.2 AND e_F430_PStotal <= 0.2 AND e_G_PStotal <= 0.2 AND e_F515_PStotal <= 0.2 AND e_R_PStotal <= 0.2 AND e_F660_PStotal <= 0.2 AND e_I_PStotal <= 0.2 AND e_F861_PStotal <= 0.2 AND e_Z_PStotal <= 0.2")
    return df

# Inplemente the query for each catalog in my databaes (IDR3.db)
df_hydra = run_query("HYDRA")
df_stripe = run_query("STRIPE82")
df_M31 = run_query("MAIN3_1")
df_M32 = run_query("MAIN3_2")
df_M33 = run_query("MAIN3_3")
df_M34 = run_query("MAIN3_4")
df_M35 = run_query("MAIN3_5")
df_M36 = run_query("MAIN3_6")
df_M37 = run_query("MAIN3_7")
df_M38 = run_query("MAIN3_8")

dfs = [df_hydra, df_stripe, df_M31, df_M32, df_M33, df_M34, df_M35, df_M36, df_M37, df_M38]

# save the dataframe
df_final = pd.concat(dfs)
df_final.to_csv("DR3_errorsall_flagallf_20r21.csv") 
