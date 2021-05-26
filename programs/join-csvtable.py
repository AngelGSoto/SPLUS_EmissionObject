import pandas as pd
from astropy.table import Table
import glob
from pathlib import Path

ROOT_PATH = Path("..") 
pattern = "result??.csv"
file_list = glob.glob(pattern)

dfs = []
for file_name in file_list:
    df = pd.read_csv(file_name)
    dfs.append(df)

# save the dataframe
df_final = pd.concat(dfs)
df_final.to_csv("table-version-incom.csv")  
