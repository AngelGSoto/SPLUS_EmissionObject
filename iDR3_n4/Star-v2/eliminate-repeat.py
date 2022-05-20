import pandas as pd
from astropy.table import Table, vstack
import numpy as np
from pathlib import Path

ROOT = Path("Final-list/") 

tab = Table.read("Halpha-DR3_PStotal-STAR_merge-clean-duplicate.ecsv", format="ascii.ecsv")

# Eliminate duplicates
sources, idu = [], []
for i, source in  enumerate(tab['ID']):
    if source in sources:
        idu.append(i)
    sources.append(source)
if idu:
    tab.remove_rows(idu)

print("The fineal number is:", len(tab))
newtabfile_ = "Halpha-DR3_PStotal-STAR_merge-clean-duplicate-unique.ecsv"
tab.write(newtabfile_, format="ascii.ecsv")

# pandas
df_tab = tab.to_pandas()
df_newtabfile_ = "Halpha-DR3_PStotal-STAR_merge-clean-duplicate-unique.csv"
df_tab.to_csv(df_newtabfile_, index = False)


