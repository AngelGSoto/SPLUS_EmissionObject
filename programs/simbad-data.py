'''
Read table with simbad infromation
'''
from __future__ import print_function
import numpy as np
from sklearn import metrics
from scipy.optimize import curve_fit
import pandas as pd
from astropy.table import Table
import seaborn as sns
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sklearn.metrics import mean_squared_error
from astropy.modeling import models, fitting
import argparse
import sys
import os
from pathlib import Path
ROOT_PATH = Path("../paper/Figs")

df = pd.read_csv("simbad.csv")
print(df.columns)

# MASKs
m1 = df['main_type'] == 'AGN_Candidate'
m2 = df['main_type'] == 'AGN'
m3 = df['main_type'] == 'EmG'
m4 = df['main_type'] == 'GinGroup'
m5 = df['main_type'] == 'Galaxy'
m6 = df['main_type'] == 'Candidate_CV*'
m7 = df['main_type'] == 'RRLyr'
m8 = df['main_type'] == 'SN'
m9 = df['main_type'] == 'HII'
m11 = df['main_type'] == 'CataclyV*'
m12 = df['main_type'] == 'FIR'
m13 = df['main_type'] == 'GinCl'
m14 = df['main_type'] == 'HII_G'
m15 = df['main_type'] == 'Seyfert_1'
m17 = df['main_type'] == 'Star'
m18 = df['main_type'] == 'PartofG'
m19 = df['main_type'] == 'RadioG'
m20 = df['main_type'] == 'IG'
m21 = df['main_type'] == 'QSO'
m22 = df['main_type'] == 'EB*'
m23 = df['main_type'] == 'Radio'
m24 = df['main_type'] == 'Seyfert_2'
m26 = df['main_type'] == 'X'
m27 = df['main_type'] == 'MolCld'
m28 = df['main_type'] == 'Cl*'
m29 = df['main_type'] == 'HMXB'
m30 = df['main_type'] == 'GinPair'
m31 = df['main_type'] == 'LSB_G'
m32 = df['main_type'] == 'WD*'
m33 = df['main_type'] == 'Candidate_RRLyr'
m34 = df['main_type'] == 'PN'
m35 = df['main_type'] == 'Blue'
m36 = df['main_type'] == 'EmObj'
m37 = df['main_type'] == 'BlueSG*'
m38 = df['main_type'] == 'StarburstG'
m39 = df['main_type'] == 'low-mass*'
m40 = df['main_type'] == 'BlueCompG'
m41 = df['main_type'] == 'UV'
m42 = df['main_type'] == 'Candidate_WD*'
m43 = df['main_type'] == 'MIR'
m44 = df['main_type'] == 'Radio(cm)'
m45 = df['main_type'] == 'Candidate_SN*'
m46 = df['main_type'] == 'QSO_Candidate'
m47 = df['main_type'] == 'BLLac'
m48 = df['main_type'] == 'PM*'
m49 = df['main_type'] == 'Possible_lensImage'
m51 = df['main_type'] == 'Nova'
m52 = df['main_type'] == 'BClG'
m53 = df['main_type'] == 'GlCl'

# Making the tables with individual object classes
df_agn = pd.concat(df[m1], df[m2])
df_EmG = pd.concat([df[m3], df[m14], df[m38], df[m40]])
GinGroup = df[m4]
df_pn = df[m34] 
df_gal = df[m5]
df_qso = pd.concat([df[m21], df[m46]])
df_cv = pd.concat([df[m6], df[m11]])
df_hii = df[m9]
df_star = df[m17]
df_sn = pd.concat(df[m8], df[m45])
df_rrly = pd.concat(df[m7], df[m33])
df_fir = df[m12]
df_GinCl = df[m13]
df_Seyfert_1 = df[m15]
df_Seyfert_2 = df[m24]
df_PartofG = df[m18]
df_RadioG = df[m19]
df_IG = df[m20]
df_EB = df[m22]
df_Radio = pd.concat(df[m23], df[m44])
df_X = df[m26]
df_MolCld = df[m27]
df_cl = df[m28]
df_HMXB = df[29]
df_GinPair = df[30]
df_LSB_G = df[m31]
df_WD = pd.concat(df[m32], df[m42])
df_Blue = df[m35]
df_EmObj = df[m36]
df_BlueSG = df[m37]
df_low_mass = df[m39]
df_uv = df[m41]
df_mier = df[m43]
df_BLLac = df[m47]
df_pm = df[m48]
df_Possible_lensImage = df[m49]
df_nova = df[m51]
df_bclg = df[m52]
df_glcl = df[m53]

# Countaining the sources of each class 
list_numbers = ["H II regions", "PN", "CV", "SN", "Nova", "BL Lac",
                "Variable Star of RR Lyr type",
                "Star", "WD", "Cluster of Stars", "HMXB", "Far-Infrared source",
                "Eclipsing binary", "Blue object", "Emission Object", "Blue supergiant star",
                "Low-mass star", "UV-emission source", "MIER", "Possible lens Image"
                "Galaxy", "Galaxy in Pair of Galaxies",
                "Emission line galaxies",
                "QSO", "AGN", "Part of a Galaxy", "X-ray source", "Molecular Cloud",
                "Galaxy in Group of Galaxies", "Radio-source", "Interacting Galaxies",
                "Low Surface Brightness Galaxy", "Radio Galaxy", "Galaxy in Cluster of Galaxies"
                "High proper-motion Star", "Seyfert 1", "Seyfert 2",
                "Brightest galaxy in a Cluster", "Globular Cluster"]

# Definition to make the colors
def colour(tab, f1, f2, f3, f4):
    xcolour = tab[f1] - tab[f2]
    ycolour = tab[f3] - tab[f4]
    return xcolour, ycolour

# Colors
cx_pn, cy_pn = colour(df_pn, "Z_PStotal", "G_PStotal", "G_PStotal", "R_PStotal")
cx_gal, cy_gal = colour(df_gal, "Z_PStotal", "G_PStotal", "G_PStotal", "R_PStotal")
cx_EmG, cy_EmG = colour(df_EmG, "Z_PStotal", "G_PStotal", "G_PStotal", "R_PStotal")
cx_qso, cy_qso = colour(df_qso, "Z_PStotal", "G_PStotal", "G_PStotal", "R_PStotal")
cx_cv, cy_cv = colour(df_cv, "Z_PStotal", "G_PStotal", "G_PStotal", "R_PStotal")
cx_hii, cy_hii = colour(df_hii, "Z_PStotal", "G_PStotal", "G_PStotal", "R_PStotal")
cx_star, cy_star = colour(df_star, "Z_PStotal", "G_PStotal", "G_PStotal", "R_PStotal")

#PLOT
# Limiting the blue and red region
x_new = np.linspace(-15.0, 1000, 200)
y = 0.45*x_new + 1.55

fig, ax = plt.subplots(figsize=(12, 12))

ax.fill_between(x_new, y, -100, color="k", alpha=0.1)
ax.plot(x_new, y, c="k", zorder=11, lw=0.5)

plt.tick_params(axis='x', labelsize=25) 
plt.tick_params(axis='y', labelsize=25)

plt.xlabel(r'$z - g$', fontsize= 25)
plt.ylabel(r'$g - r$', fontsize= 25)

ax.scatter(
        cx_pn, cy_pn,
        marker="o",
        c=sns.xkcd_rgb["cerulean"],
        label="PN",
        edgecolors="w", alpha=0.7, zorder=4
    )

ax.scatter(
        cx_gal, cy_gal,
        marker="o",
        c=sns.xkcd_rgb["dark pink"],
        label="Galaxy",
        edgecolors="w", zorder=3
    )

ax.scatter(
        cx_EmG, cy_EmG,
        marker="o",
        c=sns.xkcd_rgb["bright blue"],
        label="EmG",
        edgecolors="w", zorder=3
    )

ax.scatter(
        cx_qso, cy_qso,
        marker="o",
        c=sns.xkcd_rgb["green"],
        label="QSO",
        edgecolors="w", zorder=3
    )

ax.scatter(
        cx_cv, cy_cv,
        marker="o",
        c=sns.xkcd_rgb["periwinkle"],
        label="CV",
        edgecolors="w", zorder=5
    )

ax.scatter(
        cx_hii, cy_hii,
        marker="o",
        c=sns.xkcd_rgb["pale yellow"],
        label="HII Region",
        edgecolors="w", zorder=4
    )

ax.scatter(
        cx_star, cy_star,
        marker="o",
        c=sns.xkcd_rgb["army green"],
        label="Star",
        edgecolors="w", zorder=6
    )


ax.legend(ncol=1, fontsize=20.0, title_fontsize=30)
ax.set(xlim=[-6.8, 2.5], ylim=[-3., 5.])#, xscale="log", yscale="log")
ax.set_aspect("equal")
#ax.set(xlabel=r"$z - g$", ylabel=r"$g - r$")

fig.savefig(ROOT_PATH / "colour-digram-simbadObj.pdf")
