import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import pandas as pd
import numpy as np
from astropy.table import Table
import seaborn as sns
import argparse
import sys
import os
import glob
import json
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from pathlib import Path
from density_scatter import density_scatter
sns.set_color_codes()
ROOT_PATH = Path("../../paper/Figs")


table = Table.read("table-version-incom-class.ecsv", format="ascii.ecsv")

# Mask
mask_qso = table["CLASS"] == 0
mask_star = table["CLASS"] == 1
mask_galaxy = table["CLASS"] == 2

table_qso = table[mask_qso]
table_star = table[mask_star]
table_galaxy = table[mask_galaxy]

# # Definition for the colors
def colour(table, f1, f2, f3, f4):
    xcolour = table[f1] - table[f2]
    ycolour = table[f3] - table[f4]
    return xcolour, ycolour

# Colors
cx_qso, cy_qso = colour(table_qso, "z_PStotal", "g_PStotal", "g_PStotal", "i_PStotal")
cx_star, cy_star = colour(table_star, "z_PStotal", "g_PStotal", "g_PStotal", "i_PStotal")
cx_galaxy, cy_galaxy = colour(table_galaxy, "z_PStotal", "g_PStotal", "g_PStotal", "i_PStotal")


fig, ax = plt.subplots(figsize=(12, 12))

# ax.fill_between(x_new, y, -100, color="k", alpha=0.1)
# ax.plot(x_new, y, c="k", zorder=11, lw=0.5)

plt.tick_params(axis='x', labelsize=25) 
plt.tick_params(axis='y', labelsize=25)

plt.xlabel(r'$z - g$', fontsize= 25)
plt.ylabel(r'$g - i$', fontsize= 25)

# ax.scatter(
#         cx_qso, cy_qso,
#         marker="o",
#         c=sns.xkcd_rgb["cerulean"],
#         label="QSO",
#         edgecolors="w", alpha=0.7, zorder=3
#     )

# ax.scatter(
#         cx_star, cy_star,
#         marker="o",
#         c=sns.xkcd_rgb["dark pink"],
#         label="Star",
#         edgecolors="w", zorder=3
#     )

ax.scatter(
        cx_galaxy, cy_galaxy,
        marker="o",
        c=sns.xkcd_rgb["pale yellow"],
        label="Galaxy",
        edgecolors="w", zorder=3
    )

ax.legend(ncol=1, fontsize=20.0, title="Group", title_fontsize=30)
ax.set(xlim=[-6.8, 2.5], ylim=[-3., 5.])#, xscale="log", yscale="log")
ax.set_aspect("equal")
#ax.set(xlabel=r"$z - g$", ylabel=r"$g - r$")

fig.savefig(ROOT_PATH / "class.pdf")
