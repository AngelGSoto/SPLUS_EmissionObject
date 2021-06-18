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
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.cluster.hierarchy as shc
sns.set_color_codes()
ROOT_PATH = Path("../paper/Figs")


table_blue = Table.read("Blue0-Good-LD-Halpha-DR3_noFlag_merge-7filter-visualCleaning-Final-takeoutrepeat-Final.ecsv", format="ascii.ecsv")
table_red = Table.read("Red1-Good-LD-Halpha-DR3_noFlag_merge-7filter-visualCleaning-Final-takeoutrepeat-Final.ecsv", format="ascii.ecsv")

# Making the colors
zg_blue = table_blue['Z_PStotal'] - table_blue['G_PStotal']
gr_blue = table_blue['G_PStotal'] - table_blue['R_PStotal']

zg_red = table_red['Z_PStotal'] - table_red['G_PStotal']
gr_red = table_red['G_PStotal'] - table_red['R_PStotal']


# Equation constructed form synthetic phometry
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
        zg_blue,
        gr_blue,
        marker="o",
        c=sns.xkcd_rgb["cerulean"],
        label="Blue",
        edgecolors="w", alpha=0.7, zorder=4
    )

ax.scatter(
        zg_red,
        gr_red,
        marker="o",
        c=sns.xkcd_rgb["dark pink"],
        label="Red",
        edgecolors="w", zorder=3
    )

sns.kdeplot(
    zg_blue,
    gr_blue,
    ax=ax,
    norm=PowerNorm(0.5), zorder=10,
        cmap="Blues",
 )

sns.kdeplot(
    zg_red,
    gr_red,
    ax=ax,
    norm=PowerNorm(0.5), zorder=3,
        cmap="Reds",
 )


ax.legend(ncol=1, fontsize=20.0, title="Group", title_fontsize=30)
ax.set(xlim=[-6.8, 2.5], ylim=[-3., 5.])#, xscale="log", yscale="log")
ax.set_aspect("equal")
#ax.set(xlabel=r"$z - g$", ylabel=r"$g - r$")

fig.savefig(ROOT_PATH / "blued-red-hierarchical.pdf")
