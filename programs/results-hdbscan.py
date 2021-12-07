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
import hdbscan
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.cluster.hierarchy as shc
sns.set_color_codes()
ROOT_PATH = Path("../paper/Figs")

# Read the file
table = Table.read("../iDR3_n4/Final-list-emitters-allparam-unique.ecsv", format="ascii.ecsv")

# Colors
m = (table["e_g_PStotal"] <= 0.2) & (table["e_i_PStotal"] <= 0.2) & (table["e_z_PStotal"] <= 0.2)
m1 =  (table["e_u_PStotal"] <= 0.2) &(table["e_g_PStotal"] <= 0.2) & (table["e_i_PStotal"] <= 0.2) 
zg = table['z_PStotal'][m] - table['g_PStotal'][m]
gr = table['g_PStotal'][m] - table['r_PStotal'][m]
ri = table['r_PStotal'][m] - table['i_PStotal'][m]
ug = table['u_PStotal'][m1] - table['g_PStotal'][m1]
gr_ = table['g_PStotal'][m1] - table['r_PStotal'][m1]

# Create an array
X = np.array(list(zip(zg, gr)))
print("Shape:", X.shape)
# Standarized the data
X_std = StandardScaler().fit_transform(X)

# Applying HDBSCAN
clusterer = hdbscan.HDBSCAN(min_samples=40, min_cluster_size=80, prediction_data=True).fit(X_std) # 40 60
labels_h = clusterer.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels_h)) - (1 if -1 in labels_h else 0)
n_cluster0 = list(labels_h).count(0)
n_cluster1 = list(labels_h).count(1)
n_cluster2 = list(labels_h).count(2)
n_noise_ = list(labels_h).count(-1)

# Print parameters
print('##########################################################')
print('Estimated number of clusters: %d' % n_clusters_)
print('Estimated number of cluster points 0: %d' % n_cluster0)
print('Estimated number of cluster points 1: %d' % n_cluster1)
print('Estimated number of cluster points 2: %d' % n_cluster2)
print('Estimated number of noise points: %d' % n_noise_)
print('##########################################################')

# Getting the probabilities
prob = clusterer.probabilities_

# Add label to the table and making the colors
table_= table[m]

table_["Label"] = labels_h

mask0 = table_["Label"] == -1
mask1 = table_["Label"] == 0
mask2 = table_["Label"] == 1

# Making the colors
zg_0 = table_['z_PStotal'][mask0] - table_['g_PStotal'][mask0]
gr_0 = table_['g_PStotal'][mask0] - table_['r_PStotal'][mask0]
zg_1 = table_['z_PStotal'][mask1] - table_['g_PStotal'][mask1]
gr_1 = table_['g_PStotal'][mask1] - table_['r_PStotal'][mask1]
zg_2 = table_['z_PStotal'][mask2] - table_['g_PStotal'][mask2]
gr_2 = table_['g_PStotal'][mask2] - table_['r_PStotal'][mask2]

# Soft clustering
soft_clusters = hdbscan.all_points_membership_vectors(clusterer)

table_["P(red)"] = soft_clusters[:,0]
table_["P(Blue)"] = soft_clusters[:,1]

#Save the table
asciifile = "../iDR3_n4/Good-LD-Halpha-DR3_noFlag_merge-takeoutbad-Final-hdbscan.ecsv" 
table_.write(asciifile, format="ascii.ecsv")

# Equation constructed form synthetic phometry
# Limiting the blue and red region
x_new = np.linspace(-15.0, 1000, 200)
y = 0.47*x_new + 1.5

#############################################################
#Plot the results  ##########################################
#############################################################

#Build the cluster hierarchy

fig, ax = plt.subplots(figsize=(10, 7))
clusterer.condensed_tree_.plot(select_clusters=True,
                               selection_palette=sns.color_palette(), colorbar=True)

fig.savefig(ROOT_PATH / "cluster-hierarchy-hdbscan.pdf")
plt.clf()
##########################################################
fig, ax1 = plt.subplots(figsize=(12, 12))

ax1.fill_between(x_new, y, -100, color="k", alpha=0.1)
ax1.plot(x_new, y, c="k", zorder=11, lw=0.5)

plt.tick_params(axis='x', labelsize=25) 
plt.tick_params(axis='y', labelsize=25)

plt.xlabel(r'$z - g$', fontsize= 25)
plt.ylabel(r'$g - r$', fontsize= 25)

ax1.scatter(
        zg_0,
        gr_0,
        marker="o",
        c=sns.xkcd_rgb["grey"],
        label="Outliers",
        edgecolors="w", alpha=0.7, zorder=3
    )

ax1.scatter(
        zg_1,
        gr_1,
        marker="o",
        c=sns.xkcd_rgb["dark pink"],
        label="Red",
        edgecolors="w", zorder=3
    )

ax1.scatter(
        zg_2,
        gr_2,
        marker="o",
        c=sns.xkcd_rgb["cerulean"],
        label="Blue",
        edgecolors="w", zorder=4
    )

sns.kdeplot(
    zg_2,
    gr_2,
    ax=ax1,
    norm=PowerNorm(0.5), zorder=5,
        cmap="Blues",
 )

sns.kdeplot(
    zg_1,
    gr_1,
    ax=ax1,
    norm=PowerNorm(0.5), zorder=3,
        cmap="Reds",
 )

ax1.legend(ncol=1, fontsize=20.0, title="Group", title_fontsize=30)
ax1.set(xlim=[-6.8, 2.5], ylim=[-3., 5.])#, xscale="log", yscale="log")
ax1.text(0.05, 1.11, "HDBSCAN", fontsize=20,
                                 bbox=dict(facecolor='gray', alpha=0.2),
                                                       transform=ax.transAxes)
ax1.set_aspect("equal")
#ax.set(xlabel=r"$z - g$", ylabel=r"$g - r$")
fig.savefig(ROOT_PATH / "blued-red-hdbscan.pdf")
plt.clf()

#################################
#Soft clusters   ################
#################################
soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
color_palette = sns.color_palette('Paired', 12)
cluster_colors = [color_palette[np.argmax(x)]
                  for x in soft_clusters]

# Mask to high probabilites to belong
mask_blue = soft_clusters[:,1] > soft_clusters[:,0]
mask_red = soft_clusters[:,0] > soft_clusters[:,1]


fig, ax2 = plt.subplots(figsize=(12, 12))
plt.tick_params(axis='x', labelsize=25) 
plt.tick_params(axis='y', labelsize=25)

plt.xlabel(r'$z - g$', fontsize= 25)
plt.ylabel(r'$g - r$', fontsize= 25)
ax2.set(xlim=[-6.8, 2.5], ylim=[-3., 5.])#, xscale="log", yscale="log")
ax2.fill_between(x_new, y, -100, color="k", alpha=0.1)
ax2.plot(x_new, y, c="k", zorder=11, lw=0.5)
#ax1.scatter(zg, gr, s=50, linewidth=0.2, c=cluster_colors, edgecolors="w", alpha=0.25)

ax2.scatter(
        zg[mask_red], gr[mask_red],
        marker="o",
        c=sns.xkcd_rgb["dark pink"],
        label="Red",
        edgecolors="w", zorder=4
    )

sns.kdeplot(
    zg[mask_red], gr[mask_red],
    ax=ax2,
    norm=PowerNorm(0.5), zorder=5,
        cmap="Reds",
 )

ax2.scatter(
        zg[mask_blue], gr[mask_blue],
        marker="o",
        c=sns.xkcd_rgb["cerulean"],
        label="Blue",
        edgecolors="w", zorder=3
    )

sns.kdeplot(
    zg[mask_blue], gr[mask_blue],
    ax=ax2,
    norm=PowerNorm(0.5), zorder=3,
        cmap="Blues",
 )

ax2.legend(ncol=1, fontsize=20.0, title="Group", title_fontsize=30)
ax2.text(0.05, 1.11, "Soft Clustering for HDBSCAN", fontsize=20,
                                 bbox=dict(facecolor='gray', alpha=0.2),
                                                       transform=ax.transAxes)
ax2.set_aspect("equal")
#ax.set(xlabel=r"$z - g$", ylabel=r"$g - r$")
fig.savefig(ROOT_PATH / "blue-red-hdbscan-soft-alternative.pdf")
plt.clf()

########################################
#Dendrograms for Hierarchical Clustering

fig, ax3 = plt.subplots(figsize=(10, 7))
#plt.figure(figsize=(10, 7))
#plt.title("Customer Dendograms")
plt.xlabel('sample index',  fontsize= 25)
plt.ylabel('distance', fontsize= 25)
plt.tick_params(axis='y', labelsize=25)
dend = shc.dendrogram(shc.linkage(X_std, method='ward'),
                      truncate_mode='lastp',
                      p=12,  # show only the last p merged clusters
                      leaf_rotation=45.,
                      leaf_font_size=18.,
                      show_contracted=True, )
#dend = shc.dendrogram(shc.linkage(X, method='ward'))

# dendrogram(
#     X,
#     truncate_mode='lastp',  # show only the last p merged clusters
#     p=12,  # show only the last p merged clusters
#     leaf_rotation=90.,
#     leaf_font_size=12.,
#     show_contracted=True,  # to get a distribution impression in truncated branches
# )
plt.tight_layout()
fig.savefig(ROOT_PATH / "Customer-Dendrograms.pdf")
plt.clf()
