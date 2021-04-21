'''
Estimate the PCs emission lines splus
'''
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table, vstack, hstack
from astropy.coordinates import SkyCoord 
import numpy as np
from pathlib import Path
from astropy.table import Column
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, precision_score, recall_score
from sklearn.pipeline import Pipeline
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy import stats
from astropy.io import fits
import argparse
import sys
import seaborn as sns
import os

parser = argparse.ArgumentParser(
    description="""Make a table from the S-PLUS catalogs """)

parser.add_argument("fileName", type=str,
                    default="teste-program",
                    help="Name of table, taken the prefix ")

parser.add_argument("fileName1", type=str,
                    default="teste-program",
                    help="Name of table, taken the prefix ")

cmd_args = parser.parse_args()
file_1 = cmd_args.fileName + ".ecsv"
file_2 = cmd_args.fileName1 + ".ecsv"

table1 = Table.read(file_1, format="ascii.ecsv")
table2 = Table.read(file_2, format="ascii.ecsv")

# Add a column with the label
n1 = len(table1)
label1 = np.linspace(0, 0, num=n1, dtype = int)
table1['Label'] = label1

n2 = len(table2)
label2 = np.linspace(1, 1, num=n2, dtype = int)
table2['Label'] = label2

# Merge the tables
table_merge = vstack([table1, table2])

# Put data in form expected by scikit-learn
X = np.array(list(zip(table_merge['U_PStotal'],
 table_merge['F378_PStotal'],
 table_merge['F395_PStotal'],
 table_merge['F410_PStotal'],
 table_merge['F430_PStotal'],
 table_merge['G_PStotal'],
 table_merge['F515_PStotal'],
 table_merge['R_PStotal'],
 table_merge['F660_PStotal'],
 table_merge['I_PStotal'],
 table_merge['F861_PStotal'],
 table_merge['Z_PStotal'])))

print("Shape of array:", X.shape)

##################################################################
#Accuracy
sc = StandardScaler()
y = table_merge["Label"]
 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)

X_train = sc.fit_transform(X_train) 
X_test = sc.transform(X_test)

lda1 = LDA(n_components=1)
lda1.fit(X_train, y_train)
#lda1.fit(X_train, y_train)

test_pred = lda1.predict(X_test)

cm = confusion_matrix(y_test, test_pred)
print('Accuracy score for Testing Dataset = ', accuracy_score(test_pred, y_test))
print('Precision score for Testing Dataset = ', precision_score(test_pred, y_test, labels=[1], average="micro"))
print('Confusion matrix = ', cm)
#############################################################################################################

# Standarized the data
X_stand = StandardScaler().fit_transform(X)

lda = LDA(n_components=1)
lda.fit(X_stand, y)

X_lda = lda.transform(X_stand)

########################################################################################################
#Predicting           ##################################################################################
########################################################################################################
#X_new = []
#File_name = input('Input file name:')
tab = Table.read("Halpha-DR3_noFlag_merge.ecsv", format="ascii.ecsv")
#tab = Table.read("TAP_DR1SPLUS_HA_r_03.tab", format="ascii.tab")
X_new = np.array(list(zip(tab['U_PStotal'],
 tab['F378_PStotal'],
 tab['F395_PStotal'],
 tab['F410_PStotal'],
 tab['F430_PStotal'],
 tab['G_PStotal'],
 tab['F515_PStotal'],
 tab['R_PStotal'],
 tab['F660_PStotal'],
 tab['I_PStotal'],
 tab['F861_PStotal'],
 tab['Z_PStotal'])))

print("Data to classify:", X_new.shape)

XX_new = sc.transform(X_new)
#XX_test = lda.transform(XX_test) #Ojo slólo uso con logistic regresion 

y_pred = lda.predict(XX_new)
y_prob = lda.predict_proba(XX_new)
y_pred_dec = lda.decision_function(XX_new)

print(len(y_pred))
# for a, b in zip(y_pred, y_prob):
#     print(a, b)

#creating table files with each class selected
label_pro=("P(GoodPho)", "P(BadPho)")

for label_Pro, label_Nu in zip(label_pro, range(2)):
    tab[label_Pro] = y_prob[:,label_Nu]

def select(clas, name_file):
    m = y_pred==clas
    return tab[m].write(name_file, format='ascii.ecsv', overwrite=True)

select(0, "Good-LD-Halpha-DR3_noFlag_merge.ecsv")
select(1, "Bad-LD-Halpha-DR3_noFlag_merge.ecsv")

############################################################################################################
############################################################################################################
# ASCII file
table_merge["LD1"] = X_lda

mask1 =  table_merge['Label'] == 0
mask2 =  table_merge['Label'] == 1

final_table1 = table_merge[mask1]
final_table2 = table_merge[mask2]

lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': None}
#sns.set(style="dark")#, context="talk")
#sns.set_style('ticks')       
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
# ax1.set_xlim(-8.2, 5.7)
# ax1.set_ylim(-2.5, 1.5)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax1.set_xlim(-8, 10)
#ax1.set_ylim(-10, 10)
#ax1.set_xlim(xmin=-2.5,xmax=2.0)
plt.tick_params(axis='x', labelsize=32) 
plt.tick_params(axis='y', labelsize=32)
plt.xlabel(r'r - i', fontsize= 35)
plt.ylabel(r'r - J0660', fontsize= 35)
#print(A1[0][1], B1[0][1])        
ax1.scatter(final_table1["r - i"],  final_table1["r - J0660"], c= "red", alpha=0.8, s=90, marker='o',  zorder=3.0, label='Good')
ax1.scatter(final_table2["r - i"],  final_table2["r - J0660"],  color= sns.xkcd_rgb["pale yellow"], alpha=0.9, s=60, marker='o', label='Bad')
# pal1 = sns.dark_palette("yellow", as_cmap=True)
# for ax, ay in zip(ld1[2], ld2[2]):
#     sns.kdeplot(ax, ay, cmap=pal1);

# plt.text(0.62, 0.78, 'Other Halpha emitters',
#          transform=ax1.transAxes, fontsize=13.8)

#ax1.legend(scatterpoints=1, ncol=2, fontsize=19.8, loc='upper left', **lgd_kws)
ax1.grid()
lgd = ax1.legend(loc='center right', bbox_to_anchor=(1.27, 0.5), fontsize=7.5, **lgd_kws)
#ax2.grid(which='minor', lw=0.5)
#sns.despine(bottom=True)
plt.tight_layout()
plt.tight_layout()
#pltfile = 'Fig1-JPLUS-PC1-PC2-veri.pdf'
pltfile = 'LDA-diagramColot.pdf'
save_path = 'Plots-splus/'
file_save = os.path.join(save_path, pltfile)
plt.savefig(pltfile)
plt.clf()
