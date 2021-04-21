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
from sklearn.decomposition import PCA
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

pca_sco = PCA(n_components = 12) 
  
X_train_pca = pca_sco.fit_transform(X_train)

X_test = pca_sco.transform(X_test) 

#classifier = LogisticRegression(solver='lbfgs', max_iter = 12000)
classifier = LogisticRegression(random_state = 0)
classifier.fit(X_train_pca, y_train)

# Predicting the test set result using  
# predict function under LogisticRegression  
test_pred = classifier.predict(X_test)

cm = confusion_matrix(y_test, test_pred)
print('Accuracy score for Testing Dataset = ', accuracy_score(test_pred, y_test))
print('Precision score for Testing Dataset  = ', precision_score(test_pred, y_test, average="micro"))
print('Confusion matrix = ', cm)

#############################################################################################################

# Standarized the data
X_stand = StandardScaler().fit_transform(X)

# Creating the PCA 
pca = PCA(n_components=10)
pca.fit(X_stand)

X_pca = pca.transform(X_stand)

# Porcentages, eige-vectors and values
porc0 = []
pc_name = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
porc = pca.explained_variance_ratio_ # porcantage ratio
for perc in porc:
    porc0.append(perc)

perc1 = Table([pc_name, porc0], names=('PCs', '%'), meta={'name': 'first table'})

###########################################################################################################
einvector = Table(pca.components_, names=('V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                                          'V7', 'V8', 'V9', 'V10', 'V11', 'V12'), meta={'name': 'first table'}) # Eigevectores

einvector["PCs"] = pc_name

new_order_eigenvetor = ['PCs', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                                          'V7', 'V8', 'V9', 'V10', 'V11', 'V12']

einvector = einvector[new_order_eigenvetor]

############################################################################################################
einvalue1 = []
einvalue = pca.explained_variance_ # Eigivalues
for einvalue0 in einvalue:
    einvalue1.append(einvalue0)

einvalue2 = Table([pc_name, einvalue1], names=('PCs', 'EigenValues'), meta={'name': 'first table'})
###########################################################################################################


print("Porcentage:", pca.explained_variance_ratio_)
print("Singular Value:", pca.singular_values_)
print("Component:", pca.components_) # eigevectors
print("Sorted components:", pca.explained_variance_) # eigenvalues

############################################################################################################
############################################################################################################
# ASCII file
pca_table = Table(X_pca, names=('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6',
                               'PC7', 'PC8', 'PC9', 'PC10'), meta={'name': 'first table'})

final_table = hstack([table_merge, pca_table])

mask1 =  final_table['Label'] == 0
mask2 =  final_table['Label'] == 1

final_table1 = final_table[mask1]
final_table2 = final_table[mask2]

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
XX_new = pca.transform(XX_new)

#Fitting Logistic Regression to the Training set
classifier_new = LogisticRegression(random_state = 0)
classifier_new.fit(X_pca, y)
y_pred = classifier_new.predict(XX_new)
y_prob = classifier_new.predict_proba(XX_new)

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

select(0, "Good-PC-Halpha-DR3_noFlag_merge.ecsv")
select(1, "Bad-PC-Halpha-DR3_noFlag_merge.ecsv")

###########################################################################################################
#PLOTS
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': None}
#sns.set(style="dark")#, context="talk")
#sns.set_style('ticks')       
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
# ax1.set_xlim(-8.2, 5.7)
# ax1.set_ylim(-2.5, 1.5)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax1.set_xlim(-10, 10)
# ax1.set_ylim(-3.3, 3.2)
#ax1.set_xlim(xmin=-2.5,xmax=2.0)
plt.tick_params(axis='x', labelsize=32) 
plt.tick_params(axis='y', labelsize=32)
plt.xlabel(r'PC1', fontsize= 35)
plt.ylabel(r'PC2', fontsize= 35)

ax1.scatter(final_table1["PC1"], final_table1["PC2"],  color= sns.xkcd_rgb["aqua"], s=130, marker='o', alpha=0.8, edgecolor='black', zorder=80.0, label='Good')
ax1.scatter(final_table2["PC1"], final_table2["PC2"],  color= sns.xkcd_rgb["dark pink"], s=130, marker='o', alpha=0.8, edgecolor='black', zorder=80.0, label='bad')

ax1.grid()
lgd = ax1.legend(loc='center right', bbox_to_anchor=(1.27, 0.5), fontsize=7.5, **lgd_kws)
#ax2.grid(which='minor', lw=0.5)
#sns.despine(bottom=True)
plt.tight_layout()
plt.tight_layout()

pltfile = 'Fig1-PC1-PC2.pdf'
save_path = ' '
file_save = os.path.join(pltfile)
plt.savefig(file_save)
plt.clf()

####################################################################
#PC1 vs PC3 ########################################################
####################################################################

lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': None}
#sns.set(style="dark")#, context="talk")
#sns.set_style('ticks')       
fig = plt.figure(figsize=(12, 8))
ax2 = fig.add_subplot(111)
# ax2.set_xlim(-10.0, 8.0)
# ax2.set_ylim(-2.0, 1.5)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax2.set_xlim(-8.2, 5.7)
# ax2.set_ylim(-1.1, 1.0)
#ax1.set_xlim(xmin=-2.5,xmax=2.0)
plt.tick_params(axis='x', labelsize=32) 
plt.tick_params(axis='y', labelsize=32)
plt.xlabel(r'PC1', fontsize= 35)
plt.ylabel(r'PC3', fontsize= 35)

ax2.scatter(final_table1["PC1"], final_table1["PC3"],  color= sns.xkcd_rgb["aqua"], s=130, marker='o', alpha=0.8, edgecolor='black', zorder=80.0, label='Good')
ax2.scatter(final_table2["PC1"], final_table2["PC3"],  color= sns.xkcd_rgb["dark pink"], s=130, marker='o', alpha=0.8, edgecolor='black', zorder=80.0, label='bad')


ax2.legend(scatterpoints=1, ncol=2, fontsize=17.8, loc='upper center', **lgd_kws)
ax2.grid()
#ax2.grid(which='minor', lw=0.5)
#sns.despine(bottom=True)
plt.tight_layout()
#pltfile = 'Fig2-JPLUS-PC1-PC3-veri.pdf'
pltfile = 'Fig2-PC1-PC3.pdf'
file_save = os.path.join(pltfile)
plt.savefig(file_save)
plt.clf()
