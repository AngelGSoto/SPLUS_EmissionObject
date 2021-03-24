'''
Author: Luis A. Gutiérrez
11/09/2020
Based on sqlite-query.py
'''
import sqlite3
import os
import csv

# Create a SQL connection to our SQLite database
try:
    conn = sqlite3.connect('IDR3.db')

except sqlite3.DatabaseError as e:

        # Confirm unsuccessful connection and quit.
        print("Database connection unsuccessful.")
        quit()
        
cur = conn.cursor()


#  query
qry = ("SELECT Field, ID, RA, DEC, FWHM, ISOarea, KRON_RADIUS, nDet_magPStotal, PhotoFlagDet, U_PStotal, F378_PStotal, F395_PStotal, F410_PStotal, F430_PStotal, G_PStotal, F515_PStotal, R_PStotal, F660_PStotal, I_PStotal, F861_PStotal, Z_PStotal, e_U_PStotal, e_F378_PStotal, e_F395_PStotal, e_F410_PStotal, e_F430_PStotal, e_G_PStotal, e_F515_PStotal, e_R_PStotal, e_F660_PStotal, e_I_PStotal, e_F861_PStotal, e_Z_PStotal FROM HYDRA NATURAL JOIN STRIPE82 limit 10 WHERE PhotoFlag_U <= 3.0 AND PhotoFlag_F378 <= 3.0 AND PhotoFlag_F395 <= 3.0 AND PhotoFlag_F410 <= 3.0 AND PhotoFlag_F430 <= 3.0 AND PhotoFlag_G <= 3.0 AND PhotoFlag_F515 <= 3.0 AND PhotoFlag_R <= 3.0 AND PhotoFlag_F660 <= 3.0 AND PhotoFlag_I <= 3.0 AND PhotoFlag_F861 <= 3.0 AND  PhotoFlag_Z <= 3.0 AND R_PStotal < 16 AND e_U_PStotal <= 0.2 AND e_F378_PStotal <= 0.2 AND e_F395_PStotal <= 0.2 AND e_F410_PStotal <= 0.2 AND e_F430_PStotal <= 0.2 AND e_G_PStotal <= 0.2 AND e_F515_PStotal <= 0.2 AND e_R_PStotal <= 0.2 AND e_F660_PStotal <= 0.2 AND e_I_PStotal <= 0.2 AND e_F861_PStotal <= 0.2 AND e_Z_PStotal <= 0.2;")
#for row in cur.execute(qry):
    #print(row)

cur.execute(qry)
data = cur.fetchall()
#print(data)

# Extract the table headers
headers = [i[0] for i in cur.description]

# Open CSV file for writing
file_name = 'Halpha_SPLUS_DR3/DR3_errorsall_flag3f_r16.csv'
csv_file = csv.writer(open(file_name, 'w', newline=''),
                             delimiter=',', lineterminator='\r\n',
                             quoting=csv.QUOTE_ALL, escapechar='\\')

# Add the headers and data to the CSV file.
csv_file.writerow(headers)
csv_file.writerows(data)

# Message stating export successful.
print("Data export successful ans writting the file: {}.".format(file_name))

# Be sure to close the connection
conn.close()


# HYDRA INNER JOIN STRIPE83 ON HYDRA.ID=STRIPE83.ID INNER JOIN MAIN3_1 ON HYDRA.ID=MAIN3_1.ID INNER JOIN MAIN3_2 ON HYDRA.ID=MAIN3_2.ID INNER JOIN MAIN3_3 ON HYDRA.ID=MAIN3_3.ID INNER JOIN MAIN3_4 ON HYDRA.ID=MAIN3_4.ID INNER JOIN MAIN3_5 ON HYDRA.ID=MAIN3_5.ID INNER JOIN MAIN3_6 ON HYDRA.ID=MAIN3_6.ID INNER JOIN MAIN3_7 ON HYDRA.ID=MAIN3_7.ID INNER JOIN MAIN3_8 ON HYDRA.ID=MAIN3_8.ID 
