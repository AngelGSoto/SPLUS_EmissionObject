'''
Author: Luis A. Gutiérrez
15/03/2021
Scrito to descargar the data form my databese using queris
'''
import sqlite3
import os
import csv
import pandas as pd

# Create a SQL connection to our SQLite database

def dr3qury(columns, catalog, condition):
    try:
        conn = sqlite3.connect('IDR3.db')

    except sqlite3.DatabaseError as e:

        # Confirm unsuccessful connection and quit.
        print("Database connection unsuccessful.")
        quit()

    # Query
    qry = ("SELECT" + " " + columns + " " + "FROM" + " " +  catalog + " " +  "WHERE" + " " +  condition + ";")
    df = pd.read_sql_query(qry, conn)

    # Verify that result of SQL query is stored in the dataframe
    print(df.head(3))

    conn.close()
    return df
