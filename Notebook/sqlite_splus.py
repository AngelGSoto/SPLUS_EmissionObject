'''
Author: Luis A. Gutiérrez
15/03/2001
Using sqlite via pandas only for dr3 splus
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
    qry = ("SELECT" + " " + columns + " " + "FROM" " " + catalog + " " + "WHERE" + " " + " " + condition + ";")
    #for row in cur.execute(qry):
    #print(row)

    # Executing the query usin pandas
    df = pd.read_sql_query(qry, conn)

    # Verify that result of SQL query is stored in the dataframe
    print(df.head())

    conn.close()
    return df
