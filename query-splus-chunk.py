from pyvo.dal import TAPService
import astropy
from astropy.table import Table, vstack
import splusdata
import sys

#tap = TAPService("http://tap.roe.ac.uk/osa/")
conn = splusdata.connect('Luis', 'plutarco*80')
conn.get_tap_tables()

maxPrevSourceID = 0
merged_table = None
merged_table_list = []
chunksize = 1000
maxcount = 1000000
count = 0

#while count < maxcount:
#query_text = f"SELECT id, ra, dec FROM idr3.detection_image WHERE ra + dec > 200"
conn("SELECT id, ra, dec FROM idr3.detection_image WHERE ra + dec > 200")

   #  res_table = results.resultstable.to_table(use_names_over_ids=True)

#     if len(res_table) == 0:
#         break

#     merged_table_list.append(res_table)
#     maxPrevSourceID = res_table[-1]["sourceID"]

#     if len(res_table) < chunksize:
#         break

#     count += 1

# merged_table = vstack (merged_table_list)

# merged_table.write("votable.xml", format="votable")
