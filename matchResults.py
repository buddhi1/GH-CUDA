# This script check if the given two files are identical or it return first mismatch
# Used to compare parallel output vs sequential to verify results
import pandas as pd

# file1 = "results/Fig-large-R.poly"
# file2 = "../2019-Clipping simple polygons with degenerate intersections/polyclip2/results/Fig-large-R.poly"
file1 = "data/readPolygon/s.txt"
file2 = "../"
df_file1 = pd.read_csv (file1,  sep=" ", header=None)
df_file2 = pd.read_csv (file2,  sep=" ", header=None)

print ("Files Matching? %s" % df_file1.equals(df_file2))
# --------------------------------------------------------
# removes comma and semicolan if needed
# df_file1 = df_file1.replace(',','', regex=True)
# df_file1 = df_file1.replace(';','', regex=True)
# print (df_file1)
# --------------------------------------------------------
for i in range(len(df_file1)):
    if(not df_file1.iloc[i].equals(df_file2.iloc[i])):
        print("Line: %d" % i)
        print ("First mistmatch: "+str(df_file1.iloc[0])+" "+str(df_file2.iloc[0]))
        break
