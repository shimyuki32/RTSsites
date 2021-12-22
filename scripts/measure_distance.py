from os import getgrouplist
import pandas as pd

df_all = pd.read_csv("overall.csv")
df_all = df_all[['ENST', 'sequence', 'start', 'end']]
genegroup = df_all.groupby('ENST')

genegroup_list = []
for name, group in genegroup:
    if len(group) != 1:
            genegroup_list.append(group)

df_multiple = pd.concat(genegroup_list)
df_multiple.to_csv("temp.csv")

print("----done----")