import pandas as pd
 
filepath1 = "K+.mRNA.1.ver2_dataframe.csv"
filepath2 = "K+.mRNA.2.ver2_dataframe.csv"
filepath3 = "K++PDS.mRNA.1.ver2_dataframe.csv"
filepath4 = "K++PDS.mRNA.2.ver2_dataframe.csv"

df1 = pd.read_csv(filepath1, index_col=0)
df2 = pd.read_csv(filepath2, index_col=0)
df3 = pd.read_csv(filepath3, index_col=0)
df4 = pd.read_csv(filepath4, index_col=0)

df1t = df1.loc[:, "ENST":"class"]
df2t = df2.loc[:, "ENST":"class"]
df3t = df3.loc[:, "ENST":"class"]
df4t = df4.loc[:, "ENST":"class"]

df_all = pd.concat([df1t, df2t, df3t, df4t])
df_all.reset_index(drop=True, inplace=True)

df_all.drop(['p.value', 'rts_value'],  inplace=True, axis=1)
df_all.drop_duplicates(inplace=True)

df_all.to_csv("overall.csv")

df_all_count = df_all.loc[:, "ENST":"sequence"]
df_all_count.drop_duplicates(inplace=True)
df_all_count.to_csv("overall_genes_dropduplicate.csv", index=False, header=False)

overall_seqcount = df_all_count['sequence'].value_counts()
overall_seqcount.to_csv("overall_seqcount.csv")
overall_count = df_all_count['ENST'].value_counts()
overall_count.to_csv("overall_countgene.csv")
