import matplotlib.pyplot as plt
from numpy import nan
import pandas as pd
import seaborn as sns

df =  pd.read_csv("output/mRNA/K+_PDS/replicate1/dataframe/K+_PDS_mRNA1_whole_reacitvity.csv", index_col=0)
number = len(df.columns)

print(number)
temp = [5,25]
list = []
for i in range(temp[0], temp[1]):
    whole = df.iloc[:, i].dropna().rename("scores")
    g4 = pd.read_csv("output/mRNA/K+_PDS/replicate1/dataframe/K++PDS.mRNA.1.dataframe.csv").iloc[i, 10:40].rename("scores")

    df_whole = pd.DataFrame(whole)
    df_whole["type"] = "whole"
    df_whole["ID"] = df.columns.tolist()[i]

    df_g4 = pd.DataFrame(g4)
    df_g4["type"] = "g4"
    df_g4["ID"] = df.columns.tolist()[i]

    df_transcript = pd.concat([df_whole, df_g4])
    list.append(df_transcript)

entire_dataframe = pd.concat(list)

print(entire_dataframe)




fig ,ax = plt.subplots()
ax.set_ylim([0, 0.1])
ax.tick_params(labelbottom=False)
ax.tick_params(bottom=False)
sns.boxplot(x='ID', y='scores', data=entire_dataframe, hue='type', ax=ax ,fliersize=0)
plt.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=1, fontsize=9)
plt.savefig("temp/temp_entire.png")


# ax.boxplot([whole, g4])
# plt.show()
# plt.savefig("temp_boxplot.png")

print("done")