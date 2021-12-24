import matplotlib.pyplot as plt
from numpy import nan
import pandas as pd
import seaborn as sns
import numpy as np

df_original_whole =  pd.read_csv("output/mRNA/K+/replicate2/dataframe/K+_mRNA2_whole_reacitvity.csv", index_col=0)
df_original_g4 = pd.read_csv("output/mRNA/K+/replicate2/dataframe/K+.mRNA.2.ver2_dataframe.csv")
number = len(df_original_whole.columns)

print(number)
temp = [5,7]   
bar_list = []
diff_list = []

for i in range(number):
    whole = df_original_whole.iloc[:, i].dropna().rename("scores")
    g4 = df_original_g4.loc[i, '10':'40'].rename("scores")

    ID = df_original_g4.at[i, 'ENST'] + ' ' + df_original_g4.at[i, 'sequence']
    g4class = df_original_g4.at[i, 'class']

    df_whole = pd.DataFrame(whole)
    df_whole["label"] = "whole"

    df_g4 = pd.DataFrame(g4)
    df_g4["label"] = "g4"

    diff = whole.median() - g4.median()

    df_transcript = pd.concat([df_whole, df_g4])

    df_transcript["ID"] = ID

    df_transcript["class"] = g4class
    bar_list.append(df_transcript)
    

    diff_ser = pd.Series(
        {
            'diff' : diff,
            'class' : g4class,
        },
        name= ID
    )
    diff_list.append(diff_ser)

entire_dataframe = pd.concat(bar_list)
diff_dataframe = pd.concat(diff_list, axis=1)
diff_dataframe = diff_dataframe.T
diff_dataframe['positive'] = diff_dataframe['diff'] > 0


print(entire_dataframe)





fig ,ax = plt.subplots()
fig.set_size_inches(10, 10)
ax.set_ylim([0, 0.1])
ax.tick_params(bottom=False)
sns.boxplot(x='ID', y='scores', data=entire_dataframe, hue='label', ax=ax ,fliersize=0)

ax.set_xticklabels(diff_dataframe['class'].tolist(),rotation=90)
plt.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=1, fontsize=9)
plt.tight_layout()
plt.savefig("temp/temp_entire.png")

fig ,ax = plt.subplots()
fig.set_size_inches(10, 10)
ax.set_ylim([-0.02, 0.02])

ax.tick_params(bottom=False)
diff_dataframe.plot.bar(y='diff', ax=ax, grid=True, color=diff_dataframe.positive.map({True: 'dodgerblue', False: 'orange'}))
ax.set_xticklabels(diff_dataframe['class'].tolist(),rotation=90)
# ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
plt.tight_layout()
plt.savefig("temp/temp_boxplot.png")
# ax.boxplot([whole, g4])
# plt.show()
# plt.savefig("temp_boxplot.png")

print("done")