import matplotlib.pyplot as plt
import pandas as pd

df1 = pd.read_csv("temp.whole.dataframe.csv")['ENST00000234875GGGAGGCAGATGCGGGTGTGGAGGTGTGGG']
df2 = pd.read_csv("temp.whole.dataframe.csv")['ENST00000462138GCCTAGGGCGGGTATGAAGTGTGGGGCGGG']
df1.plot()
df2.plot()
plt.savefig("temp.png")