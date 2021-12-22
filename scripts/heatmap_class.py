from pandas.io.parsers import read_csv
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse 
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="updated selected reactivity score file")
parser.add_argument("output_prefix", type=str, help ="output file prefix")
args = parser.parse_args()

filepath = "/home/gplab/temporary/backup02/shimizu/RTSsites/K+.mRNA.1.updated_dataframe.csv"
df = read_csv(args.input_file)

groups = df.groupby('class', sort=False)

fig, axes = plt.subplots(1, 2, figsize=(12,6))
axes[1].axis('off')
sns.heatmap(df.loc[:, '0':'49'], ax=axes[1], cmap='Blues')
axes[0].axis('off')
table = axes[0].table(df.loc[:, ['class','transcript', 'annotation']].astype(str).values.tolist(), loc='center')
for pos, cell in table.get_celld().items():
    cell.set_height(1/len(df))
table.scale(0.5, 1)
table.auto_set_font_size(False)
table.set_fontsize(3)

fig.suptitle("{}".format(args.output_prefix))

plt.savefig("/home/gplab/temporary/backup02/shimizu/RTSsites/{}.png".format(args.output_prefix), dpi=300)

