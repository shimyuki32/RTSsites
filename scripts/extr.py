#This script contains a bug due to the lack of one-to-one correspondence between NMID and ENSTID
import re
import pandas as pd
import argparse 

def conv_ENSTtoNM(series: pd.core.series.Series) -> pd.core.series.Series:
    fnc_filepath = args.idconversion_file
    df_conv = pd.read_csv(fnc_filepath)
    converted_list =[]
    for i in range(len(series)):
        element = series[i]
        converted_line = df_conv.query('converted_alias==@element')
        converted_list.append(converted_line.iat[0,0])
    return pd.Series(converted_list)
    # df_temp = df_conv[df_conv['converted_alias'].isin(ser.values.tolist())]
    # return df_temp['initial_alias']

parser = argparse.ArgumentParser()
parser.add_argument("input_csv", type=str, help="csv file of reactivity score for heatmap")
parser.add_argument("region_file", type=str, help="K+ or K++PDS region info file")
parser.add_argument("idconversion_file", type=str, help="gProfile id conversion file")
parser.add_argument("output_prefix", type=str, help="output file prefix")
args = parser.parse_args()

df_mRNA_1 = pd.read_csv(args.input_csv)

#split dataframe index into geneID and sequences
temp = re.compile('(ENST[0-9]+)([A-Z]+)')
ser = df_mRNA_1.iloc[:,0].str.extract(temp, expand=True)
#make new dataframe after split
df_rs = df_mRNA_1.iloc[:, 1:]
df_sp = pd.concat([ser, df_rs], axis=1)
df_sp.rename(columns={0:'ENST',1:'sequence'}, inplace=True)

#convert Ensemble gene ID to NCBI ID
list_mRNA1 = conv_ENSTtoNM(ser[0])


#make datafame for corresponding ENSTID & NMID
df_cor = pd.concat([ser[0].rename('ENST'), list_mRNA1.rename('transcript')], axis=1)

#extract applicable rows from the original datasets
df_RTSsites = pd.read_csv(args.region_file)
df_list_type = df_RTSsites[df_RTSsites['transcript'].isin(list_mRNA1.values.tolist())]
#sort df by NMID and reset index
df_list_type = df_list_type.sort_values('transcript')
df_list_type = df_list_type.reset_index(drop=True)

#add ENST ID to the RTSsites data
df_temp = pd.merge(df_list_type, df_cor, on='transcript', how='inner')

#make complete dataframe from df_sp and df_temp 
df_pf = pd.merge(df_temp, df_sp, on=['ENST','sequence'], how='outer')

#sort by classes of rG4
class_order =['PQs_7','long_loop', 'bulge', '2-quartet','G_40', 'Others' ]
df_pf['order'] = df_pf['class'].apply(lambda x: class_order.index(x) if x in class_order else -1)
df_sorted = df_pf.sort_values('order').reset_index(drop=True).drop(columns='order')

#avoid bug (check later why it happens)
df_sorted.dropna(subset=['0'], inplace=True)
df_sorted.drop_duplicates(inplace=True)

groups = df_sorted.groupby('class')
df_PQs_7 = groups.get_group('PQs_7')

df_sorted.to_csv('/mnt/backup02/0/RTSsites/{}.updated_dataframe.csv'.format(args.output_prefix),index=False)


