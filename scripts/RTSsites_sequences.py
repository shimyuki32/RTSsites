import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
import re
import argparse

COV_CUTOFF = 200


def extract_ref_seq(file):
    cmd = 'cat {} | awk -F"," \'BEGIN{{printf "#contig,start_pos,ref"}}\
            {{gsub("\\"","",$0);if(match($0,/Gene/)){{next}}\
            else if(!a[$1]++){{printf "\\n"$1","$2","$3}}\
            else{{printf $3}}}}\'> transcripts'.format(
        file
    )
    sp.run(cmd, shell=True)
    out = pd.read_csv("transcripts")
    return out


def search_region(seq, ref):
    pattern = re.compile(seq)
    iterator = pattern.finditer(ref)
    output = []
    for match in iterator:
        output.append(match.span())
    return output


parser = argparse.ArgumentParser()
parser.add_argument(
    "reactivity_file", type=str, help="porecupine position-wise reactivity score file"
)
parser.add_argument("region_file", type=str, help="K+ or K++PDS region info file")
parser.add_argument("idconversion_file", type=str, help="gProfile id conversion file")
parser.add_argument("output_prefix", type=str, help="output file prefix")
parser.add_argument("previous_file", type=str, help="previous dataframe file")
args = parser.parse_args()

# Load gProfiler_hsapiens_~~~.csv to pd.DataFrame

print("loading {}...".format(args.idconversion_file), end="")
df_convert = pd.read_csv(args.idconversion_file)
df_convert = df_convert[["initial_alias", "converted_alias"]]
df_convert = df_convert.drop_duplicates()
df_convert = df_convert.set_index("initial_alias")
print("done")

# Load RTS sites under ~~~ condition.csv to pd.DataFrame
print("loading {}...".format(args.region_file), end="")
df_K = pd.read_csv(args.region_file)
df_K = df_K[["sequence", "transcript"]]
df_K = df_K.set_index("transcript")
df_K = df_K.join(df_convert)
df_K.to_csv("Ks")
df_K = df_K.set_index("converted_alias")
print("done")

# Load GSE133361_H9_~~~.csv to pd.DataFrame
print("loading {}...".format(args.reactivity_file), end="")
df_react = pd.read_csv(args.reactivity_file)
df_react = df_react[["Gene", "Position", "Mod_percentage", "Mod_strands"]]
df_react = df_react[df_react["Mod_strands"] >= COV_CUTOFF]
df_react = df_react[["Gene", "Position", "Mod_percentage"]]
df_react["Gene"] = df_react["Gene"].apply(lambda x: x.split(".")[0])
df_react = df_react.set_index("Gene")
# Extract reference sequence for convenience
df_refs = extract_ref_seq(args.reactivity_file)
df_refs["#contig"] = df_refs["#contig"].apply(lambda x: x.split(".")[0])
df_refs = df_refs.set_index("#contig")
print("done")

df_around_seq = pd.DataFrame(columns=['ENST', 'sequence', 'around_seq'])
up_length = 100
down_length = 100
i = 1
for name, seq in df_K.iterrows():
    try:
        print("\rsearching for {}th seq".format(i), end="")
        i += 1
        seq = seq["sequence"]
        df_r_pert = df_react.loc[name]
        refs = df_refs.loc[name]
        # if type(refs) == pd.core.series.Series:
        gap = refs["start_pos"]
        regions = search_region(seq, refs["ref"])
        for start, end in regions:
            print(": match {} at {}({},{})".format(seq, name, gap + start, gap + end))
            df = df_r_pert[
                (df_r_pert["Position"] >= gap + start - 10)
                & (df_r_pert["Position"] < gap + end + 10)
            ]["Mod_percentage"]
            if df.size != len(seq) + 20:
                print(
                    "{}/{} not enough confident, discard region".format(
                        df.size, len(seq) + 20
                    )
                )
                continue
            if start-up_length < 0:
                around_seq = refs['ref'][0:end+down_length]
            else:
                around_seq = refs['ref'][start-up_length:end+down_length]
            df_seq = pd.Series(data=[name, seq, around_seq], index=['ENST', 'sequence', 'around_seq'], name='{}'.format(i))
            df_around_seq = df_around_seq.append(df_seq)
            
            # else:
            #     for _, row in refs:
            #         gap = row["start_pos"]
            #         regions = search_region(seq, row["ref"])
            #         for start, end in regions:
            #             print(
            #                 ": match {} at {}({},{})".format(
            #                     seq, name, gap + start, gap + end
            #                 )
            #             )
            #             df = df_r_pert[
            #                 (df_r_pert["Position"] >= gap + start - 10)
            #                 & (df_r_pert["Position"] < gap + end + 10)
            #             ]["Mod_percentage"]
            #             if df.size != len(seq) + 20:
            #                 print(
            #                     "{}/{} not enough confident, discard region".format(
            #                         df.size, len(seq) + 20
            #                     )
            #                 )
            #                 continue
            #             df = df.reset_index()["Mod_percentage"]
            #             df = df.rename("{}{}".format(name, seq))
            #             df_map_regions = df_map_regions.append(df)
    except:
        continue

# Save dataframe
df_around_seq.reset_index(drop=True)

df_pre = pd.read_csv(args.previous_file)
df_new = pd.merge(df_around_seq, df_pre, how ='inner',on=['sequence','ENST'])

df_new.to_csv("{}.ver2_dataframe.csv".format(args.output_prefix))


print("\n=============== done ==============")
