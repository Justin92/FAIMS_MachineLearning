
import pandas as pd
import glob

df_list = []
for tsv in glob.glob("*.tsv"):
    df = pd.read_csv(tsv, sep = "\t")
    df_list.append(df)

combined_df = pd.concat(df_list, ignore_index=True, sort=True)
combined_df.to_csv("model_performance.tsv", sep = "\t")
