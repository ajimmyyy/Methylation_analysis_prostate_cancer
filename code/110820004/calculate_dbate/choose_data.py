import csv
import pandas as pd
import numpy as np
from tqdm import tqdm

fn_dmp = "Data/Data-origin/DMP_result_TN.csv"
fn_dbate = "all_beta_delta.csvData/Data-dbate/all_beta_delta.csv"
fn_o = "Data/Data-dbate/DMP_per_gene.csv.csv"

data_dmp_df = pd.read_csv(fn_dmp)
data_dbate_df = pd.read_csv(fn_dbate)
finded_gene = []
del_col_n = []

def find_dbate(row):
    filtered_rows = data_dbate_df[data_dbate_df[data_dbate_df.columns[0]] == row.iat[0]]
    return filtered_rows['dbate'].values[0]

def find_gene(row):
    if not row["gene"] in finded_gene:
        finded_gene.append(row["gene"])
    else:
        del_col_n.append(row.name)

tqdm.pandas(desc="find dbate")
data_dmp_df["dbate"] = data_dmp_df.progress_apply(find_dbate, axis = 1)
data_dmp_df = data_dmp_df.iloc[data_dmp_df['dbate'].abs().argsort()]

tqdm.pandas(desc="find gene")
data_dmp_df.progress_apply(find_gene, axis = 1)

for i in del_col_n:
    data_dmp_df.drop(i, inplace=True)

data_dmp_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
