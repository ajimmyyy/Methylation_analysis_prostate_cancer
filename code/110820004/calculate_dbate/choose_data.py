import csv
import pandas as pd
import numpy as np
from tqdm import tqdm
from config_loader import config

fn_dmp = config["DNP_RESULT_TRAIN_PATH"]
fn_dbate = config["ALL_BETA_APPEND_DBATE_PATH"]
# fn_o = config["DMP_PER_GENE_PATH"]
fn_o = "C:/Users/acer/Desktop/tmp/tmp.csv"

def choose_dbate_per_gene(dmp_result, all_dbeta):
    unique_genes = set()

    def find_dbate(row):
        filtered_rows = all_dbeta[all_dbeta[all_dbeta.columns[0]] == row.iat[0]]
        dbate_value = filtered_rows['dbate'].values[0]
        return dbate_value, row

    def find_gene(row):
        gene = row["gene"]
        if gene in unique_genes:
            return None
        unique_genes.add(gene)
        dbate_value, row = find_dbate(row)
        row["dbate"] = dbate_value
        return row

    tqdm.pandas(desc="find gene")
    filtered_data = dmp_result.progress_apply(find_gene, axis=1).dropna()

    return filtered_data

data_dmp_df = pd.read_csv(fn_dmp)
data_dbate_df = pd.read_csv(fn_dbate)

data_dmp_df = choose_dbate_per_gene(data_dmp_df, data_dbate_df)
data_dmp_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
