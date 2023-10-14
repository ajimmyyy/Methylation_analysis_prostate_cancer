import pandas as pd
import numpy as np
from tqdm import tqdm

fn_cutpoint = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint.csv"
fn_beta = "C:/Users/acer/Desktop/Data-origin/validation/all_beta_normalized_validation.csv"
fn_o = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint_validate.csv"
fn_o_keep = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint_validate_keep.csv"

normal_num = 20
half_total_num = 110
half_normal_num = 10
threshold_validate_MAPE = 0.8
threshold_validate_dif = 0.05

data_beta_df = pd.read_csv(fn_beta)
data_cutpoint_df = pd.read_csv(fn_cutpoint)

def check_cutpoint(row):
    cutpoint = row["cutpoint"]
    filtered_row = data_beta_df[data_beta_df[data_beta_df.columns[0]] == row["CpG"]]

    transform_row = filtered_row.to_numpy()[0]
    normal_beta = transform_row[1:normal_num + 1:2]
    tumor_beta = transform_row[normal_num + 1::2]

    FP = np.sum(normal_beta > cutpoint)
    FN = np.sum(tumor_beta < cutpoint)
    TP = np.sum(tumor_beta > cutpoint)

    sensitivity = TP/(TP+FN)
    precision = TP/(TP+FP) if TP != 0 else 0

    F1 = 2 * sensitivity * precision / (sensitivity + precision)

    if row["DNAm"] == "hypo":
        F1 = 1 - F1

    return pd.Series({"F1_validate": F1})

tqdm.pandas(desc="find cutpoint")
data_cutpoint_df["F1_validate"] = data_cutpoint_df.progress_apply(check_cutpoint, axis = 1)
data_cutpoint_df["F1_dif"] = data_cutpoint_df["F1_validate"] - data_cutpoint_df["F1"]

data_cutpoint_df = data_cutpoint_df[data_cutpoint_df["F1_validate"] > threshold_validate_MAPE]

data_cutpoint_df.to_csv(fn_o_keep, sep=',', encoding='utf-8', index=False)

data_cutpoint_df = data_cutpoint_df[abs(data_cutpoint_df["F1_dif"]) < threshold_validate_dif]

data_cutpoint_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)



    