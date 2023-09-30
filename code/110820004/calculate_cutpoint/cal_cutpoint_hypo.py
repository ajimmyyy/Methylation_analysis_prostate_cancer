import csv
import pandas as pd
import numpy as np
from tqdm import tqdm

fn_beta = "Data/Data-origin/all_beta_normalized.csv"
fn_dmp_hypo = "Data/Data-volcano/DMP_hypo.csv"
fn_o = "Data/Data-cutpoint/hypo_with_cutpoint.csv"
normal_num = 50

data_beta_df = pd.read_csv(fn_beta)

data_dmp_df = pd.read_csv(fn_dmp_hypo)
DMP_list = data_dmp_df.iloc[:, 0].tolist()
dmp_beta_df = data_beta_df.loc[data_beta_df[data_beta_df.columns[0]].isin(DMP_list)]

def cal_cutpoint(row):
    max_total = 0.0
    max_cutpoint = 0.0
    max_F1 = 0.0

    for cutpoint in np.arange(0.01, 1, 0.01):
        transform_row = row.to_numpy()
        normal_beta = transform_row[1:normal_num + 1:2]
        tumor_beta = transform_row[normal_num + 1::2]
        TN = np.sum(normal_beta > cutpoint)
        FP = np.sum(normal_beta < cutpoint)
        FN = np.sum(tumor_beta > cutpoint)
        TP = np.sum(tumor_beta < cutpoint)

        sensitivity = TP/(TP+FN)
        specificity = TN/(FP+TN)
        precision = TP/(TP+FP)

        if max_total < sensitivity + specificity:
            max_total = sensitivity + specificity
            max_cutpoint = cutpoint
            max_F1 = 2 * sensitivity * precision / (sensitivity + precision)

    return pd.Series({"cutpoint": round(max_cutpoint, 2), "F1": max_F1})

data_out = dmp_beta_df.iloc[:, 0].to_frame()
tqdm.pandas(desc="find cutpoint")
data_out[["cutpoint", "F1"]] = dmp_beta_df.progress_apply(cal_cutpoint, axis = 1)

data_out.columns.values[0] = 'CpG'
data_dmp_df.columns.values[0] = 'CpG'

data_out = pd.merge(data_dmp_df, data_out, on = ["CpG"], how = "inner")

data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)