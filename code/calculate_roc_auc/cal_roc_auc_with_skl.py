import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, auc

fn_bate = "Data/Data-origin/all_beta_normalized.csv"
fn_dmp = "Data/Data-GDC/DMP_filter_with_GDC_group.csv"
# fn_dmp = "Data/Data-cutpoint/dmp_with_cutpoint.csv"
fn_o = "Data/Data-ROC_AUC/DMP_with_auc.csv"
normal_num = 50

data_bate_df = pd.read_csv(fn_bate)

data_dmp_df = pd.read_csv(fn_dmp)
DMP_list = data_dmp_df.iloc[:, 0].tolist()
dmp_bate_df = data_bate_df.loc[data_bate_df[data_bate_df.columns[0]].isin(DMP_list)]

def cal_cutpoint(row):
    transform_row = row.to_numpy()
    TPRs = []
    FPRs = []

    predict_beta = transform_row[1::2]
    actual_beta = np.ones(278)
    actual_beta[:25] = 0

    fpr, tpr, threshold = roc_curve(actual_beta, predict_beta)

    auc1 = auc(fpr, tpr)

    return pd.Series({"auc": auc1})

data_out = dmp_bate_df.iloc[:, 0].to_frame()

tqdm.pandas(desc="find auc")
data_out["auc"] = dmp_bate_df.progress_apply(cal_cutpoint, axis = 1)

data_out.columns.values[0] = 'CpG'
data_out = pd.merge(data_dmp_df, data_out, on = ["CpG"], how = "inner")
data_out.loc[data_out["DNAm"] == "hypo", 'auc'] = 1 - data_out.loc[data_out["DNAm"] == "hypo", 'auc']

data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)