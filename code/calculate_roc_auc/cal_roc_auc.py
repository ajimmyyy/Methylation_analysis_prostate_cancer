import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

fn_bate = "Data/Data-origin/all_beta_normalized.csv"
fn_dmp = "Data/Data-GDC/DMP_filter_with_GDC.csv"
fn_o = "Data/Data-ROC_AUC/DMP_with_auc.csv"
normal_num = 50

data_bate_df = pd.read_csv(fn_bate)

data_dmp_df = pd.read_csv(fn_dmp)
DMP_list = data_dmp_df.iloc[:, 0].tolist()
dmp_bate_df = data_bate_df.loc[data_bate_df[data_bate_df.columns[0]].isin(DMP_list)]

def cal_cutpoint(row):
    auc = 0.0
    roc_point = []

    for cutpoint in np.arange(0.001, 1, 0.001):
        transform_row = row.to_numpy()
        normal_beta = transform_row[1:normal_num + 1:2]
        tumor_beta = transform_row[normal_num + 1::2]
        TN = np.sum(normal_beta < cutpoint)
        FP = np.sum(normal_beta > cutpoint)
        FN = np.sum(tumor_beta < cutpoint)
        TP = np.sum(tumor_beta > cutpoint)

        TPR = TP / (TP + FN)
        FPR = FP / (TN + FP)

        roc_point.append([FPR, TPR])
    
    roc_point.sort(key=lambda x: x[0])
    # plt.plot(np.array(roc_point)[1:, 0], np.array(roc_point)[1: ,1], c = "red", alpha = .5)

    lastx = 0.
    for x,y in roc_point:
        auc += (x - lastx) * y
        lastx = x

    return pd.Series({"auc": auc})

data_out = dmp_bate_df.iloc[:, 0].to_frame()

tqdm.pandas(desc="find roc/auc")
data_out["auc"] = dmp_bate_df.progress_apply(cal_cutpoint, axis = 1)

# plt.plot([0, 1], [0, 1],color="green")
# plt.xlabel("FPR")
# plt.ylabel("TPR")
# plt.show()

data_out.columns.values[0] = 'CpG'
data_out = pd.merge(data_dmp_df, data_out, on = ["CpG"], how = "inner")
data_out.loc[data_out["DNAm"] == "hypo", 'auc'] = 1 - data_out.loc[data_out["DNAm"] == "hypo", 'auc']

data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)