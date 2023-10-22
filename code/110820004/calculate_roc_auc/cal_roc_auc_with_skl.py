import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, auc

def cal_auc(row, total_num, normal_num):
    transform_row = row.to_numpy()
    TPRs = []
    FPRs = []

    predict_beta = transform_row[1::2]
    actual_beta = np.ones(total_num)
    actual_beta[:normal_num] = 0

    fpr, tpr, threshold = roc_curve(actual_beta, predict_beta)

    auc_score = auc(fpr, tpr)

    return pd.Series({"auc": auc_score})

if __name__ == "__main__":
    fn_bate = "C:/Users/acer/Desktop/Data-origin/train/all_beta_normalized.csv"
    # fn_dmp = "Data/110820004/Data-Comorbidity/HyperHypo_filtered_comorbidity_single.csv"
    fn_dmp = "Data/110820004/Data-Comorbidity/HyperHypo_filtered_comorbidity_group.csv"
    # fn_o = "Data/110820004/Data-ROC_AUC/comorbidity_single_auc.csv"
    fn_o = "Data/110820004/Data-ROC_AUC/comorbidity_group_auc.csv"
    # fn_o_pic = "Data/110820004/Data-ROC_AUC/ROC_curve_single.png"
    fn_o_pic = "Data/110820004/Data-ROC_AUC/ROC_curve_group.png"

    half_total_num = 278
    half_normal_num = 25

    data_bate_df = pd.read_csv(fn_bate)
    data_dmp_df = pd.read_csv(fn_dmp)
    DMP_list = data_dmp_df.iloc[:, 0].tolist()
    dmp_bate_df = data_bate_df.loc[data_bate_df[data_bate_df.columns[0]].isin(DMP_list)]

    data_out = dmp_bate_df.iloc[:, 0].to_frame()

    tqdm.pandas(desc="find auc")
    data_out["auc"] = dmp_bate_df.progress_apply(cal_auc, axis = 1, args = (half_total_num, half_normal_num,))

    data_out.columns.values[0] = 'CpG'
    data_out = pd.merge(data_dmp_df, data_out, on = ["CpG"], how = "inner")
    data_out.loc[data_out["DNAm"] == "hypo", 'auc'] = 1 - data_out.loc[data_out["DNAm"] == "hypo", 'auc']

    data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)