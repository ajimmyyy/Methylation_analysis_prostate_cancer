import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, auc

fn_bate = "C:/Users/acer/Desktop/Data-origin/train/all_beta_normalized.csv"
fn_dmp = "Data/110820004/Data-Comorbidity/HyperHypo_filtered_comorbidity_single.csv"
# fn_dmp = "Data/110820004/Data-Comorbidity/HyperHypo_filtered_comorbidity_group.csv"
fn_o = "Data/110820004/Data-ROC_AUC/comorbidity_single_auc.csv"
# fn_o = "Data/110820004/Data-ROC_AUC/comorbidity_group_auc.csv"
fn_o_pic = "Data/110820004/Data-ROC_AUC/ROC_curve_single.png"
# fn_o_pic = "Data/110820004/Data-ROC_AUC/ROC_curve_group.png"

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

    if auc1 > 0.5:
        plt.plot(fpr, tpr, color = 'gainsboro', alpha = .04)
    else:
        plt.plot(fpr, tpr, color = 'blue', label = 'AUC = %0.2f' % auc1)

    if transform_row[0] == "cg06263461":
        plt.plot(fpr, tpr, color = 'red', label = 'Max AUC = %0.2f' % auc1)
    if transform_row[0] == "cg21634602":
        plt.plot(fpr, tpr, color = 'orange', label = 'Min AUC = %0.2f' % auc1)

    return pd.Series({"auc": auc1})

data_out = dmp_bate_df.iloc[:, 0].to_frame()

tqdm.pandas(desc="find auc")
data_out["auc"] = dmp_bate_df.progress_apply(cal_cutpoint, axis = 1)

plt.legend(loc = 'lower right')
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.savefig(fn_o_pic)
plt.show()

data_out.columns.values[0] = 'CpG'
data_out = pd.merge(data_dmp_df, data_out, on = ["CpG"], how = "inner")
data_out.loc[data_out["DNAm"] == "hypo", 'auc'] = 1 - data_out.loc[data_out["DNAm"] == "hypo", 'auc']

data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
