import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, auc
from config_loader import config

fn_bate = config["ALL_BETA_NORMALIZED_TRAIN_PATH"]
# fn_dmp = config["HYPER_HYPO_COMORBIDITY_SINGLE_PATH"]
fn_dmp = config["HYPER_HYPO_COMORBIDITY_GROUP_PATH"]
# fn_o = config["COMORBIDITY_SINGLE_AUC_PATH"]
fn_o = config["COMORBIDITY_GROUP_AUC_PATH"]
# fn_o_pic = config["COMORBIDITY_SINGLE_ROC_PICTURE_PATH"]
fn_o_pic = config["COMORBIDITY_GROUP_ROC_PICTURE_PATH"]

normal_num = 50
half_total_num = 278
half_normal_num = 25

data_bate_df = pd.read_csv(fn_bate)
data_dmp_df = pd.read_csv(fn_dmp)
DMP_list = data_dmp_df.iloc[:, 0].tolist()
dmp_bate_df = data_bate_df.loc[data_bate_df[data_bate_df.columns[0]].isin(DMP_list)]

def cal_cutpoint(row):
    transform_row = row.to_numpy()
    TPRs = []
    FPRs = []

    predict_beta = transform_row[1::2]
    actual_beta = np.ones(half_total_num)
    actual_beta[:half_normal_num] = 0

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
