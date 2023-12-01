import csv
import pandas as pd
import numpy as np
from tqdm import tqdm
from calculate_dbate.cal_dbate import cal_dbate
import os

# 
# 計算癌症前後期dbate
# *由於調用上級package，無法直接執行
# *在上級新增main.py，import此檔案，並執行
# 
count_of_normal = 50
fn_i = "C:/Users/acer/Desktop/Data-origin/train/all_beta_normalized.csv"
fn_early = "Data/110820004/Data-cancer_stage/early_patient_list.txt"
fn_later = "Data/110820004/Data-cancer_stage/later_patient_list.txt"
fn_o = "Data/110820004/Data-cancer_stage/cancer_stage_with_dbate.csv"

data_df = pd.read_csv(fn_i)
with open(fn_early, "r") as file:
    early_list = [int(line.strip()) for line in file]
with open(fn_later, "r") as file:
    later_list = [int(line.strip()) for line in file]

early_list = [i for i in range(51)] + early_list
later_list = [i for i in range(51)] + later_list

early_df = data_df.iloc[:, early_list]
later_df = data_df.iloc[:, later_list]

early_dbate_df = early_df.iloc[:, 0].to_frame()
tqdm.pandas(desc="find dbeta")
reault = early_df.progress_apply(cal_dbate, axis = 1, args = (count_of_normal,))
early_dbate_df["early_dbeta"] = reault["dbeta"]
early_dbate_df.columns.values[0] = 'CpG'

later_dbate_df = later_df.iloc[:, 0].to_frame()
tqdm.pandas(desc="find dbeta")
reault = later_df.progress_apply(cal_dbate, axis = 1, args = (count_of_normal,))
later_dbate_df["later_dbeta"] = reault["dbeta"]
later_dbate_df.columns.values[0] = 'CpG'

data_out = pd.merge(early_dbate_df, later_dbate_df, on = ["CpG"], how = "inner")
data_out["dbate_dif"] = data_out["early_dbeta"] - data_out["later_dbeta"]
data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)