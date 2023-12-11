import csv
import pandas as pd
import numpy as np
from tqdm import tqdm
def cal_accuracy(df_beta, df_cpg, normal_num):
    test_beta = df_beta[df_beta['CpG'].isin(df_cpg['CpG'])]

    df = pd.merge(df_cpg, test_beta, on = ['CpG'], how = 'inner')

    beta_columns = df.columns[2:]

    label_list = []
    for beta_column in beta_columns:
        is_tumor = df[beta_column].gt(df['cutpoint'], axis=0)
        label_sum = is_tumor.astype(int).sum(axis=0)
        label_list.append(label_sum > int(len(df.index) / 2))

    normal_count = 0
    tumor_count = 0
    for i in range(len(label_list)):
        if i <= normal_num and (not label_list[i]):
            normal_count += 1
        elif i > normal_num and label_list[i]:
            tumor_count += 1
    
    print("normal: ", normal_count / normal_num)
    print("tumor: ", tumor_count / (len(label_list) - normal_num))
    print("total: ", (normal_count +  tumor_count)/ len(label_list))

normal_num = 10
fn_beta = "C:/Users/acer/Desktop/Data-origin/validation/all_beta_normalized_validation.csv"
df_beta = pd.read_csv(fn_beta)
df_beta.columns.values[0] = 'CpG'
df_beta = df_beta.iloc[:, 0::2]

# 準備資料
fn_test_cpg = "C:/Users/acer/Desktop/result/goterm_group_max_f1.csv"
df_cpg = pd.read_csv(fn_test_cpg, usecols=["CpG", "cutpoint"])

# 測試準確度
cal_accuracy(df_beta, df_cpg, normal_num)