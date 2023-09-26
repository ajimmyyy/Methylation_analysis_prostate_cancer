import csv
import pandas as pd
import numpy as np
from tqdm import tqdm

normal_num = 50
fn_i = "Data/Data-origin/all_beta_normalized.csv"
fn_o = "Data/Data-dbatea/ll_beta_delta.csv"

def cal_dbate(row):
    transform_row = row.to_numpy()
    normal_beta = transform_row[1:normal_num + 1:2]
    tumor_beta = transform_row[normal_num + 1::2]

    IQR_n = np.percentile(normal_beta,75) - np.percentile(normal_beta,25)
    Q3_n = np.percentile(normal_beta, 75)
    Q1_n = np.percentile(normal_beta, 25)
    normal_beta_filter = normal_beta[(normal_beta <= Q3_n + 1.5 * IQR_n) & (normal_beta >= Q1_n - 1.5 * IQR_n)]

    IQR_t = np.percentile(tumor_beta,75) - np.percentile(tumor_beta,25)
    Q3_t = np.percentile(tumor_beta, 75)
    Q1_t = np.percentile(tumor_beta, 25)
    tumor_beta_filter = tumor_beta[(tumor_beta <= Q3_t + 1.5 * IQR_t) & (tumor_beta >= Q1_t - 1.5 * IQR_t)]

    avg = np.mean(normal_beta_filter)
    tumor_beta_minus_avg = tumor_beta_filter - avg

    Q3_t_avg = np.percentile(tumor_beta_minus_avg, 75)
    Q1_t_avg = np.percentile(tumor_beta_minus_avg, 25)
    IQR_t_avg = Q3_t_avg - Q1_t_avg
    tumor_beta_filter_f = tumor_beta_minus_avg[(tumor_beta_minus_avg <= Q3_t_avg + 1.5 * IQR_t_avg) & (tumor_beta_minus_avg >= Q1_t_avg - 1.5 * IQR_t_avg)]

    dbate = np.mean(tumor_beta_filter_f)

    return pd.Series({"dbeta": dbate, "avg": avg})

data_df = pd.read_csv(fn_i)
data_out = data_df.iloc[:, 0].to_frame()
tqdm.pandas(desc="find dbeta")
data_out[["dbeta", "avg"]] = data_df.progress_apply(cal_dbate, axis = 1)

data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)