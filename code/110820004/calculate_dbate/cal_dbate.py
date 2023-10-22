import csv
import pandas as pd
import numpy as np
from tqdm import tqdm

# cal_dbate(row, normal_num)
# 計算各位點的delta beta值
# Parameters:
# row: DataFrame輸入列，透過apply輸入
# normal_num: Int，normal資料數量
# Return:
# dbate, avg: DataFrame行資料
# *請透過apply or progress_apply使用此function
# *請確保資料擁有各病人之beta值
def cal_dbate(row, normal_num):
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

if __name__ == "__main__":
    count_of_normal = 50
    fn_i = "C:/Users/acer/Desktop/Data-origin/train/all_beta_normalized.csv"
    fn_o = "Data/110820004/Data-dbate/all_beta_delta.csv"

    data_df = pd.read_csv(fn_i)
    data_out = data_df.iloc[:, 0].to_frame()
    tqdm.pandas(desc="find dbeta")
    data_out[["dbeta", "avg"]] = data_df.progress_apply(cal_dbate, axis = 1, args = (count_of_normal,))

    data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)