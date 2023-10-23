import pandas as pd
import numpy as np
from tqdm import tqdm

# check_cutpoint(row, normal_num):
# 驗證每個位點的切點
# Parameters:
# row: DataFrame輸入列，透過apply輸入
# normal_num: Int，normal資料數量
# Return:
# F1_validate: DataFrame行資料
# *請透過apply or progress_apply使用此function
# *請確保資料擁有各病人之beta值
def check_cutpoint(row, normal_num):
    cutpoint = row["cutpoint"]
    filtered_row = data_beta_df[data_beta_df[data_beta_df.columns[0]] == row["CpG"]]

    transform_row = filtered_row.to_numpy()[0]
    normal_beta = transform_row[1:normal_num + 1:2]
    tumor_beta = transform_row[normal_num + 1::2]

    FP = np.sum(normal_beta > cutpoint)
    FN = np.sum(tumor_beta < cutpoint)
    TP = np.sum(tumor_beta > cutpoint)

    sensitivity = TP/(TP+FN)
    precision = TP/(TP+FP) if TP != 0 else 0

    F1 = 2 * sensitivity * precision / (sensitivity + precision)

    if row["DNAm"] == "hypo":
        F1 = 1 - F1

    return pd.Series({"F1_validate": F1})

if __name__ == "__main__":
    fn_cutpoint = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint.csv"
    fn_beta = "C:/Users/acer/Desktop/Data-origin/validation/all_beta_normalized_validation.csv"
    fn_o = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint_validate.csv"
    count_of_normal = 20
    threshold_validate_MAPE = 0.8

    data_beta_df = pd.read_csv(fn_beta)
    data_cutpoint_df = pd.read_csv(fn_cutpoint)

    tqdm.pandas(desc="find cutpoint")
    data_cutpoint_df["F1_validate"] = data_cutpoint_df.progress_apply(check_cutpoint, axis = 1, args = (count_of_normal,))
    data_cutpoint_df["F1_dif"] = data_cutpoint_df["F1_validate"] - data_cutpoint_df["F1"]

    data_cutpoint_df = data_cutpoint_df[data_cutpoint_df["F1_validate"] > threshold_validate_MAPE]

    data_cutpoint_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)



    