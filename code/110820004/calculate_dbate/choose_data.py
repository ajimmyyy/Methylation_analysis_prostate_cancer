import csv
import pandas as pd
import numpy as np
import swifter
from tqdm import tqdm

# choose_dbate_per_gene(dmp_result, all_dbate)
# 選擇個基因裡dbate值最大的
# Parameters:
# dmp_result: DataFrame，DMP資料
# all_dbate: DataFrame，擁有dbate行的資料
# Return:
# dmp_result: DataFrame，dmp_result每個gene留下的資料
def choose_dbate_per_gene(dmp_result, all_dbate):
    finded_gene = []
    del_col_n = []

    def find_dbate(row):
        filtered_rows = all_dbate[all_dbate[all_dbate.columns[0]] == row.iat[0]]
        return filtered_rows['dbeta'].values[0]

    dmp_result["dbeta"] = dmp_result.swifter.apply(find_dbate, axis = 1)

    filtered_dmp = dmp_result[dmp_result['feature'].isin(['TSS200', 'TSS1500'])]
    max_indices = filtered_dmp.groupby('gene', group_keys=False)['dbeta'].apply(lambda x: x.abs().idxmax())
    dmp_result = filtered_dmp.loc[max_indices]

    return dmp_result

if __name__ == "__main__":
    fn_dmp = "C:/Users/acer/Desktop/Data-origin/train/DMP_result_TN.csv"
    fn_dbate = "Data/110820004/Data-dbate/all_beta_delta.csv"
    fn_o = "Data/110820004/Data-dbate/DMP_per_gene.csv"

    data_dmp_df = pd.read_csv(fn_dmp)
    data_dbate_df = pd.read_csv(fn_dbate)

    data_out = choose_dbate_per_gene(data_dmp_df, data_dbate_df)

    data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
