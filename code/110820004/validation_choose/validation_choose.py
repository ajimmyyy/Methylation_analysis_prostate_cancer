import pandas as pd
import numpy as np
from tqdm import tqdm

def random_choose_validation(fn_origin, normal_total, normal_num, tumor_total, tumor_num):
    data_origin_df = pd.read_csv(fn_origin)

    random_normal = np.random.choice(normal_total//2, normal_num // 2, replace=False) * 2
    random_tumor = np.random.choice(tumor_total // 2, tumor_num // 2, replace=False) * 2

    selected_columns = np.sort(np.concatenate([[0], random_normal + 1, random_normal + 2, random_tumor + normal_total + 1, random_tumor + normal_total + 2]))

    data = data_origin_df.iloc[:, selected_columns.tolist()]

    return data, selected_columns

if __name__ == "__main__":
    fn_origin = "C:/Users/acer/Desktop/Data-origin/test/all_beta_normalized_test.csv"
    fn_validation = "C:/Users/acer/Desktop/Data-origin/validation/all_beta_normalized_validation.csv"
    fn_selected_columns = "C:/Users/acer/Desktop/Data-origin/validation/validation_choose.txt"

    normal_total = 50
    normal_num = 20
    tumor_total = 500
    tumor_num = 200
    
    data, selected_columns = random_choose_validation(fn_origin, normal_total, normal_num, tumor_total, tumor_num)

    np.savetxt(fn_selected_columns, selected_columns, fmt='%d')
    data.to_csv(fn_validation, sep=',', encoding='utf-8', index=False)

