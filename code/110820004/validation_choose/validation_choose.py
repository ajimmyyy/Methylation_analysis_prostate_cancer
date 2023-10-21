import pandas as pd
import numpy as np
from tqdm import tqdm
from config_loader import config

fn_origin = config["ALL_BETA_NORMALIZED_TEST_PATH"]
fn_validation = config["ALL_BETA_NORMALIZED_VALIDATION_PATH"]
fn_selected_columns = config["VALIDATION_CHOOSE_PATH"]

normal_total = 50
normal_num = 20
tumor_total = 500
tumor_num = 200

data_origin_df = pd.read_csv(fn_origin)

random_normal = np.random.choice(normal_total//2, normal_num // 2, replace=False) * 2
random_tumor = np.random.choice(tumor_total // 2, tumor_num // 2, replace=False) * 2

print(random_normal, len(random_normal))
print(random_tumor, len(random_tumor))

selected_columns = np.sort(np.concatenate([[0], random_normal + 1, random_normal + 2, random_tumor + normal_total + 1, random_tumor + normal_total + 2]))

print(selected_columns, len(selected_columns))

data = data_origin_df.iloc[:, selected_columns.tolist()]

np.savetxt(fn_selected_columns, selected_columns, fmt='%d')
data.to_csv(fn_validation, sep=',', encoding='utf-8', index=False)

