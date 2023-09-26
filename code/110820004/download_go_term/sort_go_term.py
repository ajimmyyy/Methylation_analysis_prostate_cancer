import pandas as pd
import numpy as np
from tqdm import tqdm

fn_dmp = "GO_term/DMP_hyper_with_GO_all.csv"
fn_o = "GO_term/DMP_hyper_with_GO_all_sort.csv"

data_bate_df = pd.read_csv(fn_dmp)

data_bate_df['Comma_Count'] = data_bate_df['BP'].str.count(',')
df_sorted = data_bate_df.sort_values(by = 'Comma_Count', ascending = False)
df_sorted.drop(columns = ['Comma_Count'], inplace = True)

df_sorted.to_csv(fn_o, sep=',', encoding='utf-8', index=False)