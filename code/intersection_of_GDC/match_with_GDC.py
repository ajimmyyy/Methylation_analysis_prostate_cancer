import pandas as pd
import numpy as np
from tqdm import tqdm

fn_dmp = "Data/Data-cutpoint/dmp_with_cutpoint.csv"
GDC_path = "Data/Data-GDC/matchgene185_single_3Y10__OR2.txt"
fn_o = "Data/Data-GDC/DMP_filter_with_GDC_single.csv"

with open(GDC_path, 'r') as file:
    lines = file.readlines()

GDC = [line.strip() for line in lines]

data_dmp_df = pd.read_csv(fn_dmp)

result_df = data_dmp_df[data_dmp_df["gene"].isin(GDC)]
result_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)



