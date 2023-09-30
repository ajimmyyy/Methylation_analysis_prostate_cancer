import pandas as pd
import numpy as np
from tqdm import tqdm

fn_dmp = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint_validate.csv"
GDC_path = "Data/110820004/Data-GDC/matchgene185_single_3Y10__OR2.txt"
# GDC_path = "Data/110820004/Data-GDC/matchgene185_group_3Y10__OR2.txt"
fn_o = "Data/110820004/Data-GDC/HyperHypo_filtered_GDC_single.csv"
# fn_o = "Data/110820004/Data-GDC/HyperHypo_filtered_GDC_group.csv"

with open(GDC_path, 'r') as file:
    lines = file.readlines()

GDC = [line.strip() for line in lines]

data_dmp_df = pd.read_csv(fn_dmp)

result_df = data_dmp_df[data_dmp_df["gene"].isin(GDC)]
result_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)



