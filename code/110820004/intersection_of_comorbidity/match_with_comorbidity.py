import pandas as pd
import numpy as np
from tqdm import tqdm

fn_dmp = "Data/110820004/Data-Cutpoint/HyperHypo_cutpoint_validate.csv"
# comorbidity_path = "Data/110820004/Data-Comorbidity/matchgene185_single_3Y10__OR2.txt"
comorbidity_path = "Data/110820004/Data-Comorbidity/matchgene185_group_3Y10__OR2.txt"
# fn_o = "Data/110820004/Data-Comorbidity/HyperHypo_filtered_comorbidity_single.csv"
fn_o = "Data/110820004/Data-Comorbidity/HyperHypo_filtered_comorbidity_group.csv"

with open(comorbidity_path, 'r') as file:
    lines = file.readlines()

comorbidity = [line.strip() for line in lines]

data_dmp_df = pd.read_csv(fn_dmp)

result_df = data_dmp_df[data_dmp_df["gene"].isin(comorbidity)]
result_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)



