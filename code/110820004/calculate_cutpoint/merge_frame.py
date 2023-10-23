import pandas as pd
import numpy as np
from tqdm import tqdm

# 
# 合併 hyper,hypo
# 
fn_dmp_hyper = "Data/110820004/Data-cutpoint/hyper_with_cutpoint.csv"
fn_dmp_hypo = "Data/110820004/Data-cutpoint/hypo_with_cutpoint.csv"
fn_o = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint.csv"

data_dmp_hyper = pd.read_csv(fn_dmp_hyper)
data_dmp_hypo = pd.read_csv(fn_dmp_hypo)

merged_df = pd.concat([data_dmp_hyper, data_dmp_hypo])
merged_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)

