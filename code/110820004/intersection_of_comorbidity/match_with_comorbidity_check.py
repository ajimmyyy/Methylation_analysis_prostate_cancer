import pandas as pd
import numpy as np
from tqdm import tqdm

# 
# 共病基因single, group差集,連集
# 
fn_dmp_group = "Data/110820004/Data-comorbidity/HyperHypo_filtered_comorbidity_group.csv"
fn_dmp_single = "Data/110820004/Data-comorbidity/HyperHypo_filtered_comorbidity_single.csv"
fn_o_single = "Data/110820004/Data-comorbidity/HyperHypo_filtered_comorbidity_single_sole.csv"
fn_o_group = "Data/110820004/Data-comorbidity/HyperHypo_filtered_comorbidity_group_sole.csv"
fn_o_outer = "Data/110820004/Data-comorbidity/HyperHypo_filtered_comorbidity_both_single_group.csv"

data_dmp_df_group = pd.read_csv(fn_dmp_group)
data_dmp_df_single = pd.read_csv(fn_dmp_single)

only_in_group_df = data_dmp_df_group[~data_dmp_df_group.loc[:, "CpG"].isin(data_dmp_df_single.loc[:, "CpG"])]
only_in_single_df = data_dmp_df_single[~data_dmp_df_single.loc[:, "CpG"].isin(data_dmp_df_group.loc[:, "CpG"])]
outer_df = pd.merge(data_dmp_df_group,data_dmp_df_single, how = "outer")

only_in_group_df.to_csv(fn_o_group, sep=',', encoding='utf-8', index=False)
only_in_single_df.to_csv(fn_o_single, sep=',', encoding='utf-8', index=False)
outer_df.to_csv(fn_o_outer, sep=',', encoding='utf-8', index=False)

