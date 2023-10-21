import pandas as pd
import numpy as np
from tqdm import tqdm
from config_loader import config

fn_dmp_group = config["HYPER_HYPO_COMORBIDITY_GROUP_PATH"]
fn_dmp_single = config["HYPER_HYPO_COMORBIDITY_SINGLE_PATH"]
fn_o_single = config["HYPER_HYPO_COMORBIDITY_GROUP_GROUP_PATH"]
fn_o_group = config["HYPER_HYPO_COMORBIDITY_GROUP_SOLO_PATH"]
fn_o_outer = config["HYPER_HYPO_COMORBIDITY_BOTH_PATH"]

data_dmp_df_group = pd.read_csv(fn_dmp_group)
data_dmp_df_single = pd.read_csv(fn_dmp_single)

only_in_group_df = data_dmp_df_group[~data_dmp_df_group.loc[:, "CpG"].isin(data_dmp_df_single.loc[:, "CpG"])]
only_in_single_df = data_dmp_df_single[~data_dmp_df_single.loc[:, "CpG"].isin(data_dmp_df_group.loc[:, "CpG"])]
outer_df = pd.merge(data_dmp_df_group,data_dmp_df_single, how = "outer")

only_in_group_df.to_csv(fn_o_group, sep=',', encoding='utf-8', index=False)
only_in_single_df.to_csv(fn_o_single, sep=',', encoding='utf-8', index=False)
outer_df.to_csv(fn_o_outer, sep=',', encoding='utf-8', index=False)

