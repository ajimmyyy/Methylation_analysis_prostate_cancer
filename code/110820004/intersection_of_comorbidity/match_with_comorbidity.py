import pandas as pd
import numpy as np
from tqdm import tqdm
from config_loader import config

fn_dmp = config["HYPER_HYPO_CUTPOINT_VALIDATION_PATH"]
comorbidity_path = config["COMORBIDITY_SINGLE_PATH"]
# comorbidity_path = config["COMORBIDITY_GROUP_PATH"]
fn_o = config["HYPER_HYPO_COMORBIDITY_SINGLE_PATH"]
# fn_o = config["HYPER_HYPO_COMORBIDITY_GROUP_PATH"]

with open(comorbidity_path, 'r') as file:
    lines = file.readlines()

comorbidity = [line.strip() for line in lines]

data_dmp_df = pd.read_csv(fn_dmp)

result_df = data_dmp_df[data_dmp_df["gene"].isin(comorbidity)]

result_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)



