import pandas as pd
import numpy as np
from tqdm import tqdm
from config_loader import config

fn_dmp = config["COMORBIDITY_GROUP_AUC_PATH"]
study_dmp = config["STUDY_COUNT_GENE_PATH"]
fn_o = config["COMORBIDITY_GROUP_AUC_STUDY_COUNT_PATH"]

data_dmp_df = pd.read_csv(fn_dmp)
data_study_df = pd.read_csv(study_dmp)

filtered_data_study = data_study_df[data_study_df['gene'].isin(data_dmp_df['gene'])]
data_out = pd.merge(data_dmp_df, filtered_data_study, on = ["gene"], how = "outer")

data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
