import pandas as pd
import numpy as np
from tqdm import tqdm

# 
# 與參考研究基因表交集
# 
fn_dmp = "Data/110820004/Data-ROC_AUC/comorbidity_group_auc.csv"
study_dmp = "Data/110820004/Data-study/study_count_gene.csv"
fn_o = "Data/110820004/Data-study/comorbidity_group_auc_study_count.csv"

data_dmp_df = pd.read_csv(fn_dmp)
data_study_df = pd.read_csv(study_dmp)

# 交集(保留完整原始位點，即使不在參考研究基因表)
filtered_data_study = data_study_df[data_study_df['gene'].isin(data_dmp_df['gene'])]
data_out = pd.merge(data_dmp_df, filtered_data_study, on = ["gene"], how = "outer")

data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
