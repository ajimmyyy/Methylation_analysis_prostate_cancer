import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

fn_choose_list = "C:/Users/acer/Desktop/tmp/tmp.csv"
fn_dbate = "Data/110820004/Data-cancer_stage/cancer_stage_with_dbate.csv"

data_list = pd.read_csv(fn_choose_list)
choose_list = data_list[data_list.columns[0]].to_numpy()
data_dbeta = pd.read_csv(fn_dbate)

df_out = data_dbeta[data_dbeta["CpG"].isin(choose_list)]
df_out.to_csv("C:/Users/acer/Desktop/tmp/tmp_1.csv", sep=',', encoding='utf-8', index=False)
