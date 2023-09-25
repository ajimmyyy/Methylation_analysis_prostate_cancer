import pandas as pd
import numpy as np
from tqdm import tqdm

fn_dmp_hyper = "Data/Data-cutpoint/hyper_with_cutpoint.csv"
fn_dmp_hypo = "Data/Data-cutpoint/hypo_with_cutpoint.csv"
fn_o = "Data/Data-cutpoint/HyperHypo_cutpoint.csv"

data_dmp_hyper = pd.read_csv(fn_dmp_hyper)
data_dmp_hypo = pd.read_csv(fn_dmp_hypo)

data_dmp_hyper.loc[:, "DNAm"] = "hyper"
data_dmp_hypo.loc[:, "DNAm"] = "hypo"

merged_df = pd.concat([data_dmp_hyper, data_dmp_hypo])
merged_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)

