import pandas as pd
import numpy as np
from tqdm import tqdm
from config_loader import config

<<<<<<< Updated upstream
fn_dmp_hyper = "Data/110820004/Data-cutpoint/hyper_with_cutpoint.csv"
fn_dmp_hypo = "Data/110820004/Data-cutpoint/hypo_with_cutpoint.csv"
fn_o = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint.csv"
=======
fn_dmp_hyper = config["HYPER_CUTPOINT_PATH"]
fn_dmp_hypo = config["HYPO_CUTPOINT_PATH"]
fn_o = config["HYPER_HYPO_CUTPOINT_PATH"]
>>>>>>> Stashed changes

data_dmp_hyper = pd.read_csv(fn_dmp_hyper)
data_dmp_hypo = pd.read_csv(fn_dmp_hypo)

data_dmp_hyper.loc[:, "DNAm"] = "hyper"
data_dmp_hypo.loc[:, "DNAm"] = "hypo"

merged_df = pd.concat([data_dmp_hyper, data_dmp_hypo])
merged_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)

