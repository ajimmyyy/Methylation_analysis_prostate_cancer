import csv
import re
import pandas as pd
import numpy as np
from tqdm import tqdm

def judge_cancer_stage(row):
    if row["Sample_Group"] == "Normal":
        return pd.Series({"cancer_stage": "None"})
    
    tumor_stage = str(row["tumor_stage"])
    node_stage = str(row["node_stage"])

    if  node_stage == "N0":
        if re.search(r"T1|T2", tumor_stage):
            return pd.Series({"cancer_stage": "early"})
        else:
            return pd.Series({"cancer_stage": "later"})
    else:
        return pd.Series({"cancer_stage": "later"})

if __name__ == "__main__":
    cancer_stage_map = "Data/110820004/Data-cancer_stage/cancer_stage.csv"
    sample_sheet_dir = "Data/110820004/Data-cancer_stage/sample_sheet.csv"
    fn_o = "Data/110820004/Data-cancer_stage/patient_cancer_stage_train.csv"

    cancer_stage_df = pd.read_csv(cancer_stage_map)
    sample_sheet = pd.read_csv(sample_sheet_dir)
    
    filtered_data = cancer_stage_df[cancer_stage_df["Sentrix_ID"].isin(sample_sheet["Sentrix_ID"]) & cancer_stage_df["Sentrix_Position"].isin(sample_sheet["Sentrix_Position"])]
    data_out = pd.merge(sample_sheet, filtered_data, on = ["Sentrix_ID", "Sentrix_Position"], how = "inner")

    tqdm.pandas(desc="find cancer stage")
    data_out["cancer_stage"] = data_out.progress_apply(judge_cancer_stage, axis = 1)

    data_out.to_csv(fn_o, sep=',', encoding='utf-8', index=False)   