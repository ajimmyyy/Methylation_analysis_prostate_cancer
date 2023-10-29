import csv
import pandas as pd
import numpy as np
from tqdm import tqdm

# 
# 標註每位病患的TCGA_barcode,SentrixID,TNM
# 
if __name__ == "__main__":
    TCGA_dir = "C:/Users/acer/Desktop/special_subject/dataset/jhu-usc.edu_PRAD.HumanMethylation450.1.13.0.sdrf.txt"
    cancer_stage_dir = "Data/110820004/Data-cancer_stage/prad_tcga_pan_can_atlas_2018_clinical_data.csv"
    fn_o_cancer_stage_map = "Data/110820004/Data-cancer_stage/cancer_stage.csv"

    TCGA_list = pd.read_csv(TCGA_dir, sep = "\t")
    cancer_stage_df = pd.read_csv(cancer_stage_dir)

    ID_map = TCGA_list[["Comment [TCGA Barcode]", "Array Data File"]]
    cancer_stage_map = cancer_stage_df[["Patient ID",
                                        "American Joint Committee on Cancer Metastasis Stage Code",
                                        "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code", 
                                        "American Joint Committee on Cancer Tumor Stage Code"]]
    
    ID_map = ID_map.rename(columns = {"Comment [TCGA Barcode]": "TCGA_barcode", "Array Data File": "file_ID"})
    cancer_stage_map = cancer_stage_map.rename(columns = {"Patient ID": "TCGA_barcode", 
                                                          "American Joint Committee on Cancer Metastasis Stage Code": "metastasis_stage", 
                                                          "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": "node_stage", 
                                                          "American Joint Committee on Cancer Tumor Stage Code": "tumor_stage"})

    ID_map["TCGA_barcode"] = ID_map["TCGA_barcode"].apply(lambda x: x[:12])
    ID_map['Sentrix_ID'] = ID_map["file_ID"].str[:10]
    ID_map['Sentrix_Position'] = ID_map["file_ID"].str[11:17]
    ID_map = ID_map.drop(columns = ["file_ID"])

    filtered_data = cancer_stage_map[cancer_stage_map['TCGA_barcode'].isin(ID_map['TCGA_barcode'])]
    data_out = pd.merge(ID_map, filtered_data, on = ["TCGA_barcode"], how = "outer")
    
    data_out.to_csv(fn_o_cancer_stage_map, sep=',', encoding='utf-8', index=False)