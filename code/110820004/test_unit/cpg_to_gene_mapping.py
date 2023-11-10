import csv
import pandas as pd
import numpy as np
from tqdm import tqdm

if __name__ == "__main__":
    # test_gene = ["GAPDH", "ACTB", "MYH14;MYH14;MYH14", "PRHOXNB"]
    test_gene = ["GSTP1"]
    CpG_site_o = []

    detail_450k_dir = "C:/Users/acer/Desktop/Data-origin/humanmethylation450.csv"
    fn_o = "C:/Users/acer/Desktop/tmp/tmp.csv"

    detail_450k_df = pd.read_csv(detail_450k_dir)
    filtered_data = detail_450k_df[detail_450k_df["UCSC_RefGene_Name"].isin(test_gene)]

    filtered_data = filtered_data[["IlmnID", "UCSC_RefGene_Name"]]

    filtered_data.to_csv(fn_o, sep=',', encoding='utf-8', index=False)