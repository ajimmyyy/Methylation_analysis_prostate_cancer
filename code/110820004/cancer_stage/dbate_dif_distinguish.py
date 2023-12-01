import csv
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def paint_distribute_picture(all_df, bin_width, boundary):
    std_value = all_df["dbate_dif"].std()
    mean_value = all_df["dbate_dif"].mean()

    hist, bins = np.histogram(all_df["dbate_dif"], bins = np.arange(-boundary, boundary + bin_width, bin_width))
    
    colors = ['gray' if (x < (mean_value - 2 * std_value) or x > (mean_value + 2 * std_value)) else 'blue' for x in np.arange(-boundary, boundary + bin_width, bin_width)]

    plt.figure(figsize=(8, 6))
    plt.bar(bins[:-1], hist, width=bin_width, align='edge', color = colors)
    plt.xlabel('dbate_dif')
    plt.ylabel('frequency')

if __name__ == "__main__":
    fn_i = "Data/110820004/Data-cancer_stage/cancer_stage_with_dbate.csv"
    fn_pos_dif = "Data/110820004/Data-cancer_stage/cancer_stage_dbate_positive_dif.csv"
    fn_neg_dif = "Data/110820004/Data-cancer_stage/cancer_stage_dbate_negative_dif.csv"
    fn_between_dif = "Data/110820004/Data-cancer_stage/cancer_stage_dbate_between_dif.csv"
    fn_pic = "Data/110820004/Data-cancer_stage/cancer_stage_dbate_dif_bar.png"
    width = 0.001
    bound = 0.1

    data_df = pd.read_csv(fn_i)

    df_std = data_df["dbate_dif"].std()
    df_mean = data_df["dbate_dif"].mean()

    paint_distribute_picture(data_df, width, bound)
    plt.savefig(fn_pic)
    plt.show()

    positive_dif_df = data_df[data_df["dbate_dif"] > df_mean + (2 * df_std)]
    negative_dif_df = data_df[data_df["dbate_dif"] < df_mean - (2 * df_std)]
    between_dif_df = data_df[(data_df["dbate_dif"] > df_mean - (2 * df_std)) & (data_df["dbate_dif"] < df_mean + (2 * df_std))]

    positive_dif_df.to_csv(fn_pos_dif, sep=',', encoding='utf-8', index=False)
    negative_dif_df.to_csv(fn_neg_dif, sep=',', encoding='utf-8', index=False)
    between_dif_df.to_csv(fn_between_dif, sep=',', encoding='utf-8', index=False)

