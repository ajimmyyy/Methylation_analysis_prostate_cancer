import csv
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from dbate_dif_distinguish import paint_distribute_picture

# paint_boxplot_picture(df_hyper, df_hypo)
# 依照hyper,hypo的dbate different和畫箱型圖
# Parameters:
# df_hyper: DataFrame
# df_hypo: DataFrame
# Return:
# null
# *圖儲存於plt畫布，可使用plt.savefig, plt.show儲存,顯示
# *請確保df_hyper, df_hypo擁有dbate_dif值
def paint_boxplot_picture(df_hyper, df_hypo):
    plt.boxplot([df_hyper["dbate_dif"], df_hypo["dbate_dif"]], labels = ["hyper", "hypo"], sym = "", vert = False, patch_artist = True)
    plt.xlabel('dbate_dif')

if __name__ == "__main__":
    fn_i = "Data/110820004/Data-cancer_stage/all_data/cancer_stage_with_dbate.csv"
    fn_hyperhypo = "Data/110820004/Data-cutpoint/HyperHypo_cutpoint_validate.csv"
    fn_o = "Data/110820004/Data-cancer_stage/hyperhypo_data_only/cancer_stage_with_HyperHypo.csv"
    fn_pic_hyper = "Data/110820004/Data-cancer_stage/hyperhypo_data_only/cancer_stage_dbate_dif_hyper_bar.png"
    fn_pic_hypo = "Data/110820004/Data-cancer_stage/hyperhypo_data_only/cancer_stage_dbate_dif_hypo_bar.png"
    fn_pic_boxplot = "Data/110820004/Data-cancer_stage/hyperhypo_data_only/cancer_stage_dbate_dif_hypo_boxplot.png"
    width = 0.001
    boundary = 0.1

    data_df = pd.read_csv(fn_i)
    data_hyperhypo = pd.read_csv(fn_hyperhypo)
    merged_df = pd.merge(data_hyperhypo, data_df, on='CpG', how='inner')

    hyper_df = merged_df[merged_df["DNAm"] == "hyper"]
    hypo_df = merged_df[merged_df["DNAm"] == "hypo"]

    paint_distribute_picture(hyper_df, width, boundary)
    plt.savefig(fn_pic_hyper)
    plt.show()

    paint_distribute_picture(hypo_df, width, boundary)
    plt.savefig(fn_pic_hypo)
    plt.show()

    paint_boxplot_picture(hyper_df, hypo_df)
    plt.savefig(fn_pic_boxplot)
    plt.show()

    merged_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
