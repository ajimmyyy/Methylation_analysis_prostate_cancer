import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# Choose_hyper_hypo(dmp_result, threshold_FC_hyper, threshold_FC_hypo, threshold_Pvalue):
# 分離資料Hyper,Hypo
# Parameters:
# dmp_result: DataFrame，DMP資料
# threshold_FC_hyper: double，Hyper部分log10(FC)閥值
# threshold_FC_hypo: double，Hypo部分log10(FC)閥值
# threshold_Pvalue: double，Pvalue閥值
# Return:
# hyper, hypo: DataFrame，hyper,hypo資料
def Choose_hyper_hypo(dmp_result, threshold_FC_hyper, threshold_FC_hypo, threshold_Pvalue):
    dmp_result = dmp_result[dmp_result["adj.P.Val"] < threshold_Pvalue]

    hyper = dmp_result[dmp_result["logFC"] > threshold_FC_hyper]
    hypo = dmp_result[dmp_result["logFC"] < threshold_FC_hypo]

    return hyper, hypo

if __name__ == "__main__":
    fn_i = "Data/110820004/Data-dbate/DMP_per_gene.csv"
    fn_o = "Data/110820004/Data-volcano/DMP_hyper.csv"
    fn_o_hypo = "Data/110820004/Data-volcano/DMP_hypo.csv"
    pic_o = "Data/110820004/Data-volcano/DMP_volcano_plot.jpg"

    threshold_FC_hyper = 0.37
    threshold_FC_hypo = -0.258
    threshold_Pvalue = -np.log10(0.05)

    data_df = pd.read_csv(fn_i)

    uper_x, lower_x = Choose_hyper_hypo(data_df, threshold_FC_hyper, threshold_FC_hypo, threshold_Pvalue)

    plt.figure(figsize=(32,18))
    plt.style.use("ggplot")
    plt.xlabel("log(FC)", fontweight = "bold")
    plt.ylabel("-log(P.adj)", fontweight = "bold")

    plt.scatter((data_df["logFC"]), -np.log10(data_df["adj.P.Val"]), c = "grey", s = 15, alpha = .5, label='other')
    plt.axvline(threshold_FC_hyper,color="black", linestyle="--")
    plt.axvline(threshold_FC_hypo,color="black", linestyle="--")
    plt.scatter((uper_x["logFC"]), -np.log10(uper_x["adj.P.Val"]), c = "red", s = 15, alpha = .5, label='hyper')
    plt.scatter((lower_x["logFC"]), -np.log10(lower_x["adj.P.Val"]), c = "blue", s = 15, alpha = .5, label='hypo')
    plt.axhline(threshold_Pvalue,color="black", linestyle="--")

    TSS_x = uper_x[uper_x["feature"].str.startswith('TSS')]
    plt.scatter((TSS_x["logFC"]), -np.log10(TSS_x["adj.P.Val"]), c = "orange", s = 15, alpha = .5, label='TSS')
    plt.legend()

    plt.savefig(pic_o)
    plt.close()

    # uper_x.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
    # lower_x.to_csv(fn_o_hypo, sep=',', encoding='utf-8', index=False)
    print(len(uper_x))
    print(len(TSS_x[TSS_x["feature"] == 'TSS1500']))
    print(len(TSS_x[TSS_x["feature"] == 'TSS200']))
