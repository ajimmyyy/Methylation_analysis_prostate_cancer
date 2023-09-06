import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

fn_i = "DMP_per_gene.csv"
fn_o = "DMP_hyper.csv"
fn_o_hypo = "DMP_hypo.csv"
pic_o = "DMP_with_volcano.jpg"

threshold_FC = 0.37
threshold_Pvalue = -np.log10(0.05)

data_df = pd.read_csv(fn_i)

plt.figure(figsize=(32,18))
plt.style.use("ggplot")
plt.xlabel("log(FC)", fontweight = "bold")
plt.ylabel("-log(P.adj)", fontweight = "bold")

plt.scatter((data_df["logFC"]), -np.log10(data_df["adj.P.Val"]), c = "grey", s = 15, alpha = .5, label='other')

uper_x = data_df[data_df["logFC"] > threshold_FC]
lower_x = data_df[data_df["logFC"] < -threshold_FC]

plt.axvline(threshold_FC,color="black", linestyle="--")
plt.axvline(-threshold_FC,color="black", linestyle="--")
plt.scatter((uper_x["logFC"]), -np.log10(uper_x["adj.P.Val"]), c = "red", s = 15, alpha = .5, label='hyper')
plt.scatter((lower_x["logFC"]), -np.log10(lower_x["adj.P.Val"]), c = "blue", s = 15, alpha = .5, label='hypo')

plt.axhline(threshold_Pvalue,color="black", linestyle="--")

TSS_x = uper_x[uper_x["feature"].str.startswith('TSS')]

plt.scatter((TSS_x["logFC"]), -np.log10(TSS_x["adj.P.Val"]), c = "orange", s = 15, alpha = .5, label='TSS')

plt.legend()

plt.savefig(pic_o)
plt.close()
uper_x.to_csv(fn_o, sep=',', encoding='utf-8', index=False)
lower_x.to_csv(fn_o_hypo, sep=',', encoding='utf-8', index=False)
print(len(uper_x))
print(len(TSS_x[TSS_x["feature"] == 'TSS1500']))
print(len(TSS_x[TSS_x["feature"] == 'TSS200']))
