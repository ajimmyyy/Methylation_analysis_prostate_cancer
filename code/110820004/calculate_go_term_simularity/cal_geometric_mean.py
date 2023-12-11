import pandas as pd
import numpy as np
from scipy.stats import gmean

fn_gene_list = "Data/110820004/Data-ROC_AUC/comorbidity_group_auc.csv"
fn_MF = "Data/110820004/Data-gosemsim/MF/mgeneSim_hyper.csv"
fn_CC = "Data/110820004/Data-gosemsim/CC/mgeneSim_hyper.csv"
fn_BP = "Data/110820004/Data-gosemsim/BP/mgeneSim_hyper.csv"
fn_o = "Data/110820004/Data-gosemsim/geometric_mean/mgeneSim_geometric_mean_hyper.csv"

df_gene = pd.read_csv(fn_gene_list, usecols=['gene', 'DNAm'])
gene_list = df_gene[df_gene["DNAm"] == "hyper"]['gene'].to_list()

df_MF = pd.read_csv(fn_MF, index_col=0)
df_CC = pd.read_csv(fn_CC, index_col=0)
df_BP = pd.read_csv(fn_BP, index_col=0)

df_MF = df_MF.reindex(index = gene_list, columns = gene_list, fill_value=0)
df_CC = df_CC.reindex(index = gene_list, columns = gene_list, fill_value=0)
df_BP = df_BP.reindex(index = gene_list, columns = gene_list, fill_value=0)

df_geometric_mean = pd.DataFrame(index=gene_list, columns=gene_list)
for i in gene_list:
    for j in gene_list:
        values = [value for value in [df_MF.loc[i, j], df_CC.loc[i, j], df_BP.loc[i, j]] if value != 0]
        if values:
            geometric_mean = gmean(values)
            df_geometric_mean.loc[i, j] = geometric_mean
        else:
            df_geometric_mean.loc[i, j] = 0

df_geometric_mean.to_csv(fn_o, sep=',', encoding='utf-8')