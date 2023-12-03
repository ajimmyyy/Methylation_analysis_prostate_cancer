import pandas as pd

fn_origin = "Data/110820004/Data-ROC_AUC/comorbidity_group_auc.csv"
fn_go = "Data/110820004/Data-gosemsim/mgeneSim_results.csv"
fn_miss = "Data/110820004/Data-gosemsim/gene_miss.txt"

df_origin = pd.read_csv(fn_origin, usecols=['gene'])
df_go = pd.read_csv(fn_go)

list_origin = df_origin['gene'].tolist()
list_go = df_go.columns.to_list()

genes_miss = [gene for gene in list_origin if gene not in list_go]

with open(fn_miss, 'w') as file:
    for gene in genes_miss:
        file.write(gene + '\n')