import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    fn_gosemsim = "Data/110820004/Data-gosemsim/mgeneSim_hypo.csv"
    fn_o = "Data/110820004/Data-gosemsim/goterm_hypo_hotmap.png"
    
    df_go = pd.read_csv(fn_gosemsim)
    attributes = df_go.columns.tolist()
    df_go.index = attributes

    hotmap = sns.clustermap(df_go, method='average', annot=True, figsize=(200, 200))
    hotmap.savefig(fn_o) 

