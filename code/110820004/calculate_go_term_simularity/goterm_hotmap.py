import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy

if __name__ == "__main__":
    fn_gosemsim = "Data/110820004/Data-gosemsim/mgeneSim_hyper.csv"
    fn_o = "Data/110820004/Data-gosemsim/goterm_hyper_hotmap.png"
    fn_cluster = "Data/110820004/Data-gosemsim/goterm_hyper_cluster.png"
    fn_group = "Data/110820004/Data-gosemsim/goterm_group_hyper.csv"
    threshold_hyper = 2.9
    threshold_hypo = 2.1
    
    df_go = pd.read_csv(fn_gosemsim)
    attributes = df_go.columns.tolist()
    df_go.index = attributes

    row_linkage = hierarchy.linkage(df_go, method='average')

    fig = plt.figure(figsize=(20, 20))
    dn = hierarchy.dendrogram(row_linkage, labels=attributes)
    plt.savefig(fn_cluster)
    plt.cla()

    clusters = hierarchy.fcluster(row_linkage, t=threshold_hyper, criterion='distance')
    data = {
        "gene": attributes,
        "cluster": clusters
    }
    df = pd.DataFrame(data)
    df.to_csv(fn_group, sep=',', encoding='utf-8', index=False)

    hotmap = sns.clustermap(df_go, method='average', annot=True, figsize=(200, 200))
    hotmap.savefig(fn_o) 

