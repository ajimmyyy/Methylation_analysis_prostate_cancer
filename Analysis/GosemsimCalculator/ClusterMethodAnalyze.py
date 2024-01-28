from configparser import ConfigParser
import pandas as pd
import numpy as np
from MakeFile.FileSaver import FileSaver
import matplotlib.pyplot as plt
import seaborn as sns

def GetClusterList(clusterDf):
    cluster = []
    for clusterValue in clusterDf['cluster'].unique():
        gene = clusterDf[clusterDf['cluster'] == clusterValue]['gene'].tolist()
        cluster.append(gene)
    return cluster

def JaccardSimilarity(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union != 0 else 0


def ReorderColumns(df):
    for i in range(len(df.index)):
        max_similarity = df.iloc[i, i]
        max_index = i
        for j in range(i + 1, len(df.index)):
            if df.iloc[i, j] > max_similarity:
                max_similarity = df.iloc[i, j]
                max_index = j
        if max_index != i:
            df.iloc[:, [i, max_index]] = df.iloc[:, [max_index, i]]

    return df

def ClusterSimilarity(clusters1, clusters2):
    data = []
    for cluster1 in clusters1:
        row = []
        for cluster2 in clusters2:
            similarity = JaccardSimilarity(cluster1, cluster2)
            row.append(similarity)
        data.append(row)

    return ReorderColumns(pd.DataFrame(data))

if __name__ == "__main__":
    _configPath = "Analysis/GosemsimCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _wardDf = pd.read_csv(_config["Paths"]["GENE_CLUSTER_MF_HYPER_HIERARCY_WARD_PATH"])
    _averageDf = pd.read_csv(_config["Paths"]["GENE_CLUSTER_MF_HYPER_HIERARCY_AVERAGE_PATH"])
    _kmedoidsDf = pd.read_csv(_config["Paths"]["GENE_CLUSTER_MF_HYPER_KMEDOIDS_PATH"])

    _wardCluster = GetClusterList(_wardDf)
    _averageCluster = GetClusterList(_averageDf)
    _kmedoidsCluster = GetClusterList(_kmedoidsDf)

    _wardKmedoidsSimilarity = ClusterSimilarity(_kmedoidsCluster, _wardCluster)
    _wardAverageSimilarity = ClusterSimilarity(_averageCluster, _wardCluster)

    ax = sns.heatmap(_wardKmedoidsSimilarity, annot=True)
    ax.set_title("Hierarchy(Ward) vs KMedoids Similarity")
    plt.xlabel('Hierarchy', fontsize = 15)
    plt.ylabel('KMedoids', fontsize = 15)
    plt.show()

    ax = sns.heatmap(_wardAverageSimilarity, annot=True)
    ax.set_title("Hierarchy(Ward) vs UPGMA Similarity")
    plt.xlabel('Hierarchy(Ward)', fontsize = 15) 
    plt.ylabel('UPGMA', fontsize = 15)
    plt.show()
