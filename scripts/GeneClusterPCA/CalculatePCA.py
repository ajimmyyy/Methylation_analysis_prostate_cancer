from configparser import ConfigParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Models.GeneFilter import GeneFilter
from MakeFile.FileSaver import FileSaver
from sklearn.decomposition import PCA

if __name__ == "__main__":
    _configPath = "Analysis/GeneClusterPCA/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"])
    _betaDataDf.columns.values[0] = "CpG"
    _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"], usecols=["CpG", "gene"])
    _geneWardCluster = pd.read_csv(_config["Paths"]["GENE_CLUSTER_MF_HYPER_WARD_PATH"])

    _geneFilter = GeneFilter()
    _geneWardCluster = _geneFilter.IntersectData(_aucDf, _geneWardCluster, "gene")
    _geneWardCluster = _geneFilter.IntersectData(_betaDataDf, _geneWardCluster, "CpG")

    FileSaver.SaveData(_geneWardCluster, "C:/Users/acer/Desktop/test/test.csv")

    _clusters = _geneWardCluster['cluster'].unique()

    for cluster in _clusters:
        clusterData = _geneWardCluster[_geneWardCluster['cluster'] == cluster].iloc[:, 1:-2].T

        pca = PCA(n_components = 3)
        transformed_data = pca.fit_transform(clusterData)

        print("cluster: ", cluster)
        print(clusterData.shape)
        print(transformed_data.shape)

        pc1 = transformed_data[:, 0]
        pc2 = transformed_data[:, 1]

        colors = np.array(['red'] * 50 + ['blue'] * (len(pc1) - 50))
        plt.figure(figsize=(8,6))
        plt.scatter(pc1, pc2,c=colors,cmap='rainbow')
        plt.xlabel('First principal component')
        plt.ylabel('Second Principal Component')
        plt.show()

