from configparser import ConfigParser
import pandas as pd
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

    FileSaver.SaveDataframe(_geneWardCluster, "C:/Users/acer/Desktop/test/test.csv")

    _clusters = _geneWardCluster['cluster'].unique()

    for cluster in _clusters:
        clusterData = _geneWardCluster[_geneWardCluster['cluster'] == cluster].iloc[:, 1:-2].T

        pca = PCA()
        transformed_data = pca.fit_transform(clusterData)

        explainedVarianceRatio = pca.explained_variance_ratio_
        original_genes = _geneWardCluster[_geneWardCluster['cluster'] == cluster]['gene'].values
        for i, gene in enumerate(original_genes):
            print(f'Cluster {cluster} - {gene}: {explainedVarianceRatio[i]*100:.2f}%')

        cumulativeVarianceRatio = explainedVarianceRatio.cumsum()
        plt.plot(range(1, len(cumulativeVarianceRatio) + 1), cumulativeVarianceRatio, marker='o')
        plt.xlabel('Number of Principal Components')
        plt.ylabel('Cumulative Variance Ratio')
        plt.title(f'Cluster {cluster} - Cumulative Variance Ratio of Principal Components')
        plt.show()

