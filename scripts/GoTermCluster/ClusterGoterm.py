from configparser import ConfigParser
from Models.GeneFilter import GeneFilter
from Models.GosemsimCalculator import GosemsimCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver
import matplotlib.pyplot as plt

if __name__ == "__main__":
    _configPath = "Analysis/GoTermCluster/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"])
    _gosemsimDf = pd.read_csv(_config["Paths"]["GOSEMSIM_MEAN_HYPER_PATH"], index_col = 0)

    _gosemsimCalculator = GosemsimCalculator()
    _geneCluster, _silhouette = _gosemsimCalculator.ClusterHierarchy(_gosemsimDf, "ward", range(2, 25))
    plt.show()
    # _geneCluster, _silhouette = _gosemsimCalculator.ClusterKMedoids(_gosemsimDf, clusterNumRange=range(2, 25))
    # plt.show()
    # _geneClusterNJ, NJtree = _gosemsimCalculator.ClusterNeighborJioning(_gosemsimDf, clusterNum = 10)
    # plt.show()

    _geneFilter = GeneFilter()
    _geneCluster = _geneFilter.IntersectData(_aucDf, _geneCluster, "gene")

    FileSaver.SaveData(_geneCluster, _config["Paths"]["GENE_CLUSTER_MEAN_HYPER_HIERARCY_WARD_PATH"])
    FileSaver.SaveData(_silhouette, _config["Paths"]["GENE_MEAN_HYPER_HIERARCY_WARD_SILHOUETTE_PATH"])