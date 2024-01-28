from configparser import ConfigParser
from Models.GosemsimCalculator import GosemsimCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver
import matplotlib.pyplot as plt

if __name__ == "__main__":
    _configPath = "Analysis/GosemsimCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _gosemsimDf = pd.read_csv(_config["Paths"]["GOSEMSIM_MF_HYPER_PATH"], index_col = 0)

    _gosemsimCalculator = GosemsimCalculator()
    # _geneClusterHierarchy, silhouetteHierarchy = _gosemsimCalculator.ClusterHierarchy(_gosemsimDf, "ward")
    # plt.show()
    _geneClusterKMedoids, silhouetteKMedoids = _gosemsimCalculator.ClusterKMedoids(_gosemsimDf)
    plt.show()
    # _geneClusterNJ, NJtree = _gosemsimCalculator.ClusterNeighborJioning(_gosemsimDf)
    # plt.show()

    # FileSaver.SaveDataframe(_geneClusterHierarchy, _config["Paths"]["GENE_CLUSTER_MF_HYPER_HIERARCY_WARD_PATH"])
    # FileSaver.SavePlot(silhouetteHierarchy, _config["Paths"]["GENE_MF_HYPER_HIERARCY_WARD_SILHOUETTE_PATH"])