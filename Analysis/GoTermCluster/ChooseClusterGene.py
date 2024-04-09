from configparser import ConfigParser
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/GoTermCluster/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _geneCluster = pd.read_csv( _config["Paths"]["GENE_CLUSTER_MF_HYPER_WARD_PATH"])

    idx_max_dbeta = _geneCluster.groupby('cluster')['dbeta'].idxmax()
    idx_max_auc = _geneCluster.groupby('cluster')['auc'].idxmax()
    idx_max_F1 = _geneCluster.groupby('cluster')['F1'].idxmax()

    _geneCluster = _geneCluster.loc[idx_max_F1]

    FileSaver.SaveDataframe(_geneCluster, _config["Paths"]["MF_HYPER_WARD_CHOOSE_PATH"])

