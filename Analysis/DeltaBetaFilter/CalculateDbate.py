from configparser import ConfigParser
from Models.DeltaBetaFilter import DeltaBetaFilter
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/DeltaBetaFilter/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"]) 
    _betaDataDf.columns.values[0] = "CpG"

    _dbateFilter = DeltaBetaFilter()
    _dfOut = _dbateFilter.CalculateDeltaBeta(_betaDataDf, 25)

    FileSaver.SaveData(_dfOut, _config["Paths"]["DBETA_DATA_PATH"])