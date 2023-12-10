from configparser import ConfigParser
from DeltaBetaFilter import DeltaBetaFilter
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Src/Models/DeltaBetaFilterModel/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"]) 
    _betaDataDf.columns.values[0] = "CpG"

    _dbateFilter = DeltaBetaFilter()
    _dfOut = _dbateFilter.CalculateDeltaBeta(_betaDataDf, 50)

    FileSaver.SaveDataframe(_dfOut, _config["Paths"]["DBETA_DATA_PATH"])