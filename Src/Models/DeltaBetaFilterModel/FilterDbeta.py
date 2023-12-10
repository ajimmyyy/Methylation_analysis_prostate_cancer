from configparser import ConfigParser
from DeltaBetaFilter import DeltaBetaFilter
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Src/Models/DeltaBetaFilterModel/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _dbetaDataDf = pd.read_csv(_config["Paths"]["DBETA_DATA_PATH"]) 
    _dmpDataDf = pd.read_csv(_config["Paths"]["DMP_DATA_PATH"]) 
    _dmpDataDf.columns.values[0] = "CpG"

    _dbateFilter = DeltaBetaFilter()
    _dfOut = _dbateFilter.FilterDeltaBeta(_dbetaDataDf, _dmpDataDf, onlyPromoter = True)

    FileSaver.SaveDataframe(_dfOut, _config["Paths"]["Filter_DBEAT_DATA_PATH"])