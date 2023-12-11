from configparser import ConfigParser
from DeltaBetaFilter import DeltaBetaFilter
import pandas as pd
import numpy as np
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Src/Models/DeltaBetaFilterModel/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    threshold_dbeta_hyper = 0.368
    threshold_dbeta_hypo = -0.238
    threshold_Pvalue = -np.log10(0.05)

    _dmpDataDf = pd.read_csv(_config["Paths"]["Filter_DBEAT_DATA_PATH"]) 

    _dbateFilter = DeltaBetaFilter()
    _hyperDf, _hypoDf = _dbateFilter.DetermineDNAm(_dmpDataDf, threshold_dbeta_hyper, threshold_dbeta_hypo, threshold_Pvalue)
    _fig = _dbateFilter.DrawVolcanoPlot(_dmpDataDf, _hyperDf, _hypoDf)

    FileSaver.SaveDataframe(_hyperDf, _config["Paths"]["DMP_HYPER_DATA_PATH"])
    FileSaver.SaveDataframe(_hypoDf, _config["Paths"]["DMP_HYPO_DATA_PATH"])
    FileSaver.SavePlot(_fig, _config["Paths"]["VOLCANO_PLOT_PATH"])