from configparser import ConfigParser
from Models.AucCalculator import AucCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/AucCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"])
    _betaDataDf.columns.values[0] = "CpG"
    _comorbidityGroupDf = pd.read_csv(_config["Paths"]["CUTPOINT_VALIDATE_GROUP_PATH"])

    _aucCalculator = AucCalculator()
    _outDf = _aucCalculator.CalculateAuc(_betaDataDf, _comorbidityGroupDf, 25)

    FileSaver.SaveData(_outDf, _config["Paths"]["AUC_GROUP_DATA_PATH"])