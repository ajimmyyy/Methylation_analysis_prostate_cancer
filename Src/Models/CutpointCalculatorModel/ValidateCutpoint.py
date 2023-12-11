from configparser import ConfigParser
from ValidateModel.ValidateData import ValidateData
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Src/Models/CutpointCalculatorModel/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["VALIDATE_DATA_PATH"])
    _betaDataDf.columns.values[0] = "CpG"
    _cutpointDataDf = pd.read_csv(_config["Paths"]["CUTPOINT_DATA_PATH"])

    _validateData = ValidateData()
    _dataOut, _dataNotFind = _validateData.ValidateCutpoint(_cutpointDataDf, _betaDataDf, 20)

    FileSaver.SaveDataframe(_dataOut, _config["Paths"]["CUTPOINT_VALIDATE_DATA_PATH"])
    