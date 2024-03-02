from configparser import ConfigParser
from ValidateModel.ValidateData import ValidateData
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/TestAccuracy/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["TEST_DATA_PATH"])
    _betaDataDf.columns.values[0] = "CpG"
    _testDf = pd.read_csv(_config["Paths"]["MF_HYPER_WARD_CHOOSE_PATH"])

    _validateData = ValidateData()
    _results = _validateData.ValidateGeneAccuracy(_testDf, _betaDataDf, 25, chooseOn = "F1")
    
    for i in range(len(_results)):
        print(_results[i])
