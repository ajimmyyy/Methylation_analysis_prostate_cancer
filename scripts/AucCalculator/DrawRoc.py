from configparser import ConfigParser
from Models.AucCalculator import AucCalculator
import pandas as pd
import matplotlib.pyplot as plt
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/AucCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"])
    _betaDataDf.columns.values[0] = "CpG"
    _comorbidityGroupDf = pd.read_csv(_config["Paths"]["CUTPOINT_VALIDATE_GROUP_PATH"])

    _aucCalculator = AucCalculator()
    _fig = _aucCalculator.DrawRoc(_betaDataDf, _comorbidityGroupDf, 25)
    
    plt.show()
    FileSaver.SaveData(_fig, _config["Paths"]["ROC_GROUP_CURVE_PATH"])