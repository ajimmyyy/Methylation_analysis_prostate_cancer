from configparser import ConfigParser
from Models.CutpointCalculator import CutpointCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/CutpointCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"])
    _betaDataDf.columns.values[0] = "CpG"
    _hyperDataDf = pd.read_csv(_config["Paths"]["DMP_HYPER_DATA_PATH"]) 
    _hypoDataDf = pd.read_csv(_config["Paths"]["DMP_HYPO_DATA_PATH"])

    _cutpointCalculator = CutpointCalculator()
    _dfOut = _cutpointCalculator.CalculateCutpoint(_betaDataDf, [_hyperDataDf, _hypoDataDf], 25, ["hyper", "hypo"], "CpG", "mid")

    FileSaver.SaveDataframe(_dfOut, _config["Paths"]["CUTPOINT_DATA_PATH"])