from configparser import ConfigParser
import pandas as pd
from Models.CancerStageCalculator import CancerStageCalculator

if __name__ == "__main__":
    _configPath = "Analysis/CancerStageCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDf = pd.read_csv(_config["Paths"]["TCGA_PARD_450K"])
    _sampleSheetDf = pd.read_csv(_config["Paths"]["TCGA_PARD_SAMPLE_SHEET"])

    _cancerStageCalculator = CancerStageCalculator()

    _earlyDf = _cancerStageCalculator.CancerStageSplit(_betaDf, _sampleSheetDf, ["T1", "T2"])
    _nromalDf = _betaDf.iloc[:, :51]

    _mergedDf = pd.concat([_nromalDf, _earlyDf], axis=1)

    _mergedDf.to_csv(_config["Paths"]["EARLY_STAGE_DATA_PATH"], index=False)