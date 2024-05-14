from configparser import ConfigParser
import pandas as pd
from Models.CancerStageCalculator import CancerStageCalculator
from ValidateModel.ValidateData import ValidateData

if __name__ == "__main__":
    _configPath = "Analysis/CancerStageCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDf = pd.read_csv(_config["Paths"]["TCGA_PARD_450K"], index_col=0)
    _sampleSheetDf = pd.read_csv(_config["Paths"]["TCGA_PARD_SAMPLE_SHEET"])
    _betaDf2 = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"], index_col=0)
    _sampleSheetDf2 = pd.read_csv(_config["Paths"]["SAMPLE_SHEET_PATH"])

    _cancerStageCalculator = CancerStageCalculator()
    _earlyDf = _cancerStageCalculator.CancerStageSplit(_betaDf, _sampleSheetDf, ["T3", "T4"])
    _nromalDf = _betaDf.iloc[:, :50]
    _earlyDf2 = _cancerStageCalculator.CancerStageSplit(_betaDf2, _sampleSheetDf2, ["T3", "T4"], stageColName="ajcc_clinical_t", sampleColName="sample_id")
    _nromalDf2 = _betaDf2.iloc[:, :25]

    _mergedDf = pd.concat([_nromalDf, _nromalDf2, _earlyDf, _earlyDf2], axis=1)
    _mergedDf.dropna(inplace=True)

    _mergedDf.to_csv(_config["Paths"]["LATE_STAGE_DATA_PATH"], index=True)