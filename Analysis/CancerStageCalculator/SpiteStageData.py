from configparser import ConfigParser
import pandas as pd
from Models.CancerStageCalculator import CancerStageCalculator

if __name__ == "__main__":
    _configPath = "Analysis/CancerStageCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDf = pd.read_csv(_config["Paths"]["TCGA_PARD_450K"], index_col=0)
    _betaDf2 = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"], index_col=0)
    _sampleSheetDf = pd.read_csv(_config["Paths"]["TCGA_PARD_SAMPLE_SHEET"])
    _sampleSheetDf2 = pd.read_csv(_config["Paths"]["SAMPLE_SHEET_PATH"])

    _cancerStageCalculator = CancerStageCalculator()
    _earlyDf = _cancerStageCalculator.CancerStageSplit(_betaDf, _sampleSheetDf, ["T1", "T2"])
    _nromalDf = _betaDf.iloc[:, :50]
    _earlyDf2 = _cancerStageCalculator.CancerStageSplit(_betaDf2, _sampleSheetDf2, ["T1", "T2"], stageColName="ajcc_pathologic_t", sampleColName="sample_id")
    _nromalDf2 = _betaDf2.iloc[:, :25]

    _normalMerge = pd.merge(_nromalDf, _nromalDf2, left_index=True, right_index=True, how="inner")
    _tumorMerge = pd.merge(_earlyDf, _earlyDf2, left_index=True, right_index=True, how="inner")
    _earlyMerge = pd.concat([_normalMerge, _tumorMerge], axis=1)

    _earlyMerge.to_csv(_config["Paths"]["EARLY_STAGE_DATA_PATH"], index=True)

    # _betaDf = pd.read_csv(_config["Paths"]["TCGA_PARD_450K"])
    # _betaDf2 = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"])
    # _betaDf2.columns.values[0] = "CpG"
    # _sampleSheetDf = pd.read_csv(_config["Paths"]["TCGA_PARD_SAMPLE_SHEET"])
    # _sampleSheetDf2 = pd.read_csv(_config["Paths"]["SAMPLE_SHEET_PATH"])

    # _name = _betaDf.iloc[:, :1]
    # _betaDf = _betaDf.iloc[:, 51:]
    # _cancerStageCalculator = CancerStageCalculator()
    # _earlyDf = _cancerStageCalculator.CancerStageSplit(_betaDf, _sampleSheetDf, ["T1", "T2"])

    # _name2 = _betaDf2.iloc[:, :1]
    # _betaDf2 = _betaDf2.iloc[:, 26:]
    # _earlyDf2 = _cancerStageCalculator.CancerStageSplit(_betaDf2, _sampleSheetDf2, ["T1", "T2"], stageColName="ajcc_pathologic_t", sampleColName="sample_id")

    # _mergedDf = pd.concat([_name, _earlyDf], axis=1)
    # _mergedDf2 = pd.concat([_name2, _earlyDf2], axis=1)

    # _earlyMerge = pd.merge(_mergedDf, _mergedDf2, on="CpG", how="outer")
    # _earlyMerge.to_csv(_config["Paths"]["EARLY_STAGE_DATA_PATH"], index=False)

    