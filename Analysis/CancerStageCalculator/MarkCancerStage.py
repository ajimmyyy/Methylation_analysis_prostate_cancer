from configparser import ConfigParser
from Models.CancerStageCalculator import CancerStageCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/CancerStageCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _cancerStageDf = pd.read_csv(_config["Paths"]["DATASET_CANCER_STAGE_PATH"])
    _sampleSheet = pd.read_csv(_config["Paths"]["SAMPLE_SHEET_PATH"])
    
    _cancerStageDf = _cancerStageDf[_cancerStageDf["Sentrix_ID"].isin(_sampleSheet["Sentrix_ID"]) & _cancerStageDf["Sentrix_Position"].isin(_sampleSheet["Sentrix_Position"])]
    _sampleSheet = pd.merge(_sampleSheet, _cancerStageDf, on = ["Sentrix_ID", "Sentrix_Position"], how = "inner")

    _cancerStageCalculator = CancerStageCalculator()
    _sampleSheet, _earlyPatient, _laterPatient = _cancerStageCalculator.MarkCancerStage(_sampleSheet)

    FileSaver.SaveDataframe(_sampleSheet, _config["Paths"]["TRAIN_CANCER_STAGE"])
    FileSaver.SaveList(_earlyPatient["Sample_Name"].to_list(), _config["Paths"]["TRAIN_EARLY_CANCER_STAGE"])
    FileSaver.SaveList(_laterPatient["Sample_Name"].to_list(), _config["Paths"]["TRAIN_LATER_CANCER_STAGE"])