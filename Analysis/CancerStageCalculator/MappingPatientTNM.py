from configparser import ConfigParser
from Models.CancerStageCalculator import CancerStageCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/CancerStageCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _TcgaList = pd.read_csv(_config["Paths"]["TCGA_ID_PATH"], sep = "\t")
    _cancerStageDf = pd.read_csv(_config["Paths"]["PAIENT_CANCER_STAGE_PATH"])

    _IdMap = _TcgaList[["Comment [TCGA Barcode]", "Array Data File"]]
    _cancerStageMap = _cancerStageDf[["Patient ID",
                                        "American Joint Committee on Cancer Metastasis Stage Code",
                                        "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code", 
                                        "American Joint Committee on Cancer Tumor Stage Code"]]
    
    _IdMap = _IdMap.rename(columns = {"Comment [TCGA Barcode]": "TCGA_barcode", "Array Data File": "file_ID"})
    _cancerStageMap = _cancerStageMap.rename(columns = {"Patient ID": "TCGA_barcode", 
                                                          "American Joint Committee on Cancer Metastasis Stage Code": "metastasis_stage", 
                                                          "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": "node_stage", 
                                                          "American Joint Committee on Cancer Tumor Stage Code": "tumor_stage"})

    _IdMap["TCGA_barcode"] = _IdMap["TCGA_barcode"].apply(lambda x: x[:12])
    _IdMap['Sentrix_ID'] = _IdMap["file_ID"].str[:10]
    _IdMap['Sentrix_Position'] = _IdMap["file_ID"].str[11:17]
    _IdMap = _IdMap.drop(columns = ["file_ID"])

    _cancerStageCalculator = CancerStageCalculator()
    _outDf = _cancerStageCalculator.MappingPatientTNM(_cancerStageMap, _IdMap)

    FileSaver.SaveDataframe(_outDf, _config["Paths"]["DATASET_CANCER_STAGE_PATH"])