from configparser import ConfigParser
from Models.GeneFilter import GeneFilter
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/GeneFilter/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _inDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"])
    _selectDf = pd.read_csv(_config["Paths"]["STUDY_DATA_PATH"])

    _geneFilter = GeneFilter()
    _studyCountDf = _geneFilter.IntersectData(_inDf, _selectDf, "gene")

    FileSaver.SaveData(_studyCountDf, _config["Paths"]["AUC_GROUP_STUDY_DATA_PATH"])