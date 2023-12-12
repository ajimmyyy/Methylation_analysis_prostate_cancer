from configparser import ConfigParser
from GeneFilter import GeneFilter
import pandas as pd
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Src/Models/GeneFilterModel/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _inDf = pd.read_csv(_config["Paths"]["CUTPOINT_VALIDATE_DATA_PATH"])

    with open(_config["Paths"]["GROUP_COMORBIDITY_PATH"], 'r') as file:
        _lines = file.readlines()
    _groupComorbidityList = [line.strip() for line in _lines]

    with open(_config["Paths"]["SINGLE_COMORBIDITY_PATH"], 'r') as file:
        _lines = file.readlines()
    _singleComorbidityList = [line.strip() for line in _lines]

    _geneFilter = GeneFilter()
    _groupDf = _geneFilter.IntersectData(_inDf, _groupComorbidityList, "gene")
    _singleDf = _geneFilter.IntersectData(_inDf, _singleComorbidityList, "gene")

    FileSaver.SaveDataframe(_groupDf, _config["Paths"]["CUTPOINT_VALIDATE_GROUP_PATH"])
    FileSaver.SaveDataframe(_singleDf, _config["Paths"]["CUTPOINT_VALIDATE_SINGLE_PATH"])