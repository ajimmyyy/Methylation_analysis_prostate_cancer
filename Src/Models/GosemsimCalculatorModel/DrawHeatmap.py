from configparser import ConfigParser
from GosemsimCalculator import GosemsimCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver
import matplotlib.pyplot as plt

if __name__ == "__main__":
    _configPath = "Src/Models/GosemsimCalculatorModel/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _gosemsimDf = pd.read_csv(_config["Paths"]["GOSEMSIM_MF_HYPER_PATH"], index_col = 0)

    _gosemsimCalculator = GosemsimCalculator()
    _fig = _gosemsimCalculator.DrawHeatmap(_gosemsimDf, "ward")

    FileSaver.SavePlot(_fig, _config["Paths"]["GENE_MF_HYPER_WARD_HEATMAP_PATH"])