from configparser import ConfigParser
from Models.GosemsimCalculator import GosemsimCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver
import matplotlib.pyplot as plt

if __name__ == "__main__":
    _configPath = "Analysis/GosemsimCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _gosemsimDf = pd.read_csv(_config["Paths"]["MF_GONTO_SIMILARITY_HYPER_PATH"], index_col = 0)

    _gosemsimCalculator = GosemsimCalculator()
    _fig = _gosemsimCalculator.DrawHeatmap(_gosemsimDf, "average")

    plt.show()

    # FileSaver.SavePlot(_fig, _config["Paths"]["GENE_MF_HYPER_AVERAGE_HEATMAP_PATH"])