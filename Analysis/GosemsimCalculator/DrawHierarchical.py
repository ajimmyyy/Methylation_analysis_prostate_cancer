from configparser import ConfigParser
from Models.GosemsimCalculator import GosemsimCalculator
import pandas as pd
from MakeFile.FileSaver import FileSaver
import matplotlib.pyplot as plt

if __name__ == "__main__":
    _configPath = "Analysis/GosemsimCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _gosemsimDf = pd.read_csv(_config["Paths"]["GOSEMSIM_MF_HYPER_PATH"], index_col = 0)

    _gosemsimCalculator = GosemsimCalculator()
    _fig = _gosemsimCalculator.DrawHierarchy(_gosemsimDf, "ward")

    plt.show()

    FileSaver.SaveData(_fig, _config["Paths"]["GENE_MF_HYPER_WARD_HIERARCHICAL_PATH"])