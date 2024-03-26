from configparser import ConfigParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from RandomForest import TransformTrainData
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _testDf = pd.read_csv(_config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)

    _testDf = _testDf.iloc[:, 1::2]
    _testDf.to_csv('C:/Users/acer/Desktop/test/test.csv')

    _trainDf = pd.read_csv(_config["Paths"]["TRAIN_BETA_DATA_PATH"], index_col=0)
    _trainDf = TransformTrainData(_trainDf, 25)
    _extroDf = pd.read_csv(_config["Paths"]["450K_DATA_PATH"], index_col=0)
    _extroDf = TransformTrainData(_extroDf, 10)
    
    _mixTrainDf = pd.merge(_trainDf, _extroDf, left_index=True, right_index=True, how='outer')
    FileSaver.SaveDataframe(_mixTrainDf, "C:/Users/acer/Desktop/test/mixTrain.csv")