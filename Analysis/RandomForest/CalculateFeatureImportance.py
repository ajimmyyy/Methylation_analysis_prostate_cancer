from configparser import ConfigParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from sklearn.ensemble import RandomForestClassifier
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/RandomForest/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _forest = joblib.load(_config["Paths"]["RANDOM_FOREST_TREE_PATH"])

    importances = _forest.feature_importances_

