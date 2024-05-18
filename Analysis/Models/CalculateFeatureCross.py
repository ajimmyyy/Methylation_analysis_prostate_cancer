from configparser import ConfigParser
import pandas as pd

if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _rfDf = pd.read_csv(_config["Paths"]["RANDOM_FOREST_FEATURES_SELECTION_PATH"], index_col=0)
    _xgbDf = pd.read_csv(_config["Paths"]["XGBOOST_FEATURES_SELECTION_PATH"], index_col=0)

    index = _rfDf.index.intersection(_xgbDf.index)
    crossDf = _rfDf.loc[index]

    crossDf.to_csv(_config["Paths"]["CROSS_FEATURE_PATH"], index=True)