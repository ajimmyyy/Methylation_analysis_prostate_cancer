from configparser import ConfigParser
from Models.DeltaBetaFilter import DeltaBetaFilter
import pandas as pd
import pickle
import tqdm
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/CancerStageCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"])
    _betaDf.columns.values[0] = 'CpG'
    with open(_config["Paths"]["TRAIN_EARLY_CANCER_STAGE"], 'rb') as f:
        _earlyList = pickle.load(f)
    with open(_config["Paths"]["TRAIN_LATER_CANCER_STAGE"], 'rb') as f:
        _laterList = pickle.load(f)

    _earlyList = [i for i in range(51)] + _earlyList
    _laterList = [i for i in range(51)] + _laterList

    _earlyDf = _betaDf.iloc[:, _earlyList]
    _laterDf = _betaDf.iloc[:, _laterList]

    _deltaBetaFilter = DeltaBetaFilter()
    _earlyDbetaDf = _deltaBetaFilter.CalculateDeltaBeta(_earlyDf, 50)
    _earlyDbetaDf = _earlyDbetaDf.drop("avg", axis=1)
    _earlyDbetaDf.rename(columns={"dbeta": "early_dbeta"}, inplace=True)
    _laterDbetaDf = _deltaBetaFilter.CalculateDeltaBeta(_laterDf, 50)
    _laterDbetaDf = _laterDbetaDf.drop("avg", axis=1)
    _laterDbetaDf.rename(columns={"dbeta": "later_dbeta"}, inplace=True)

    outDf = pd.merge(_earlyDbetaDf, _laterDbetaDf, on = ["CpG"], how = "inner")
    outDf["dbate_dif"] = outDf["early_dbeta"] - outDf["later_dbeta"]
    outDf.to_csv(_config["Paths"]["CANCER_STAGE_DBETA"], sep=',', encoding='utf-8', index=False)