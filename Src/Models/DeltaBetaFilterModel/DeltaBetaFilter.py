import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from configparser import ConfigParser
from MakeFile.FileSaver import FileSaver

class DeltaBetaFilter:
    def CalculateDeltaBeta(self, betaDf, normalCount):
        tqdm.pandas(desc="find dbeta")
        betaDf[["dbeta", "avg"]] = betaDf.progress_apply(self.__CalculateRowDeltaBeta, axis = 1, args = (normalCount,))

        return betaDf.loc[:, ["CpG", "dbeta", "avg"]]

    def __CalculateRowDeltaBeta(self, row, normalCount):
        rowList = row.to_numpy()

        normalBeta = rowList[1:normalCount + 1:2]
        tumorBeta = rowList[normalCount + 1::2]

        normalBeta = self.__RemoveOutliers(normalBeta)
        tumorBeta = self.__RemoveOutliers(tumorBeta)

        avg = np.mean(normalBeta)
        tumorBeta = tumorBeta - avg

        tumorBeta = self.__RemoveOutliers(tumorBeta)
        dbate = np.mean(tumorBeta)

        return pd.Series({"dbeta": dbate, "avg": avg})
    
    def __RemoveOutliers(self, row):
        q1 = np.percentile(row, 25)
        q3 = np.percentile(row, 75)
        iqr = q3 - q1

        return [x for x in row if (x >= q1 - 1.5 * iqr) and (x <= q3 + 1.5 * iqr)]

if __name__ == "__main__":
    _configPath = "Src/Models/DeltaBetaFilterModel/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)
    
    _dbateFilter = DeltaBetaFilter()

    _betaDataDf = pd.read_csv(_config["Input"]["BETA_DATA_PATH"]) 
    _betaDataDf.columns.values[0] = "CpG"

    _dfOut = _dbateFilter.CalculateDeltaBeta(_betaDataDf, 50)

    FileSaver.SaveDataframe(_dfOut, _config["Output"]["DBETA_DATA_PATH"])

    
