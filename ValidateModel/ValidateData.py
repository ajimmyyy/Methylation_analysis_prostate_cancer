import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

class ValidateData:
    NOT_FIND = -1    

    # ValidateCutpoint(self, cutpointDf, betaDf, normalCount, thresholdF1 = 0.8)
    # 驗證每個位點的切點
    # Parameters:
    # cutpointDf: DataFrame，擁有cutpoint的資料
    # betaDf: DataFrame，擁有beta的資料
    # normalCount: Int，normal資料數量
    # thresholdF1: Float，F1閥值，預設0.8
    # Return:
    # cutpointDf, notfindDf: DataFrame，驗證後的資料和找不到的資料
    # *請確保資料擁有各病人之beta值
    def ValidateCutpoint(self, cutpointDf, betaDf, normalCount, thresholdF1 = 0.8):
        tqdm.pandas(desc="find cutpoint")
        cutpointDf["F1_validate"] = cutpointDf.progress_apply(self.__CalculateF1, axis = 1, args = (betaDf, normalCount,))
        
        notfindDf = cutpointDf[cutpointDf["F1_validate"] != self.NOT_FIND]
        cutpointDf = cutpointDf[cutpointDf["F1"] > thresholdF1]
        cutpointDf = cutpointDf[cutpointDf["F1_validate"] > thresholdF1]
        
        return cutpointDf, notfindDf
    
    def __CalculateF1(self, row, betaDf, normalCount):
        cutpoint = row["cutpoint"]
        betaDf = betaDf[betaDf["CpG"] == row["CpG"]]
        try:
            betaDf = betaDf.to_numpy()[0]
        except:
            return pd.Series({"F1_validate": self.NOT_FIND})
        normalBeta = betaDf[1:normalCount + 1:2]
        tumorBeta = betaDf[normalCount + 1::2]

        if row["DNAm"] == "hyper":
            fp = np.sum(normalBeta > cutpoint)
            fn = np.sum(tumorBeta < cutpoint)
            tp = np.sum(tumorBeta > cutpoint)
        if row["DNAm"] == "hypo":
            fp = np.sum(normalBeta < cutpoint)
            fn = np.sum(tumorBeta > cutpoint)
            tp = np.sum(tumorBeta < cutpoint)

        precision = tp / (tp + fp) if (tp + fp) != 0 else 0
        sensitivity  = tp / (tp + fn) if (tp + fn) != 0 else 0

        f1 = 2 * sensitivity * precision / (sensitivity + precision) if (precision + sensitivity) != 0 else 0

        return pd.Series({"F1_validate": f1})