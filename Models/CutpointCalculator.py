import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from multipledispatch import dispatch

class CutpointCalculator:
    # CalculateCutpoint(self, betaDf, selectDf, normalCount, type, selectLine = "CpG", cutpointType = "mid")
    # 計算各位點的切點
    # Parameters:
    # betaDf: DataFrame，擁有dbate的資料
    # selectDf: DataFrame，要選取的資料(將回傳此DF和其cutpoint)
    # normalCount: Int，normal資料數量
    # type: String，"Hyper"或"Hypo"
    # selectLine: String，betaDf和selectDf的交集行名，預設"CpG"
    # cutpointType: String，若遇到相同F1的cutpoint選擇，"mid"取中間值,"min"取最小值,"max"取最大值
    # Return:
    # DataFrame，selectDf各selectLine(CpG)位點的cutpoint
    # *請確保資料擁有各病人之beta值
    @dispatch(pd.DataFrame, pd.DataFrame, int, str, str, str)
    def CalculateCutpoint(self, betaDf, selectDf, normalCount, type, selectLine = "CpG", cutpointType = "mid"):
        betaDf = betaDf[betaDf[selectLine].isin(selectDf[selectLine])]

        tqdm.pandas(desc="find cutpoint")
        betaDf.loc[:, ("DNAm","cutpoint", "F1")] = betaDf.progress_apply(self.__CalculateRowCutpoint, axis = 1, args = (normalCount, type, cutpointType,))
        betaDf = betaDf[["CpG", "DNAm","cutpoint", "F1"]]

        assert betaDf.shape[0] == selectDf.shape[0], "size don't match"

        mergeDf = pd.merge(selectDf, betaDf, on = "CpG", how = "inner")
        return mergeDf
    
    # CalculateCutpoint(self, betaDf, selectDfList, normalCount, typeList, selectLine = "CpG", cutpointType = "mid")
    # 計算各位點的切點
    # Parameters:
    # betaDf: DataFrame，擁有dbate的資料
    # selectDfList: List[DataFrame]，要選取的資料(將回傳此DF和其cutpoint)
    # normalCount: Int，normal資料數量
    # typeList: List[String]，"Hyper"或"Hypo"
    # selectLine: String，betaDf和selectDf的交集行名，預設"CpG"
    # cutpointType: String，若遇到相同F1的cutpoint選擇，"mid"取中間值,"min"取最小值,"max"取最大值
    # Return:
    # DataFrame，合併資料，selectDfList各selectLine(CpG)位點的cutpoint
    # *請確保資料擁有各病人之beta值
    @dispatch(pd.DataFrame, list, int, list, str, str)
    def CalculateCutpoint(self, betaDf, selectDfList, normalCount, typeList, selectLine = "CpG", cutpointType = "mid"):
        assert len(selectDfList) == len(typeList), "size don't match"

        dfList = []
        for i in range(len(selectDfList)):
            dfList.append(self.CalculateCutpoint(betaDf, selectDfList[i], normalCount, typeList[i], selectLine, cutpointType))

        mergeDf = pd.concat(dfList)
        return mergeDf
    
    def __CalculateRowCutpoint(self, row, normalCount, type, cutpointType):
        maxSS = 0.0
        minCutpoint = 0.0
        maxCutpoint = 0.0
        maxF1 = 0.0

        for cutpoint in np.arange(0.01, 1, 0.01):
            rowList = row.to_numpy()
            normalBeta = rowList[1:normalCount + 1]
            tumorBeta = rowList[normalCount + 1:]

            if type == "hyper":
                tn = np.sum(normalBeta < cutpoint)
                fp = np.sum(normalBeta > cutpoint)
                fn = np.sum(tumorBeta < cutpoint)
                tp = np.sum(tumorBeta > cutpoint)
            elif type == "hypo":
                tn = np.sum(normalBeta > cutpoint)
                fp = np.sum(normalBeta < cutpoint)
                fn = np.sum(tumorBeta > cutpoint)
                tp = np.sum(tumorBeta < cutpoint)

            precision = tp / (tp + fp) if (tp + fp) != 0 else 0
            sensitivity  = tp / (tp + fn) if (tp + fn) != 0 else 0
            specificity = tn / (fp + tn) if (fp + tn) != 0 else 0

            if maxSS < sensitivity + specificity:
                maxSS = sensitivity + specificity
                minCutpoint = cutpoint
                maxCutpoint = cutpoint
                maxF1 = 2 * sensitivity * precision / (sensitivity + precision) if (precision + sensitivity) != 0 else 0
            elif maxSS == sensitivity + specificity:
                maxCutpoint = cutpoint

        if cutpointType == "min":
            maxCutpoint = minCutpoint
        elif cutpointType == "mid":
            maxCutpoint = (minCutpoint + maxCutpoint) / 2
        elif cutpointType == "max":
            pass
        
        return pd.Series({"DNAm": type,"cutpoint": round(maxCutpoint, 2), "F1": maxF1})
