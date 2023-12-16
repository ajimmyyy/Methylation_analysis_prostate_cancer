import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from itertools import combinations

class ValidateData:
    NOT_FIND = -1

    # randomChooseValidation(betaDf, normalCount, normalChoose, tumorCount, tumorChoose)
    # 隨機選取驗證集
    # Parameters:
    # betaDf: DataFrame，資料集
    # normalCount: Int，normal資料數量
    # normalChoose: Int，要選取的normal資料數量
    # tumorCount: Int，tumor資料數量
    # tumorChoose: Int，要選取的tumor資料數量
    # Return:
    # data: DataFrame，驗證集
    # selectColumns: List，選取表
    def randomChooseValidation(betaDf, normalCount, normalChoose, tumorCount, tumorChoose):
        randomNormal = np.random.choice(normalCount//2, normalChoose // 2, replace=False) * 2
        randomTumor = np.random.choice(tumorCount // 2, tumorChoose // 2, replace=False) * 2

        selectColumns = np.sort(np.concatenate([[0], randomNormal + 1, randomNormal + 2, randomTumor + normalCount + 1, randomTumor + normalCount + 2]))

        data = betaDf.iloc[:, selectColumns.tolist()]

        return data, selectColumns

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
    
    def ValidateGeneAccuracy(self, cutpointDf, betaDfs, normalCount, testCount = 3, method = "count", chooseOn = "F1"):
        allCutpointDf = cutpointDf[["CpG", "cutpoint"]]
        combinationsList = list(combinations(allCutpointDf.index, testCount))

        cutpointDfList = []
        for combination in combinationsList:
            selectedColumns = allCutpointDf.loc[list(combination)]
            cutpointDfList.append(selectedColumns)
        
        if method == "count":
            result = self.__PredictResultCount(cutpointDfList, betaDfs, normalCount, "CpG", chooseOn)

        return result

    
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
    
    def __PredictResultCount(self, cutpointDfList, betaDf, normalCount, mergeOn, chooseOn):
        maxScore = 0
        result = []
        for cutpointDf in cutpointDfList:
            testBetaDf = betaDf[betaDf[mergeOn].isin(cutpointDf[mergeOn])]
            df = pd.merge(cutpointDf, testBetaDf, on = [mergeOn], how = 'inner')
            betaColumns = df.columns[2:]

            predictList = []
            for betaColumn in betaColumns:
                isTumor = df[betaColumn].gt(df['cutpoint'], axis=0)
                labelSum = isTumor.astype(int).sum(axis=0)
                predictList.append(labelSum > int(len(df.index) / 2))
            
            f1, accuracy = self.__CalculateAccuracy(predictList, normalCount)
            if chooseOn == "F1":
                if f1 > maxScore:
                    maxScore = f1
                    result.clear()
                    result.append(cutpointDf)
                elif f1 == maxScore:
                    result.append(cutpointDf)
            elif chooseOn == "accuracy":
                if accuracy > maxScore:
                    maxScore = accuracy
                    result.clear()
                    result.append(cutpointDf)
                elif accuracy == maxScore:
                    result.append(cutpointDf)
                    
        print("max score:", maxScore)
        return result
        
    def __CalculateAccuracy(self, predictList, normalCount):
        tn = predictList[:normalCount].count(False)
        fp = predictList[:normalCount].count(True)
        fn = predictList[normalCount + 1:].count(False)
        tp = predictList[normalCount + 1:].count(True)

        precision = tp / (tp + fp) if (tp + fp) != 0 else 0
        sensitivity  = tp / (tp + fn) if (tp + fn) != 0 else 0

        f1 = 2 * sensitivity * precision / (sensitivity + precision) if (precision + sensitivity) != 0 else 0
        accuracy = (tp + tn) / len(predictList)

        return f1, accuracy


