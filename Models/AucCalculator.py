import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, RocCurveDisplay, auc

class AucCalculator:
    # CalculateAuc(self, betaDf, originDf, normalCount, selectLine = "CpG")
    # 計算每個位點的AUC值
    # Parameters:
    # betaDf: DataFrame，擁有dbate的資料
    # originDf: DataFrame，要篩選的資料
    # normalCount: Int，normal資料數量
    # selectLine: String，要篩選的行
    # Return:
    # DataFrame
    # *請確保資料擁有各病人之beta值
    def CalculateAuc(self, betaDf, originDf, normalCount, selectLine = "CpG"):
        betaDf = betaDf[betaDf[selectLine].isin(originDf[selectLine])]

        tqdm.pandas(desc="find auc")
        betaDf.loc[:, ["auc"]] = betaDf.progress_apply(self.__CalculateRowAuc, axis = 1, args = (normalCount,))

        betaDf = betaDf[["CpG", "auc"]]
        mergeDf = pd.merge(originDf, betaDf, on = "CpG", how = "inner")
        mergeDf.loc[mergeDf["DNAm"] == "hypo", 'auc'] = 1 - mergeDf.loc[mergeDf["DNAm"] == "hypo", 'auc']

        return mergeDf
    
    # DrawRoc(self, betaDf, originDf, normalCount, selectLine = "CpG")
    # 畫出ROC曲線
    # Parameters:
    # betaDf: DataFrame，擁有dbate的資料
    # originDf: DataFrame，要篩選的資料
    # normalCount: Int，normal資料數量
    # selectLine: String，要篩選的行
    # Return:
    # Figure，ROC曲線圖
    def DrawRoc(self, betaDf, originDf, normalCount, selectLine = "CpG"):
        betaDf = betaDf[betaDf[selectLine].isin(originDf[selectLine])]
        fig = plt.figure(figsize=(32,18))

        tqdm.pandas(desc="draw roc")
        betaDf.progress_apply(self.__DrawRowRoc, axis = 1, args = (normalCount,))

        return fig

    def __CalculateRowAuc(self, row, normalCount):
        fpr, tpr = self.__CalculateRowRoc(row, normalCount)
        auc_score = auc(fpr, tpr)
        
        return pd.Series({"auc": auc_score})
    
    def __CalculateRowRoc(self, row, normalCount):
        rowList = row.to_numpy()

        predictBeta = rowList[1:]
        actual_beta = np.ones(int(len(predictBeta)))
        actual_beta[:normalCount] = 0

        fpr, tpr, threshold = roc_curve(actual_beta, predictBeta)
        return fpr, tpr
    
    def __DrawRowRoc(self, row, normalCount):
        fpr, tpr = self.__CalculateRowRoc(row, normalCount)
        plt.plot(fpr,tpr,marker = 'o')
