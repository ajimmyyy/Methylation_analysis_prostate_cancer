import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# 以Deltabeta進行操作
class DeltaBetaFilter:

    # CalculateDeltaBeta(self, betaDf, normalCount)
    # 計算各位點的delta beta值
    # Parameters:
    # betaDf: DataFrame，包含各病人,CpG位點的beta值
    # normalCount: Int，normal資料數量
    # Return:
    # DataFrame，各CpG位點的dbeta,avg
    # *請確保資料擁有各病人之beta值
    def CalculateDeltaBeta(self, betaDf, normalCount):
        tqdm.pandas(desc="find dbeta")
        betaDf.loc[:, ["dbeta", "avg"]] = betaDf.progress_apply(self.__CalculateRowDeltaBeta, axis = 1, args = (normalCount,))

        return betaDf.loc[:, ["CpG", "dbeta", "avg"]]
    
    # FilterDeltaBeta(self, dbetaDf, dmpDf, onlyPromoter = True)
    # 篩選個基因裡dbate
    # Parameters:
    # dbetaDf: DataFrame，擁有dbate的資料
    # dmpDf: DataFrame，DMP資料
    # onlyPromoter: bool，是否只篩選promoter，預設true
    # Return:
    # DataFrame，合併資料
    def FilterDeltaBeta(self, dbetaDf, dmpDf, onlyPromoter = True):
        dbetaDf = dbetaDf[["CpG", "dbeta"]]
        merged_df = pd.merge(dmpDf, dbetaDf, on='CpG')

        if onlyPromoter:
            merged_df = merged_df[merged_df['feature'].isin(['TSS200', 'TSS1500'])]

        merged_df['abs_dbeta'] = merged_df['dbeta'].abs()
        max_dbeta_index = merged_df.groupby('gene')['abs_dbeta'].idxmax()
        merged_df = merged_df.loc[max_dbeta_index]
        merged_df = merged_df.drop(columns=['abs_dbeta'], axis=1)

        return merged_df
    
    # DetermineDNAm(self, dmpDf, thresholdDbetaHyper, thresholdDbetaHypo, thresholdPvalue, dbetaRow = "dbeta")
    # 分離資料Hyper,Hypo
    # Parameters:
    # dmpDf: DataFrame，DMP資料
    # thresholdDbetaHyper: double，Hyper部分dbeta閥值
    # thresholdDbetaHypo: double，Hypo部分dbeta閥值
    # thresholdPvalue: double，Pvalue閥值
    # dbetaRow: string，dbeta行名，預設"deltaBeta"
    # Return:
    # hyper, hypo: DataFrame，hyper,hypo資料
    def DetermineDNAm(self, dmpDf, thresholdDbetaHyper, thresholdDbetaHypo, thresholdPvalue, dbetaRow = "dbeta"):
        dmpDf = dmpDf[dmpDf["adj.P.Val"] < thresholdPvalue]

        hyper = dmpDf[dmpDf[dbetaRow] > thresholdDbetaHyper]
        hypo = dmpDf[dmpDf[dbetaRow] < thresholdDbetaHypo]

        return hyper, hypo
    
    # DrawVolcanoPlot(self, dmpDf, hyperDf, hypoDf, dbetaRow = "dbeta")
    # 畫出火山圖
    # Parameters:
    # dmpDf: DataFrame，DMP資料
    # hyperDf: DataFrame，hyper資料
    # hypoDf: DataFrame，hypo資料
    # dbetaRow: string，dbeta行名，預設"dbeta"
    # Return:
    # Figure，火山圖
    def DrawVolcanoPlot(self, dmpDf, hyperDf, hypoDf, dbetaRow = "dbeta"):
        fig = plt.figure(figsize=(32,18))
        plt.style.use("ggplot")
        plt.xlabel("dbeta", fontweight = "bold")
        plt.ylabel("-log(P.adj)", fontweight = "bold")

        plt.scatter((dmpDf[dbetaRow]), -np.log10(dmpDf["adj.P.Val"]), c = "grey", s = 15, alpha = .5, label='other')
        plt.scatter((hyperDf[dbetaRow]), -np.log10(hyperDf["adj.P.Val"]), c = "red", s = 15, alpha = .5, label='hyper')
        plt.scatter((hypoDf[dbetaRow]), -np.log10(hypoDf["adj.P.Val"]), c = "blue", s = 15, alpha = .5, label='hypo')

        plt.legend()

        return fig

    def __CalculateRowDeltaBeta(self, row, normalCount):
        rowList = row.to_numpy()

        normalBeta = rowList[1:normalCount + 1]
        tumorBeta = rowList[normalCount + 1:]

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