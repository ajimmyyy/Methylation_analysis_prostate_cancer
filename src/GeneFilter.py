import pandas as pd
import numpy as np
from tqdm import tqdm
from multipledispatch import dispatch

class GeneFilter:
    # IntersectData(self, originDf, selectList, selectType = "gene")
    # 取得交集資料
    # Parameters:
    # originDf: DataFrame，擁有selectType的資料
    # selectList: List，要選取的資料表
    # selectType: String，要篩選的行
    # Return:
    # DataFrame，originDf與selectList交集資料
    @dispatch(pd.DataFrame, list, str)
    def IntersectData(self, originDf, selectList, selectType = "gene"):
        intersectDf = originDf[originDf[selectType].isin(selectList)]
        return intersectDf
    
    # IntersectData(self, originDf, selectDf, selectType = "gene")
    # 取得交集資料
    # Parameters:
    # originDf: DataFrame，擁有selectType的資料
    # selectDf: DataFrame，擁有selectType的要選取的資料
    # selectType: String，要篩選的行
    # Return:
    # DataFrame，originDf與selectDf交集資料
    @dispatch(pd.DataFrame, pd.DataFrame, str)
    def IntersectData(self, originDf, selectDf, selectType = "gene"):
        intersectDf = pd.merge(originDf, selectDf, on = [selectType], how = "inner")
        return intersectDf
    
    # unionData(self, dataFrame1, dataFrame2, selectType = "gene")
    # 取得聯集資料
    # Parameters:
    # dataFrame1: DataFrame，擁有selectType的資料
    # dataFrame2: DataFrame，擁有selectType的資料
    # selectType: String，要篩選的行
    # Return:
    # DataFrame，dataFrame1與dataFrame2聯集資料
    def unionData(self, dataFrame1, dataFrame2, selectType = "gene"):
        unionDf = pd.merge(dataFrame1, dataFrame2, on = [selectType], how = "outer")
        return unionDf
    
    # differnceData(self, dataFrame1, dataFrame2, selectType = "gene")
    # 取得差集資料
    # Parameters:
    # dataFrame1: DataFrame，擁有selectType的資料
    # dataFrame2: DataFrame，擁有selectType的資料
    # selectType: String，要篩選的行
    # Return:
    # DataFrame，dataFrame1對dataFrame2的差集資料
    def differnceData(self, dataFrame1, dataFrame2, selectType = "gene"):
        differnceDf = dataFrame1[~dataFrame1.loc[:, selectType].isin(dataFrame2.loc[:, selectType])]
        return differnceDf