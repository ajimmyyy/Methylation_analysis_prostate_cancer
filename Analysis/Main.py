from configparser import ConfigParser
import pandas as pd
import numpy as np

from Models.DeltaBetaFilter import DeltaBetaFilter
from Models.CutpointCalculator import CutpointCalculator
from Models.GeneFilter import GeneFilter
from Models.AucCalculator import AucCalculator
from Models.GosemsimCalculator import GosemsimCalculator
from Models.CancerStageCalculator import CancerStageCalculator

from MakeFile.FileSaver import FileSaver
from ValidateModel.ValidateData import ValidateData

if __name__ == "__main__":
    _configPath = "Analysis/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    # 載入資料
    NORMAL_NUM = 50
    VALIDATE_NORMAL_NUM = 20
    _hyperThreshold = 0.368
    _hypoThreshold = -0.238
    _pValueThreshold = -np.log10(0.05)

    _betaDataDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"])
    _betaDataDf.columns.values[0] = "CpG"
    _dmpDataDf = pd.read_csv(_config["Paths"]["DMP_DATA_PATH"]) 
    _dmpDataDf.columns.values[0] = "CpG"

    _betaValidateDataDf = pd.read_csv(_config["Paths"]["VALIDATE_DATA_PATH"])
    _betaValidateDataDf.columns.values[0] = "CpG"

    _gosemsimDf = pd.read_csv(_config["Paths"]["GOSEMSIM_MF_HYPER_PATH"], index_col = 0)

    with open(_config["Paths"]["GROUP_COMORBIDITY_PATH"], 'r') as file:
        _lines = file.readlines()
    _groupComorbidityList = [line.strip() for line in _lines]

    _deltaBetaFilter = DeltaBetaFilter()
    # 計算dbeta
    _dbetaDf = _deltaBetaFilter.CalculateDeltaBeta(_betaDataDf, NORMAL_NUM)
    # 篩選dbeta
    _dbetaDf = _deltaBetaFilter.FilterDeltaBeta(_dbetaDf, _dmpDataDf, onlyPromoter = True)
    # 分hyper,hypo
    _hyperDf, _hypoDf = _deltaBetaFilter.DetermineDNAm(_dmpDataDf, _hyperThreshold, _hypoThreshold, _pValueThreshold)

    _validateData = ValidateData()
    _cutpointCalculator = CutpointCalculator()
    # 計算切點
    _cutpointDf = _cutpointCalculator.CalculateCutpoint(_betaDataDf, [_hyperDf, _hypoDf], NORMAL_NUM, ["hyper", "hypo"], "CpG", "mid")
    # 驗證切點
    _cutpointDf = _validateData.ValidateCutpoint(_cutpointDf, _betaValidateDataDf, VALIDATE_NORMAL_NUM)

    _geneFilter = GeneFilter()
    # 與共病基因交集
    _groupDf = _geneFilter.IntersectData(_cutpointDf, _groupComorbidityList, "gene")

    _aucCalculator = AucCalculator()
    # 計算AUC
    _aucDf = _aucCalculator.CalculateAuc(_betaDataDf, _groupDf, NORMAL_NUM)

    # 
    # 計算基因相似度
    # 使用R語言gosemsim計算
    # 

    _gosemsimCalculator = GosemsimCalculator()
    # 分群
    _geneCluster, silhouetteFig = _gosemsimCalculator.ClusterGoterm(_gosemsimDf, "ward", range(5, 51))
    # 合併資料
    _geneCluster = _geneFilter.IntersectData(_aucDf, _geneCluster, "gene")

    # 挑選資料
    idx_max_dbeta = _geneCluster.groupby('cluster')['dbeta'].idxmax()
    idx_max_auc = _geneCluster.groupby('cluster')['auc'].idxmax()
    idx_max_F1 = _geneCluster.groupby('cluster')['F1'].idxmax()

    print("Dbeta Max:")
    print(_geneCluster.loc[idx_max_dbeta])
    print("AUC Max:")
    print(_geneCluster.loc[idx_max_auc])
    print("F1 Max:")
    print(_geneCluster.loc[idx_max_F1])


    



