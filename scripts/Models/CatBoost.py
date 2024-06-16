from sklearn.feature_selection import RFECV
from catboost import CatBoostRegressor
from imblearn.over_sampling import SVMSMOTE
import utils
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    # filter out the CpG sites
    # _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"])
    _aucDf = pd.read_csv("Data/Processed/Models/XGBoost/xgboost_feature_selection.csv")
    _aucDf = _aucDf[_aucDf['DNAm'] == "hyper"]
    keepFeature = _aucDf["CpG"].tolist()
    
    # read the training and testing data
    _trainDf = pd.read_csv(_config["Paths"]["TRAIN_BETA_DATA_PATH"], index_col=0)
    _testDf = pd.read_csv(_config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)
    
    # split the training, testing data into X and Y
    _trainX, _trainY, _testX, _testY = utils.SpliteTrainTest(_trainDf, _testDf, 25, keepFeature)

    # oversample the training data
    _trainX, _trainY = SVMSMOTE().fit_resample(_trainX, _trainY) 
    print(_trainY.value_counts())

    model = CatBoostRegressor(
        loss_function='RMSE',
        eval_metric='RMSE',
        use_best_model=True)

    model.fit(_trainX, _trainY, eval_set=(_testX, _testY), verbose=0, plot=True)