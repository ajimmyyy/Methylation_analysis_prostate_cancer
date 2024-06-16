from configparser import ConfigParser
import pandas as pd
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import StackingClassifier
from sklearn.feature_selection import RFECV, RFE
from imblearn.over_sampling import SVMSMOTE
import utils
from MakeFile.FileSaver import FileSaver


if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    # filter out the CpG sites
    _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"])
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

    features = _trainX.columns
    # filter the features
    rf = RandomForestClassifier(n_estimators=300, n_jobs=-1)
    selector = RFE(rf, n_features_to_select=50, step=1, verbose=1)
    selector.fit(_trainX, _trainY)
    features = [features[i] for i in range(len(features)) if selector.support_[i]]

    xgb = XGBClassifier(n_estimators=300, max_depth=4, min_child_weight=2, subsample=0.9, colsample_bytree=0.8, reg_alpha=0.1)
    selector = RFE(xgb, n_features_to_select=50, step=1, verbose=1)
    selector.fit(_trainX, _trainY)
    features = [features[i] for i in range(len(features)) if selector.support_[i]]

    print("Selected features num: ", len(features))

    # train the model
    _trainX, _trainY, _testX, _testY = utils.SpliteTrainTest(_trainDf, _testDf, 25, features)
    lr = RandomForestClassifier(n_estimators=100, n_jobs=-1)
    lr.fit(_trainX, _trainY)
    # test the model
    print("RandomForest")
    utils.TestModelPerformance(lr, _trainX, _trainY)
    utils.TestModelPerformance(lr, _testX, _testY)

    # train the model
    estimators = [
        ('rf', RandomForestClassifier(n_estimators=100, n_jobs=-1)),
        ('lr', LogisticRegression(penalty='l1', C=1.0, solver='liblinear'))
    ]
    clf = StackingClassifier(
        estimators=estimators, final_estimator=RandomForestClassifier(n_estimators=100, n_jobs=-1)
    )
    clf.fit(_trainX, _trainY)
    # test the model
    print("Stacking Classifier")
    utils.TestModelPerformance(clf, _trainX, _trainY)
    utils.TestModelPerformance(clf, _testX, _testY)